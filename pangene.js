#!/usr/bin/env k8

/**************
 * From k8.js *
 **************/

Array.prototype.delete_at = function(i) {
	for (let j = i; j < this.length - 1; ++j)
		this[j] = this[j + 1];
	--this.length;
}

function* getopt(argv, ostr, longopts) {
	if (argv.length == 0) return;
	let pos = 0, cur = 0;
	while (cur < argv.length) {
		let lopt = "", opt = "?", arg = "";
		while (cur < argv.length) { // skip non-option arguments
			if (argv[cur][0] == "-" && argv[cur].length > 1) {
				if (argv[cur] == "--") cur = argv.length;
				break;
			} else ++cur;
		}
		if (cur == argv.length) break;
		let a = argv[cur];
		if (a[0] == "-" && a[1] == "-") { // a long option
			pos = -1;
			let c = 0, k = -1, tmp = "", o;
			const pos_eq = a.indexOf("=");
			if (pos_eq > 0) {
				o = a.substring(2, pos_eq);
				arg = a.substring(pos_eq + 1);
			} else o = a.substring(2);
			for (let i = 0; i < longopts.length; ++i) {
				let y = longopts[i];
				if (y[y.length - 1] == "=") y = y.substring(0, y.length - 1);
				if (o.length <= y.length && o == y.substring(0, o.length)) {
					k = i, tmp = y;
					++c; // c is the number of matches
					if (o == y) { // exact match
						c = 1;
						break;
					}
				}
			}
			if (c == 1) { // find a unique match
				lopt = tmp;
				if (pos_eq < 0 && longopts[k][longopts[k].length-1] == "=" && cur + 1 < argv.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				}
			}
		} else { // a short option
			if (pos == 0) pos = 1;
			opt = a[pos++];
			let k = ostr.indexOf(opt);
			if (k < 0) {
				opt = "?";
			} else if (k + 1 < ostr.length && ostr[k+1] == ":") { // requiring an argument
				if (pos >= a.length) {
					arg = argv[cur+1];
					argv.delete_at(cur + 1);
				} else arg = a.substring(pos);
				pos = -1;
			}
		}
		if (pos < 0 || pos >= argv[cur].length) {
			argv.delete_at(cur);
			pos = 0;
		}
		if (lopt != "") yield { opt: `--${lopt}`, arg: arg };
		else if (opt != "?") yield { opt: `-${opt}`, arg: arg };
		else yield { opt: "?", arg: "" };
	}
}

function* k8_readline(fn) {
	let buf = new Bytes();
	let file = new File(fn);
	while (file.readline(buf) >= 0) {
		yield buf.toString();
	}
	file.close();
	buf.destroy();
}

/*******
 * GFA *
 *******/

class GFA {
	constructor() {
		this.seg = [], this.arc = [], this.segname = {}, this.idx = [], this.walk = [], this.err = 0;
	}
	#seg_add(name) {
		if (name in this.segname) {
			return this.segname[name];
		} else {
			const sid = this.seg.length;
			this.segname[name] = sid;
			this.seg.push({ name:name, len:-1, sname:null, soff:-1, rank:-1 });
			return sid;
		}
	}
	#index() {
		const n_vtx = this.seg.length * 2;
		for (let v = 0; v < n_vtx; ++v)
			this.idx[v] = { o:0, n:0 };
		this.arc = this.arc.sort(function(a,b) { return a.v - b.v });
		for (let i = 1, st = 0; i <= this.arc.length; ++i)
			if (i == this.arc.length || this.arc[i].v != this.arc[st].v)
				this.idx[this.arc[st].v] = { o:st, n:i-st }, st = i;
		// reorder such that rank==0 is the first
		for (let v = 0; v < n_vtx; ++v) {
			const ov = this.idx[v].o;
			const nv = this.idx[v].n;
			let i0 = -1, n0 = 0;
			for (let i = 0; i < nv; ++i)
				if (this.arc[ov + i].rank == 0)
					++n0, i0 = i;
			if (n0 > 1) this.err |= 2;
			if (i0 > 0) { // then swap [0] and [i0]
				const tmp = this.arc[ov];
				this.arc[ov] = this.arc[ov + i0];
				this.arc[ov + i0] = tmp;
			}
		}
	}
	#parse_S(line) {
		const t = line.split("\t");
		if (t.length < 3) return;
		const sid = this.#seg_add(t[1]);
		let s = this.seg[sid];
		if (t[2] != "*") s.len = t[2].length;
		for (let j = 3; j < t.length; ++j) {
			let m;
			if ((m = /^(LN:i|SN:Z|SO:i|SR:i):(\S+)/.exec(t[j])) == null) continue; // TODO: parse other tags
			if (m[1] == "LN:i") s.len = parseInt(m[2]);
			else if (m[1] == "SN:Z") s.sname = m[2];
			else if (m[1] == "SO:i") s.soff = parseInt(m[2]);
			else if (m[1] == "SR:i") s.rank = parseInt(m[2]);
		}
	}
	#parse_L(line) {
		const t = line.split("\t");
		if (t.length < 5) return;
		if (t[2] != '+' && t[2] != '-') return;
		if (t[4] != '+' && t[4] != '-') return;
		const sid1 = this.#seg_add(t[1]);
		const sid2 = this.#seg_add(t[3]);
		const v = sid1 * 2 | (t[2] == '+'? 0 : 1);
		const w = sid2 * 2 | (t[4] == '+'? 0 : 1);
		let m, ov = 0, ow = 0, rank = -1;
		if (t.length >= 6) {
			const re_cigar = /(\d+)([MIDSN])/g;
			while ((m = re_cigar.exec(t[5])) != null) {
				if (m[2] == 'M' || m[2] == 'D' || m[2] == 'N') ov += parseInt(m[1]);
				if (m[2] == 'M' || m[2] == 'I' || m[2] == 'S') ow += parseInt(m[1]);
			}
			for (let j = 6; j < t.length; ++j)
				if ((m = /^(SR:i):(\S+)/.exec(t[j])) != null) // TODO: parse other tags
					rank = parseInt(m[2]);
		}
		this.arc.push({ v:v,   w:w,   ov:ov, ow:ow, rank:rank, ori:true });
		this.arc.push({ v:w^1, w:v^1, ov:ow, ow:ov, rank:rank, ori:false });
	}
	#parse_W(line) {
		const t = line.split("\t");
		if (t.length < 7) return;
		const re_walk = /([><])([^\s><]+)/g;
		let m, walk = { asm:t[1]+"#"+t[2], sample:t[1], hap:parseInt(t[2]), sname:t[3], st:-1, en:-1, v:[], lof:[] };
		if (t[4] != "*") walk.st = parseInt(t[4]);
		if (t[5] != "*") walk.st = parseInt(t[5]);
		while ((m = re_walk.exec(t[6])) != null) {
			if (this.segname[m[2]] != null) {
				const sid = this.segname[m[2]];
				const v = sid * 2 | (m[1] == '>'? 0 : 1);
				walk.v.push(v);
			}
		}
		for (let k = 7; k < t.length; ++k) {
			if (/^lf:B:i/.test(t[k])) {
				let s = t[k].substr(7).split(",");
				for (let j = 0; j < s.length; ++j)
					s[j] = parseInt(s[j]);
				walk.lof = s;
			}
		}
		this.walk.push(walk);
	}
	#parse_line(line) {
		if (line[0] == 'S') this.#parse_S(line);
		else if (line[0] == 'L') this.#parse_L(line);
		else if (line[0] == 'W') this.#parse_W(line);
	}
	toString() {
		let lines = [];
		for (let i = 0; i < this.seg.length; ++i) {
			const s = this.seg[i];
			let t = ['S', s.name, '*'];
			if (s.len >= 0) t.push(`LN:i:${s.len}`);
			if (s.sname != null && s.soff >= 0)
				t.push(`SN:Z:${s.sname}`, `SO:i:${s.soff}`);
			if (s.rank >= 0) t.push(`SR:i:${s.rank}`);
			lines.push(t.join("\t"));
		}
		for (let i = 0; i < this.arc.length; ++i) {
			const a = this.arc[i];
			if (!a.ori) continue;
			let t = ['L', this.seg[a.v>>1].name, a.v&1? '-' : '+', this.seg[a.w>>1].name, a.w&1? '-' : '+'];
			if (a.ov == 0 && a.ow == 0) t.push('0M');
			else t.push(a.ov + ':' + a.ow);
			if (a.rank >= 0) t.push(`SR:i:${a.rank}`);
			lines.push(t.join("\t"));
		}
		// TODO: output W-lines
		return lines.join("\n");
	}
	from_string(str) {
		for (const line of str.split("\n"))
			this.#parse_line(line);
		this.#index();
	}
	from_file(fn) {
		for (const line of k8_readline(fn))
			this.#parse_line(line);
		this.#index();
	}
}

/********************************
 * Intrusive double linked list *
 ********************************/

class LinkedList {
	constructor() {
		this.size = 0;
		this.head = null;
		this.tail = null;
	}
	push(node) { // node MUST have .prev and .next
		if (this.head == null && this.tail == null) {
			this.head = this.tail = node;
		} else {
			this.tail.next = node;
			node.prev = this.tail;
			this.tail = node;
		}
		++this.size;
	}
	push_list(list) {
		if (list.head == null && list.tail == null) return;
		if (this.head == null && this.tail == null) {
			this.head = list.head;
			this.tail = list.tail;
		} else {
			this.tail.next = list.head;
			list.head.prev = this.tail;
			this.tail = list.tail;
		}
		this.size += list.size;
	}
	delete(node) {
		if (this.head == node && this.tail == node) {
			this.head = this.tail = null;
		} else if (this.tail == node) {
			this.tail = node.prev, this.tail.next = null;
		} else if (this.head == node) {
			this.head = node.next, this.head.prev = null;
		} else {
			node.prev.next = node.next;
			node.next.prev = node.prev;
		}
		--this.size;
	}
}

/**************************
 * Program Structure Tree *
 **************************/

class BackEdgeNode {
	constructor(a) {
		this.a = a;
		this.recent_size = -1;
		this.recent_cec = -1;
		this.prev = null;
		this.next = null;
	}
}

class SegEdgeGraph {
	constructor() {
		this.n_node = 0;
		this.n_seg = 0;
		this.end_cat = [];
		this.arc = [];
		this.idx = [];
		this.dfs_dis = [];
		this.dfs_fin = [];
		this.dfs_par = [];
	}
	from_gfa(g) {
		this.n_seg = g.seg.length;
		const n_vtx = g.seg.length * 2;
		// collect tips
		let tip = [];
		for (let v = 0; v < n_vtx; ++v)
			if (g.idx[v].n == 0)
				tip.push(v^1);
		// collect "link" edges
		let a = [];
		for (let v = 0; v < n_vtx; ++v) {
			const n = g.idx[v].n, off = g.idx[v].o;
			for (let i = 0; i < n; ++i)
				a.push([v^1, g.arc[off + i].w]);
		}
		a.sort(function(x, y) { return x[0] - y[0] });
		// index a[]
		let idx = [];
		for (let v = 0; v < n_vtx; ++v)
			idx[v] = { o:0, n:0 };
		for (let i = 1, i0 = 0; i <= a.length; ++i)
			if (i == a.length || a[i0][0] != a[i][0])
				idx[a[i0][0]] = { o:i0, n:i-i0 }, i0 = i;
		// connected components from a[]
		let x = 0;
		for (let v = 0; v < n_vtx; ++v)
			this.end_cat[v] = -1;
		for (let v = 0; v < n_vtx; ++v) {
			if (this.end_cat[v] >= 0) continue;
			let stack = [v];
			while (stack.length > 0) { // a DFS
				const w = stack.pop();
				this.end_cat[w] = x;
				const n = idx[w].n, off = idx[w].o;
				for (let i = 0; i < n; ++i) {
					const u = a[off + i][1];
					if (this.end_cat[u] < 0) {
						this.end_cat[u] = x;
						stack.push(u);
					} else if (this.end_cat[u] != x) {
						throw Error("Wrong!");
					}
				}
			}
			++x;
		}
		this.n_node = x;
		// generate the graph
		this.arc = [];
		for (let i = 0; i < g.seg.length; ++i) {
			this.arc.push({ v:this.end_cat[i*2],   w:this.end_cat[i*2|1], seg:i, ori:1,  pair:-1, cec:-1, dfs_type:0 });
			this.arc.push({ v:this.end_cat[i*2|1], w:this.end_cat[i*2],   seg:i, ori:-1, pair:-1, cec:-1, dfs_type:0 });
		}
		// add super node
		if (tip.length > 0) {
			const super_node = this.n_node++;
			let seg_id = g.seg.length;
			for (const v of tip) {
				this.arc.push({ v:super_node, w:this.end_cat[v], seg:seg_id, ori:1,  pair:-1, cec:-1, dfs_type:0 });
				this.arc.push({ v:this.end_cat[v], w:super_node, seg:seg_id, ori:-1, pair:-1, cec:-1, dfs_type:0 });
				++seg_id;
			}
		}
		// index arc[]
		for (let i = 0; i < this.n_node; ++i)
			this.idx[i] = { n:0, o:0 };
		this.arc.sort(function(x, y) { return x.v - y.v });
		for (let i = 1, i0 = 0; i <= this.arc.length; ++i)
			if (i == this.arc.length || this.arc[i0].v != this.arc[i].v)
				this.idx[this.arc[i0].v] = { o:i0, n:i-i0 }, i0 = i;
		// populate arc[].pair
		let vtx2arc = [];
		for (let v = 0; v < n_vtx; ++v)
			vtx2arc[v] = -1;
		for (let a = 0; a < this.arc.length; ++a) {
			const arc = this.arc[a];
			if (arc.ori > 0) vtx2arc[arc.seg*2] = a;
			else vtx2arc[arc.seg*2+1] = a;
		}
		for (let a = 0; a < this.arc.length; ++a) {
			const arc = this.arc[a];
			if (arc.ori > 0) arc.pair = vtx2arc[arc.seg*2+1];
			else arc.pair = vtx2arc[arc.seg*2];
		}
	}
	print_graph(g) {
		for (let i = 0; i < this.arc.length; ++i) {
			const x = this.arc[i];
			print('E', i, `${x.v} -- ${x.w}`, x.seg < g.seg.length? g.seg[x.seg].name : '*', x.pair);
		}
	}
	dfs_traverse1(v, t, state) {
		if (state[v] != 0) return;
		this.dfs_dis[v] = t.dis++;
		state[v] = 2; // in stack
		let stack = [[v, 0]];
		while (stack.length > 0) {
			const [w, i] = stack.pop();
			const n = this.idx[w].n, off = this.idx[w].o;
			if (i < n) {
				let a = this.arc[off + i];
				stack.push([w, i + 1]); // repush to the stack
				if (a.dfs_type == 3) continue;
				const u = a.w;
				if (state[u] == 0) { // not visited before
					state[u] = 2; // in stack
					this.dfs_dis[u] = t.dis++;
					this.dfs_par[u] = w;
					stack.push([u, 0]);
					a.dfs_type = 1; // a tree edge
					this.arc[a.pair].dfs_type = 3; // wont' traverse this edge
				} else if (state[u] == 2) {
					a.dfs_type = 2; // a back edge
					this.arc[a.pair].dfs_type = 3;
				}
			} else {
				state[w] = 1; // out of stack
				this.dfs_fin[w] = t.fin++;
			}
		}
	}
	dfs_traverse() {
		for (let v = 0; v < this.n_node; ++v)
			this.dfs_dis[v] = this.dfs_fin[v] = this.dfs_par[v] = -1;
		let t = { dis:0, fin:0 }, state = [];
		for (let v = 0; v < this.n_node; ++v) state[v] = 0; // not visited
		this.dfs_traverse1(this.n_node - 1, t, state); // we can traverse every node due to super node
		for (let v = 0; v < this.n_node; ++v)
			if (state[v] == 0)
				this.dfs_traverse1(v, t, state);
		if (t.dis != this.n_node || t.fin != this.n_node)
			throw Error("DFS bug");
	}
	dfs_debug(g) {
		this.dfs_traverse();
		let v_dis = [];
		for (let v = 0; v < this.dfs_dis.length; ++v)
			v_dis[this.dfs_dis[v]] = v;
		for (let j = 0; j < v_dis.length; ++j) {
			const v = v_dis[j];
			const n = this.idx[v].n, off = this.idx[v].o;
			for (let i = 0; i < n; ++i) {
				const a = this.arc[off + i];
				if (a.dfs_type == 1 || a.dfs_type == 2)
					print(`A${a.dfs_type} ${v} -> ${a.w}`);
			}
		}
	}
	dfs_pst1(v, visited, cec_entry, sese) {
		if (visited[v] != 0) return;
		visited[v] = 1;
		let stack = [[v, 0, -1]];
		while (stack.length > 0) {
			const [w, i, b] = stack.pop();
			const n = this.idx[w].n, off = this.idx[w].o;
			if (i == n) continue;
			stack.push([w, i + 1, b]); // repush to the stack
			let a = this.arc[off + i];
			if (a.dfs_type == 3) continue;
			const u = a.w;
			let b2 = b;
			if (a.cec >= 0) {
				if (cec_entry[a.cec] != -1) // if there is a start, close it
					sese[cec_entry[a.cec]].en = off + i;
				sese.push({ st:off+i, en:-1, par:b, unflt:-1, i:-1 }); // create a new entry
				b2 = cec_entry[a.cec] = sese.length - 1;
			}
			if (visited[u] != 0) continue;
			visited[u] = 1;
			stack.push([u, 0, b2]);
		}
	}
	pst(g) {
		this.dfs_traverse();

		// put vertices in the order of their discovery time
		let v_dis = [];
		for (let v = 0; v < this.dfs_dis.length; ++v)
			v_dis[this.dfs_dis[v]] = v;

		// vs[v] keeps track of earliest node, bracket list and back edges ending at v
		let vs = [];
		for (let v = 0; v < this.n_node; ++v)
			vs[v] = { hi:this.n_node, blist:null, be_end:[], be_end_cap:[] };

		// find cycle equivalent class
		let cec = 0; // cycle equivalent class
		for (let t = v_dis.length - 1; t >= 0; --t) {
			const v = v_dis[t];
			const n = this.idx[v].n, off = this.idx[v].o;

			// compute hi0, the earliest discovery time among back edges
			let hi0 = this.n_node;
			for (let i = 0; i < n; ++i) { // traverse back edges
				if (this.arc[off + i].dfs_type !== 2) continue;
				const w = this.arc[off + i].w;
				hi0 = hi0 < this.dfs_dis[w]? hi0 : this.dfs_dis[w];
			}

			// compute hi1 and hi2, the earliest and the second earliest time among descendants
			let hi1 = this.n_node, hi2 = this.n_node;
			let blist = new LinkedList(); // initial bracket list
			for (let i = 0; i < n; ++i) { // traverse tree edges
				if (this.arc[off + i].dfs_type !== 1) continue;
				const w = this.arc[off + i].w;
				if (hi1 > vs[w].hi) hi2 = hi1, hi1 = vs[w].hi;
				else if (hi2 > vs[w].hi) hi2 = vs[w].hi;
				blist.push_list(vs[w].blist); // merge blists from v's children
			}
			vs[v].hi = hi0 < hi1? hi0 : hi1;

			// compute the final bracket list
			for (const b of vs[v].be_end_cap) // delete capping back edges ending at v
				blist.delete(b);
			for (const b of vs[v].be_end) { // delete (normal) back edges ending at v
				blist.delete(b);
				if (this.arc[b.a].cec < 0)
					this.arc[b.a].cec = cec++;
			}
			for (let i = 0; i < n; ++i) { // traverse back edges starting at v
				if (this.arc[off + i].dfs_type != 2) continue;
				const w = this.arc[off + i].w;
				if (w === v) continue; // no loop!
				const e = new BackEdgeNode(off + i);
				blist.push(e);
				vs[w].be_end.push(e);
			}
			if (hi2 < hi0) { // then create a capping back edge
				const w = v_dis[hi2];
				if (w != v) { // no loop!
					const d = new BackEdgeNode(-1); // capping back edge
					blist.push(d);
					vs[w].be_end_cap.push(d);
				}
			}
			vs[v].blist = blist;

			if (0) {
				let l = [];
				for (let p = blist.head; p != null; p = p.next) {
					if (p.a < 0) l.push('*');
					else l.push(`${this.arc[p.a].v}:${this.arc[p.a].w}`);
				}
				print('X', v, `hi0=${hi0},hi1=${hi1},hi2=${hi2}`, l.join(","));
			}

			// determine the category for tree edge (parent(v),v)
			if (this.dfs_par[v] >= 0 && blist.size > 0) { // not a root (there may be multiple roots if the graph is disconnected)
				const u = this.dfs_par[v]; // v's parent
				const n = this.idx[u].n, off = this.idx[u].o;
				let e = -1; // the tree edge from u to v
				for (let i = 0; i < n; ++i)
					if (this.arc[off + i].w === v && this.arc[off + i].dfs_type === 1)
						e = off + i;
				if (e < 0) throw Error(`Bug: failed to find tree edge ${u}->${v}`);
				const b = blist.tail;
				if (b.recent_size !== blist.size) {
					b.recent_size = blist.size;
					b.recent_cec = cec++;
				}
				if (b.recent_cec < 0) throw Error(`Bug: recent_cec not set when processing edge ${e}`);
				this.arc[e].cec = b.recent_cec;
				if (b.recent_size === 1 && b.a >= 0) // the tree edge e and back edge b.a are equivalent
					this.arc[b.a].cec = this.arc[e].cec;
			}
		}
		if (0) {
			for (let i = 0; i < this.arc.length; ++i) {
				const a = this.arc[i];
				if (a.dfs_type == 1 || a.dfs_type == 2)
					print(`${a.v} -> ${a.w}`, a.seg < g.seg.length? g.seg[a.seg].name : "*", a.dfs_type, a.cec);
			}
		}

		// construct initial PST
		let state = [], sese = [], cec_entry = [];
		for (let v = 0; v < this.n_node; ++v) state[v] = 0; // not visited
		for (let c = 0; c < cec; ++c) cec_entry[c] = -1;
		for (let t = 0; t < v_dis.length; ++t) {
			const v = v_dis[t];
			if (state[v] == 0)
				this.dfs_pst1(v, state, cec_entry, sese);
		}

		// filter out open bubbles and point bubbles
		let sese_flt = [];
		for (let i = 0; i < sese.length; ++i) {
			let b = sese[i], flt = false;
			if (b.en < 0) flt = true; // an open bubble
			else if (this.arc[b.st].seg >= this.n_seg || this.arc[b.en].seg >= this.n_seg) flt = true; // involving the dummy node
			else if (this.arc[b.st].w == this.arc[b.en].v && this.idx[this.arc[b.en].v].n == 2) flt = true; // a point bubble
			if (flt) {
				if (b.par >= 0) b.unflt = sese[b.par].unflt;
				else b.unflt = -1;
			} else {
				b.unflt = i, b.par = sese[b.par].unflt;
				b.i = sese_flt.length;
				const par = b.par < 0? -1 : sese[b.par].i;
				sese_flt.push({ st:b.st, en:b.en, par:par });
			}
		}
		return sese_flt;
	}
	print_pst(sese, g) {
		for (let i = 0; i < sese.length; ++i) {
			const st = this.arc[sese[i].st].seg * 2 + (this.arc[sese[i].st].ori > 0? 0 : 1);
			const en = sese[i].en < 0? "*" : this.arc[sese[i].en].seg * 2 + (this.arc[sese[i].en].ori > 0? 0 : 1);
			print(i, sese[i].par, "><"[st&1] + g.seg[st>>1].name, "><"[en&1] + g.seg[en>>1].name);
		}
	}
	print_bandage_csv(g) {
		print("segment,label");
		for (let i = 0; i < this.arc.length; ++i) {
			const a = this.arc[i];
			if (a.seg < g.seg.length && (a.dfs_type == 1 || a.dfs_type == 2))
				print(`${g.seg[a.seg].name},${a.cec}`);
		}
	}
}

/***************
 * Subcommands *
 ***************/

function pg_cmd_parse_gfa(args) {
	if (args.length == 0) {
		print("Usage: pangene.js parse-gfa <in.gfa>");
		return;
	}
	let g = new GFA();
	g.from_file(args[0]);
	print(g);
}

function pg_cmd_test(args) {
	if (args.length == 0) {
		print("Usage: pangene.js parse-gfa <in.gfa>");
		return;
	}
	let g = new GFA();
	g.from_file(args[0]);
	let e = new SegEdgeGraph();
	e.from_gfa(g);
	//e.print_graph(g);
	//e.dfs_debug();
	const sese = e.pst(g);
	e.print_pst(sese, g);
}

/*****************
 * Main function *
 *****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: pangene.js <command> [arguments]");
		print("Commands:");
		print("  parse-gfa     parse a GFA file (for debugging only)");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'parse-gfa') pg_cmd_parse_gfa(args);
	else if (cmd == 'test') pg_cmd_test(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
