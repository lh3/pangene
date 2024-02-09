#!/usr/bin/env k8

const pg_version = "r194-dirty";

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
			this.seg.push({ name:name, len:-1, sname:null, soff:-1, rank:-1, cec:-1 });
			return sid;
		}
	}
	#index() {
		const n_vtx = this.seg.length * 2;
		for (let v = 0; v < n_vtx; ++v)
			this.idx[v] = { o:0, n:0 };
		this.arc.sort(function(a,b) { return a.v - b.v });
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
		//this.arc.push({ v:w^1, w:v^1, ov:ow, ow:ov, rank:rank, ori:false });
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
	static int_hash(x) {
		x = Math.imul((x >> 16) ^ x, 0x45d9f3b) & 0xffffffff;
		x = Math.imul((x >> 16) ^ x, 0x45d9f3b) & 0xffffffff;
		return (x >> 16) ^ x;
	}
	#traverse_bubble(vs, ve, flag, f, max_n) {
		let stack = [vs], list = [];
		flag[vs] = f;
		while (stack.length) {
			const v = stack.pop();
			const off = this.idx[v].o, n = this.idx[v].n;
			for (let i = 0; i < n; ++i) {
				const a = this.arc[off + i], w = a.w;
				if (w == (vs^1)) continue; // don't pass the starting vertex on either strand
				if (w == (ve^1)) return []; // if reaching the reverse complement of ve, there is no bubble
				if (flag[w] != f) {
					flag[w] = f;
					if (w == ve) continue; // reaching the end vertex; flag it but don't add to list[] or stack[]
					if (flag[w^1] != f) list.push(w>>1);
					stack.push(w);
				}
			}
			if (list.length > max_n) break;
		}
		return list.length > max_n? [] : list;
	}
	get_bubble_id(vs, ve, flag, f, max_n) {
		let f_for = f, f_rev = f + this.seg.length * 2;
		let list_for = this.#traverse_bubble(vs,   ve,   flag, f_for, max_n);
		let list_rev = this.#traverse_bubble(ve^1, vs^1, flag, f_rev, max_n);
		if (list_for.length != list_rev.length) return [];
		let n_in = 0;
		for (let i = 0; i < list_for.length; ++i)
			if (flag[list_for[i]<<1|0] == f_rev || flag[list_for[i]<<1|1] == f_rev)
				++n_in;
		if (n_in != list_for.length) return [];
		for (let i = 0; i < list_for.length; ++i) {
			for (let rev = 0; rev < 2; ++rev) {
				const v = list_for[i]<<1 | rev;
				const off = this.idx[v].o, n = this.idx[v].n;
				for (let j = 0; j < n; ++j) {
					const a = this.arc[off + j];
					if (flag[a.w] != f_for && flag[a.w] != f_rev) // if reaching other vertices not reachable from start or end, there is no bubble
						return [];
				}
			}
		}
		return list_for;
	}
	get_bubble(vs, ve, flag, f, max_n) {
		let b = [], a = this.get_bubble_id(vs, ve, flag, f, max_n);
		for (let i = 0; i < a.length; ++i)
			b[i] = this.seg[a[i]].name;
		return b;
	}
	get_bubble_all(max_ext) {
		// initialize cec_par[]
		let max_cec = -1, cec_par = []; // parent cec
		for (let i = 0; i < this.seg.length; ++i)
			max_cec = max_cec > this.seg[i].cec? max_cec : this.seg[i].cec;
		for (let i = 0; i <= max_cec; ++i)
			cec_par[i] = -1;

		// initialize flag arrays
		const n_vtx = this.seg.length * 2;
		let f1 = 0, f2 = 0;
		let flag1 = [], flag2 = []; // flag1 for end finding; flag2 for get_bubble()
		for (let v = 0; v < n_vtx; ++v)
			flag1[v] = flag2[v] = -1;

		// look for bubbles
		let bb = [];
		for (let vs = 0; vs < n_vtx; ++vs) {
			// test if vs needs to be checked
			const cec = this.seg[vs>>1].cec;
			if (cec < 0 || this.idx[vs].n == 0) continue;
			if (this.idx[vs].n == 1) {
				const w = this.arc[this.idx[vs].o].w ^ 1;
				if (this.idx[w].n < 2) continue;
			}

			// BFS starting from vs
			let queue = [vs], ve = [], ext = 0;
			flag1[vs] = f1;
			while (queue.length) { // this loop is similar to the loop in #traverse_bubble()
				const v = queue.shift();
				const off = this.idx[v].o, n = this.idx[v].n;
				for (let i = 0; i < n; ++i) {
					const a = this.arc[off + i], w = a.w;
					if (flag1[w] != f1) {
						if (flag1[w^1] != f1) ++ext;
						if (w == (vs^1)) continue;
						flag1[w] = f1;
						if (this.seg[w>>1].cec == cec) {
							ve.push(w);
							continue;
						}
						queue.push(w);
					}
				}
				if (ext > max_ext) break;
			}

			// test potential bubbles
			for (let i = 0; i < ve.length; ++i) {
				let r = this.get_bubble_id(vs, ve[i], flag2, f2, max_ext);
				if (r.length > 0 && vs < ve[i]) {
					let name_list = [];
					for (let j = 0; j < r.length; ++j) {
						if (this.seg[r[j]].cec >= 0)
							cec_par[this.seg[r[j]].cec] = cec;
						name_list.push(this.seg[r[j]].name);
					}
					bb.push({ cec:cec, cec_par:-1, vs:vs, ve:ve[i], list:name_list });
				}
				++f2;
			}
			++f1;
		}

		// update parent cec
		for (let i = 0; i < bb.length; ++i)
			bb[i].cec_par = cec_par[bb[i].cec];
		for (let i = 0; i < bb.length; ++i)
			print("BB", bb[i].cec, bb[i].cec_par, bb[i].cec, "><"[bb[i].vs&1] + this.seg[bb[i].vs>>1].name, "><"[bb[i].ve&1] + this.seg[bb[i].ve>>1].name, bb[i].list.length, bb[i].list.join(","));
		return bb;
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

class NetGraph {
	constructor(g, ref) {
		this.n_node = 0;
		this.end_cat = [];
		this.arc = [];
		this.idx = [];
		this.dfs_dis = [];
		this.dfs_fin = [];
		this.dfs_par = [];
		this.gfa = g;
		this.ref = ref;
		this.#convert_gfa();
	}
	#convert_gfa() {
		const g = this.gfa;
		const n_vtx = g.seg.length * 2;
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
		// collect tips
		let tip = [];
		for (let v = 0; v < n_vtx; ++v)
			if (g.idx[v].n == 0)
				tip.push(v^1);
		if (this.ref && g.walk.length > 0) {
			let f = [];
			for (let v = 0; v < n_vtx; ++v) f[v] = 0;
			for (let i = 0; i < tip.length; ++i) f[tip[i]] = 1;
			for (const w of g.walk) {
				if (w.asm != this.ref || w.v.length < 2) continue;
				let t1 = w.v[0], t2 = w.v[w.v.length-1]^1;
				if (f[t1] == 0) f[t1] = 2;
				if (f[t2] == 0) f[t2] = 2;
			}
			for (let v = 0; v < n_vtx; ++v)
				if (f[v] == 2)
					tip.push(v);
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
	dfs_pst1(v, visited, cec_entry, sese) { // compute SESEs and their hierarchy; this is different from Johnson's PhD thesis
		if (visited[v] != 0) return;
		visited[v] = 1;
		let stack = [[v, 0, -1]];
		while (stack.length > 0) {
			const [w, i, b] = stack.pop(); // b: the index of the bubble that leads to w
			const n = this.idx[w].n, off = this.idx[w].o;
			if (i == n) continue;
			stack.push([w, i + 1, b]); // repush to the stack
			let a = this.arc[off + i];
			if (a.dfs_type == 3) continue; // blocked edge
			const u = a.w;
			let b2 = b;
			if (a.cec >= 0) {
				let par = b;
				if (cec_entry[a.cec] != -1) // if there is a start, close it
					sese[cec_entry[a.cec]].en = off + i, par = sese[cec_entry[a.cec]].par;
				sese.push({ st:off+i, en:-1, par:par, unflt:-1, i:-1 }); // create a new entry with the same parrent
				b2 = cec_entry[a.cec] = sese.length - 1;
			}
			if (visited[u] != 0) continue;
			visited[u] = 1;
			stack.push([u, 0, b2]);
		}
	}
	#dbg_blist(blist) {
		let l = [];
		for (let p = blist.head; p != null; p = p.next) {
			if (p.a < 0) l.push('*');
			else l.push(`${this.arc[p.a].v},${this.arc[p.a].w}`);
		}
		print('X', v, `hi0=${hi0},hi1=${hi1},hi2=${hi2}`, l.join(";"));
	}
	mark_cec() {
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
		let cec = 1; // cycle equivalent class; class 0 is reserved for tree edges not in cycles
		for (let t = v_dis.length - 1; t >= 0; --t) {
			const v = v_dis[t];
			const n = this.idx[v].n, off = this.idx[v].o;

			// compute hi0, the earliest discovery time among back edges
			let hi0 = this.n_node;
			for (let i = 0; i < n; ++i) { // traverse back edges
				if (this.arc[off + i].dfs_type !== 2) continue;
				const w = this.arc[off + i].w;
				if (v === w) continue;
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
			if (hi2 < hi0 && hi2 < t) { // then create a capping back edge; this line is different from Johnson et al
				const w = v_dis[hi2];
				const d = new BackEdgeNode(-1); // capping back edge
				blist.push(d);
				vs[w].be_end_cap.push(d);
			}
			vs[v].blist = blist;
			if (0) this.#dbg_blist(blist);

			// determine the category for tree edge (parent(v),v)
			if (this.dfs_par[v] >= 0) { // not a root (there may be multiple roots if the graph is disconnected)
				const u = this.dfs_par[v]; // v's parent
				const n = this.idx[u].n, off = this.idx[u].o;
				let e = -1; // the tree edge from u to v
				for (let i = 0; i < n; ++i)
					if (this.arc[off + i].w === v && this.arc[off + i].dfs_type === 1)
						e = off + i;
				if (e < 0) throw Error(`Bug: failed to find tree edge ${u}->${v}`);
				if (blist.size > 0) {
					const b = blist.tail;
					if (b.recent_size !== blist.size) {
						b.recent_size = blist.size;
						b.recent_cec = cec++;
					}
					if (b.recent_cec < 0) throw Error(`Bug: recent_cec not set when processing edge ${e}`);
					this.arc[e].cec = b.recent_cec;
					if (b.recent_size === 1 && b.a >= 0) // the tree edge e and back edge b.a are equivalent
						this.arc[b.a].cec = this.arc[e].cec;
				} else this.arc[e].cec = 0; // we won't come here given a control flow graph
			}
		}

		// mark cycle equivalence classes on the original GFA
		for (let e = 0; e < this.arc.length; ++e) {
			const a = this.arc[e];
			if (a.seg < this.gfa.seg.length && (a.dfs_type == 1 || a.dfs_type == 2))
				this.gfa.seg[a.seg].cec = a.cec;
		}
		return cec;
	}
	pst() {
		const cec = this.mark_cec();
		let v_dis = [];
		for (let v = 0; v < this.dfs_dis.length; ++v)
			v_dis[this.dfs_dis[v]] = v;

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
			else if (this.arc[b.st].seg >= this.gfa.seg.length || this.arc[b.en].seg >= this.gfa.seg.length) flt = true; // involving the dummy node
			else if (this.arc[b.st].w == this.arc[b.en].v && this.idx[this.arc[b.en].v].n == 2) flt = true; // a point bubble
			if (flt) {
				if (b.par >= 0) b.unflt = sese[b.par].unflt;
				else b.unflt = -1;
			} else {
				b.unflt = i;
				if (b.par >= 0) b.par = sese[b.par].unflt;
				b.i = sese_flt.length;
				const par = b.par < 0? -1 : sese[b.par].i;
				sese_flt.push({ st:b.st, en:b.en, par:par, vs:-1, ve:-1 });
			}
		}
		this.#cal_vs_ve(sese_flt);
		return sese_flt;
	}
	#cal_vs_ve(sese) {
		for (let i = 0; i < sese.length; ++i) {
			if (sese[i].en < 0) continue;
			sese[i].vs = this.arc[sese[i].st].seg * 2 + (this.arc[sese[i].st].ori > 0? 0 : 1);
			sese[i].ve = this.arc[sese[i].en].seg * 2 + (this.arc[sese[i].en].ori > 0? 0 : 1);
		}
	}
	print_pst(sese, max_ext) {
		const g = this.gfa;
		let flag = [];
		for (let v = 0; v < g.seg.length; ++v) flag[v] = -1;
		for (let i = 0; i < sese.length; ++i) {
			const st = sese[i].vs, en = sese[i].ve;
			const b = this.gfa.get_bubble(st, en, flag, i, max_ext);
			if (b.length == 0) {
				print('FB', i, sese[i].par, this.arc[sese[i].st].cec, "><"[st&1] + g.seg[st>>1].name, "><"[en&1] + g.seg[en>>1].name);
			} else {
				print('BB', i, sese[i].par, this.arc[sese[i].st].cec, "><"[st&1] + g.seg[st>>1].name, "><"[en&1] + g.seg[en>>1].name, 0, b.length, b.join(","));
			}
		}
	}
	print_bandage_csv() {
		const g = this.gfa;
		print("segment,label");
		for (let i = 0; i < this.arc.length; ++i) {
			const a = this.arc[i];
			if (a.seg < g.seg.length && (a.dfs_type == 1 || a.dfs_type == 2) && a.cec >= 0)
				print(`${g.seg[a.seg].name},${a.cec}`);
		}
	}
	print_dfs() { // for debugging only
		const g = this.gfa;
		if (this.dfs_dis.length == 0) this.dfs_traverse();
		let v_dis = [];
		for (let v = 0; v < this.dfs_dis.length; ++v)
			v_dis[this.dfs_dis[v]] = v;
		for (let j = 0; j < v_dis.length; ++j) {
			const v = v_dis[j];
			const n = this.idx[v].n, off = this.idx[v].o;
			for (let i = 0; i < n; ++i) {
				const a = this.arc[off + i];
				if (a.dfs_type == 1 || a.dfs_type == 2)
					print('DF', ["tree", "back"][a.dfs_type-1], `${v},${a.w}`, (a.seg < g.seg.length? "><"[a.ori>0?0:1] + g.seg[a.seg].name : "*"));
			}
		}
	}
	print_cycle_equiv() { // for debugging only
		const g = this.gfa;
		for (let i = 0; i < this.arc.length; ++i) {
			const a = this.arc[i];
			if (a.dfs_type == 1 || a.dfs_type == 2)
				print('EC', a.cec, ["tree", "back"][a.dfs_type-1], `${a.v},${a.w}`, (a.seg < g.seg.length? "><"[a.ori>0?0:1] + g.seg[a.seg].name : "*"));
		}
	}
	walk_ht(sese) {
		const g = this.gfa;
		let st = [], en = [], ht = [];
		for (let v = 0; v < g.seg.length * 2; ++v)
			st[v] = [], en[v] = { walk:-1, a:[] };
		for (let i = 0; i < sese.length; ++i) {
			if (sese[i].en < 0) continue;
			ht[i] = [];
			st[sese[i].vs].push({ en:sese[i].ve, bid:i, ori:1 });
			st[sese[i].ve^1].push({ en:sese[i].vs^1, bid:i, ori:-1 });
		}
		for (let j = 0; j < g.walk.length; ++j) {
			const vtx = g.walk[j].v;
			for (let i = 0; i < vtx.length; ++i) {
				const v = vtx[i];
				for (let k = 0; k < st[v].length; ++k) {
					let e = en[st[v][k].en];
					if (e.walk != j)
						e.walk = j, e.a = [];
					e.a.push({ st_off:i, bid:st[v][k].bid, ori:st[v][k].ori });
				}
				if (en[v].walk != j) continue;
				for (let k = 0; k < en[v].a.length; ++k) {
					const x = en[v].a[k];
					ht[x.bid].push({ walk:j, st_off:x.st_off, en_off:i, bid:x.bid, ori:x.ori });
					//print(j, x.st_off, i, x.bid, x.ori);
				}
			}
		}
		return ht;
	}
	count_allele(sese, ht, max_ext) {
		const g = this.gfa;
		for (let i = 0; i < sese.length; ++i) {
			let gene_hash = {}, gene_list = [];
			for (let j = 0; j < ht[i].length; ++j) { // get the list of genes
				const x = ht[i][j];
				const w = g.walk[x.walk];
				for (let k = x.st_off + 1; k < x.en_off; ++k) {
					const v = w.v[k];
					if (gene_hash[v>>1] == null)
						gene_hash[v>>1] = 1, gene_list.push(g.seg[v>>1].name);
				}
			}
			sese[i].n_gene = gene_list.length;
			sese[i].gene = [];
			sese[i].al = [];
			if (gene_list.length > max_ext) continue;
			sese[i].gene = gene_list;
			let al = {};
			for (let j = 0; j < ht[i].length; ++j) { // get alleles
				const x = ht[i][j];
				const w = g.walk[x.walk];
				let a = [];
				if (x.ori > 0) {
					for (let k = x.st_off; k <= x.en_off; ++k)
						a.push(w.v[k]);
				} else if (x.ori < 0) {
					for (let k = x.en_off; k >= x.st_off; --k)
						a.push(w.v[k]^1);
				}
				const s = a.join(",");
				if (al[s] == null) al[s] = { a:a.slice(0), n:0 };
				++al[s].n;
			}
			for (const key in al)
				sese[i].al.push({ n:al[key].n, a:al[key].a });
			sese[i].al.sort(function(a,b) { return b.n - a.n });
		}
	}
	print_pst_walk(sese, min_n_allele, max_ext) {
		const g = this.gfa;
		let flag = [];
		for (let v = 0; v < g.seg.length; ++v) flag[v] = -1;
		for (let i = 0; i < sese.length; ++i) {
			const vs = sese[i].vs, ve = sese[i].ve;
			const b = this.gfa.get_bubble(vs, ve, flag, i, max_ext);
			if (b.length == 0) {
				print('FB', i, sese[i].par, this.arc[sese[i].st].cec, "><"[vs&1] + g.seg[vs>>1].name, "><"[ve&1] + g.seg[ve>>1].name);
				print('//');
				continue;
			}
			const gene = sese[i].gene;
			const gene_list = gene.length == 0? sese[i].n_gene : `${gene.length}\t${gene.join(",")}`;
			if (sese[i].al.length < min_n_allele) continue;
			print('BB', i, sese[i].par, this.arc[sese[i].st].cec, "><"[vs&1] + g.seg[vs>>1].name, "><"[ve&1] + g.seg[ve>>1].name, sese[i].al.length, gene_list);
			for (let j = 0; j < sese[i].al.length; ++j) {
				let a = [];
				for (let k = 0; k < sese[i].al[j].a.length; ++k) {
					const v = sese[i].al[j].a[k];
					a.push("><"[v&1], g.seg[v>>1].name);
				}
				print('AL', sese[i].al[j].n, a.join(""));
			}
			print('//');
		}
	}
}

/***************
 * Subcommands *
 ***************/

function pg_cmd_call(args) {
	let opt = { print_pst:true, print_bandage:false, print_cec:false, print_dfs:false, max_ext:100, ignore_walk:false, use_algo2:false, min_n_allele:2, ref:null };
	for (const o of getopt(args, "bedamwc:r:", [])) {
		if (o.opt == "-b") opt.print_bandage = true, opt.print_pst = false;
		else if (o.opt == "-e") opt.print_cec = true, opt.print_pst = false;
		else if (o.opt == "-d") opt.print_dfs = true, opt.print_pst = false;
		else if (o.opt == "-m") opt.max_ext = parseInt(o.arg);
		else if (o.opt == "-w") opt.ignore_walk = true;
		else if (o.opt == "-c") opt.min_n_allele = parseInt(o.arg);
		else if (o.opt == "-r") opt.ref = o.arg;
		else if (o.opt == "-a") opt.use_algo2 = true;
	}
	if (args.length == 0) {
		print("Usage: pangene.js call [options] <in.gfa>");
		print("Options:");
		print("  General:");
		print(`    -m INT   don't output gene lists longer than INT [${opt.max_ext}]`);
		print(`    -c INT   min number of alleles [${opt.min_n_allele}]`);
		print("    -b       output equivalent classes for Bandage visualization");
		print("    -r INT   reference assembly []");
		print("    -w       ignore walks");
		print("  Debugging:");
		print("    -a       alternative algorithm for bubble finding (no alleles)");
		print("    -d       output DFS traversal");
		print("    -e       output cycle equivalent class");
		return;
	}
	let g = new GFA();
	g.from_file(args[0]);
	let e = new NetGraph(g, opt.ref);
	const sese = e.pst();
	if (opt.print_dfs) e.print_dfs();
	if (opt.print_bandage) e.print_bandage_csv();
	if (opt.print_cec) e.print_cycle_equiv();
	if (opt.print_pst) {
		print("CC", "FB  bbID  parID  side1  side2");
		print("CC", "BB  bbID  parID  side1  side2  #alleles  #genes  geneList");
		print("CC", "AL  #hap  walk");
		print("CC");
		if (opt.use_algo2) {
			g.get_bubble_all(opt.max_ext);
		} else if (!opt.ignore_walk && g.walk.length > 0) {
			let ht = e.walk_ht(sese);
			e.count_allele(sese, ht, opt.max_ext);
			e.print_pst_walk(sese, opt.min_n_allele, opt.max_ext);
		} else {
			e.print_pst(sese, opt.max_ext);
		}
	}
}

function pg_cmd_call2html(args) {
	let endpoint = "/view", graph = null;
	for (const o of getopt(args, "e:g:", [])) {
		if (o.opt == "-e") endpoint = o.arg;
		else if (o.opt == "-g") graph = o.arg;
	}
	if (args.length == 0) {
		print("Usage: pangene.js call2html [options] <pangene-call.out>");
		print("Options:");
		print(`  -e STR     endpoint [${endpoint}]`);
		print(`  -g STR     graph name []`);
		return;
	}
	print(`<head>`);
	print(`<title>List of variants</title>`);
	print(`<style type="text/css">`);
	print(`  table { font-family: "helvetica neue", helvetica, arial, sans-serif; font-size: 0.8em; text-align: left; }`);
	print(`  th, td { padding: 2px; }`);
	print(`  a { text-decoration: none; color: blue; }`);
	print(`</style>`);
	print(`</head>`);
	print(`<body>`);
	print(`<table border="1" style="border-collapse: collapse; max-width: 1024px; width: 100%;">`);
	print(`<tr><th>VarID<th>Parent<th>#alleles<th>End genes<th>Genes</tr>`);
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		if (t[0] != "BB") continue;
		if (t.length < 9) continue;
		const st = (t[4][0] == ">"? "&gt;" : "&lt;") + t[4].substr(1);
		const en = (t[5][0] == ">"? "&gt;" : "&lt;") + t[5].substr(1);
		const genes = [t[4].substr(1), t[8], t[5].substr(1)].join(",");
		let link = `${endpoint}?`;
		if (graph != null) link += `graph=${graph}&`;
		link += `gene=${genes}&step=0&ori=` + t[4].substr(1);
		const gene_space = t[8].replace(/,/g, ", ");
		let out = `<tr><td style="text-align: right;">${t[1]}<td style="text-align: right;">${t[2]}<td style="text-align: right;">${t[6]}`;
		out += `<td style="white-space: nowrap;"><a href="${link}" target="_blank">${st} &mdash; ${en}</a><td>${gene_space}</tr>`;
		print(out);
	}
	print(`</table>`);
	print(`</body>`);
}

function pg_cmd_calldiff(args) {
	for (const o of getopt(args, "", [])) {
	}
	if (args.length < 2) {
		print("Usage: pangene.js calldiff <call1.out> <call2.out>");
		return;
	}
	let h = {};
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		if (t[0] != "BB") continue;
		const g1 = t[4].substr(1), g2 = t[5].substr(1);
		const key = g1 < g2? `${g1}\t${g2}` : `${g2}\t${g1}`;
		h[key] = [false, t.slice(1).join("\t")];
	}
	for (const line of k8_readline(args[1])) {
		let t = line.split("\t");
		if (t[0] != "BB") continue;
		const g1 = t[4].substr(1), g2 = t[5].substr(1);
		const key = g1 < g2? `${g1}\t${g2}` : `${g2}\t${g1}`;
		if (key in h)
			h[key][0] = true;
		else
			print("B2", t.slice(1).join("\t"));
	}
	for (const key in h)
		if (h[key][0] == false)
			print("B1", h[key][1]);
}

function pg_cmd_getaa(args) {
	let species = null, excl_decay = false, keep_thru = false, canon_only = false;
	for (const o of getopt(args, "s:erc", [])) {
		if (o.opt == "-s") species = o.arg;
		else if (o.opt == "-e") excl_decay = true;
		else if (o.opt == "-r") keep_thru = true;
		else if (o.opt == "-c") canon_only = true;
	}
	if (args.length < 2) {
		print("Usage: pangene.js getaa [options] <anno.gtf> <proteins.faa>");
		print("Options:");
		print("  -s STR     species name []");
		print("  -e         exclude transcripts that are not protein_coding");
		print("  -r         keep readthrough transcripts");
		print("  -c         canonical only (for GenCode)");
		return;
	}
	const re = /([^\s"]+) "([^\s"]+)"/g;
	let h = {};
	for (const line of k8_readline(args[0])) { // read symbols
		if (line[0] == '#') continue;
		let m, t = line.split("\t");
		if (t[2] !== "CDS") continue;
		if (t[0] === "MT" || t[0] === "chrM" || t[0] === "chrMT") continue;
		let gid = null, gname = null, pid = null, pver = null, ttype = null, gtype = null, thru = false, canon = false;
		while ((m = re.exec(t[8])) != null) {
			if (m[1] == "gene_id") {
				gid = m[2];
			} else if (m[1] == "protein_id") {
				pid = m[2];
			} else if (m[1] == "protein_version") {
				pver = m[2];
			} else if (m[1] == "gene_name") {
				gname = m[2];
			} else if (m[1] == "transcript_biotype" || m[1] == "transcript_type") {
				ttype = m[2];
			} else if (m[1] == "gene_biotype" || m[1] == "gene_type") {
				gtype = m[2];
			} else if (m[1] === "tag" && m[2] === "readthrough_transcript") {
				thru = true;
			} else if (m[1] === "tag" && m[2] === "Ensembl_canonical") {
				canon = true;
			}
		}
		if (gtype !== "protein_coding") continue;
		if (canon_only && !canon) continue;
		if (excl_decay && ttype !== "protein_coding") continue;
		if (!keep_thru && thru) continue;
		let gene = gname != null? gname : gid;
		if (gene == null) throw Error("failed to parse the gene name");
		if (species != null) gene = `${gene}_${species}`;
		let prot = pver != null? `${pid}.${pver}` : pid;
		h[prot] = `${gene}:${prot} ${t[0]}`;
	}
	let skip = false;
	for (const line of k8_readline(args[1])) {
		let m;
		if ((m = /^>([^\s\|]+)/.exec(line)) != null) {
			const pid = m[1];
			if (h[pid] != null) {
				print(`>${h[m[1]]}`);
				skip = false;
			} else {
				warn(`WARNING: skip "${m[1]}"`);
				skip = true;
			}
		} else if (!skip) {
			print(line);
		}
	}
}

function pg_cmd_gfa2matrix(args) {
	let copy_number = false;
	for (const o of getopt(args, "c", []))
		if (o.opt == "-c") copy_number = true;
	if (args.length == 0) {
		print("Usage: pangene.js gfa2matrix [-c] <in.gfa>");
		return;
	}
	let g = new GFA();
	g.from_file(args[0]);
	let mat = [], asm_h = {}, asm_a = [];
	for (const w of g.walk) {
		if (asm_h[w.asm] == null) {
			asm_h[w.asm] = asm_a.length;
			asm_a.push(w.asm);
		}
	}
	for (let i = 0; i < g.seg.length; ++i)
		mat[i] = Array(asm_a.length).fill(0);
	for (const w of g.walk) {
		const asm_id = asm_h[w.asm];
		for (const v of w.v)
			++mat[v>>1][asm_id];
	}
	if (copy_number == false) {
		for (let i = 0; i < mat.length; ++i)
			for (let j = 0; j < mat[i].length; ++j)
				if (mat[i][j] > 1) mat[i][j] = 1;
	}
	print('Gene', asm_a.join("\t"));
	for (let i = 0; i < mat.length; ++i)
		print(g.seg[i].name, mat[i].join("\t"));
}

function pg_cmd_flt_mmseqs(args) {
	let sim = 0.9, qonly = false, rev = false;
	for (const o of getopt(args, "s:qv", [])) {
		if (o.opt == "-s") sim = parseFloat(o.arg);
		else if (o.opt == "-q") qonly = true;
		else if (o.opt == "-v") rev = true;
	}
	if (args.length == 0) {
		print("Usage: pangene.js flt-mmseqs [-v] [-q] [-s 0.9] <mmseqs.2.txt> | cut -f1 | uniq > filtered.txt");
		return;
	}
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t");
		let flt = false;
		if (qonly) {
			t[2] = parseFloat(t[2]);
			if (t[2] < sim) continue;
			const qal = parseInt(t[7]) - parseInt(t[6]) + 1;
			const qlen = parseInt(t[12]);
			if (qal < qlen * sim) flt = true;
		} else {
			const qlen = parseInt(t[12]);
			const tlen = parseInt(t[13]);
			const qs = parseInt(t[6]) - 1, qe = parseInt(t[7]);
			const ts = parseInt(t[8]) - 1, te = parseInt(t[9]);
			const l0 = qs < ts? qs : ts;
			const l1 = tlen - te < qlen - qe? tlen - te : qlen - qe;
			const n_iden = parseInt(t[3]) - parseInt(t[4]) - parseInt(t[5]);
			if (n_iden < (l0 + l1 + parseInt(t[3])) * sim) flt = true;
		}
		if (rev) {
			if (flt) print(line);
		} else {
			if (!flt) print(line);
		}
	}
}

function pg_cmd_bed2paf(args) {
	for (const o of getopt(args, "", [])) {
	}
	if (args.length < 3) {
		print("Usage: pangene.js bed2paf <in.bed> <seq1.fai> <seq2.fai>");
		return;
	}
	let s1 = {};
	for (const line of k8_readline(args[1])) {
		const t = line.split("\t", 2);
		s1[t[0]] = parseInt(t[1]);
	}
	let s2 = {};
	for (const line of k8_readline(args[2])) {
		const t = line.split("\t", 2);
		s2[t[0]] = parseInt(t[1]);
	}
	let h = {};
	for (const line of k8_readline(args[0])) {
		let t = line.split("\t", 6);
		let k = (t[0] in s1)? 0 : (t[0] in s2)? 1 : -1;
		if (k < 0) continue;
		t.push(k);
		let g = t[3].split(":")[0];
		t[4] = parseInt(t[4]);
		if (h[g] == null) h[g] = [];
		h[g].push(t);
	}
	for (const g in h) {
		if (h[g].length != 2) continue;
		const a = h[g];
		let n = [0, 0];
		for (let i = 0; i < a.length; ++i)
			++n[a[i][6]];
		if (n[0] != 1 || n[1] != 1) continue;
		const k = a[0][6] == 0? 0 : 1, l = 1 - k;
		print(a[k][0], s1[a[k][0]], a[k][1], a[k][2], a[k][5] == a[l][5]? '+' : '-',
			a[l][0], s2[a[l][0]], a[l][1], a[l][2],
			a[k][4] < a[l][4]? a[k][4] : a[l][4],
			a[k][4] < a[l][4]? a[l][4] : a[k][4], 60, `pn:Z:${g}`);
	}
}

/*****************
 * Main function *
 *****************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: pangene.js <command> [arguments]");
		print("Commands:");
		print("  call           call variants from a pangene graph");
		print("  call2html      generate a HTML page from call output");
		print("  calldiff       compare two call files");
		print("  bed2paf        generate PAF from a pair of samples");
		print("  gfa2matrix     generate gene_presence_absence.Rtab from pangene GFA");
		print("  getaa          generate protein files from Ensembl or GenCode annotations");
		print("  version        print version number");
		print("  flt-mmseqs     compare pangene sets in mmseqs format 2");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'call') pg_cmd_call(args);
	else if (cmd == 'call2html') pg_cmd_call2html(args);
	else if (cmd == 'calldiff') pg_cmd_calldiff(args);
	else if (cmd == 'getaa') pg_cmd_getaa(args);
	else if (cmd == 'bed2paf') pg_cmd_bed2paf(args);
	else if (cmd == 'gfa2matrix') pg_cmd_gfa2matrix(args);
	else if (cmd == 'flt-mmseqs') pg_cmd_flt_mmseqs(args);
	else if (cmd == 'version') {
		print(pg_version);
	} else throw Error("unrecognized command: " + cmd);
}

main(arguments);
