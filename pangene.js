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
	#test_inv(v0, mark, max_ext) {
		let q = [];
		q.push(v0);
		while (q.length) {
			const v = q.shift();
			mark[v] = v0;
			const o = this.idx[v].o, n = this.idx[v].n;
			for (let i = 0; i < n; ++i) {
			}
		}
	}
	collect_inv(max_ext) {
		const n_vtx = this.seg.length * 2;
		let mark = [];
		for (let i = 0; i < n_vtx; ++i)
			mark[i] = -1;
		for (let v = 0; v < n_vtx; ++v) {
			const n = this.idx[v].n;
			if (n < 2) continue;
		}
	}
}

/**************************
 * Program Structure Tree *
 **************************/

class LinkedList { // double linked list
	constructor() {
		this.size = 0;
		this.p = null;
		this.q = null;
	}
	push(data) {
		let ptr = { data:data, p:null, q:null };
		if (this.p == null && this.q == null) {
			this.p = this.q = ptr;
		} else {
			this.q.q = ptr;
			ptr.p = this.q;
			this.q = ptr;
		}
		++this.size;
		return this.q;
	}
	pop() {
		let ptr;
		if (this.p == null && this.q == null) return null;
		else if (this.p == this.q) ptr = this.p, this.p = this.q = null;
		else ptr = this.q, this.q = ptr.p, this.q.q = null;
		--this.size;
		return ptr.data;
	}
	shift() {
		let ptr;
		if (this.p == null && this.q == null) return null;
		else if (this.p == this.q) ptr = this.p, this.p = this.q = null;
		else ptr = this.p, this.p = ptr.q, this.p.p = null;
		--this.size;
		return ptr.data;
	}
	top() {
		return this.q != null? this.q.data : null;
	}
	delete(ptr) {
		ptr.p.q = ptr.q;
		ptr.q.p = ptr.p;
	}
	concat(bl) {
		for (let p = bl.p; p != null; p = p.q)
			this.push(p.data);
	}
}

class SegEdgeGraph {
	constructor() {
		this.cat = [];
		this.n_node = 0;
		this.arc = [];
		this.idx = [];
		this.dfs_time = [];
	}
	from_gfa(g) {
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
		for (let v = 0; v < n_vtx; ++v) this.cat[v] = -1;
		for (let v = 0; v < n_vtx; ++v) {
			if (this.cat[v] >= 0) continue;
			let stack = [v];
			while (stack.length > 0) { // a DFS
				const w = stack.pop();
				this.cat[w] = x;
				const n = idx[w].n, off = idx[w].o;
				for (let i = 0; i < n; ++i) {
					const u = a[off + i][1];
					if (this.cat[u] < 0) {
						this.cat[u] = x;
						stack.push(u);
					} else if (this.cat[u] != x) {
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
			this.arc.push({ v:this.cat[i*2],   w:this.cat[i*2|1], seg:i, ori:1,  dfs_type:0 });
			this.arc.push({ v:this.cat[i*2|1], w:this.cat[i*2],   seg:i, ori:-1, dfs_type:0 });
		}
		// index arc[]
		for (let i = 0; i < this.n_node; ++i)
			this.idx[i] = { n:0, o:0 };
		this.arc.sort(function(x, y) { return x.v - y.v });
		for (let i = 1, i0 = 0; i <= this.arc.length; ++i)
			if (i == this.arc.length || this.arc[i0].v != this.arc[i].v)
				this.idx[this.arc[i0].v] = { o:i0, n:i-i0 }, i0 = i;
	}
	dfs_traverse() {
		for (let v = 0; v < this.n_node; ++v)
			this.dfs_time[v] = [-1, -1]; // discovery time and finishing time
		let t_dis = 0, t_fin = 0, state = [];
		for (let v = 0; v < this.n_node; ++v)
			state[v] = 0; // not visited
		for (let v = 0; v < this.n_node; ++v) {
			if (state[v] != 0) continue; // visited before
			this.dfs_time[v][0] = t_dis++;
			state[v] = 2; // in stack
			let stack = [[v, 0]];
			while (stack.length > 0) {
				const [w, i] = stack.pop();
				const n = this.idx[w].n, off = this.idx[w].o;
				if (i < n) {
					stack.push([w, i + 1]); // repush to the stack
					const u = this.arc[off + i].w;
					if (state[u] == 0) { // not visited before
						state[u] = 2; // in stack
						this.dfs_time[u][0] = t_dis++;
						stack.push([u, 0]);
						this.arc[off + i].dfs_type = 1; // a tree edge
					} else if (state[u] == 2) {
						this.arc[off + i].dfs_type = 2; // a back edge
					}
				} else {
					state[w] = 1; // out of stack
					this.dfs_time[w][1] = t_fin++;
				}
			}
		}
		if (t_dis != this.n_node || t_fin != this.n_node)
			throw Error("DFS bug");
	}
	cycle_equiv() {
		this.dfs_traverse();
		// set discover time
		/*
		let td = [];
		for (let j = 0; j < dfs.length; ++j)
			td[dfs[j]] = j;
		for (let j = dfs.length - 1; j >= 0; --j) {
			const v = dfs[j];
			const n = this.idx[v].n, off = this.idx[v].o;
			let hi0 = this.n_node;
			for (let i = 0; i < n; ++i) {
				const w = this.arc[off + i].w;
				hi0 = hi0 < td[w]? hi0 : td[w];
			}
		}
		*/
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
	e.cycle_equiv();
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
