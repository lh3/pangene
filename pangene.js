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
		const n_vtx = this.selength * 2;
		for (let v = 0; v < n_vtx; ++v)
			this.idx[v] = { o:0, n:0 };
		this.arc = this.arc.sort(function(a,b) { return a.v - b.v });
		let st = 0;
		for (let i = 1; i <= this.arc.length; ++i)
			if (i == this.arc.length || this.arc[i].v != this.arc[st].v)
				this.idx[this.arc[st].v] = { o:st, n:i-st }, st = i;
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
	}
	from_file(fn) {
		for (const line of k8_readline(fn))
			this.#parse_line(line);
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
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
