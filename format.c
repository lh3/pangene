#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include "pgpriv.h"

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = PG_REALLOC(char, s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void pg_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'l' && *(p+1) == 'd') {
				int c, i, l = 0;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

static inline void pg_write_bed_hit(kstring_t *out, const pg_data_t *d, int32_t aid, const pg_hit_t *a)
{
	const pg_genome_t *g = &d->genome[aid];
	int32_t i;
	char buf[16];
	pg_sprintf_lite(out, "%s\t%ld\t%ld\t%s\t%d\t%c\t", g->ctg[a->cid].name, a->cs, a->ce, d->prot[a->pid].name, a->score, "+-"[a->rev]);
	pg_sprintf_lite(out, "%ld\t%ld\t0\t%d\t", a->cs, a->ce, a->n_exon);
	for (i = 0; i < a->n_exon; ++i)
		pg_sprintf_lite(out, "%d,", g->exon[a->off_exon + i].oe - g->exon[a->off_exon + i].os);
	pg_sprintf_lite(out, "\t");
	for (i = 0; i < a->n_exon; ++i)
		pg_sprintf_lite(out, "%d,", g->exon[a->off_exon + i].os);
	snprintf(buf, 15, "%.4f", (double)a->mlen / a->blen);
	pg_sprintf_lite(out, "\tft:i:%d\tio:i:%d\tis:i:%d\t\trk:i:%d\trp:i:%d\tsd:i:%d\tvt:i:%d\tps:i:%d\tbr:i:%d\tcm:i:%ld\ts2:i:%d\tid:f:%s\tdm:Z:%s\n",
		a->flt, a->flt_iso_ov, a->flt_iso_scat, a->rank, a->rep, a->shadow, a->vtx, a->pseudo, a->weak_br, a->cm, a->score2, buf, a->pid_dom < 0? "*" : d->prot[a->pid_dom].name);
}

static void pg_write_bed_genome(const pg_data_t *d, int32_t aid, int32_t is_walk)
{
	int32_t i;
	const pg_genome_t *g;
	kstring_t out = {0,0,0};
	assert(aid < d->n_genome);
	g = &d->genome[aid];
	for (i = 0; i < g->n_hit; ++i) {
		const pg_hit_t *a = &g->hit[i];
		if (is_walk && a->flt) continue;
		out.l = 0;
		pg_write_bed_hit(&out, d, aid, a);
		fwrite(out.s, 1, out.l, stdout);
	}
	free(out.s);
}

void pg_write_bed(const pg_data_t *d, int32_t is_walk)
{
	int32_t j;
	for (j = 0; j < d->n_genome; ++j)
		pg_write_bed_genome(d, j, is_walk);
}

static void pg_write_seg(const pg_graph_t *g)
{
	const pg_data_t *d = g->d;
	kstring_t out = {0,0,0};
	int32_t i;
	for (i = 0; i < g->n_seg; ++i) {
		const pg_seg_t *s = &g->seg[i];
		int32_t gid = s->gid;
		int32_t pid = g->d->gene[gid].rep_pid;
		out.l = 0;
		pg_sprintf_lite(&out, "S\t%s\t*\tLN:i:%d\tng:i:%d\tnc:i:%d\tc1:i:%d\tc2:i:%d\tpp:Z:%s\n",
			d->gene[gid].name, d->prot[pid].len, s->n_genome, s->tot_cnt, s->n_dom, s->n_sub, d->prot[pid].name);
		fwrite(out.s, 1, out.l, stdout);
	}
	free(out.s);
}

static void pg_write_arc(const pg_graph_t *g)
{
	const pg_data_t *d = g->d;
	kstring_t out = {0,0,0};
	int32_t i;
	for (i = 0; i < g->n_arc; ++i) {
		pg_arc_t *a = &g->arc[i];
		uint32_t v = a->x>>32, w = (uint32_t)a->x;
		out.l = 0;
		pg_sprintf_lite(&out, "L\t%s\t%c\t%s\t%c\t0M\t", d->gene[g->seg[v>>1].gid].name, "+-"[v&1], d->gene[g->seg[w>>1].gid].name, "+-"[w&1]);
		pg_sprintf_lite(&out, "ng:i:%d\tnc:i:%d\tad:i:%d\ts1:i:%d\ts2:i:%d\n", a->n_genome, a->tot_cnt, a->avg_dist, a->s1, a->s2);
		fwrite(out.s, 1, out.l, stdout);
	}
	free(out.s);
}

void pg_write_graph(const pg_graph_t *g)
{
	pg_write_seg(g);
	pg_write_arc(g);
}

static int32_t pg_parse_sample(kstring_t *buf, const char *name) // parse "sample#hap#name"
{
	const char *p, *q;
	int32_t i, hap = 0;
	buf->l = 0;
	for (p = q = name, i = 0;; ++p) {
		if (*p == 0 || *p == '#') {
			if (i == 0) {
				if (p - q == 0) return -1;
				str_copy(buf, q, p);
				buf->s[p - q] = 0;
			} else if (i == 1) {
				char *r;
				hap = strtol(q, &r, 10);
				if (r == p && hap >= 0) return hap;
				else return -1;
			}
			q = p + 1, ++i;
			if (*p == 0) break;
		}
	}
	return i == 3? hap : -1;
}

void pg_write_walk(pg_graph_t *q)
{
	int32_t i, i0, j;
	kstring_t out = {0,0,0}, buf = {0,0,0};
	pg_data_t *d = q->d;
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		pg_hit_sort(&d->genome[j], 1);
		for (i0 = 0, i = 1; i <= g->n_hit; ++i) {
			if (i == g->n_hit || g->hit[i].cid != g->hit[i0].cid) {
				int32_t k, n, hap, cid = g->hit[i0].cid;
				hap = pg_parse_sample(&buf, g->ctg[cid].name);
				out.l = 0;
				if (hap >= 0)
					pg_sprintf_lite(&out, "W\t%s\t%d", buf.s, hap);
				else if (g->label)
					pg_sprintf_lite(&out, "W\t%s\t0", g->label);
				else
					pg_sprintf_lite(&out, "W\t%d\t0", j);
				pg_sprintf_lite(&out, "\t%s\t*\t*\t", g->ctg[cid].name);
				for (k = i0, n = 0; k < i; ++k) {
					const pg_hit_t *a = &g->hit[k];
					if (a->flt) continue;
					pg_sprintf_lite(&out, "%c%s", "><"[a->rev], d->gene[d->prot[a->pid].gid].name);
					++n;
				}
				if (n > 0) {
					pg_sprintf_lite(&out, "\tlf:B:i");
					for (k = i0; k < i; ++k) {
						const pg_hit_t *a = &g->hit[k];
						if (a->flt) continue;
						pg_sprintf_lite(&out, ",%d", a->lof);
					}
					puts(out.s);
				}
				i0 = i;
			}
		}
		pg_hit_sort(&d->genome[j], 0);
	}
	free(buf.s);
	free(out.s);
}
