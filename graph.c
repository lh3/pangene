#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "pgpriv.h"
#include "ksort.h"

void pg_post_process(const pg_opt_t *opt, pg_data_t *d)
{
	int32_t i;
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %d genes and %d proteins\n", __func__, pg_timestamp(), d->n_gene, d->n_prot);
	pg_flag_primary(d);
	for (i = 0; i < d->n_genome; ++i) {
		pg_genome_t *g = &d->genome[i];
		int32_t n_pseudo, n_shadow;
		n_pseudo = pg_flag_pseudo(d->prot, g);
		pg_hit_sort(g, 0);
		n_shadow = pg_flag_shadow(opt, d->prot, g, 0, 0);
		if (pg_verbose >= 3)
			fprintf(stderr, "[M::%s::%s] genome %d: %d pseudo, %d shadow\n", __func__, pg_timestamp(), i, n_pseudo, n_shadow);
	}
}

pg_graph_t *pg_graph_init(pg_data_t *d)
{
	pg_graph_t *g;
	g = PG_CALLOC(pg_graph_t, 1);
	g->d = d;
	g->seg = PG_CALLOC(pg_seg_t, d->n_gene);
	return g;
}

void pg_graph_destroy(pg_graph_t *q)
{
	free(q->g2s); free(q->seg); free(q->arc); free(q->idx);
	free(q);
}

void pg_gen_g2s(pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int32_t i;
	if (q->g2s) free(q->g2s);
	q->g2s = PG_MALLOC(int32_t, d->n_gene);
	for (i = 0; i < d->n_gene; ++i)
		q->g2s[i] = -1;
	for (i = 0; i < q->n_seg; ++i)
		q->g2s[q->seg[i].gid] = i;
}

void pg_graph_flag_vtx(pg_graph_t *q)
{
	int32_t j, i;
	for (j = 0; j < q->d->n_genome; ++j) {
		pg_genome_t *g = &q->d->genome[j];
		for (i = 0; i < g->n_hit; ++i)
			g->hit[i].vtx = (q->g2s[q->d->prot[g->hit[i].pid].gid] >= 0);
	}
}

typedef struct {
	uint64_t x;
	int32_t n, dist;
	int32_t s1, s2;
} pg_tmparc_t;

#define sort_key_x(a) ((a).x)
KRADIX_SORT_INIT(pg_tmparc, pg_tmparc_t, sort_key_x, 8)

void pg_gen_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	int32_t j, i, i0, *seg_cnt = 0;
	int64_t n_arc = 0, m_arc = 0, n_arc1 = 0, m_arc1 = 0;
	pg_tmparc_t *p, *arc = 0, *arc1 = 0;

	q->n_arc = 0;
	seg_cnt = PG_MALLOC(int32_t, q->n_seg);
	for (i = 0; i < q->n_seg; ++i)
		q->seg[i].n_genome = q->seg[i].tot_cnt = 0;
	for (j = 0; j < q->d->n_genome; ++j) {
		pg_genome_t *g = &q->d->genome[j];
		uint32_t w, v = (uint32_t)-1;
		int64_t vpos = -1;
		int32_t vcid = -1, si = -1;
		pg_flag_shadow(opt, q->d->prot, g, 1, 1); // this requires sorting by pg_hit_t::cs
		pg_hit_sort(g, 1); // sort by pg_hit_t::cm
		n_arc1 = 0;
		memset(seg_cnt, 0, q->n_seg * sizeof(int32_t));
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *a = &g->hit[i];
			uint32_t sid;
			if (!pg_hit_arc(a)) continue;
			sid = q->g2s[q->d->prot[a->pid].gid];
			w = (uint32_t)sid<<1 | a->rev;
			++seg_cnt[sid];
			if (a->cid != vcid) v = (uint32_t)-1, vpos = -1;
			if (v != (uint32_t)-1) {
				PG_EXTEND0(pg_tmparc_t, arc1, n_arc1, m_arc1);
				p = &arc1[n_arc1++], p->x = (uint64_t)v<<32|w, p->dist = a->cm - vpos, p->s1 = si, p->s2 = a->score;
				PG_EXTEND0(pg_tmparc_t, arc1, n_arc1, m_arc1);
				p = &arc1[n_arc1++], p->x = (uint64_t)(w^1)<<32|(v^1), p->dist = a->cm - vpos, p->s1 = a->score, p->s2 = si;
			}
			v = w, vpos = a->cm, vcid = a->cid, si = a->score;
		}
		pg_hit_sort(g, 0); // sort by pg_hit_t::cs
		assert(n_arc1 <= INT32_MAX);
		for (i = 0; i < q->n_seg; ++i)
			q->seg[i].n_genome += (seg_cnt[i] > 0), q->seg[i].tot_cnt += seg_cnt[i];
		radix_sort_pg_tmparc(arc1, arc1 + n_arc1);
		for (i = 1, i0 = 0; i <= n_arc1; ++i) {
			if (i == n_arc1 || arc1[i0].x != arc1[i].x) {
				int32_t j, max_s1 = 0, max_s2 = 0;
				uint64_t dist = 0;
				for (j = i0; j < i; ++j) {
					dist += arc1[j].dist;
					max_s1 = max_s1 > arc1[j].s1? max_s1 : arc1[j].s1;
					max_s2 = max_s2 > arc1[j].s2? max_s2 : arc1[j].s2;
				}
				PG_EXTEND0(pg_tmparc_t, arc, n_arc, m_arc);
				p = &arc[n_arc++];
				p->x = arc1[i0].x;
				p->n = i - i0;
				p->dist = (int32_t)((double)dist / (i - i0) + .499);
				p->s1 = max_s1, p->s2 = max_s2;
				i0 = i;
			}
		}
	}
	free(arc1);
	free(seg_cnt);

	assert(n_arc <= INT32_MAX);
	radix_sort_pg_tmparc(arc, arc + n_arc);
	q->n_arc = 0;
	for (i0 = 0, i = 1; i <= n_arc; ++i) {
		if (i == n_arc || arc[i].x != arc[i0].x) {
			int64_t dist = 0, s1 = 0, s2 = 0;
			int32_t n = 0;
			pg_arc_t *p;
			for (j = i0; j < i; ++j) {
				n += arc[j].n;
				dist += (uint64_t)arc[j].dist * arc[j].n;
				s1 += arc[j].s1;
				s2 += arc[j].s2;
			}
			PG_EXTEND(pg_arc_t, q->arc, q->n_arc, q->m_arc);
			p = &q->arc[q->n_arc++];
			memset(p, 0, sizeof(*p));
			p->x = arc[i0].x;
			p->n_genome = i - i0;
			p->tot_cnt = n;
			p->avg_dist = (int64_t)((double)dist / n + .499);
			p->s1 = (int32_t)((double)s1 / (i - i0) + .499);
			p->s2 = (int32_t)((double)s2 / (i - i0) + .499);
			i0 = i;
		}
	}
	free(arc);
}

void pg_graph_rm_del(pg_graph_t *q)
{
	int32_t i, k;
	for (i = 0, k = 0; i < q->n_arc; ++i) {
		pg_arc_t *a = &q->arc[i];
		uint32_t v = a->x>>32, w = (uint32_t)a->x;
		if (!(a->del || q->seg[v>>1].del || q->seg[w>>1].del))
			q->arc[k++] = q->arc[i];
	}
	q->n_arc = k;
}

static void pg_graph_cut_low_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	int32_t i, n_aflt = 0;
	for (i = 0; i < q->n_arc; ++i)
		if (q->arc[i].n_genome < opt->min_arc_cnt)
			q->arc[i].del = 1, ++n_aflt;
	pg_graph_rm_del(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] filtered %d low-occurrence arcs\n", __func__, pg_timestamp(), n_aflt);
}

static uint64_t *pg_arc_index_core(const pg_graph_t *q)
{
	uint64_t *idx;
	int32_t i, i0;
	idx = PG_CALLOC(uint64_t, q->n_seg * 2);
	for (i0 = 0, i = 1; i <= q->n_arc; ++i)
		if (i == q->n_arc || q->arc[i].x>>32 != q->arc[i0].x>>32)
			idx[q->arc[i0].x>>32] = (uint64_t)i0<<32 | (i - i0), i0 = i;
	return idx;
}

static void pg_arc_index(pg_graph_t *q)
{
	if (q->idx) free(q->idx);
	q->idx = pg_arc_index_core(q);
}

static void pg_graph_flt_high_occ(const pg_opt_t *opt, pg_graph_t *q)
{
	int32_t i, i0, k, n_high_occ = 0, n_high_deg = 0;
	for (i = 0; i < q->n_seg; ++i) // filter segments occurring too many times
		if (q->seg[i].tot_cnt > opt->max_avg_occ * q->seg[i].n_genome)
			q->seg[i].del = 1, ++n_high_occ;
	for (i0 = 0, i = 1; i <= q->n_arc; ++i) {
		if (i == q->n_arc || q->arc[i].x>>32 != q->arc[i0].x>>32) {
			int32_t sid = q->arc[i0].x>>32>>1;
			if (i - i0 > opt->max_degree && !q->seg[sid].del)
				q->seg[sid].del = 1, ++n_high_deg;
			i0 = i;
		}
	}
	if (pg_verbose >= 3) {
		fprintf(stderr, "[M::%s::%s] %d high-occurrence segments\n", __func__, pg_timestamp(), n_high_occ);
		fprintf(stderr, "[M::%s::%s] %d high-degree segments additionally\n", __func__, pg_timestamp(), n_high_deg);
	}
	for (i = k = 0; i < q->n_seg; ++i) {
		if (!q->seg[i].del)
			q->seg[k++] = q->seg[i];
		else if (pg_verbose >= 3)
			fprintf(stderr, "[M::%s] dropped %s\n", __func__, q->d->gene[q->seg[i].gid].name);
	}
	q->n_seg = k;
	pg_gen_g2s(q);
}

static int32_t pg_mark_branch_flt_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	uint32_t v, n_vtx = q->n_seg * 2, n_flt = 0;
	for (v = 0; v < n_vtx; ++v) {
		int32_t i, max_s1 = 0, off = q->idx[v]>>32, n = (int32_t)q->idx[v];
		for (i = 0; i < n; ++i)
			max_s1 = max_s1 > q->arc[off + i].s1? max_s1 : q->arc[off + i].s1;
		for (i = 0; i < n; ++i)
			if (q->arc[off + i].s1 < max_s1 * (1.0 - opt->branch_diff))
				q->arc[off + i].branch_flt = 1, ++n_flt;
	}
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] marked %d diverged branches\n", __func__, pg_timestamp(), n_flt);
	return n_flt;
}

static inline const pg_arc_t *pg_get_arc(const pg_graph_t *q, uint32_t v, uint32_t w)
{
	int32_t i, n = (int32_t)q->idx[v];
	const pg_arc_t *a = &q->arc[q->idx[v]>>32];
	for (i = 0; i < n; ++i)
		if ((uint32_t)a[i].x == w)
			return &a[i];
	return 0;
}

static int32_t pg_mark_branch_flt_hit(const pg_opt_t *opt, pg_graph_t *q) // call after pg_mark_branch_flt_arc()
{
	pg_data_t *d = q->d;
	int32_t i, j, n_flt = 0;
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		uint32_t v = (uint32_t)-1;
		int32_t vi = -1;
		pg_hit_sort(g, 1); // sort by pg_hit_t::cm
		for (i = 0; i < g->n_hit; ++i) {
			pg_hit_t *a = &g->hit[i];
			const pg_arc_t *e;
			uint32_t w;
			int32_t sid;
			if (!pg_hit_arc(a)) continue;
			sid = q->g2s[q->d->prot[a->pid].gid];
			if (vi >= 0 && a->cid != g->hit[vi].cid) v = (uint32_t)-1;
			w = (uint32_t)sid<<1 | a->rev;
			if (v != (uint32_t)-1) {
				e = pg_get_arc(q, v, w);
				if (e && e->branch_flt) g->hit[vi].branch_flt = 1;
				e = pg_get_arc(q, w^1, v^1);
				if (e && e->branch_flt) a->branch_flt = 1;
			}
			v = w, vi = i;
		}
		for (i = 0; i < g->n_hit; ++i)
			if (g->hit[i].branch_flt)
				++n_flt;
		pg_hit_sort(g, 0); // sort by pg_hit_t::cs
	}
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] marked %d diverged hits\n", __func__, pg_timestamp(), n_flt);
	return n_flt;
}

void pg_debug_gene(const pg_graph_t *q, const char *name)
{
	int32_t gid, sid;
	uint32_t j;
	gid = pg_dict_get(q->d->d_gene, name);
	assert(gid >= 0);
	sid = q->g2s[gid];
	assert(sid >= 0);
	for (j = 0; j < q->n_arc; ++j) {
		const pg_arc_t *a = &q->arc[j];
		if (a->x>>32>>1 == sid)
			fprintf(stderr, "Z\t%c%s\t%c%s\t%d\t%d\n", "><"[a->x>>32&1], q->d->gene[q->seg[a->x>>32>>1].gid].name, "><"[a->x&1], q->d->gene[q->seg[(uint32_t)a->x>>1].gid].name, a->n_genome, a->branch_flt);
	}
}

void pg_graph_gen(const pg_opt_t *opt, pg_graph_t *q)
{
	int32_t i;

	// graph 1: initial vertices
	pg_gen_vtx(opt, q);
	pg_graph_flag_vtx(q);
	pg_gen_arc(opt, q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] round-1 graph: %d genes and %d arcs\n", __func__, pg_timestamp(), q->n_seg, q->n_arc);

	// graph 2: after removing high-occurrence vertices
	pg_graph_flt_high_occ(opt, q);
	pg_graph_flag_vtx(q);
	pg_gen_arc(opt, q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] round-2 graph: %d genes and %d arcs\n", __func__, pg_timestamp(), q->n_seg, q->n_arc);

	// graph 3: branching filtering; vertices not changed
	for (i = 0; i < opt->n_branch_flt; ++i) {
		pg_arc_index(q);
		pg_mark_branch_flt_arc(opt, q);
		pg_mark_branch_flt_hit(opt, q);
		pg_gen_arc(opt, q);
	}
	if (opt->min_arc_cnt > 1)
		pg_graph_cut_low_arc(opt, q);
	pg_arc_index(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] round-3 graph: %d genes and %d arcs\n", __func__, pg_timestamp(), q->n_seg, q->n_arc);
}
