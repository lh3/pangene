#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "pgpriv.h"
#include "ksort.h"

void pg_post_process(const pg_opt_t *opt, pg_data_t *d)
{
	int32_t i, j;
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %d genes and %d proteins\n", __func__, pg_timestamp(), d->n_gene, d->n_prot);
	pg_cap_score_dom(d);
	pg_flag_representative(d);
	if (!(opt->flag & PG_F_NO_JOINT_PSEUDO)) {
		int32_t n;
		n = pg_flag_pseudo_joint(opt, d);
		if (pg_verbose >= 3)
			fprintf(stderr, "[M::%s::%s] %d pseudogene hits identified jointly\n", __func__, pg_timestamp(), n);
	}
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		int32_t n_shadow, tot;
		for (i = 0, tot = 0; i < g->n_hit; ++i)
			if (!g->hit[i].flt) ++tot;
		n_shadow = pg_shadow(opt, d, j, 0);
		fprintf(stderr, "[M::%s::%s] genome[%d]: %s; %d hits remain, of which %d are shadowed\n",
				__func__, pg_timestamp(), j, g->label, tot, n_shadow);
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

#define get_score(a) ((a)->score_ori > (a)->score_dom? (a)->score_ori : (a)->score_dom)

static inline int32_t pg_get_score(const pg_graph_t *q, const pg_hit_t *a, int32_t ori)
{
	return ori || a->score_ori > a->score_dom || a->pid_dom0 < 0 || q->g2s[q->d->prot[a->pid_dom0].gid] >= 0? a->score_ori : a->score_dom;
}

void pg_gen_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	int32_t j, i, i0, *seg_cnt = 0, use_ori = !!(opt->flag&PG_F_ORI_FOR_BRANCH);
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
		pg_shadow(opt, q->d, j, 0);
		pg_hit_sort(g, 1); // sort by pg_hit_t::cm
		n_arc1 = 0;
		memset(seg_cnt, 0, q->n_seg * sizeof(int32_t));
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *a = &g->hit[i];
			uint32_t sid;
			if (a->flt || a->shadow) continue;
			sid = q->g2s[q->d->prot[a->pid].gid];
			assert(sid != (uint32_t)-1);
			w = (uint32_t)sid<<1 | a->rev;
			++seg_cnt[sid];
			if (a->cid != vcid) v = (uint32_t)-1, vpos = -1;
			if (v != (uint32_t)-1) {
				PG_GROW0(pg_tmparc_t, arc1, n_arc1, m_arc1);
				p = &arc1[n_arc1++], p->x = (uint64_t)v<<32|w, p->dist = a->cm - vpos, p->s1 = si, p->s2 = pg_get_score(q, a, use_ori);
				PG_GROW0(pg_tmparc_t, arc1, n_arc1, m_arc1);
				p = &arc1[n_arc1++], p->x = (uint64_t)(w^1)<<32|(v^1), p->dist = a->cm - vpos, p->s1 = pg_get_score(q, a, use_ori), p->s2 = si;
			}
			v = w, vpos = a->cm, vcid = a->cid, si = pg_get_score(q, a, use_ori);
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
				PG_GROW0(pg_tmparc_t, arc, n_arc, m_arc);
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
			PG_GROW(pg_arc_t, q->arc, q->n_arc, q->m_arc);
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

static void pg_hard_delete(pg_graph_t *q)
{
	int32_t i, k;
	for (i = k = 0; i < q->n_seg; ++i) {
		pg_seg_t *s = &q->seg[i];
		if (!s->del) {
			q->seg[k++] = *s;
		} else if (pg_verbose >= 3) {
			if (q->idx)
				fprintf(stderr, "#del\t%s\tavg_occ=%.1f\tdeg=%d,%d\tdist_deg=%d,%d\n", q->d->gene[s->gid].name, (double)s->tot_cnt / q->d->n_genome,
						(uint32_t)q->idx[(uint32_t)i<<1], (uint32_t)q->idx[(uint32_t)i<<1|1], s->n_dist_loci[0], s->n_dist_loci[1]);
			else
				fprintf(stderr, "#del\t%s\tavg_occ=%.1f\tdeg=*,*\tdist_deg=%d,%d\n", q->d->gene[s->gid].name, (double)s->tot_cnt / q->d->n_genome, s->n_dist_loci[0], s->n_dist_loci[1]);
		}
	}
	q->n_seg = k;
}

static void pg_flt_high_occ(int32_t max_avg_occ, int32_t max_degree, int32_t max_dist_loci, pg_graph_t *q)
{
	int32_t i, i0, n_high_occ = 0, n_high_deg = 0, n_high_loci = 0;
	for (i = 0; i < q->n_seg; ++i) // filter segments occurring too many times
		if (q->seg[i].tot_cnt > max_avg_occ * q->d->n_genome)
			q->seg[i].del = 1, ++n_high_occ;
	for (i0 = 0, i = 1; i <= q->n_arc; ++i) {
		if (i == q->n_arc || q->arc[i].x>>32 != q->arc[i0].x>>32) {
			int32_t sid = q->arc[i0].x>>32>>1;
			if (i - i0 > max_degree && !q->seg[sid].del)
				q->seg[sid].del = 1, ++n_high_deg;
			i0 = i;
		}
	}
	for (i = 0; i < q->n_seg; ++i) { // filter segments occurring too many times
		pg_seg_t *s = &q->seg[i];
		int32_t m = s->n_dist_loci[0] > s->n_dist_loci[1]? s->n_dist_loci[0] : s->n_dist_loci[1];
		if (m > max_dist_loci && !s->del)
			q->seg[i].del = 1, ++n_high_loci;
	}
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] filtered %d high-occurrence segments, %d high-degree segments and %d segments connecting distant loci\n",
				__func__, pg_timestamp(), n_high_occ, n_high_deg, n_high_loci);
	pg_hard_delete(q);
	pg_gen_g2s(q);
	pg_graph_flag_vtx(q);
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
			fprintf(stderr, "Z\t%c%s\t%c%s\t%d\t%d\n", "><"[a->x>>32&1], q->d->gene[q->seg[a->x>>32>>1].gid].name, "><"[a->x&1], q->d->gene[q->seg[(uint32_t)a->x>>1].gid].name, a->n_genome, a->weak_br);
	}
}

void pg_graph_gen(const pg_opt_t *opt, pg_graph_t *q)
{
	int32_t i;

	// graph 1: initial vertices
	PG_SET_FILTER(q->d, pseudo == 1);
	pg_gen_vtx(opt, q);
	pg_graph_flag_vtx(q);
	PG_SET_FILTER(q->d, vtx == 0);
	pg_gen_arc(opt, q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] round-1 graph: %d genes and %d arcs\n", __func__, pg_timestamp(), q->n_seg, q->n_arc);

	// graph 2: after removing high-occurrence vertices
	pg_flt_high_occ(opt->max_avg_occ * 2, opt->max_degree * 2, opt->max_dist_loci, q); // max_dist_loci is not applied as the relevant info not calculated at this point
	PG_SET_FILTER(q->d, vtx == 0);
	pg_gen_arc(opt, q); // don't apply "PG_SET_FILTER(q->d, shadow == 1)" as this would filter out CYP2D7
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] round-2 graph: %d genes and %d arcs\n", __func__, pg_timestamp(), q->n_seg, q->n_arc);

	// graph 3: branching filtering; vertices not changed
	for (i = 0; i < opt->n_branch_flt; ++i) {
		double r = 1.0 + (double)(opt->n_branch_flt - 1 - i) / opt->n_branch_flt;
		int32_t max_avg_occ = (int32_t)(opt->max_avg_occ * r + .499);
		int32_t max_degree = (int32_t)(opt->max_degree * r + .499);
		int32_t max_dist_loci = (int32_t)(opt->max_dist_loci * r + .499);
		pg_arc_index(q);
		pg_mark_branch_flt_arc(opt, q);
		pg_mark_branch_flt_hit(opt, q);
		PG_SET_FILTER(q->d, weak_br == 2);
		if (i > 0) {
			pg_flt_high_occ(max_avg_occ, max_degree, max_dist_loci, q);
			PG_SET_FILTER(q->d, vtx == 0);
		}
		pg_gen_arc(opt, q);
	}
	PG_SET_FILTER(q->d, shadow == 1);
	if (opt->min_arc_cnt > 1)
		pg_graph_cut_low_arc(opt, q);
	pg_arc_index(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] round-3 graph: %d genes and %d arcs\n", __func__, pg_timestamp(), q->n_seg, q->n_arc);
}
