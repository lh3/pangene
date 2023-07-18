#include <stdio.h>
#include <assert.h>
#include "pgpriv.h"

static pg128_t **pg_gen_rep_pos(const pg_data_t *d)
{
	int32_t i, j;
	pg128_t **a;
	a = PG_MALLOC(pg128_t*, d->n_genome);
	for (j = 0; j < d->n_genome; ++j) {
		const pg_genome_t *g = &d->genome[j];
		pg128_t *aj;
		aj = a[j] = PG_CALLOC(pg128_t, d->n_gene);
		for (i = 0; i < d->n_gene; ++i)
			aj[i].x = (uint64_t)-1;
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *b = &g->hit[i];
			if (b->rep && b->rank == 0) {
				int32_t gid = d->prot[b->pid].gid;
				aj[gid].x = b->cid, aj[gid].y = b->cm;
			}
		}
	}
	return a;
}

static int32_t pg_n_close(int32_t close_thres, int32_t n_genome, pg128_t *const*a, int32_t g1, int32_t g2)
{
	int32_t j, n_close = 0;
	for (j = 0; j < n_genome; ++j) {
		const pg128_t *a1 = &a[j][g1], *a2 = &a[j][g2];
		int64_t d;
		if (a1->x == (uint64_t)-1 || a2->x == (uint64_t)-1) continue;
		if (a1->x != a2->x) continue;
		d = (int64_t)a1->y - (int64_t)a2->y;
		if (d >= -close_thres && d <= close_thres) ++n_close;
	}
	return n_close;
}

int32_t pg_mark_branch_flt_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	uint32_t v, n_vtx = q->n_seg * 2, n_flt1 = 0, n_flt2 = 0;
	int32_t j, *max_gid;
	pg128_t **pos;
	max_gid = PG_MALLOC(int32_t, q->n_seg * 2); // we don't need an array this large, but the program wastes a lot more memory elsewhere
	pos = pg_gen_rep_pos(q->d);
	for (v = 0; v < n_vtx; ++v) {
		int32_t i, n_max, max_s1 = 0, n = (int32_t)q->idx[v];
		pg_arc_t *a = &q->arc[q->idx[v]>>32];
		if (n < 2) continue;
		for (i = 0; i < n; ++i)
			max_s1 = max_s1 > a[i].s1? max_s1 : a[i].s1;
		for (i = n_max = 0; i < n; ++i)
			if (a[i].s1 == max_s1)
				max_gid[n_max++] = q->seg[(uint32_t)a[i].x>>1].gid;
		assert(n_max > 0);
		for (i = 0; i < n; ++i) {
			if (a[i].s1 < max_s1 * (1.0 - opt->branch_diff)) {
				int32_t n_close = 0;
				for (j = 0; j < n_max; ++j)
					n_close += pg_n_close(opt->close_thres, q->d->n_genome, pos, max_gid[j], q->seg[(uint32_t)a[i].x>>1].gid);
				if (n_close == 0) a[i].weak_br = 2, ++n_flt2;
				else a[i].weak_br = 1, ++n_flt1;
			}
		}
	}
	for (j = 0; j < q->d->n_genome; ++j) free(pos[j]);
	free(pos);
	free(max_gid);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] marked %d closely diverged branches and %d distantly diverged branches\n", __func__, pg_timestamp(), n_flt1, n_flt2);
	return n_flt1 + n_flt2;
}

int32_t pg_mark_branch_flt_hit(const pg_opt_t *opt, pg_graph_t *q) // call after pg_mark_branch_flt_arc()
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
			if (a->flt) continue;
			sid = q->g2s[q->d->prot[a->pid].gid];
			if (vi >= 0 && a->cid != g->hit[vi].cid) v = (uint32_t)-1;
			w = (uint32_t)sid<<1 | a->rev;
			if (v != (uint32_t)-1) {
				e = pg_get_arc(q, v, w);
				if (e && e->weak_br) g->hit[vi].weak_br = e->weak_br;
				e = pg_get_arc(q, w^1, v^1);
				if (e && e->weak_br) a->weak_br = e->weak_br;
			}
			v = w, vi = i;
		}
		for (i = 0; i < g->n_hit; ++i)
			if (g->hit[i].weak_br)
				++n_flt;
		pg_hit_sort(g, 0); // sort by pg_hit_t::cs
	}
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] marked %d diverged hits\n", __func__, pg_timestamp(), n_flt);
	return n_flt;
}
