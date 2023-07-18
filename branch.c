#include <stdio.h>
#include "pgpriv.h"

int32_t pg_mark_branch_flt_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	uint32_t v, n_vtx = q->n_seg * 2, n_flt = 0;
	for (v = 0; v < n_vtx; ++v) {
		int32_t i, max_s1 = 0, off = q->idx[v]>>32, n = (int32_t)q->idx[v];
		for (i = 0; i < n; ++i)
			max_s1 = max_s1 > q->arc[off + i].s1? max_s1 : q->arc[off + i].s1;
		for (i = 0; i < n; ++i)
			if (q->arc[off + i].s1 < max_s1 * (1.0 - opt->branch_diff))
				q->arc[off + i].weak_br = 1, ++n_flt;
	}
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] marked %d diverged branches\n", __func__, pg_timestamp(), n_flt);
	return n_flt;
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
			if (!pg_hit_arc(a)) continue;
			sid = q->g2s[q->d->prot[a->pid].gid];
			if (vi >= 0 && a->cid != g->hit[vi].cid) v = (uint32_t)-1;
			w = (uint32_t)sid<<1 | a->rev;
			if (v != (uint32_t)-1) {
				e = pg_get_arc(q, v, w);
				if (e && e->weak_br) g->hit[vi].weak_br = 1;
				e = pg_get_arc(q, w^1, v^1);
				if (e && e->weak_br) a->weak_br = 1;
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

