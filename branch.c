#include <stdio.h>
#include <assert.h>
#include "pgpriv.h"

static pg128_t **pg_gen_rep_pos(const pg_data_t *d) // only work if sorted
{
	int32_t i, j;
	pg128_t **a;
	a = PG_MALLOC(pg128_t*, d->n_genome);
	for (j = 0; j < d->n_genome; ++j) {
		const pg_genome_t *g = &d->genome[j];
		pg128_t *aj;
		int32_t r = 0;
		aj = a[j] = PG_CALLOC(pg128_t, d->n_gene);
		for (i = 0; i < d->n_gene; ++i)
			aj[i].x = (uint64_t)-1;
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *b = &g->hit[i];
			if (!b->shadow && !b->flt) {
				int32_t gid = d->prot[b->pid].gid;
				aj[gid].x = (uint64_t)b->cid<<32 | r;
				aj[gid].y = b->cm;
				++r;
			}
		}
	}
	return a;
}

static int32_t pg_n_local(int32_t local_dist, int32_t local_count, int32_t n_genome, pg128_t *const*a, int32_t g1, int32_t g2)
{
	int32_t j, n_local = 0;
	for (j = 0; j < n_genome; ++j) {
		const pg128_t *a1 = &a[j][g1], *a2 = &a[j][g2];
		int64_t d;
		int32_t c;
		if (a1->x == (uint64_t)-1 || a2->x == (uint64_t)-1) continue;
		if (a1->x>>32 != a2->x>>32) continue;
		d = (int64_t)a1->y - (int64_t)a2->y;
		c = (int32_t)a1->x - (int32_t)a2->x;
		if ((d >= -local_dist && d <= local_dist) || (c >= -local_count && c <= local_count))
			++n_local;
	}
	return n_local;
}

int32_t pg_mark_branch_flt_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	uint32_t v, n_vtx = q->n_seg * 2, n_flt1 = 0, n_flt2 = 0;
	int32_t j, *max_gid;
	pg128_t **pos;
	max_gid = PG_MALLOC(int32_t, q->n_seg * 2); // we don't need an array this large, but the program uses a lot more memory elsewhere
	pos = pg_gen_rep_pos(q->d);
	for (v = 0; v < n_vtx; ++v) {
		int32_t i, n_max, max_s1 = 0, n = (int32_t)q->idx[v];
		pg_arc_t *a = &q->arc[q->idx[v]>>32];
		if (n < 2) continue;
		for (i = 0; i < n; ++i) // find max score
			max_s1 = max_s1 > a[i].s1? max_s1 : a[i].s1;
		for (i = n_max = 0; i < n; ++i)
			if (a[i].s1 == max_s1) // find arcs with the max score
				max_gid[n_max++] = q->seg[(uint32_t)a[i].x>>1].gid;
		assert(n_max > 0);
		for (i = 0; i < n; ++i) {
			if (a[i].s1 < max_s1 * (1.0 - opt->branch_diff)) {
				int32_t n_local = 0, gid = q->seg[(uint32_t)a[i].x>>1].gid;
				for (j = 0; j < n_max; ++j)
					n_local += pg_n_local(opt->local_dist, opt->local_count, q->d->n_genome, pos, max_gid[j], gid);
				if (n_local == 0 || a[i].s1 < max_s1 * (1.0 - opt->branch_diff_cut)) a[i].weak_br = 2, ++n_flt2;
				else a[i].weak_br = 1, ++n_flt1;
				//fprintf(stderr, "B\t%s\t%s\t%s\t%.4f\t%d\n", q->d->gene[q->seg[a[i].x>>33].gid].name, q->d->gene[max_gid[0]].name, q->d->gene[gid].name, (double)a[i].s1 / max_s1, n_local);
			}
		}
	}
	for (j = 0; j < q->d->n_genome; ++j) free(pos[j]);
	free(pos);
	free(max_gid);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] marked %d locally diverged branches and %d distantly diverged branches\n", __func__, pg_timestamp(), n_flt1, n_flt2);
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
			if (a->flt || a->shadow) continue;
			sid = q->g2s[q->d->prot[a->pid].gid];
			if (vi >= 0 && a->cid != g->hit[vi].cid) v = (uint32_t)-1;
			w = (uint32_t)sid<<1 | a->rev;
			if (v != (uint32_t)-1) {
				//fprintf(stderr, "L\t%s\t%s\t%s\n", g->ctg[a->cid].name, d->gene[q->seg[v>>1].gid].name, d->gene[q->seg[w>>1].gid].name);
				e = pg_get_arc(q, v, w);
				if (e && e->weak_br)
					g->hit[vi].weak_br = g->hit[vi].weak_br > e->weak_br? g->hit[vi].weak_br : e->weak_br;
				e = pg_get_arc(q, w^1, v^1);
				if (e && e->weak_br)
					a->weak_br = a->weak_br > e->weak_br? a->weak_br : e->weak_br;
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
