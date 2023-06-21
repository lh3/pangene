#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kalloc.h"

pg_graph_t *pg_graph_init(pg_data_t *d)
{
	pg_graph_t *g;
	g = PG_CALLOC(pg_graph_t, 1);
	g->d = d;
	g->v = PG_CALLOC(pg_vertex_t, d->n_gene);
	return g;
}

void pg_graph_destroy(pg_graph_t *q)
{
	free(q->g2v); free(q->v); free(q->a);
	free(q);
}

void pg_graph_flag_vtx(pg_graph_t *q)
{
	int32_t j, i;
	for (j = 0; j < q->d->n_genome; ++j) {
		pg_genome_t *g = &q->d->genome[j];
		for (i = 0; i < g->n_hit; ++i)
			if (q->g2v[q->d->prot[g->hit[i].pid].gid] >= 0)
				g->hit[i].vtx = 1;
	}
}

void pg_graph_flag_shadow(void *km, const pg_opt_t *opt, const pg_data_t *d, pg_genome_t *g) // similar to pg_gs_overlap() in geneset.c
{
	int32_t i0, i, j;
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		if (!g->hit[i].vtx || g->hit[i].pseudo) continue;
		while (i0 < i) {
			if (g->hit[i0].cid == g->hit[i].cid && g->hit[i0].ce > g->hit[i].cs)
				break;
			++i0;
		}
		for (j = i0; j < i; ++j) {
			uint64_t x;
			int32_t li, lj, gi, gj;
			double cov_short;
			if (g->hit[j].ce <= g->hit[i].cs) continue;
			if (!g->hit[j].vtx || g->hit[j].pseudo) continue;
			gi = d->prot[g->hit[j].pid].gid;
			gj = d->prot[g->hit[i].pid].gid;
			if (gi == gj) continue; // same gene
			x = pg_hit_overlap(g, &g->hit[j], &g->hit[i]);
			lj = pg_cds_len(&g->hit[j], g->exon);
			li = pg_cds_len(&g->hit[i], g->exon);
			cov_short = (double)(x>>32) / (li < lj? li : lj);
			assert(cov_short <= 1.0);
			if (cov_short < opt->min_ov_ratio) continue; // overlap too short
			if (g->hit[i].score2 < opt->max_score2_ratio * g->hit[j].score2) {
				g->hit[i].shadow = 1;
			} else if (g->hit[j].score2 < opt->max_score2_ratio * g->hit[i].score2) {
				g->hit[j].shadow = 1;
			}
		}
	}
}

void pg_graph_add_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int32_t i;
}
