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

void pg_graph_add_arc1(void *km, const pg_opt_t *opt, pg_graph_t *q, int32_t aid)
{
	const pg_data_t *d = q->d;
	const pg_genome_t *g = &d->genome[aid];
	int32_t i, n_hit = 0, m_hit = 0;
	pg_hit_t *hit = 0;
	for (i = 0; i < g->n_hit; ++i) {
		if (q->g2v[d->prot[g->hit[i].pid].gid] >= 0) {
			Kexpand(km, pg_hit_t, hit, n_hit, m_hit);
			hit[n_hit++] = g->hit[i];
		}
	}
	free(hit);
}

void pg_graph_add_arc(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int32_t i;
	for (i = 0; i < d->n_genome; ++i)
		pg_graph_add_arc1(0, opt, q, i);
}
