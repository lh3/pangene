#include "pgpriv.h"

pg_graph_t *pg_graph_init(const pg_data_t *d)
{
	pg_graph_t *g;
	g = PG_CALLOC(pg_graph_t, 1);
	g->d = d;
	g->v = PG_CALLOC(pg_vertex_t, d->n_gene);
	return g;
}

void pg_graph_destroy(pg_graph_t *g)
{
	free(g->v);
	free(g);
}
