#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kalloc.h"

void pg_post_process(const pg_opt_t *opt, pg_data_t *d)
{
	int32_t i;
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] %d genes and %d proteins\n", __func__, pg_timestamp(), d->n_gene, d->n_prot);
	for (i = 0; i < d->n_genome; ++i) {
		pg_genome_t *g = &d->genome[i];
		int32_t n_pseudo, n_shadow;
		n_pseudo = pg_flag_pseudo(0, d->prot, g);
		n_shadow = pg_flag_shadow(opt, d->prot, g, 0);
		pg_hit_sort(0, g, 0);
		if (pg_verbose >= 3)
			fprintf(stderr, "[M::%s::%s] genome %d: %d pseudo, %d shadow\n", __func__, pg_timestamp(), i, n_pseudo, n_shadow);
	}
}

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

static void pg_gen_vertex_aux(void *km, const pg_opt_t *opt, const pg_data_t *d, pg_genome_t *g, uint64_t *cnt)
{
	int32_t i;
	int8_t *flag;
	flag = Kcalloc(km, int8_t, d->n_gene);
	for (i = 0; i < g->n_hit; ++i) {
		const pg_hit_t *a = &g->hit[i];
		int32_t gid;
		if (a->rank != 0) continue;
		gid = d->prot[a->pid].gid;
		if (a->shadow == 0) flag[gid] |= 1;
		else flag[gid] |= 2;
	}
	for (i = 0; i < d->n_gene; ++i) {
		if (flag[i]&1) cnt[i] += 1ULL<<32;
		else if (flag[i]&2) cnt[i] += 1ULL;
	}
	kfree(km, flag);
}

void pg_gen_vertex(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int32_t i;
	uint64_t *cnt;
	cnt = PG_CALLOC(uint64_t, d->n_gene);
	for (i = 0; i < d->n_genome; ++i)
		pg_gen_vertex_aux(0, opt, d, &d->genome[i], cnt);
	for (i = 0; i < d->n_gene; ++i) {
		int32_t pri = cnt[i]>>32, sec = (int32_t)cnt[i];
		if (pri >= d->n_genome * opt->min_vertex_ratio) {
			pg_vertex_t *p;
			PG_EXTEND(pg_vertex_t, q->v, q->n_v, q->m_v);
			p = &q->v[q->n_v++];
			p->gid = i, p->pri = pri, p->sec = sec;
		}
	}
	free(cnt);
	q->g2v = PG_MALLOC(int32_t, d->n_gene);
	for (i = 0; i < d->n_gene; ++i)
		q->g2v[i] = -1;
	for (i = 0; i < q->n_v; ++i)
		q->g2v[q->v[i].gid] = i;
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] selected %d vertices out of %d genes\n", __func__, pg_timestamp(), q->n_v, q->d->n_gene);
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

void pg_graph_add_arc(const pg_opt_t *opt, pg_graph_t *q)
{
}
