#include <string.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kalloc.h"

static void pg_gen_vtx_algo1_aux(void *km, const pg_opt_t *opt, const pg_data_t *d, pg_genome_t *g, uint64_t *cnt)
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

static void pg_gen_vtx_algo1(const pg_opt_t *opt, pg_graph_t *q)
{
	pg_data_t *d = q->d;
	int32_t i;
	uint64_t *cnt;
	cnt = PG_CALLOC(uint64_t, d->n_gene);
	for (i = 0; i < d->n_genome; ++i) {
		pg_flag_shadow(0, opt, d, &d->genome[i], 0, 0, 0);
		pg_gen_vtx_algo1_aux(0, opt, d, &d->genome[i], cnt);
	}
	for (i = 0; i < d->n_gene; ++i) {
		int32_t pri = cnt[i]>>32, sec = (int32_t)cnt[i];
		if (pri >= d->n_genome * opt->min_vertex_ratio) {
			pg_seg_t *p;
			PG_EXTEND0(pg_seg_t, q->seg, q->n_seg, q->m_seg);
			p = &q->seg[q->n_seg++];
			p->gid = i, p->cnt_dom = pri, p->cnt_sub = sec;
		}
	}
	free(cnt);
	pg_gen_g2s(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] selected %d vertices out of %d genes\n", __func__, pg_timestamp(), q->n_seg, q->d->n_gene);
}

void pg_gen_vtx(const pg_opt_t *opt, pg_graph_t *q)
{
	pg_gen_vtx_algo1(opt, q);
}
