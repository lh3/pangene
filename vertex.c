#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"

void pg_gen_vtx(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int8_t *flag;
	int32_t i, j;
	uint32_t **aux = 0;
	pg128_t *cnt;
	pg_seg_t *tmp;
	uint64_t *srt;

	flag = PG_CALLOC(int8_t, d->n_gene);
	cnt = PG_MALLOC(pg128_t, d->n_gene);
	for (i = 0; i < d->n_gene; ++i)
		cnt[i].x = d->prot[d->gene[i].rep_pid].avg_score, cnt[i].y = i;

	aux = PG_CALLOC(uint32_t*, d->n_genome);
	for (j = 0; j < d->n_genome; ++j) {
		aux[j] = PG_MALLOC(uint32_t, d->n_gene);
		for (i = 0; i < d->n_gene; ++i)
			aux[j][i] = (d->n_gene + 1) << 1;
	}

	for (j = 0; j < d->n_genome; ++j) {
		const pg_genome_t *g = &d->genome[j];
		memset(flag, 0, sizeof(*flag) * d->n_gene);
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *a = &g->hit[i];
			int32_t gid;
			if (a->rank != 0 || a->flt) continue;
			gid = d->prot[a->pid].gid;
			//if (strcmp(d->gene[gid].name, "KIR3DS1") == 0) fprintf(stderr, "X\t%s\t%s\t%d\t%d\t%s\n", g->label, d->gene[gid].name, a->score, a->shadow, a->pid_dom >= 0? d->prot[a->pid_dom].name : "*");
			if (a->shadow) {
				assert(a->pid_dom >= 0);
				flag[gid] |= 2;
				if (aux[j][gid] == (uint32_t)(d->n_gene + 1) << 1)
					aux[j][gid] = (uint32_t)d->prot[a->pid_dom].gid<<1;
			} else {
				flag[gid] |= 1;
				aux[j][gid] = (uint32_t)d->n_gene<<1;
			}
		}
		for (i = 0; i < d->n_gene; ++i) {
			if (flag[i]&1) cnt[i].x += 1ULL<<32; // it is possible that flag[i]==3; in this case, we count the non-shadowed hit first
			else if (flag[i]&2) cnt[i].y += 1ULL<<32;
		}
	}
	free(flag);

	for (j = 0; j < d->n_gene; ++j)
		if (d->gene[j].preferred)
			cnt[j].x |= 1ULL<<63;

	// generate segments
	radix_sort_pg128x(cnt, cnt + d->n_gene);
	for (i = d->n_gene - 1; i >= 0; --i) {
		int32_t n_dom = cnt[i].x<<1>>33, n_sub = cnt[i].y>>32;
		int32_t gid = (int32_t)cnt[i].y, x = 0, y = 0;
		for (j = 0; j < d->n_genome; ++j)
			if (aux[j][gid]>>1 == d->n_gene)
				++x, y += (aux[j][gid]&1);
		if (opt->flag & PG_F_WRITE_VTX_SEL)
			printf("g\t%s\t%d\t%d\t%d\t%d\t%c\n", d->gene[gid].name, (int32_t)cnt[i].x, x, y, n_sub, "NY"[d->gene[gid].included]);
		if (d->gene[gid].included || (n_dom >= d->n_genome * opt->min_vertex_ratio && y < x)) {
			pg_seg_t *p;
			PG_GROW0(pg_seg_t, q->seg, q->n_seg, q->m_seg);
			p = &q->seg[q->n_seg++];
			p->gid = gid, p->n_dom = n_dom, p->n_sub = n_sub;
			if (x > 0) {
				for (j = 0; j < d->n_genome; ++j)
					if (aux[j][gid]>>1 < d->n_gene)
						aux[j][aux[j][gid]>>1] |= 1;
			}
		}
		//else fprintf(stderr, "Y\t%s\tn_dom=%d\n", d->gene[gid].name, n_dom);
	}
	for (j = 0; j < d->n_genome; ++j) free(aux[j]);
	free(aux);

	// sort by pg_seg_t::gid
	srt = PG_CALLOC(uint64_t, q->n_seg);
	for (i = 0; i < q->n_seg; ++i)
		srt[i] = (uint64_t)q->seg[i].gid<<32 | i;
	radix_sort_pg64(srt, srt + q->n_seg);
	tmp = PG_CALLOC(pg_seg_t, q->n_seg);
	for (i = 0; i < q->n_seg; ++i)
		tmp[i] = q->seg[(int32_t)srt[i]];
	free(srt);
	free(q->seg);
	q->seg = tmp, q->m_seg = q->n_seg;

	free(cnt);
	pg_gen_g2s(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] selected %d vertices out of %d genes\n", __func__, pg_timestamp(), q->n_seg, q->d->n_gene);
}
