#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kalloc.h"

void pg_gen_vtx(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int8_t *flag;
	int32_t i, j;
	uint32_t **aux = 0;
	pg128_t *cnt;

	flag = PG_CALLOC(int8_t, d->n_gene);
	cnt = PG_MALLOC(pg128_t, d->n_gene);
	for (i = 0; i < d->n_gene; ++i)
		cnt[i].x = d->prot[d->gene[i].pri_pid].avg_score2, cnt[i].y = i;

	if (opt->flag & PG_F_MERGE_ORTHO) {
		aux = PG_CALLOC(uint32_t*, d->n_genome);
		for (j = 0; j < d->n_genome; ++j) {
			aux[j] = PG_MALLOC(uint32_t, d->n_gene);
			for (i = 0; i < d->n_gene; ++i)
				aux[j][i] = (d->n_gene + 1) << 1;
		}
	}

	for (j = 0; j < d->n_genome; ++j) {
		const pg_genome_t *g = &d->genome[j];
		memset(flag, 0, sizeof(*flag) * d->n_gene);
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *a = &g->hit[i];
			int32_t gid;
			if (a->rank != 0 || a->pri == 0) continue;
			gid = d->prot[a->pid].gid;
			if (a->shadow) {
				assert(a->pid_dom >= 0);
				flag[gid] |= 2;
				if (aux) aux[j][gid] = (uint32_t)d->prot[a->pid_dom].gid<<1;
			} else {
				flag[gid] |= 1;
				if (aux) aux[j][gid] = (uint32_t)d->n_gene<<1;
			}
		}
		for (i = 0; i < d->n_gene; ++i) {
			assert(flag[i] != 3);
			if (flag[i]&1) cnt[i].x += 1ULL<<32;
			else if (flag[i]&2) cnt[i].y += 1ULL<<32;
		}
	}
	free(flag);

	if (opt->flag & PG_F_MERGE_ORTHO) {
		pg_seg_t *tmp;
		uint64_t *srt;

		// generate segments
		radix_sort_pg128x(cnt, cnt + d->n_gene);
		for (i = d->n_gene - 1; i >= 0; --i) {
			int32_t n_dom = cnt[i].x>>32, n_sub = cnt[i].y>>32;
			int32_t gid = (int32_t)cnt[i].y, x = 0, y = 0;
			for (j = 0; j < d->n_genome; ++j)
				if (aux[j][gid]>>1 == d->n_gene)
					++x, y += (aux[j][gid]&1);
			//printf("X\t%s\t%d\t%d\t%d\t%d\t%d\n", d->gene[gid].name, n_dom, n_sub, (int32_t)cnt[i].x, x, y);
			if (n_dom >= d->n_genome * opt->min_vertex_ratio && y < x) {
				pg_seg_t *p;
				PG_EXTEND0(pg_seg_t, q->seg, q->n_seg, q->m_seg);
				p = &q->seg[q->n_seg++];
				p->gid = gid, p->n_dom = n_dom, p->n_sub = n_sub;
				if (x > 0) {
					for (j = 0; j < d->n_genome; ++j)
						if (aux[j][gid]>>1 < d->n_gene)
							aux[j][aux[j][gid]>>1] |= 1;
				}
			}
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
	} else {
		for (i = 0; i < d->n_gene; ++i) {
			int32_t n_dom = cnt[i].x>>32, n_sub = cnt[i].y>>32;
			if (n_dom >= d->n_genome * opt->min_vertex_ratio) {
				pg_seg_t *p;
				PG_EXTEND0(pg_seg_t, q->seg, q->n_seg, q->m_seg);
				p = &q->seg[q->n_seg++];
				p->gid = i, p->n_dom = n_dom, p->n_sub = n_sub;
			}
		}
	}

	free(cnt);
	pg_gen_g2s(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] selected %d vertices out of %d genes\n", __func__, pg_timestamp(), q->n_seg, q->d->n_gene);
}
