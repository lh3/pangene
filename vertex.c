#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kalloc.h"

void pg_gen_vtx(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int8_t *flag;
	int32_t i, j, n_pair = 0, m_pair = 0;
	pg128_t *cnt, *pair = 0, *p;

	flag = PG_CALLOC(int8_t, d->n_gene);
	cnt = PG_MALLOC(pg128_t, d->n_gene);
	for (i = 0; i < d->n_gene; ++i)
		cnt[i].x = pg_hash_uint32(i), cnt[i].y = i;

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
				if (opt->flag & PG_F_MERGE_ORTHO) {
					PG_EXTEND(pg128_t, pair, n_pair, m_pair);
					//if (strcmp(d->gene[gid].name, "0_0_3569") == 0) printf("Y\t%d\t%s\n", j, d->gene[d->prot[a->pid_dom].gid].name);
					p = &pair[n_pair++];
					p->x = (uint64_t)gid<<32 | d->prot[a->pid_dom].gid;
					p->y = 0;
				}
			} else {
				flag[gid] |= 1;
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
		uint64_t *idx;
		int32_t i0, k, *masked;
		pg_seg_t *tmp;

		// collapse identical pairs
		radix_sort_pg128x(pair, pair + n_pair);
		for (i0 = 0, i = 1, k = 0; i <= n_pair; ++i)
			if (i == n_pair || pair[i0].x != pair[i].x)
				pair[k].x = pair[i0].x, pair[k++].y = i - i0, i0 = i;
		n_pair = k;

		// index pairs
		idx = PG_CALLOC(uint64_t, d->n_gene);
		for (i0 = 0, i = 1; i <= n_pair; ++i)
			if (i == n_pair || pair[i0].x>>32 != pair[i].x>>32)
				idx[pair[i0].x>>32] = (uint64_t)i0<<32 | (i - i0), i0 = i;

		// generate segments
		masked = PG_CALLOC(int32_t, d->n_gene);
		radix_sort_pg128x(cnt, cnt + d->n_gene);
		for (i = d->n_gene - 1; i >= 0; --i) {
			int32_t n_dom = cnt[i].x>>32, n_sub = cnt[i].y>>32;
			int32_t gid = (int32_t)cnt[i].y;
			//printf("X\t%s\t%d\t%d\t%d\n", d->gene[gid].name, n_dom, n_sub, masked[gid]);
			if (n_dom >= d->n_genome * opt->min_vertex_ratio && n_dom > masked[gid]) {
				pg_seg_t *p;
				int32_t off = idx[gid]>>32, n = (int32_t)idx[gid];
				PG_EXTEND0(pg_seg_t, q->seg, q->n_seg, q->m_seg);
				p = &q->seg[q->n_seg++];
				p->gid = gid, p->n_dom = n_dom, p->n_sub = n_sub;
				for (j = 0; j < n; ++j)
					masked[(int32_t)pair[off + j].x] += pair[off + j].y;
			}
		}
		free(masked);

		// sort by pg_seg_t::gid
		for (i = 0; i < q->n_seg; ++i)
			idx[i] = (uint64_t)q->seg[i].gid<<32 | i; // reuse idx[] for sorting
		radix_sort_pg64(idx, idx + q->n_seg);
		tmp = PG_CALLOC(pg_seg_t, q->n_seg);
		for (i = 0; i < q->n_seg; ++i)
			tmp[i] = q->seg[(int32_t)idx[i]];
		free(idx);
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

	free(pair);
	free(cnt);
	pg_gen_g2s(q);
	if (pg_verbose >= 3)
		fprintf(stderr, "[M::%s::%s] selected %d vertices out of %d genes\n", __func__, pg_timestamp(), q->n_seg, q->d->n_gene);
}
