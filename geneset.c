#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "pgpriv.h"
#include "kalloc.h"

static inline int32_t pg_cds_len(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len = 0;
	for (i = 0; i < a->n_exon; ++i)
		len += e[a->off_exon + i].oe - e[a->off_exon + i].os;
	return len;
}

uint64_t *pg_gs_overlap(void *km, const pg_opt_t *opt, const pg_data_t *d, int32_t aid, int32_t *n_ov)
{
	const pg_genome_t *g = &d->genome[aid];
	int32_t i0, i, j;
	int32_t np = 0, mp = 0, k;
	pg128_t *pairs = 0;
	uint64_t *ret;
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		while (i0 < i) {
			if (g->hit[i0].cid == g->hit[i].cid && g->hit[i0].ce > g->hit[i].cs)
				break;
			++i0;
		}
		for (j = i0; j < i; ++j) {
			uint64_t x;
			int32_t li, lj, gi, gj;
			double cov_short;
			pg128_t *p;
			if (g->hit[i].rank > 0 || g->hit[j].rank > 0) continue; // only consider the top hit
			if (g->hit[j].ce <= g->hit[i].cs) continue;
			gi = d->prot[g->hit[j].pid].gid;
			gj = d->prot[g->hit[i].pid].gid;
			if (gi == gj) continue; // same gene
			x = pg_hit_overlap(g, &g->hit[j], &g->hit[i]);
			lj = pg_cds_len(&g->hit[j], g->exon);
			li = pg_cds_len(&g->hit[i], g->exon);
			cov_short = (double)(x>>32) / (li < lj? li : lj);
			assert(cov_short <= 1.0);
			if (cov_short < opt->min_ov_ratio) continue; // overlap too short
			Kexpand(km, pg128_t, pairs, np, mp);
			p = &pairs[np++];
			p->x = (uint64_t)(gi < gj? gi : gj) << 32 | (gi < gj? gj : gi);
			p->y = x;
		}
	}
	radix_sort_pg128x(pairs, pairs + np);
	for (i = 1, i0 = 0, k = 0; i <= np; ++i) {
		if (i == np || pairs[i].x != pairs[i0].x) {
			uint64_t max = 0;
			for (j = i0; j < i; ++j)
				max = max > pairs[j].y? max : pairs[j].y;
			pairs[k].x = pairs[i0].x;
			pairs[k++].y = max;
			i0 = i;
		}
	}
	ret = Kmalloc(km, uint64_t, k * 2);
	for (i = 0, j = 0; i < k; ++i) {
		ret[j++] = pairs[i].x;
		ret[j++] = pairs[i].x<<32 | pairs[i].x>>32;
	}
	kfree(km, pairs);
	radix_sort_pg64(ret, ret + k * 2);
	*n_ov = k * 2;
	return ret;
}

uint64_t *pg_gs_gene_score(void *km, const pg_data_t *d, int32_t aid)
{
	const pg_genome_t *g = &d->genome[aid];
	int32_t i;
	uint64_t *sc;
	sc = Kcalloc(km, uint64_t, d->n_gene);
	for (i = 0; i < g->n_hit; ++i) {
		const pg_hit_t *a = &g->hit[i];
		int32_t score, gid;
		score = (int32_t)(pow(a->score, (double)a->mlen / a->blen) + 1.0); // then score is at least 1
		gid = d->prot[a->pid].gid;
		sc[gid] = sc[gid] > score? sc[gid] : score;
	}
	return sc;
}

void pg_gs_choose1(void *km, const pg_opt_t *opt, const pg_data_t *d, int32_t aid, uint64_t *cnt)
{
	int32_t i0, i, j, off, n_ov;
	uint64_t *ov, *idx, *sc;
	int8_t *flag;
	const pg_genome_t *g = &d->genome[aid];

	ov = pg_gs_overlap(km, opt, d, aid, &n_ov);
	idx = Kcalloc(km, uint64_t, d->n_gene);
	for (i = 1, i0 = 0; i <= n_ov; ++i) {
		if (i == n_ov || ov[i]>>32 != ov[i0]>>32) {
			assert(ov[i0]>>32 < d->n_gene && (int32_t)ov[i0] < d->n_gene);
			idx[ov[i0]>>32] = i - i0;
			i0 = i;
		}
	}
	for (i = 0, off = 0; i < d->n_gene; ++i) {
		idx[i] = (uint64_t)off<<32 | idx[i];
		off += (int32_t)idx[i];
	}

	sc = pg_gs_gene_score(km, d, aid);
	for (i = 0; i < d->n_gene; ++i)
		sc[i] = sc[i]<<32 | i;
	radix_sort_pg64(sc, sc + d->n_gene);
	flag = Kcalloc(km, int8_t, d->n_gene);
	for (i = 0; i < g->n_hit; ++i)
		flag[d->prot[g->hit[i].pid].gid] |= 1;
	for (i = d->n_gene - 1; i >= 0; --i) {
		int32_t gid = (int32_t)sc[i];
		if (!(flag[gid] & 1)) continue; // no alignment; then don't count
		if (!(flag[gid] & 2)) {
			int32_t off = idx[gid]>>32, c = (int32_t)idx[gid];
			cnt[gid] += 1ULL<<32;
			for (j = 0; j < c; ++j)
				flag[(int32_t)ov[off + j]] |= 2;
		} else cnt[gid] += 1ULL;
	}
	kfree(km, flag);

	kfree(km, sc);
	kfree(km, idx);
	kfree(km, ov);
}

void pg_gen_vertex(const pg_opt_t *opt, pg_graph_t *q)
{
	const pg_data_t *d = q->d;
	int32_t i;
	uint64_t *cnt;
	cnt = PG_CALLOC(uint64_t, d->n_gene);
	for (i = 0; i < d->n_genome; ++i)
		pg_gs_choose1(0, opt, d, i, cnt);
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
		fprintf(stderr, "[M::%s::%.3f*%.2f] selected %d vertices out of %d genes\n", __func__, pg_realtime(), pg_percent_cpu(), q->n_v, q->d->n_gene);
}
