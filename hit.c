#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kalloc.h"
#include "ksort.h"

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(pg128x, pg128_t, sort_key_128x, 8) 

#define sort_key_generic(a) (a)
KRADIX_SORT_INIT(pg64, uint64_t, sort_key_generic, 8) 

int64_t pg_hit_cal_cm(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len, half;
	for (i = 0, len = 0; i < a->n_exon; ++i)
		len += e[i].oe - e[i].os;
	half = len>>1;
	for (i = 0, len = 0; i < a->n_exon; ++i) {
		if (len <= half && half < len + e[i].oe - e[i].os)
			return a->cs + e[i].os + half - len;
		len += e[i].oe - e[i].os;
	}
	abort();
	return -1;
}

void pg_hit_sort(void *km, pg_genome_t *g, int32_t by_cm)
{
	int32_t i, *n, *off;
	pg_hit_t *a;
	pg128_t *tmp;

	n = Kcalloc(km, int32_t, g->n_ctg);
	for (i = 0; i < g->n_hit; ++i)
		++n[g->hit[i].cid];
	off = Kcalloc(km, int32_t, g->n_ctg + 1);
	for (i = 1; i <= g->n_ctg; ++i)
		off[i] = off[i - 1] + n[i - 1];

	tmp = Kmalloc(km, pg128_t, g->n_hit);
	for (i = 0; i < g->n_hit; ++i) {
		pg128_t t;
		t.x = by_cm? g->hit[i].cm : g->hit[i].cs;
		t.y = i;
		tmp[off[g->hit[i].cid]++] = t;
	}

	for (i = 1, off[0] = 0; i <= g->n_ctg; ++i)
		off[i] = off[i - 1] + n[i - 1];
	for (i = 0; i < g->n_ctg; ++i)
		radix_sort_pg128x(&tmp[off[i]], &tmp[off[i+1]]);
	kfree(km, off);
	kfree(km, n);

	a = PG_MALLOC(pg_hit_t, g->n_hit);
	for (i = 0; i < g->n_hit; ++i)
		a[i] = g->hit[tmp[i].y];

	kfree(km, tmp);
	free(g->hit);
	g->hit = a;
}

uint64_t pg_hit_overlap(const pg_genome_t *g, const pg_hit_t *aa, const pg_hit_t *ab)
{
	int32_t l_inter = 0, l_union = 0, x;
	int64_t z = 0;
	const pg_hit_t *a[2];
	const pg_exon_t *e[2], *ee[2];
	if (aa->cid != ab->cid || !(aa->cs < ab->ce && aa->ce > ab->cs)) return 0;
	assert(aa->n_exon > 0 && ab->n_exon > 0);
	a[0] = aa, a[1] = ab;
	e[0] = &g->exon[a[0]->off_exon], e[1] = &g->exon[a[1]->off_exon];
	ee[0] = e[0] + a[0]->n_exon, ee[1] = e[1] + a[1]->n_exon;
	while (e[0] < ee[0] && e[1] < ee[1]) {
		int32_t x = a[0]->cs + e[0]->os < a[1]->cs + e[1]->os? 0 : 1, y = !x;
		z = z > a[x]->cs + e[x]->os? z : a[x]->cs + e[x]->os;
		if (a[x]->cs + e[x]->oe < a[y]->cs + e[y]->oe) { // x ends earlier
			int64_t o = (a[x]->cs + e[x]->oe) - (a[y]->cs + e[y]->os);
			l_inter += o > 0? o : 0;
			l_union += a[x]->cs + e[x]->oe - z;
			z = a[x]->cs + e[x]->oe;
			++e[x];
		} else { // y contained in x
			l_inter += e[y]->oe - e[y]->os;
			l_union += a[y]->cs + e[y]->oe - z;
			z = a[y]->cs + e[y]->oe;
			++e[y];
		}
	}
	x = e[0] < ee[0]? 0 : 1;
	while (e[x] < ee[x]) {
		z = z > a[x]->cs + e[x]->os? z : a[x]->cs + e[x]->os;
		l_union += a[x]->cs + e[x]->oe - z;
		++e[x];
	}
	assert(l_inter <= l_union);
	//fprintf(stderr, "%d,%c%c - %d:[%ld,%ld) <=> %d:[%ld,%ld): %d,%d\n", aa->cid, "+-"[aa->rev], "+-"[ab->rev], aa->pid, (long)aa->cs, (long)aa->ce, ab->pid, (long)ab->cs, (long)ab->ce, l_inter, l_union);
	return (uint64_t)l_inter<<32 | l_union;
}

int32_t pg_flag_pseudo(void *km, const pg_prot_t *prot, pg_genome_t *g)
{
	int32_t i, i0, j, n_pseudo = 0;
	pg128_t *a;
	a = Kmalloc(km, pg128_t, g->n_hit);
	for (i = 0; i < g->n_hit; ++i) {
		a[i].x = (uint64_t)g->hit[i].pid<<32 | g->hit[i].rank;
		a[i].y = i;
	}
	radix_sort_pg128x(a, a + g->n_hit);
	for (i = 1, i0 = 0; i <= g->n_hit; ++i) {
		if (i == g->n_hit || a[i].x>>32 != a[i0].x>>32) {
			int32_t n_multi = 0;
			if (i - i0 > 1) {
				for (j = i0; j < i; ++j)
					if (g->hit[a[j].y].n_exon > 1)
						++n_multi;
			}
			if (n_multi > 0 && i - i0 - n_multi > 0) {
				int32_t j1 = -1;
				for (j = i0; j < i; ++j) {
					if (g->hit[a[j].y].n_exon == 1)
						g->hit[a[j].y].pseudo = 1, ++n_pseudo;
					else if (j1 < 0)
						j1 = j;
				}
				assert(j1 >= 0);
				if (g->hit[a[j1].y].rank > 0) { // promote the first multi-exon hit to rank 0
					for (j = i0; j < j1; ++j)
						g->hit[a[j].y].rank++;
					g->hit[a[j1].y].rank = 0;
				}
			}
			i0 = i;
		}
	}
	kfree(km, a);
	return n_pseudo;
}

static inline uint32_t kh_hash_uint32(uint32_t key)
{
	key += ~(key << 15);
	key ^=  (key >> 10);
	key +=  (key << 3);
	key ^=  (key >> 6);
	key += ~(key << 11);
	key ^=  (key >> 16);
	return key;
}

static inline int32_t pg_cds_len(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len = 0;
	for (i = 0; i < a->n_exon; ++i)
		len += e[a->off_exon + i].oe - e[a->off_exon + i].os;
	return len;
}

int32_t pg_flag_shadow(const pg_opt_t *opt, const pg_prot_t *prot, pg_genome_t *g, int32_t check_vtx, int32_t check_pri)
{
	int32_t i, i0, n_shadow = 0;
	for (i = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		if (check_vtx && ai->vtx == 0) continue;
		if (check_pri && ai->pri == 0) continue;
		ai->overlap = ai->shadow = 0;
	}
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		int32_t j, li, gi;
		uint32_t hi;
		if (check_vtx && ai->vtx == 0) continue;
		if (check_pri && ai->pri == 0) continue;
		while (i0 < i && !(g->hit[i0].cid == ai->cid && g->hit[i0].ce > ai->cs)) // update i0
			++i0;
		gi = prot[ai->pid].gid;
		hi = kh_hash_uint32(gi);
		li = pg_cds_len(ai, g->exon);
		for (j = i0; j < i; ++j) {
			uint64_t x;
			int32_t lj, gj;
			double cov_short;
			uint32_t hj;
			uint64_t si, sj;
			pg_hit_t *aj = &g->hit[j];
			if (check_vtx && aj->vtx == 0) continue;
			if (aj->ce <= ai->cs) continue; // no overlap
			gj = prot[aj->pid].gid;
			hj = kh_hash_uint32(gj);
			if (gi == gj) continue; // ignore iso-forms of the same gene
			x = pg_hit_overlap(g, aj, ai);
			lj = pg_cds_len(aj, g->exon);
			cov_short = (double)(x>>32) / (li < lj? li : lj);
			assert(cov_short <= 1.0);
			if (cov_short < opt->min_ov_ratio) continue; // overlap too short
			si = (uint64_t)ai->score2<<32 | hi;
			sj = (uint64_t)aj->score2<<32 | hj;
			ai->overlap = aj->overlap = 1;
			if (si < sj) ai->shadow = 1;
			else aj->shadow = 1;
		}
	}
	for (i = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		if (check_vtx && ai->vtx == 0) continue;
		if (check_pri && ai->pri == 0) continue;
		if (ai->shadow) ++n_shadow;
	}
	return n_shadow;
}

void pg_flag_primary(pg_data_t *d)
{
	pg128_t *z;
	int32_t i, j;
	z = PG_CALLOC(pg128_t, d->n_prot);
	for (i = 0; i < d->n_gene; ++i) d->gene[i].pri_pid = -1;
	for (i = 0; i < d->n_prot; ++i) z[i].y = i, d->prot[i].pri = 0;
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i) {
			if (g->hit[i].rank == 0)
				z[g->hit[i].pid].x += 1ULL<<32 | g->hit[i].score2; // NB: assuming each protein has only one rank=0 hit
			g->hit[i].pri = 0;
		}
	}
	radix_sort_pg128x(z, z + d->n_prot);
	for (i = d->n_prot - 1; i >= 0; --i) {
		int32_t pid = z[i].y;
		int32_t gid = d->prot[pid].gid;
		if (d->gene[gid].pri_pid < 0) {
			d->gene[gid].pri_pid = pid;
			d->prot[pid].pri = 1;
		}
	}
	free(z);
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i)
			if (d->prot[g->hit[i].pid].pri)
				g->hit[i].pri = 1;
	}
}
