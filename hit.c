#include "pgpriv.h"
#include "kalloc.h"
#include "ksort.h"

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(pg128x, pg128_t, sort_key_128x, 8) 

int64_t pg_hit_cal_cm(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len, half;
	for (i = 0, len = 0; i < a->n_exon; ++i)
		len += e[i].oen - e[i].ost;
	half = len>>1;
	for (i = 0, len = 0; i < a->n_exon; ++i) {
		if (len <= half && half < len + e[i].oen - e[i].ost)
			return a->cs + e[i].ost + half - len;
		len += e[i].oen - e[i].ost;
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
