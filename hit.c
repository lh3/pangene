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
	//fprintf(stderr, "%d,%c%c: [%ld,%ld) <=> [%ld,%ld): %d,%d\n", aa->cid, "+-"[aa->rev], "+-"[ab->rev], (long)aa->cs, (long)aa->ce, (long)ab->cs, (long)ab->ce, l_inter, l_union);
	return (uint64_t)l_inter<<32 | l_union;
}
