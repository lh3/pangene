#include <assert.h>
#include <stdio.h>
#include "pgpriv.h"
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

void pg_hit_sort(pg_genome_t *g, int32_t by_cm)
{
	int32_t i, *n, *off;
	pg_hit_t *a;
	pg128_t *tmp;

	n = PG_CALLOC(int32_t, g->n_ctg);
	for (i = 0; i < g->n_hit; ++i)
		++n[g->hit[i].cid];
	off = PG_CALLOC(int32_t, g->n_ctg + 1);
	for (i = 1; i <= g->n_ctg; ++i)
		off[i] = off[i - 1] + n[i - 1];

	tmp = PG_MALLOC(pg128_t, g->n_hit);
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
	free(off);
	free(n);

	a = PG_MALLOC(pg_hit_t, g->n_hit);
	for (i = 0; i < g->n_hit; ++i)
		a[i] = g->hit[tmp[i].y];

	free(tmp);
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

int32_t pg_flag_pseudo(const pg_prot_t *prot, pg_genome_t *g)
{
	int32_t i, i0, j, n_pseudo = 0;
	pg128_t *a;
	a = PG_MALLOC(pg128_t, g->n_hit);
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
	free(a);
	return n_pseudo;
}

typedef struct {
	int32_t c[2];
	int64_t s[2];
} pseudo_joint_aux_t;

int32_t pg_flag_pseudo_joint(const pg_opt_t *opt, pg_data_t *d) // call after pg_flag_pseudo()
{
	int32_t i, j, n_pseudo = 0;
	pseudo_joint_aux_t *aux;
	aux = PG_CALLOC(pseudo_joint_aux_t, d->n_prot);
	for (j = 0; j < d->n_genome; ++j) {
		const pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *a = &g->hit[i];
			if (a->rank == 0) {
				int32_t w = a->n_exon == 1? 0 : 1;
				aux[a->pid].c[w]++;
				aux[a->pid].s[w] += a->score;
			}
		}
	}
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i) {
			pg_hit_t *a = &g->hit[i];
			pseudo_joint_aux_t *p = &aux[a->pid];
			if (a->n_exon == 1 && p->c[1] >= d->n_genome * opt->min_vertex_ratio && !a->pseudo
				&& p->s[1] * p->c[0] > p->s[0] * p->c[1]) // multi-exon should have higher score
			{
				a->pseudo = 1, ++n_pseudo;
			}
		}
	}
	free(aux);
	return n_pseudo;
}

static inline int32_t pg_cds_len(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len = 0;
	for (i = 0; i < a->n_exon; ++i)
		len += e[a->off_exon + i].oe - e[a->off_exon + i].os;
	return len;
}

int32_t pg_flag_shadow(const pg_opt_t *opt, const pg_prot_t *prot, pg_genome_t *g, int32_t cmp_iso)
{
	int32_t i, i0, n_shadow = 0;
	pg128_t *tmp;
	tmp = PG_CALLOC(pg128_t, g->n_hit);
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		int32_t j, li, gi;
		uint32_t hi;
		if (ai->flt) continue;
		ai->overlap = ai->shadow = 0;
		while (i0 < i && !(g->hit[i0].cid == ai->cid && g->hit[i0].ce > ai->cs)) // update i0
			++i0;
		gi = prot[ai->pid].gid;
		hi = pg_hash_uint32(gi);
		li = pg_cds_len(ai, g->exon);
		for (j = i0; j < i; ++j) {
			uint64_t x;
			int32_t lj, gj, shadow = -1;
			double cov_short;
			uint32_t hj;
			uint64_t si, sj;
			pg_hit_t *aj = &g->hit[j];
			if (aj->ce <= ai->cs) continue; // no overlap
			if (aj->flt) continue;
			gj = prot[aj->pid].gid;
			hj = pg_hash_uint32(gj);
			if (!cmp_iso && gi == gj && ai->pid != aj->pid) continue; // ignore iso-forms of the same gene
			x = pg_hit_overlap(g, aj, ai);
			lj = pg_cds_len(aj, g->exon);
			cov_short = (double)(x>>32) / (li < lj? li : lj);
			assert(cov_short <= 1.0);
			if ((ai->pid != aj->pid && cov_short < opt->min_ov_ratio) || (ai->pid == aj->pid && x>>32 == 0)) continue; // overlap too short
			si = (uint64_t)ai->score2<<32 | hi;
			sj = (uint64_t)aj->score2<<32 | hj;
			ai->overlap = aj->overlap = 1;
			if (ai->weak_br == aj->weak_br) {
				shadow = (si < sj || (si == sj && ai->pid > aj->pid))? 0 : 1; // 0 for i and 1 for j
			} else if (ai->weak_br > aj->weak_br) { // i is worse
				shadow = 0;
			} else { // j is worse
				shadow = 1;
			}
			if (shadow == 0) { // mark i
				ai->shadow = 1;
				if (tmp[i].y < sj) tmp[i].y = sj, tmp[i].x = aj->pid;
			} else { // mark j
				aj->shadow = 1;
				if (tmp[j].y < si) tmp[j].y = si, tmp[j].x = ai->pid;
			}
		}
	}
	for (i = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		if (ai->flt) continue;
		ai->pid_dom = tmp[i].y == 0? -1 : tmp[i].x;
		if (ai->shadow) ++n_shadow;
	}
	free(tmp);
	return n_shadow;
}

void pg_flag_representative(pg_data_t *d) // flag representative isoform
{
	pg128_t *z;
	int32_t i, j;
	z = PG_CALLOC(pg128_t, d->n_prot);
	for (i = 0; i < d->n_gene; ++i) d->gene[i].rep_pid = -1;
	for (i = 0; i < d->n_prot; ++i) z[i].y = i, d->prot[i].rep = 0;
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i) {
			pg_hit_t *a = &g->hit[i];
			if (a->rank == 0 && a->flt == 0)
				z[a->pid].x += (uint64_t)a->score2<<32 | 1; // NB: assuming each protein has only one rank=0 hit
			a->rep = 0;
		}
	}
	for (i = 0; i < d->n_prot; ++i) {
		d->prot[i].n = (uint32_t)z[i].x;
		d->prot[i].avg_score2 = d->prot[i].n? (int32_t)((double)(z[i].x>>32) / d->prot[i].n + .499) : 0;
	}
	radix_sort_pg128x(z, z + d->n_prot);
	for (i = d->n_prot - 1; i >= 0; --i) {
		int32_t pid = z[i].y;
		int32_t gid = d->prot[pid].gid;
		if (d->gene[gid].rep_pid < 0) {
			d->gene[gid].rep_pid = pid;
			d->prot[pid].rep = 1;
		}
	}
	free(z);
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i)
			if (d->prot[g->hit[i].pid].rep)
				g->hit[i].rep = 1;
	}
}
