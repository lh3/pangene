#include <stdio.h>
#include <assert.h>
#include "pgpriv.h"

// Compute the CDS overlap length and union length. The union length is actually not used at the moment.
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

// Compute the CDS length of an alignment
static inline int32_t pg_cds_len(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len = 0;
	for (i = 0; i < a->n_exon; ++i)
		len += e[a->off_exon + i].oe - e[a->off_exon + i].os;
	return len;
}

/*
 * All rountines below follow a similar logic and have some code duplications.
 */

// select among overlapping isoforms of the same gene
int32_t pg_flt_ov_isoform(const pg_opt_t *opt, pg_data_t *d, int32_t id)
{
	const pg_prot_t *prot = d->prot;
	pg_genome_t *g = &d->genome[id];
	int32_t i, j, i0, n_flt = 0, gi, gj;
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		uint32_t hi;
		if (ai->flt) continue;
		while (i0 < i && !(g->hit[i0].cid == ai->cid && g->hit[i0].ce > ai->cs)) // update i0
			++i0;
		gi = prot[ai->pid].gid;
		hi = pg_hash_uint32(ai->pid);
		for (j = i0; j < i; ++j) {
			uint32_t hj;
			uint64_t x, si, sj;
			pg_hit_t *aj = &g->hit[j];
			if (aj->flt || aj->ce <= ai->cs) continue; // no overlap
			gj = prot[aj->pid].gid;
			if (gi != gj) continue; // ignore isoforms from different genes
			hj = pg_hash_uint32(aj->pid);
			x = pg_hit_overlap(g, aj, ai);
			if (x>>32 == 0) continue; // no overlap on CDS
			si = (uint64_t)ai->score_adj<<33 | (uint64_t)d->gene[gi].preferred<<32 | hi;
			sj = (uint64_t)aj->score_adj<<33 | (uint64_t)d->gene[gj].preferred<<32 | hj;
			if (si < sj || (si == sj && ai->rank > aj->rank))
				ai->flt_iso_ov = 1;
			else aj->flt_iso_ov = 1;
		}
	}
	for (i = 0; i < g->n_hit; ++i)
		if (g->hit[i].flt_iso_ov)
			g->hit[i].flt = 1, ++n_flt;
	return n_flt;
}

// test overlap between same or different genes
int32_t pg_shadow(const pg_opt_t *opt, pg_data_t *d, int32_t id)
{
	const pg_prot_t *prot = d->prot;
	pg_genome_t *g = &d->genome[id];
	int32_t i, i0, n_shadow = 0;
	pg128_t *tmp;
	tmp = PG_CALLOC(pg128_t, g->n_hit);
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		pg_hit_t *ai = &g->hit[i];
		int32_t j, li, gi;
		uint32_t hi;
		if (ai->flt) continue;
		ai->shadow = 0;
		while (i0 < i && !(g->hit[i0].cid == ai->cid && g->hit[i0].ce > ai->cs)) // update i0
			++i0;
		gi = prot[ai->pid].gid;
		hi = pg_hash_uint32(ai->pid);
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
			hj = pg_hash_uint32(aj->pid);
			x = pg_hit_overlap(g, aj, ai);
			if (x>>32 == 0) continue; // no overlap on CDS
			lj = pg_cds_len(aj, g->exon);
			cov_short = (double)(x>>32) / (li < lj? li : lj);
			assert(cov_short <= 1.0);
			if (gi != gj && cov_short < opt->min_ov_ratio) continue; // overlap too short
			si = (uint64_t)ai->score_adj<<33 | (uint64_t)d->gene[gi].preferred<<32 | hi;
			sj = (uint64_t)aj->score_adj<<33 | (uint64_t)d->gene[gj].preferred<<32 | hj;
			if (gi == gj) { // don't consider weak_br for different isoforms of the same gene
				shadow = (si < sj || (si == sj && ai->rank > aj->rank))? 0 : 1; // 0 for i and 1 for j
			} else if (ai->weak_br == aj->weak_br) {
				shadow = (si < sj || (si == sj && ai->rank > aj->rank))? 0 : 1; // 0 for i and 1 for j
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
