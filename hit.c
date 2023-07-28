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

int32_t pg_flt_subopt_isoform(const pg_prot_t *prot, int32_t n_gene, pg_genome_t *g)
{
	uint64_t *best;
	int32_t i, gid, n_flt = 0;
	best = PG_CALLOC(uint64_t, n_gene);
	for (i = 0; i < g->n_hit; ++i) {
		const pg_hit_t *a = &g->hit[i];
		if (a->flt || a->rank > 0) continue;
		gid = prot[a->pid].gid;
		if (a->score > best[gid]>>32)
			best[gid] = (uint64_t)a->score << 32 | a->pid;
	}
	for (i = 0; i < g->n_hit; ++i) {
		pg_hit_t *a = &g->hit[i];
		if (a->flt || a->rank > 0) continue;
		gid = prot[a->pid].gid;
		if (a->pid != (int32_t)best[gid])
			a->flt = a->flt_iso_sub_self = 1, ++n_flt;
	}
	free(best);
	return n_flt;
}

int32_t pg_flt_subopt_joint(const pg_opt_t *opt, pg_data_t *d) // call after pg_flt_ov_isoform()
{
	int32_t i, j, gid, n_flt = 0, *best;
	best = PG_CALLOC(int32_t, d->n_gene);
	for (j = 0; j < d->n_genome; ++j) {
		const pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i) {
			const pg_hit_t *a = &g->hit[i];
			if (a->flt) continue;
			gid = d->prot[a->pid].gid;
			if (a->score > best[gid])
				best[gid] = a->score;
		}
	}
	for (j = 0; j < d->n_genome; ++j) {
		pg_genome_t *g = &d->genome[j];
		for (i = 0; i < g->n_hit; ++i) {
			pg_hit_t *a = &g->hit[i];
			if (a->flt) continue;
			gid = d->prot[a->pid].gid;
			if (a->score < (double)best[gid] * (1.0 - opt->max_div))
				a->flt = a->flt_iso_sub_joint = 1, ++n_flt;
		}
	}
	free(best);
	return n_flt;
}

int32_t pg_flt_chain_shadow(const pg_prot_t *prot, int32_t n_prot, pg_genome_t *g)
{
	int32_t i, n_flt = 0;
	int8_t *flag;
	flag = PG_CALLOC(int8_t, n_prot);
	for (i = 0; i < n_prot; ++i) flag[i] = 1;
	for (i = 0; i < g->n_hit; ++i)
		if (!g->hit[i].flt_iso_ov)
			flag[g->hit[i].pid] = 0;
	for (i = 0; i < g->n_hit; ++i) {
		pg_hit_t *a = &g->hit[i];
		if (a->pid_dom0 >= 0 && flag[a->pid_dom0])
			a->flt = a->flt_chain = 1, ++n_flt;
	}
	free(flag);
	return n_flt;
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
				z[a->pid].x += (uint64_t)a->score<<32 | 1; // NB: assuming each protein has only one rank=0 hit
			a->rep = 0;
		}
	}
	for (i = 0; i < d->n_prot; ++i) {
		d->prot[i].n = (uint32_t)z[i].x;
		d->prot[i].avg_score = d->prot[i].n? (int32_t)((double)(z[i].x>>32) / d->prot[i].n + .499) : 0;
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
