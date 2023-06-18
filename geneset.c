#include <assert.h>
#include "pgpriv.h"
#include "kalloc.h"

static inline int32_t pg_cds_len(const pg_hit_t *a, const pg_exon_t *e)
{
	int32_t i, len = 0;
	for (i = 0; i < a->n_exon; ++i)
		len += e[a->off_exon + i].oe - e[a->off_exon + i].os;
	return len;
}

void pg_gs_build1(void *km, const pg_opt_t *opt, const pg_data_t *d, int32_t aid)
{
	const pg_genome_t *g = &d->genome[aid];
	int32_t i0, i, j;
	int32_t np = 0, mp = 0;
	pg128_t *pairs = 0;
	for (i = 1, i0 = 0; i < g->n_hit; ++i) {
		while (i0 < i) {
			if (g->hit[i0].ce > g->hit[i].cs)
				break;
			++i0;
		}
		for (j = i0; j < i; ++j) {
			uint64_t x;
			int32_t li, lj, gi, gj;
			double cov_short;
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
		}
	}
}
