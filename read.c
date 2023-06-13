#include <stdlib.h>
#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include "pgpriv.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

int pg_verbose = 3;

pg_data_t *pg_data_init(void)
{
	pg_data_t *d;
	d = PG_CALLOC(pg_data_t, 1);
	d->d_ctg = pg_dict_init();
	d->d_gene = pg_dict_init();
	d->d_prot = pg_dict_init();
	return d;
}

void pg_data_destroy(pg_data_t *d)
{
	int32_t i;
	for (i = 0; i < d->n_genome; ++i) {
		pg_genome_t *g = &d->genome[i];
		free(g->ctg); free(g->hit); free(g->exon);
	}
	free(d->genome); free(d->prot);
	pg_dict_destroy_copy(d->d_ctg);
	pg_dict_destroy_copy(d->d_gene);
	pg_dict_destroy_copy(d->d_prot);
}

typedef struct {
	int32_t n_exon, m_exon;
	pg_exon_t *exon;
} pg_exons_t;

static inline void pg_add_exon(pg_exons_t *tmp, int32_t st)
{
	pg_exon_t *p;
	PG_EXTEND(pg_exon_t, tmp->exon, tmp->n_exon, tmp->m_exon);
	p = &tmp->exon[tmp->n_exon++];
	p->ost = p->oen = st;
}

static void pg_parse_cigar(pg_data_t *d, pg_genome_t *g, pg_hit_t *hit, pg_exons_t *tmp, const char *cg)
{
	const char *p = cg;
	char *r;
	int64_t x = 0;
	int32_t i, n_fs = 0;
	pg_exon_t *t;
	tmp->n_exon = 0;
	pg_add_exon(tmp, 0);
	while (*p) {
		int64_t l;
		l = strtol(p, &r, 10);
		if (*r == 'N' || *r == 'U' || *r == 'V') {
			int64_t st, en;
			if (*r == 'N') st = x, en = x + l;
			else if (*r == 'U') st = x + 1, en = x + l - 2;
			else st = x + 2, en = x + l - 1;
			tmp->exon[tmp->n_exon - 1].oen = st;
			pg_add_exon(tmp, en);
			x += l;
		} else if (*r == 'M' || *r == 'X' || *r == '=' || *r == 'D') {
			x += l * 3;
		} else if (*r == 'F' || *r == 'G') {
			x += l, ++n_fs;
		}
		p = r + 1;
	}
	tmp->exon[tmp->n_exon - 1].oen = x;
	assert(x == hit->ce - hit->cs);
	PG_EXTEND(pg_exon_t, g->exon, g->n_exon + tmp->n_exon - 1, g->m_exon);
	t = &g->exon[g->n_exon];
	if (!hit->rev) {
		memcpy(t, tmp->exon, tmp->n_exon * sizeof(pg_exon_t));
	} else {
		for (i = tmp->n_exon - 1; i >= 0; --i, ++t) {
			t->ost = x - tmp->exon[i].oen;
			t->oen = x - tmp->exon[i].ost;
		}
	}
	hit->n_exon = tmp->n_exon;
	hit->off_exon = g->n_exon;
	hit->fs = n_fs;
	g->n_exon += tmp->n_exon;
}

int32_t pg_read_paf(pg_data_t *d, const char *fn, int32_t sep)
{
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t dret, absent;
	void *d_ctg, *hit_rank;
	pg_genome_t *g;
	pg_exons_t buf = {0,0,0};

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return -1;

	hit_rank = pg_dict_init();
	d_ctg = pg_dict_init();
	PG_EXTEND(pg_genome_t, d->genome, d->n_genome, d->m_genome);
	g = &d->genome[d->n_genome++];
	memset(g, 0, sizeof(*g));

	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p, *q, *r;
		int32_t i;
		pg_hit_t hit;
		hit.pid = hit.cid = hit.off_exon = hit.n_exon = -1;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == '\t' || *p == 0) {
				int32_t c = *p;
				*p = 0;
				if (i == 0) { // query name
					int32_t gid, pid;
					const char *tmp;
					for (r = q; r < p && *r != sep; ++r) {}
					if (*r == sep) {
						*r = 0;
						tmp = *pg_dict_put(d->d_gene, q, pg_dict_size(d->d_gene), 1, &gid, 0);
						*r = sep;
					} else {
						tmp = *pg_dict_put(d->d_gene, q, pg_dict_size(d->d_gene), 1, &gid, 0);
					}
					tmp = *pg_dict_put(d->d_prot, q, pg_dict_size(d->d_prot), 1, &pid, &absent);
					if (absent) { // protein is new
						d->n_prot++;
						PG_EXTEND(pg_prot_t, d->prot, pid, d->m_prot);
					}
					d->prot[pid].name = tmp;
					d->prot[pid].gid = gid;
					hit.pid = pid;
					hit.rank = pg_dict_inc(hit_rank, d->prot[pid].name, 0);
				} else if (i == 1) {
					d->prot[hit.pid].len = strtol(q, &r, 10);
				} else if (i == 2) {
					hit.qs = strtol(q, &r, 10);
				} else if (i == 3) {
					hit.qe = strtol(q, &r, 10);
				} else if (i == 4) {
					if (*q != '+' && *q != '-') break;
					hit.rev = *q == '+'? 0 : 1;
				} else if (i == 5) {
					int32_t cid;
					const char **ret;
					ret = pg_dict_put(d_ctg, q, pg_dict_size(d_ctg), 0, &cid, &absent);
					if (absent) { // a new contig not seen in this PAF file
						const char *name;
						PG_EXTEND(pg_ctg_t, g->ctg, g->n_ctg, g->m_ctg);
						name = *pg_dict_put(d->d_ctg, q, pg_dict_size(d->d_ctg), 1, 0, 0);
						g->ctg[g->n_ctg++].name = *ret = name;
					}
					assert(cid < g->m_ctg);
					hit.cid = cid;
				} else if (i == 6) {
					g->ctg[hit.cid].len = strtol(q, &r, 10);
				} else if (i == 7) {
					hit.cs = strtol(q, &r, 10);
				} else if (i == 8) {
					hit.ce = strtol(q, &r, 10);
				} else if (i == 9)  {
					hit.mlen = strtol(q, &r, 10);
				} else if (i == 10) {
					hit.blen = strtol(q, &r, 10);
				} else if (i >= 12) {
					if (strncmp(q, "ms:i:", 5) == 0) {
						hit.score = strtol(q + 5, &r, 10);
					} else if (strncmp(q, "cg:Z:", 5) == 0) {
						pg_parse_cigar(d, g, &hit, &buf, q + 5);
					}
				}
				q = p + 1, ++i;
				if (c == 0) break;
			}
		}
		if (hit.n_exon >= 1) {
			PG_EXTEND(pg_hit_t, g->hit, g->n_hit, g->m_hit);
			g->hit[g->n_hit++] = hit;
		}
	}
	pg_dict_destroy(d_ctg);
	pg_dict_destroy(hit_rank);
	ks_destroy(ks);
	gzclose(fp);
	return 0;
}
