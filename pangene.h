#ifndef PANGENE_H
#define PANGENE_H

#include <stdint.h>

#define PG_VERSION "0.0-r18-dirty"

typedef struct {
	double min_prot_ratio; // filter out a protein if less than 50% of proteins are aligned
	double min_ov_ratio; // consider two proteins of different genes overlap if 50% of the short protein overlap
	double min_vertex_ratio; // a gene is considered as a vertex if it is primary in 33% of the assemblies
} pg_opt_t;

typedef struct {
	int32_t os, oe;
} pg_exon_t;

typedef struct {
	const char *name;
	int32_t len;
	int32_t gid;
} pg_prot_t;

typedef struct {
	const char *name;
} pg_gene_t;

typedef struct {
	uint32_t pid:31, rev:1; // protein ID
	int32_t qs, qe;
	int32_t cid; // contig ID
	int32_t mlen, blen, fs;
	int32_t score, rank;
	int32_t n_exon, off_exon;
	int64_t cs, cm, ce;
} pg_hit_t;

typedef struct {
	const char *name;
	int64_t len;
} pg_ctg_t;

typedef struct {
	int32_t n_ctg, m_ctg;
	pg_ctg_t *ctg;
	int32_t n_hit, m_hit;
	pg_hit_t *hit;
	int32_t n_exon, m_exon;
	pg_exon_t *exon;
} pg_genome_t;

typedef struct {
	void *d_ctg, *d_gene, *d_prot;
	int32_t n_genome, m_genome;
	pg_genome_t *genome;
	int32_t n_gene, m_gene;
	pg_gene_t *gene;
	int32_t n_prot, m_prot;
	pg_prot_t *prot;
} pg_data_t;

typedef struct {
	int32_t gid, pri, sec;
} pg_vertex_t;

typedef struct {
	const pg_data_t *d;
	int32_t n_v, m_v;
	pg_vertex_t *v;
} pg_graph_t;

extern int pg_verbose;

void pg_opt_init(pg_opt_t *opt);

pg_data_t *pg_data_init(void);
void pg_data_destroy(pg_data_t *d);
int32_t pg_read_paf(const pg_opt_t *opt, pg_data_t *d, const char *fn, int32_t gene_sep);

pg_graph_t *pg_graph_init(const pg_data_t *d);
void pg_graph_destroy(pg_graph_t *g);

void pg_write_bed(const pg_data_t *d, int32_t aid);

#endif
