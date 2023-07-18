#ifndef PANGENE_H
#define PANGENE_H

#include <stdint.h>

#define PG_VERSION "0.0-r100-dirty"

#define PG_F_WRITE_BED_RAW      0x1
#define PG_F_WRITE_BED_WALK     0x2
#define PG_F_WRITE_BED_FLAG     0x4
#define PG_F_WRITE_NO_WALK      0x8
#define PG_F_NO_MERGE_ORTHO     0x10
#define PG_F_WRITE_VTX_SEL      0x20

typedef struct {
	uint64_t x, y;
} pg128_t;

typedef struct {
	uint32_t flag;
	int32_t gene_delim;
	double min_prot_ratio; // filter out a protein if less than 50% of the protein is aligned
	double min_prot_iden; // filter out a protein if identity below 20%
	double min_ov_ratio; // consider two proteins of different genes overlap if 50% of the short protein overlap
	double min_vertex_ratio; // a gene is considered as a vertex if it is primary in 33% of the assemblies
	double branch_diff;
	int32_t max_avg_occ;
	int32_t max_degree;
	int32_t n_branch_flt;
	int32_t min_arc_cnt;
	int32_t close_thres;
} pg_opt_t;

typedef struct {
	int32_t os, oe; // start and end relative to pg_hit_t::cs
} pg_exon_t;

typedef struct {
	const char *name;
	int32_t len;
	uint32_t gid:31, rep:1; // rep: if this protein is the representative isoform of the gene
	int32_t n, avg_score2;
} pg_prot_t;

typedef struct {
	const char *name;
	int32_t len; // longest ORF
	int32_t rep_pid; // pid of the representative isoform
} pg_gene_t;

typedef struct {
	int32_t pid; // protein ID
	int32_t qs, qe;
	int32_t cid; // contig ID
	int32_t mlen, blen, fs;
	int32_t rank;
	int32_t score, score2;
	int32_t n_exon, off_exon;
	int32_t pid_dom;
	uint32_t rev:1, flt:1, pseudo:1, vtx:1, overlap:1, shadow:1, rep:1, weak_br:2, dummy:23;
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
	int32_t gid, n_dom, n_sub;
	int32_t n_genome;
	int32_t tot_cnt;
	uint32_t del:1, dummy:31;
} pg_seg_t;

typedef struct {
	uint64_t x; // v<<32|w
	int32_t n_genome;
	int32_t tot_cnt;
	int32_t avg_dist;
	int32_t s1, s2;
	uint32_t del:1, weak_br:2, dummy:29;
} pg_arc_t;

typedef struct {
	pg_data_t *d;
	int32_t *g2s;
	int32_t n_seg, m_seg;
	pg_seg_t *seg;
	int32_t n_arc, m_arc;
	pg_arc_t *arc;
	uint64_t *idx;
} pg_graph_t;

extern int pg_verbose;

void pg_opt_init(pg_opt_t *opt);

pg_data_t *pg_data_init(void);
void pg_data_destroy(pg_data_t *d);
int32_t pg_read_paf(const pg_opt_t *opt, pg_data_t *d, const char *fn);
void pg_post_process(const pg_opt_t *opt, pg_data_t *d);

pg_graph_t *pg_graph_init(pg_data_t *d);
void pg_graph_gen(const pg_opt_t *opt, pg_graph_t *q);
void pg_graph_destroy(pg_graph_t *g);

void pg_write_bed(const pg_data_t *d, int32_t is_walk);
void pg_write_graph(const pg_graph_t *g);
void pg_write_walk(pg_graph_t *g);

#endif
