#ifndef PANGENE_H
#define PANGENE_H

#include <stdint.h>

struct pg_dict_s;
typedef struct pg_dict_s pg_dict_t;

typedef struct {
	int32_t ost, oen;
} pg_exon_t;

typedef struct {
	const char *name;
	int32_t len;
	int32_t gid;
} pg_prot_t;

typedef struct {
	uint32_t pid:31, rev:1; // protein ID
	int32_t qs, qe;
	int32_t cid; // contig ID
	int32_t mlen, blen;
	int32_t score;
	int32_t n_exon;
	int64_t cs, ce;
	int64_t off_exon;
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
	pg_dict_t *d_ctg, *d_gene, *d_prot;
	int32_t n_genome, m_genome;
	pg_genome_t *genome;
	int32_t n_prot, m_prot;
	pg_prot_t *prot;
} pg_data_t;

#endif
