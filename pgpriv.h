#ifndef PGPRIV_H
#define PGPRIV_H

#include <stddef.h>
#include <stdlib.h>
#include "pangene.h"

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#define PG_MALLOC(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define PG_CALLOC(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define PG_REALLOC(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define PG_GROW(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = PG_REALLOC(type, ptr, (__m)); \
		} \
	} while (0)

#define PG_GROW0(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			size_t old_m = (__m); \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = PG_REALLOC(type, ptr, (__m)); \
			memset((ptr) + old_m, 0, ((__m) - old_m) * sizeof(type)); \
		} \
	} while (0)

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

void radix_sort_pg128x(pg128_t *st, pg128_t *en);
void radix_sort_pg64(uint64_t *st, uint64_t *en);

void *pg_dict_init(int32_t do_copy);
void pg_dict_destroy(void *h_);
int32_t pg_dict_size(const void *h_);
const char **pg_dict_put(void *d_, const char *s, int32_t v0, int32_t *v1, int32_t *absent_);
int32_t pg_dict_get(const void *d_, const char *s);
int32_t pg_dict_inc(void *h_, const char *s, int32_t v0);

void pg_sprintf_lite(kstring_t *s, const char *fmt, ...);

double pg_cputime(void);
long pg_peakrss(void);
double pg_realtime(void);
double pg_percent_cpu(void);
const char *pg_timestamp(void);

int64_t pg_hit_cal_cm(const pg_hit_t *a, const pg_exon_t *e);
void pg_hit_sort(pg_genome_t *g, int32_t by_cm);
uint64_t pg_hit_overlap(const pg_genome_t *g, const pg_hit_t *aa, const pg_hit_t *ab);
int32_t pg_flag_pseudo(const pg_prot_t *prot, pg_genome_t *g);
int32_t pg_flag_shadow(const pg_opt_t *opt, const pg_prot_t *prot, pg_genome_t *g, int32_t check_vtx, int32_t check_pri);
void pg_flag_representative(pg_data_t *d);

void pg_gen_g2s(pg_graph_t *q);
void pg_graph_flag_vtx(pg_graph_t *q);
void pg_gen_vtx(const pg_opt_t *opt, pg_graph_t *q);

int32_t pg_mark_branch_flt_arc(const pg_opt_t *opt, pg_graph_t *q);
int32_t pg_mark_branch_flt_hit(const pg_opt_t *opt, pg_graph_t *q);

static inline uint32_t pg_hash_uint32(uint32_t key)
{
	key += ~(key << 15);
	key ^=  (key >> 10);
	key +=  (key << 3);
	key ^=  (key >> 6);
	key += ~(key << 11);
	key ^=  (key >> 16);
	return key;
}

static inline int32_t pg_hit_arc(const pg_hit_t *a)
{
	return (a->rep && a->vtx && !a->shadow && !a->pseudo);
}

static inline const pg_arc_t *pg_get_arc(const pg_graph_t *q, uint32_t v, uint32_t w)
{
	int32_t i, n = (int32_t)q->idx[v];
	const pg_arc_t *a = &q->arc[q->idx[v]>>32];
	for (i = 0; i < n; ++i)
		if ((uint32_t)a[i].x == w)
			return &a[i];
	return 0;
}

#endif
