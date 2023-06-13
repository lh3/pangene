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

#define PG_EXTEND(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = PG_REALLOC(type, ptr, (__m)); \
		} \
	} while (0)

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

void *pg_dict_init(int32_t do_copy);
void pg_dict_destroy(void *h_);
int32_t pg_dict_size(const void *h_);
const char **pg_dict_put(void *d_, const char *s, int32_t v0, int32_t *v1, int32_t *absent_);
int32_t pg_dict_inc(void *h_, const char *s, int32_t v0);

void pg_sprintf_lite(kstring_t *s, const char *fmt, ...);

double pg_cputime(void);
long pg_peakrss(void);
double pg_realtime(void);
double pg_percent_cpu(void);

#endif
