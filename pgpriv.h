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

pg_dict_t *pg_dict_init(void);
const char *pg_dict_put(pg_dict_t *d, const char *s, int32_t *id, int32_t *absent);
const char *pg_dict_get(const pg_dict_t *d, const char *s, int32_t *id);
int32_t pg_dict_size(const pg_dict_t *d);
void pg_dict_destroy(pg_dict_t *d);

void pg_sprintf_lite(kstring_t *s, const char *fmt, ...);

double pg_cputime(void);
long pg_peakrss(void);
double pg_realtime(void);

#endif
