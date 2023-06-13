#include "pgpriv.h"
#include "khashl.h"
KHASHL_CMAP_INIT(KH_LOCAL, pg_strhash_t, pg_sh, const char*, int32_t, kh_hash_str, kh_eq_str)

typedef struct {
	pg_strhash_t *h;
	int32_t do_copy;
} pg_dict_t;

char *pg_strdup(const char *src)
{
	int32_t len;
	char *dst;
	len = strlen(src);
	dst = PG_MALLOC(char, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

char *pg_strndup(const char *src, size_t n)
{
	char *dst;
	dst = PG_MALLOC(char, n + 1);
	strncpy(dst, src, n);
	dst[n] = 0;
	return dst;
}

void *pg_dict_init(int32_t do_copy)
{
	pg_dict_t *d;
	d = PG_CALLOC(pg_dict_t, 1);
	d->do_copy = !!do_copy;
	d->h = pg_sh_init();
	return d;
}

void pg_dict_destroy(void *h_)
{
	pg_dict_t *d = (pg_dict_t*)h_;
	if (d->do_copy) {
		khint_t k;
		for (k = 0; k < kh_end(d->h); ++k)
			if (kh_exist(d->h, k))
				free((char*)kh_key(d->h, k));
	}
	pg_sh_destroy(d->h);
	free(d);
}

int32_t pg_dict_size(const void *d)
{
	return kh_size(((pg_dict_t*)d)->h);
}

const char **pg_dict_put(void *d_, const char *s, int32_t v0, int32_t *v1, int32_t *absent_)
{
	pg_dict_t *d = (pg_dict_t*)d_;
	int absent;
	khint_t k;
	k = pg_sh_put(d->h, s, &absent);
	if (absent) {
		kh_key(d->h, k) = d->do_copy? pg_strdup(s) : s;
		kh_val(d->h, k) = v0;
	}
	if (absent_) *absent_ = absent;
	if (v1) *v1 = kh_val(d->h, k);
	return &kh_key(d->h, k);
}

int32_t pg_dict_inc(void *d_, const char *s, int32_t v0)
{
	pg_dict_t *d = (pg_dict_t*)d_;
	int absent;
	khint_t k;
	k = pg_sh_put(d->h, s, &absent);
	if (absent) {
		kh_key(d->h, k) = d->do_copy? pg_strdup(s) : s;
		kh_val(d->h, k) = v0;
	}
	else ++kh_val(d->h, k);
	return kh_val(d->h, k);
}
