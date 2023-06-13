#include "pgpriv.h"
#include "khashl.h"
KHASHL_CMAP_INIT(KH_LOCAL, pg_strhash_t, pg_sh, const char*, int32_t, kh_hash_str, kh_eq_str)

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

void *pg_dict_init(void)
{
	return pg_sh_init();
}

void pg_dict_destroy(void *h_)
{
	pg_sh_destroy((pg_strhash_t*)h_);
}

void pg_dict_destroy_copy(void *h_)
{
	pg_strhash_t *h = (pg_strhash_t*)h_;
	khint_t k;
	for (k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k))
			free((char*)kh_key(h, k));
	pg_sh_destroy(h);
}

int32_t pg_dict_size(const void *h_)
{
	return kh_size((pg_strhash_t*)h_);
}

const char **pg_dict_put(void *h_, const char *s, int32_t v0, int32_t do_copy, int32_t *v1, int32_t *absent_)
{
	pg_strhash_t *h = (pg_strhash_t*)h_;
	int absent;
	khint_t k;
	k = pg_sh_put(h, s, &absent);
	if (absent) {
		kh_key(h, k) = do_copy? pg_strdup(s) : s;
		kh_val(h, k) = v0;
	}
	if (absent_) *absent_ = absent;
	if (v1) *v1 = kh_val(h, k);
	return &kh_key(h, k);
}

int32_t pg_dict_inc(void *h_, const char *s, int32_t v0)
{
	pg_strhash_t *h = (pg_strhash_t*)h_;
	int absent;
	khint_t k;
	k = pg_sh_put(h, s, &absent);
	if (absent) kh_key(h, k) = s, kh_val(h, k) = v0;
	else ++kh_val(h, k);
	return kh_val(h, k);
}
