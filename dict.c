#include "pgpriv.h"
#include "khashl.h"
KHASHL_CMAP_INIT(KH_LOCAL, pg_strhash_t, pg_sh, const char*, int32_t, kh_hash_str, kh_eq_str)

struct pg_dict_s {
	int32_t n_str, m_str;
	char **str;
	pg_strhash_t *h;
};

pg_dict_t *pg_dict_init(void)
{
	pg_dict_t *d;
	d = PG_CALLOC(pg_dict_t, 1);
	d->h = pg_sh_init();
	return d;
}

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

const char *pg_dict_put(pg_dict_t *d, const char *s, int32_t *id, int32_t *absent_)
{
	int absent;
	khint_t k;
	k = pg_sh_put(d->h, s, &absent);
	if (absent) {
		if (d->n_str == d->m_str) {
			d->m_str += (d->m_str>>1) + 16;
			d->str = PG_REALLOC(char*, d->str, d->m_str);
		}
		kh_key(d->h, k) = d->str[d->n_str] = pg_strdup(s);
		kh_val(d->h, k) = d->n_str++;
	}
	if (absent_) *absent_ = absent;
	if (id) *id = kh_val(d->h, k);
	return kh_key(d->h, k);
}

const char *pg_dict_get(const pg_dict_t *d, const char *s, int32_t *id)
{
	khint_t k;
	k = pg_sh_get(d->h, s);
	if (id) *id = k == kh_end(d->h)? -1 : kh_val(d->h, k);
	return k == kh_end(d->h)? 0 : kh_key(d->h, k);
}

int32_t pg_dict_size(const pg_dict_t *d)
{
	return d->n_str;
}

void pg_dict_destroy(pg_dict_t *d)
{
	int32_t i;
	if (d == 0) return;
	pg_sh_destroy(d->h);
	for (i = 0; i < d->n_str; ++i)
		free(d->str[i]);
	free(d->str);
	free(d);
}
