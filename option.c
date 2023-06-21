#include <string.h>
#include "pangene.h"

int pg_verbose = 3;

void pg_opt_init(pg_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->min_prot_ratio = 0.5;
	opt->min_ov_ratio = 0.5;
	opt->min_vertex_ratio = 0.33;
}
