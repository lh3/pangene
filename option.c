#include <string.h>
#include "pangene.h"

int pg_verbose = 3;

void pg_opt_init(pg_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->flag = 0;
	opt->gene_delim = ':';
	opt->min_prot_iden = 0.5;
	opt->min_prot_ratio = 0.5;
	opt->min_ov_ratio = 0.5;
	opt->min_vertex_ratio = 0.05;
	opt->max_avg_occ = 10;
	opt->max_degree = 10;
	opt->n_branch_flt = 3;
	opt->min_arc_cnt = 1;
	opt->local_dist = 2000000;
	opt->local_count = 20;
	opt->branch_diff = 0.03;
}
