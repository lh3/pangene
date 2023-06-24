#include <stdio.h>
#include "pgpriv.h"
#include "ketopt.h"

static ko_longopt_t long_options[] = {
	{ "bed",             ko_no_argument,       301 },
	{ "version",         ko_no_argument,       401 },
	{ 0, 0, 0 }
};

static int32_t pg_usage(FILE *fp, const pg_opt_t *opt)
{
	fprintf(fp, "Usage: pangene [options] <in.paf> [...]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  -d CHAR       gene name delimiter [%c]\n", opt->gene_delim);
	fprintf(fp, "  -e FLOAT      min protein identity [%g]\n", opt->min_prot_iden);
	fprintf(fp, "  -l FLOAT      min protein alignment fraction [%g]\n", opt->min_prot_ratio);
	fprintf(fp, "  -f FLOAT      min overlap fraction [%g]\n", opt->min_ov_ratio);
	fprintf(fp, "  -p FLOAT      min primary ratio to select a gene [%g]\n", opt->min_vertex_ratio);
	fprintf(fp, "  -c            max number of average occurrence [%d]\n", opt->max_avg_occ);
	fprintf(fp, "  --bed         output BED12 (mainly for debugging)\n");
	fprintf(fp, "  --version     print version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c, bed_out = 0;
	pg_opt_t opt;
	pg_data_t *d;
	pg_graph_t *g;

	pg_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "d:e:l:f:p:v:", long_options)) >= 0) {
		if (c == 'd') opt.gene_delim = *o.arg;
		else if (c == 'e') opt.min_prot_iden = atof(o.arg);
		else if (c == 'l') opt.min_prot_ratio = atof(o.arg);
		else if (c == 'p') opt.min_vertex_ratio = atof(o.arg);
		else if (c == 'f') opt.min_ov_ratio = atof(o.arg);
		else if (c == 'v') pg_verbose = atoi(o.arg);
		else if (c == 301) bed_out = 1;
		else if (c == 401) {
			puts(PG_VERSION);
			return 0;
		}
	}
	if (argc - o.ind < 1)
		return pg_usage(stderr, &opt);

	pg_realtime();
	d = pg_data_init();
	for (i = o.ind; i < argc; ++i)
		pg_read_paf(&opt, d, argv[i]);
	pg_post_process(&opt, d);
	if (bed_out) {
		for (i = 0; i < d->n_genome; ++i)
			pg_write_bed(d, i);
		return 0; // TODO: deallocate
	}
	g = pg_graph_init(d);
	pg_graph_gen(&opt, g);
	pg_graph_write(g);
	pg_graph_destroy(g);
	pg_data_destroy(d);

	if (pg_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PG_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, pg_realtime(), pg_cputime(), pg_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
