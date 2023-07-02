#include <stdio.h>
#include "pgpriv.h"
#include "ketopt.h"

static ko_longopt_t long_options[] = {
	{ "bed",             ko_optional_argument, 301 },
	{ "version",         ko_no_argument,       401 },
	{ 0, 0, 0 }
};

static int32_t pg_usage(FILE *fp, const pg_opt_t *opt)
{
	fprintf(fp, "Usage: pangene [options] <in.paf> [...]\n");
	fprintf(fp, "Options:\n");
	fprintf(fp, "  Input preprocessing:\n");
	fprintf(fp, "    -d CHAR       gene-protein delimiter [%c]\n", opt->gene_delim);
	fprintf(fp, "    -e FLOAT      drop an alignment if its identity <FLOAT [%g]\n", opt->min_prot_iden);
	fprintf(fp, "    -l FLOAT      drop an alignment if <FLOAT fraction of the protein aligned [%g]\n", opt->min_prot_ratio);
	fprintf(fp, "  Graph construction:\n");
	fprintf(fp, "    -f FLOAT      min overlap fraction [%g]\n", opt->min_ov_ratio);
	fprintf(fp, "    -p FLOAT      gene considered if dominant in FLOAT fraction of genes [%g]\n", opt->min_vertex_ratio);
	fprintf(fp, "    -a INT        prune an arc if it is supported by <INT genomes [%d]\n", opt->min_arc_cnt);
	fprintf(fp, "    -b FLOAT      drop a branch if weaker than the best by FLOAT [%g]\n", opt->branch_diff);
	fprintf(fp, "    -c INT        drop a gene if average occurrence is >INT [%d]\n", opt->max_avg_occ);
	fprintf(fp, "  Output:\n");
	fprintf(fp, "    -w            Suppress walk lines (W-lines)\n");
	fprintf(fp, "    --bed[=STR]   output 12-column BED where STR is walk, raw or flag [walk]\n");
	fprintf(fp, "    --version     print version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c;
	pg_opt_t opt;
	pg_data_t *d;

	pg_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "d:e:l:f:p:b:c:a:wv:", long_options)) >= 0) {
		if (c == 'd') opt.gene_delim = *o.arg;
		else if (c == 'e') opt.min_prot_iden = atof(o.arg);
		else if (c == 'l') opt.min_prot_ratio = atof(o.arg);
		else if (c == 'f') opt.min_ov_ratio = atof(o.arg);
		else if (c == 'p') opt.min_vertex_ratio = atof(o.arg);
		else if (c == 'b') opt.branch_diff = atof(o.arg);
		else if (c == 'c') opt.max_avg_occ = atoi(o.arg);
		else if (c == 'a') opt.min_arc_cnt = atoi(o.arg);
		else if (c == 'w') opt.flag |= PG_F_WRITE_NO_WALK;
		else if (c == 'v') pg_verbose = atoi(o.arg);
		else if (c == 301) {
			if (o.arg == 0 || strcmp(o.arg, "walk") == 0) opt.flag |= PG_F_WRITE_BED_WALK;
			else if (strcmp(o.arg, "raw") == 0) opt.flag |= PG_F_WRITE_BED_RAW;
			else if (strcmp(o.arg, "flag") == 0) opt.flag |= PG_F_WRITE_BED_FLAG;
			else {
				fprintf(stderr, "ERROR: unrecognized --bed argument. Should be 'raw' or 'walk'.\n");
				return 1;
			}
		} else if (c == 401) {
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
	if (opt.flag & PG_F_WRITE_BED_RAW) {
		pg_write_bed(d, 0);
	} else {
		pg_graph_t *g;
		g = pg_graph_init(d);
		pg_graph_gen(&opt, g);
		if (opt.flag & PG_F_WRITE_BED_WALK) {
			pg_write_bed(d, 1);
		} else if (opt.flag & PG_F_WRITE_BED_FLAG) {
			pg_write_bed(d, 0);
		} else {
			pg_write_graph(g);
			if (!(opt.flag & PG_F_WRITE_NO_WALK))
				pg_write_walk(g);
		}
		pg_graph_destroy(g);
	}
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
