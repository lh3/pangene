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
	fprintf(fp, "    -X STR/@FILE  exclude genes in STR list or in @FILE []\n");
	fprintf(fp, "    -I STR/@FILE  attempt to include genes in the output graph []\n");
	fprintf(fp, "    -e FLOAT      drop an alignment if its identity <FLOAT [%g]\n", opt->min_prot_iden);
	fprintf(fp, "    -l FLOAT      drop an alignment if <FLOAT fraction of the protein aligned [%g]\n", opt->min_prot_ratio);
	fprintf(fp, "    -m FLOAT      score adjustment coefficient [%g]\n", opt->score_adj_coef);
	fprintf(fp, "  Graph construction:\n");
	fprintf(fp, "    -f FLOAT      min overlap fraction [%g]\n", opt->min_ov_ratio);
	fprintf(fp, "    -J            skip joint pseudogene filtering\n");
	fprintf(fp, "    -p FLOAT      gene considered if dominant in FLOAT fraction of genes [%g]\n", opt->min_vertex_ratio);
	fprintf(fp, "    -c INT        drop a gene if average occurrence is >INT [%d]\n", opt->max_avg_occ);
	fprintf(fp, "    -g INT        drop a gene if its in- or out-degree >INT [%d]\n", opt->max_degree);
	fprintf(fp, "    -r INT        drop a gene if it connects >INT distant loci [%d]\n", opt->max_dist_loci);
	fprintf(fp, "    -b FLOAT      demote a branching arc if weaker than the best by FLOAT [%g]\n", opt->branch_diff);
	fprintf(fp, "    -B FLOAT      cut a branching arc if weaker by FLOAT [%g]\n", opt->branch_diff_cut);
	fprintf(fp, "    -y FLOAT      cut a distant branching arc if weaker by FLOAT [%g]\n", opt->branch_diff_dist);
	fprintf(fp, "    -T INT        apply branch cutting for INT times [%d]\n", opt->n_branch_flt);
	fprintf(fp, "    -F            don't consider genes on different contigs as distant\n");
	fprintf(fp, "    -a INT        prune an arc if it is supported by <INT genomes [%d]\n", opt->min_arc_cnt);
	fprintf(fp, "  Output:\n");
	fprintf(fp, "    -w            Suppress walk lines (W-lines)\n");
	fprintf(fp, "    --bed[=STR]   output 12-column BED where STR is walk, raw or flag [walk]\n");
	fprintf(fp, "    --version     print version number\n");
	return fp == stdout? 0 : 1;
}

static inline int64_t pg_parse_num2(const char *str, char **q)
{
	double x;
	char *p;
	x = strtod(str, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9, ++p;
	else if (*p == 'M' || *p == 'm') x *= 1e6, ++p;
	else if (*p == 'K' || *p == 'k') x *= 1e3, ++p;
	if (q) *q = p;
	return (int64_t)(x + .499);
}

static inline int64_t pg_parse_num(const char *str)
{
	return pg_parse_num2(str, 0);
}

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c;
	pg_opt_t opt;
	pg_data_t *d;

	pg_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "d:e:l:f:g:p:b:B:y:Fr:c:a:wv:GD:C:T:X:I:P:m:JO", long_options)) >= 0) {
		if (c == 'd') opt.gene_delim = *o.arg;
		else if (c == 'e') opt.min_prot_iden = atof(o.arg);
		else if (c == 'l') opt.min_prot_ratio = atof(o.arg);
		else if (c == 'f') opt.min_ov_ratio = atof(o.arg);
		else if (c == 'm') opt.score_adj_coef = atof(o.arg);
		else if (c == 'p') opt.min_vertex_ratio = atof(o.arg);
		else if (c == 'J') opt.flag |= PG_F_NO_JOINT_PSEUDO;
		else if (c == 'b') opt.branch_diff = atof(o.arg);
		else if (c == 'B') opt.branch_diff_cut = atof(o.arg);
		else if (c == 'O') opt.flag |= PG_F_ORI_FOR_BRANCH;
		else if (c == 'y') opt.branch_diff_dist = atof(o.arg);
		else if (c == 'r') opt.max_dist_loci = atoi(o.arg);
		else if (c == 'F') opt.flag |= PG_F_FRAG_MODE;
		else if (c == 'T') opt.n_branch_flt = atof(o.arg);
		else if (c == 'c') opt.max_avg_occ = atoi(o.arg);
		else if (c == 'g') opt.max_degree = atoi(o.arg);
		else if (c == 'a') opt.min_arc_cnt = atoi(o.arg);
		else if (c == 'w') opt.flag |= PG_F_WRITE_NO_WALK;
		else if (c == 'G') opt.flag |= PG_F_WRITE_VTX_SEL;
		else if (c == 'D') opt.local_dist = pg_parse_num(o.arg);
		else if (c == 'C') opt.local_count = atoi(o.arg);
		else if (c == 'X') opt.excl = pg_read_list_dict(o.arg);
		else if (c == 'I') opt.incl = pg_read_list_dict(o.arg);
		else if (c == 'P') opt.preferred = pg_read_list_dict(o.arg);
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
	if (opt.excl) pg_dict_destroy(opt.excl);
	if (opt.incl) pg_dict_destroy(opt.incl);
	if (opt.preferred) pg_dict_destroy(opt.preferred);

	if (pg_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, PG_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB\n", __func__, pg_realtime(), pg_cputime(), pg_peakrss() / 1024.0 / 1024.0 / 1024.0);
	}
	return 0;
}
