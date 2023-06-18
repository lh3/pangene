#include <stdio.h>
#include "pgpriv.h"
#include "ketopt.h"

static ko_longopt_t long_options[] = {
	{ "bed",             ko_no_argument,       301 },
	{ 0, 0, 0 }
};

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t i, c, gene_sep = ':', bed_out = 0;
	pg_opt_t opt;
	pg_data_t *d;

	pg_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "d:", long_options)) >= 0) {
		if (c == 301) bed_out = 1;
	}
	if (argc - o.ind < 1) {
		fprintf(stderr, "Usage: pangene [options] <in.paf> [...]\n");
		return 1;
	}

	pg_realtime();
	d = pg_data_init();
	for (i = o.ind; i < argc; ++i) {
		if (pg_verbose >= 3)
			fprintf(stderr, "[M::%s] reading file '%s'\n", __func__, argv[i]);
		pg_read_paf(&opt, d, argv[i], gene_sep);
	}
	if (bed_out) {
		for (i = 0; i < d->n_genome; ++i)
			pg_write_bed(d, i);
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
