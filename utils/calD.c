#include <zlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };

typedef struct {
	int n_called, n_ABBA, n_BABA;
} stat_t;

double n_unequal_jk(double theta0, int g, const int *m, const double *theta, double *var)
{
	double *h, est;
	int j;
	int64_t n = 0;
	h = calloc(g, sizeof(double));
	for (j = 0; j < g; ++j) n += m[j];
	for (j = 0; j < g; ++j) h[j] = (double)n / m[j];
	for (j = 0, est = 0.0; j < g; ++j) est += (1. - (double)m[j] / n) * theta[j];
	est = g * theta0 - est;
	*var = 0;
	for (j = 0; j < g; ++j) {
		double x = (h[j] * theta0 - (h[j] - 1) * theta[j] - est);
		*var += 1. / (h[j] - 1) * x * x;
	}
	*var /= g;
	free(h);
	return est;
}

int main(int argc, char *argv[])
{
	stat_t *a = 0;
	int i, c, win_size = 1000000, max_a = 0, shift = 0;
	gzFile fp[4];
	kseq_t *ks[4];

	while ((c = getopt(argc, argv, "w:")) >= 0) {
	}
	if (optind + 4 > argc) {
		fprintf(stderr, "Usage: calD <1.fa> <2.fa> <test.fa> <outgroup.fa>\n");
		return 1;
	}

	srand48(11);
	for (i = 0; i < 4; ++i) {
		fp[i] = gzopen(argv[optind + i], "rb");
		if (fp[i] == 0) {
			fprintf(stderr, "(EE) fail to open file\n");
			return 2;
		}
		ks[i] = kseq_init(fp[i]);
	}
	while (kseq_read(ks[0]) >= 0) {
		int l, min_l = 1<<30;
		for (i = 1; i < 4; ++i)
			if (kseq_read(ks[i]) < 0) break;
		if (i != 4) break; // finish
		for (i = 0; i < 4; ++i)
			if (min_l > ks[i]->seq.l)
				min_l = ks[i]->seq.l;
		for (l = 0; l < min_l; ++l) {
			int c[4], b[4], win = shift + l / win_size;
			if (win >= max_a) {
				int old_max = max_a;
				max_a = win + 1;
				kroundup32(max_a);
				a = realloc(a, max_a * sizeof(stat_t));
				memset(&a[old_max], 0, (max_a - old_max) * sizeof(stat_t));
			}
			for (i = 0; i < 4; ++i) {
				if (islower(ks[i]->seq.s[l])) break;
				c[i] = seq_nt16_table[(uint8_t)ks[i]->seq.s[l]];
				b[i] = bitcnt_table[c[i]];
				if (b[i] == 0 || b[i] > 2) break;
				if (b[i] == 2) { // pick a random allele
					int y, z, x = drand48() < .5? 0 : 1;
					for (y = 8, z = 0; y; y >>= 1)
						if ((c[i]&y) && z++ == x) break;
					b[i] = c[i] & y;
				} else b[i] = c[i];
			}
			if (i < 4 || bitcnt_table[c[0]|c[1]|c[2]|c[3]] > 2) continue; // filtered or not biallelic
			++a[win].n_called;
			if (b[2] == b[3] || b[0] == b[1]) continue; // not ABBA or BABA
			if (b[0] == b[3]) ++a[win].n_ABBA;
			else ++a[win].n_BABA;
		}
		shift += min_l / win_size + 1;
	}
	for (i = 0; i < 4; ++i) {
		kseq_destroy(ks[i]);
		gzclose(fp[i]);
	}

	{ // block jack-knife
		long s_ABBA = 0, s_BABA = 0;
		int k, *m;
		double theta0, *theta, e, var;
		for (i = 0; i < shift; ++i) s_ABBA += a[i].n_ABBA, s_BABA += a[i].n_BABA;
		theta0 = (double)(s_ABBA - s_BABA) / (s_ABBA + s_BABA);
		m = malloc(shift * sizeof(int));
		theta = malloc(shift * sizeof(double));
		for (i = k = 0; i < shift; ++i) {
			if (a[i].n_ABBA + a[i].n_BABA == 0) continue;
			m[k] = a[i].n_ABBA + a[i].n_BABA;
			theta[k++] = (double)((s_ABBA - a[i].n_ABBA) - (s_BABA - a[i].n_BABA)) / ((s_ABBA - a[i].n_ABBA) + (s_BABA - a[i].n_BABA));
		}
		e = n_unequal_jk(theta0, k, m, theta, &var);
		free(m); free(theta);
		printf("R\t%f\t%f\t%f\n", e, sqrt(var), e / sqrt(var));
	}

	for (i = 0; i < shift; ++i) printf("B\t%d\t%d\t%d\n", a[i].n_called, a[i].n_ABBA, a[i].n_BABA);

	free(a);
	return 0;
}
