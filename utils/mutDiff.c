#include <zlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <assert.h>

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
	int n_called, n_ABA, n_BAA;
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
	int i, c, win_size = 1000000, max_a = 0, shift = 0, list_only = 0, binary = 0, step = 100;
	int32_t n_counts = 3;
	gzFile fp[3];
	kseq_t *ks[3];

	while ((c = getopt(argc, argv, "blw:s:")) >= 0) {
		if (c == 'l') list_only = 1;
		else if (c == 's') step = atoi(optarg);
		else if (c == 'w') step = atoi(optarg);
		else if (c == 'b') binary = list_only = 1;
	}
	if (optind + 3 > argc) {
		fprintf(stderr, "Usage: mutDiff [-bl] [-s %d] <1.fa> <2.fa> <outgroup.fa>\n", step);
		return 1;
	}
	srand48(11);
	for (i = 0; i < 3; ++i) {
		fp[i] = gzopen(argv[optind + i], "rb");
		if (fp[i] == 0) {
			fprintf(stderr, "(EE) fail to open file\n");
			return 2;
		}
		ks[i] = kseq_init(fp[i]);
	}
	if (binary) fwrite(&n_counts, 4, 1, stdout);
	while (kseq_read(ks[0]) >= 0) {
		int32_t z[3], l, min_l = 1<<30;
		for (i = 1; i < 3; ++i)
			if (kseq_read(ks[i]) < 0) break;
		if (i != 3) break; // finish
		for (i = 0; i < 3; ++i)
			if (min_l > ks[i]->seq.l)
				min_l = ks[i]->seq.l;
		z[0] = z[1] = z[2] = 0;
		l = (min_l + step - 1) / step;
		if (binary) fwrite(&l, 4, 1, stdout); // write length
		for (l = 0; l < min_l; ++l) {
			int c[3], b[3], win = shift + l / win_size;
			int is_bi_called, is_good_diff;
			if (win >= max_a) {
				int old_max = max_a;
				max_a = win + 1;
				kroundup32(max_a);
				a = realloc(a, max_a * sizeof(stat_t));
				memset(&a[old_max], 0, (max_a - old_max) * sizeof(stat_t));
			}
			for (i = 0; i < 3; ++i) {
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
			is_bi_called = i < 3 || bitcnt_table[c[0]|c[1]|c[2]] > 2? 0 : 1;
			a[win].n_called += is_bi_called;
			is_good_diff = is_bi_called && b[0] != b[1]? 1 : 0;
			if (is_good_diff) {
				assert(b[0] == b[2] || b[1] == b[2]);
				if (b[0] == b[2]) ++a[win].n_ABA;
				else ++a[win].n_BAA;
			}
			if (binary) {
				if (l && l%step == 0) {
					fwrite(z, 4, 3, stdout);
					z[0] = z[1] = z[2] = 0;
				}
				z[0] += is_bi_called;
				z[1] += (is_good_diff && b[0] == b[2]);
				z[2] += (is_good_diff && b[1] == b[2]);
			} else if (list_only && is_good_diff) {
				printf("%s\t%d\t%d\t%d\n", ks[0]->name.s, l + 1, (b[0] == b[2]), (b[1] == b[2]));
			}
		}
		if (binary) fwrite(z, 4, 3, stdout);
		shift += min_l / win_size + 1;
	}
	for (i = 0; i < 3; ++i) {
		kseq_destroy(ks[i]);
		gzclose(fp[i]);
	}

	if (!list_only) { // block jack-knife
		long s_ABA = 0, s_BAA = 0;
		int k, *m;
		double theta0, *theta, e, var;
		for (i = 0; i < shift; ++i) s_ABA += a[i].n_ABA, s_BAA += a[i].n_BAA;
		theta0 = (double)(s_ABA - s_BAA) / (s_ABA + s_BAA);
		m = malloc(shift * sizeof(int));
		theta = malloc(shift * sizeof(double));
		for (i = k = 0; i < shift; ++i) {
			if (a[i].n_ABA + a[i].n_BAA == 0) continue;
			m[k] = a[i].n_ABA + a[i].n_BAA;
			theta[k++] = (double)((s_ABA - a[i].n_ABA) - (s_BAA - a[i].n_BAA)) / ((s_ABA - a[i].n_ABA) + (s_BAA - a[i].n_BAA));
		}
		e = n_unequal_jk(theta0, k, m, theta, &var);
		free(m); free(theta);
		printf("R\t%f\t%f\t%f\n", e, sqrt(var), e / sqrt(var));
		for (i = 0; i < shift; ++i)
			printf("B\t%d\t%d\t%d\n", a[i].n_called, a[i].n_ABA, a[i].n_BAA);
	}

	free(a);
	return 0;
}
