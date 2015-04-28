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

/****************
 * From seqtk.c *
 ****************/

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
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)

typedef kh_reg_t reghash_t;

reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int beg = -1, end = -1;
		reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s);
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s);
						if (end < 0) end = -1;
					}
				}
			}
		}
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg < 0) beg = 0, end = INT_MAX;
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return h;
}

void stk_reg_destroy(reghash_t *h)
{
	khint_t k;
	if (h == 0) return;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

void stk_mask(kseq_t *seq, const khash_t(reg) *h, int is_complement, int mask_chr)
{
	unsigned i, j;
	khiter_t k;
	k = kh_get(reg, h, seq->name.s);
	if (k == kh_end(h)) { // not found in the hash table
		if (is_complement) {
			if (mask_chr) {
				for (j = 0; j < seq->seq.l; ++j)
					seq->seq.s[j] = mask_chr;
			} else {
				for (j = 0; j < seq->seq.l; ++j)
					seq->seq.s[j] = tolower(seq->seq.s[j]);
			}
		}
	} else {
		reglist_t *p = &kh_val(h, k);
		if (!is_complement) {
			for (i = 0; i < p->n; ++i) {
				unsigned beg = p->a[i]>>32, end = p->a[i];
				if (beg >= seq->seq.l) continue;
				if (end > seq->seq.l) end = seq->seq.l;
				if (!mask_chr) for (j = beg; j < end; ++j) seq->seq.s[j] = tolower(seq->seq.s[j]);
				else for (j = beg; j < end; ++j) seq->seq.s[j] = mask_chr;
			}
		} else {
			int8_t *mask = calloc(seq->seq.l, 1);
			for (i = 0; i < p->n; ++i) {
				unsigned beg = p->a[i]>>32, end = p->a[i];
				if (end >= seq->seq.l) end = seq->seq.l;
				for (j = beg; j < end; ++j) mask[j] = 1;
			}
			if (mask_chr) {
				for (j = 0; j < seq->seq.l; ++j)
					if (mask[j] == 0) seq->seq.s[j] = mask_chr;
			} else {
				for (j = 0; j < seq->seq.l; ++j)
					if (mask[j] == 0) seq->seq.s[j] = tolower(seq->seq.s[j]);
			}
			free(mask);
		}
	}
}

/*******************
 * Block JackKnife *
 *******************/

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

/*****************
 * Main function *
 *****************/

typedef struct {
	int n_called, n_ABBA, n_BABA, tid, pos;
	int n_hom, n_het;
} stat_t;

int main(int argc, char *argv[])
{
	stat_t *a = 0;
	int i, c, win_size = 1000000, max_a = 0, shift = 0, n_seqs = 0, is_3 = 0, q_thres = 0, n_files, step = 100, is_binary = 0, is_binning = 0;
	int32_t tmp;
	reghash_t *h = 0;
	gzFile fp[4];
	kseq_t *ks[4];

	while ((c = getopt(argc, argv, "w:3q:M:s:bB")) >= 0) {
		if (c == 'w') win_size = atoi(optarg);
		else if (c == '3') is_3 = 1;
		else if (c == 'q') q_thres = atoi(optarg);
		else if (c == 'M') h = stk_reg_read(optarg);
		else if (c == 'b') is_binary = 1;
		else if (c == 'B') is_binning = 1;
		else if (c == 's') step = atoi(optarg);
	}
	n_files = argc - optind;
	if (n_files < 4) {
		fprintf(stderr, "Usage: calD [options] <U.fa> <V.fa> <X.fa> <Y.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -w INT     window size for block Jack-Knife [%d]\n", win_size);
		fprintf(stderr, "  -q INT     quality threshold [%d]\n", q_thres);
		fprintf(stderr, "  -M FILE    mask regions in BED FILE [null]\n");
		fprintf(stderr, "  -s INT     step size (effective with -b) [%d]\n", step);
		fprintf(stderr, "  -b         binary output\n");
		return 1;
	}

	srand48(11);
	for (i = 0; i < n_files; ++i) {
		fp[i] = gzopen(argv[optind + i], "rb");
		if (fp[i] == 0) {
			fprintf(stderr, "(EE) fail to open file\n");
			return 2;
		}
		ks[i] = kseq_init(fp[i]);
	}
	if (is_binary) {
		tmp = 2;
		fwrite(&tmp, 4, 1, stdout);
	}
	while (kseq_read(ks[0]) >= 0) {
		int l, min_l = 1<<30, z[2];
		for (i = 1; i < n_files; ++i)
			if (kseq_read(ks[i]) < 0) break;
		if (i != n_files) break; // finish
		// apply quality and/or regional mask
		for (i = 0; i < n_files; ++i) {
			kseq_t *ksi = ks[i];
			if (ks[i]->qual.l && q_thres > 0) { // apply quality threshold
				for (l = 0; l < ksi->seq.l; ++l)
					if (ksi->qual.s[l] - 33 < q_thres)
						ksi->seq.s[l] = tolower(ksi->seq.s[l]);
			}
		}
		if (h) stk_mask(ks[0], h, 0, 'X'); // region masking
		// find the smallest chromosome length
		for (i = 0; i < n_files; ++i)
			if (min_l > ks[i]->seq.l)
				min_l = ks[i]->seq.l;
		if (is_binary) {
			tmp = (min_l + step - 1) / step;
			fwrite(&tmp, 4, 1, stdout); // write length
		}
		// core loop
		z[0] = z[1] = 0;
		for (l = 0; l < min_l; ++l) {
			int c[4], b[4], win = shift + l / win_size;
			if (l && l % step == 0) {
				if (is_binning) z[0] = !!z[0], z[1] = !!z[1];
				if (is_binary) fwrite(z, 4, 2, stdout);
				z[0] = z[1] = 0;
			}
			if (win >= max_a) {
				int old_max = max_a;
				max_a = win + 1;
				kroundup32(max_a);
				a = realloc(a, max_a * sizeof(stat_t));
				memset(&a[old_max], 0, (max_a - old_max) * sizeof(stat_t));
			}
			a[win].tid = n_seqs; a[win].pos = l / win_size;
			for (i = 0; i < n_files; ++i) {
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
			if (b[0] == b[1] && b[1] != b[2] && b[2] == b[3]) ++a[win].n_hom;
			else if (b[0] != b[1] && b[2] != b[3]) ++a[win].n_het;
			if (b[2] == b[3] || b[0] == b[1]) continue; // not ABBA or BABA
			if (b[0] == b[3]) ++a[win].n_ABBA, ++z[0];
			else ++a[win].n_BABA, ++z[1];
		}
		if (is_binning) z[0] = !!z[0], z[1] = !!z[1];
		if (is_binary) fwrite(z, 4, 2, stdout);
		shift += min_l / win_size + 1;
		++n_seqs;
	}
	for (i = 0; i < n_files; ++i) {
		kseq_destroy(ks[i]);
		gzclose(fp[i]);
	}
	if (h) stk_reg_destroy(h);

	if (!is_3) { // block jack-knife for D
		long s_ABBA = 0, s_BABA = 0;
		int k, *m;
		double theta0, *theta, e, var;
		for (i = 0; i < shift; ++i) s_ABBA += a[i].n_ABBA, s_BABA += a[i].n_BABA;
		theta0 = (double)(s_BABA - s_ABBA) / (s_ABBA + s_BABA);
		m = malloc(shift * sizeof(int));
		theta = malloc(shift * sizeof(double));
		for (i = k = 0; i < shift; ++i) {
			if (a[i].n_ABBA + a[i].n_BABA == 0) continue;
			m[k] = a[i].n_ABBA + a[i].n_BABA;
			theta[k++] = (double)((s_BABA - a[i].n_BABA) - (s_ABBA - a[i].n_ABBA)) / ((s_ABBA - a[i].n_ABBA) + (s_BABA - a[i].n_BABA));
		}
		e = n_unequal_jk(theta0, k, m, theta, &var);
		free(m); free(theta);
		printf("D4\t%f\t%f\t%f\n", e, sqrt(var), e / sqrt(var));
	} else { // F3 statistics; not quite finished yet
		long s_hom = 0, s_het = 0;
		int k, *m;
		double theta0, *theta, e, var;
		for (i = 0; i < shift; ++i) s_hom += a[i].n_hom, s_het += a[i].n_het;
		theta0 = (double)(2 * s_hom - s_het) / (2 * s_hom + s_het);
		m = malloc(shift * sizeof(int));
		theta = malloc(shift * sizeof(double));
		for (i = k = 0; i < shift; ++i) {
			if (a[i].n_hom + a[i].n_het == 0) continue;
			m[k] = 2 * a[i].n_hom + a[i].n_het;
			theta[k++] = (double)(2 * (s_hom - a[i].n_hom) - (s_het - a[i].n_het)) / (2 * (s_hom - a[i].n_hom) + (s_het - a[i].n_het));
		}
		e = n_unequal_jk(theta0, k, m, theta, &var);
		free(m); free(theta);
		printf("D3\t%f\t%f\t%f\n", e, sqrt(var), e / sqrt(var));
	}

	free(a);
	return 0;
}
