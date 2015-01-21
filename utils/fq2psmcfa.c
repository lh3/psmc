#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static double N_RATIO = 0.9;
static int BLOCK_LEN = 100;
static int par1_b = 1, par1_e = 2709520; // this is the b36 coordinate
static int par2_b = 154584237, par2_e = 154913754;

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

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int c, len, n_min_good = 10000, min_qual = 10, mask_pseudo = 0, tv_only = 0, ts_only = 0, cpg_only = 0, cpg_excl = 0;
	while ((c = getopt(argc, argv, "q:xg:s:vncC")) >= 0) {
		switch (c) {
		case 'q': min_qual = atoi(optarg); break;
		case 'x': mask_pseudo = 1; break;
		case 'v': tv_only = 1; break;
		case 'n': ts_only = 1; break;
		case 'c': cpg_only = 1; break;
		case 'C': cpg_excl = 1; break;
		case 'g': n_min_good = atoi(optarg); break;
		case 's': BLOCK_LEN = atoi(optarg); break;
		}
	}
	if (tv_only + ts_only + cpg_only + cpg_excl > 1) {
		fprintf(stderr, "[E::%s] only one of the options -c, -n, -v and -C can be applied\n", __func__);
		return 2;
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fq2psmcfa [-cnvx] [-q %d] [-g %d] [-s %d] <in.fq>\n", min_qual, n_min_good, BLOCK_LEN);
		return 1;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((len = kseq_read(seq)) >= 0) {
		int i, l = 0, nN = 0, is_hetb = 0, n_good_bases = 0;
		char *ss;
		// filtering
		if (mask_pseudo && (strcmp(seq->name.s, "X") == 0 || strcmp(seq->name.s, "chrX") == 0)) {
			for (i = par1_b - 1; i < par1_e && i < seq->seq.l; ++i) seq->seq.s[i] = tolower(seq->seq.s[i]);
			for (i = par2_b - 1; i < par2_e && i < seq->seq.l; ++i) seq->seq.s[i] = tolower(seq->seq.s[i]);
		}
		if (seq->qual.l) {
			for (i = 0; i < seq->seq.l; ++i)
				if (seq->qual.s[i] - 33 < min_qual) seq->seq.s[i] = tolower(seq->seq.s[i]);
		}
		if (tv_only) {
			for (i = 0; i < seq->seq.l; ++i) {
				int c = seq_nt16_table[(int)seq->seq.s[i]];
				if (c == 5 || c == 10) seq->seq.s[i] = tolower(seq->seq.s[i]);
			}
		} else if (ts_only) {
			int pre = -1;
			for (i = 0; i < seq->seq.l; ++i) {
				int c = seq_nt16_table[(int)seq->seq.s[i]];
				if (i > 0 && (c == 4 || c == 5) && (pre == 2 || pre == 10)) {
					seq->seq.s[i] = tolower(seq->seq.s[i]);
					seq->seq.s[i-1] = tolower(seq->seq.s[i-1]);
				} else if (c == 3 || c == 9 || c == 6 || c == 12)
					seq->seq.s[i] = tolower(seq->seq.s[i]);
				pre = c;
			}
		} else if (cpg_only) {
			int pre = -1;
			for (i = 0; i < seq->seq.l; ++i) {
				int c = seq_nt16_table[(int)seq->seq.s[i]];
				if (c == 3 || c == 9 || c == 6 || c == 12) {
					seq->seq.s[i] = tolower(seq->seq.s[i]);
				} else if (i > 0 && pre == 10 && c != 4 && c != 5) {
					seq->seq.s[i-1] = tolower(seq->seq.s[i-1]);
				} else if (i > 0 && c == 5 && pre != 2 && pre != 10) {
					seq->seq.s[i] = tolower(seq->seq.s[i]);
				}
				pre = c;
			}
		} else if (cpg_excl) {
			int pre = -1;
			for (i = 0; i < seq->seq.l; ++i) {
				int c = seq_nt16_table[(int)seq->seq.s[i]];
				if (i > 0 && (c == 4 || c == 5) && (pre == 2 || pre == 10)) {
					seq->seq.s[i] = tolower(seq->seq.s[i]);
					seq->seq.s[i-1] = tolower(seq->seq.s[i-1]);
				}
				pre = c;
			}
		}
		//
		ss = (char*)calloc((len + BLOCK_LEN - 1) / BLOCK_LEN + 2, 1);
        for (i = 0; i != len; ++i) {
			int is_N, is_het = 0, c = islower(seq->seq.s[i])? 15 : seq_nt16_table[(int)seq->seq.s[i]];
			if (i && i%BLOCK_LEN == 0) {
				ss[l++] = (float)nN/BLOCK_LEN > N_RATIO? 'N' : (is_hetb? 'K' : 'T');
				is_hetb = 0; nN = 0;
			}
			is_N = (c == 15)? 1 : 0;
			if (is_N) ++nN;
			else if (bitcnt_table[c] == 2) is_het = 1;
			if (is_het) is_hetb = 1;
			if (!is_N) ++n_good_bases;
        }
		ss[l++] = (float)nN/BLOCK_LEN > N_RATIO? 'N' : (is_hetb? 'K' : 'T');
		//fprintf(stderr, "[fq2psmcfa] %s: %d, %d\n", seq->name.s, n_good_bases, len);
		if ((double)n_good_bases / len >= 0.2 && n_good_bases >= n_min_good) {
			printf(">%s", seq->name.s);
			for (i = 0; i < l ; ++i) {
				if (i%60 == 0) putchar('\n');
				putchar(ss[i]);
			}
			putchar('\n');
			fflush(stdout);
		}
		free(ss);
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
