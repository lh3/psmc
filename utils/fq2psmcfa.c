#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static double N_RATIO = 0.9;
static int BLOCK_LEN = 100;
static int par1_b = 1, par1_e = 2709520;
static int par2_b = 154584237, par2_e = 154913754;

unsigned char aln_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,16 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
	15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
	15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
	15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

static int count_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int c, len, n_min_good = 10000, min_qual = 10, mask_pseudo = 0;
	while ((c = getopt(argc, argv, "q:xg:")) >= 0) {
		switch (c) {
		case 'q': min_qual = atoi(optarg); break;
		case 'x': mask_pseudo = 1; break;
		case 'g': n_min_good = atoi(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: fq2psmcfa [-x] [-q %d] [-g %d] <in.fq>\n", min_qual, n_min_good);
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
		} else fprintf(stderr, "[fq2psmcfa] the input is not FASTQ. Quality filter disabled.\n");
		//
		ss = (char*)calloc((len + BLOCK_LEN - 1) / BLOCK_LEN + 2, 1);
        for (i = 0; i != len; ++i) {
			int is_N, is_het = 0, c = islower(seq->seq.s[i])? 15 : aln_nt16_table[(int)seq->seq.s[i]];
			if (i && i%BLOCK_LEN == 0) {
				ss[l++] = (float)nN/BLOCK_LEN > N_RATIO? 'N' : (is_hetb? 'K' : 'T');
				is_hetb = 0; nN = 0;
			}
			is_N = (c == 15)? 1 : 0;
			if (is_N) ++nN;
			else if (count_table[c] == 2) is_het = 1;
			if (is_het) is_hetb = 1;
			if (!is_N) ++n_good_bases;
        }
		ss[l++] = (float)nN/BLOCK_LEN > N_RATIO? 'N' : (is_hetb? 'K' : 'T');
		fprintf(stderr, "[fq2psmcfa] %s: %d, %d\n", seq->name.s, n_good_bases, len);
		if ((double)n_good_bases / len >= 0.33 && n_good_bases >= n_min_good) {
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
