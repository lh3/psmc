#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <zlib.h>
#include <stdint.h>
#include <unistd.h>

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

char *seq_nt16_rev_table = "XACMGRSVTWYHKDBN";
uint8_t seq_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int32_t l, c, step = 100, bin = 0;
	double n_ratio = 0.9;
	while ((c = getopt(argc, argv, "s:n:b")) >= 0) {
		switch (c) {
		case 'b': bin = 1; break;
		case 's': step = atoi(optarg); break;
		case 'n': n_ratio = atof(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: cntcpg [-b] [-s 100] [-n 0.9] <in.fa>\n\n");
		fprintf(stderr, "Output: 5 numbers: #CpG #CpG-ts #nonCpG #nonCpG-ts+tv #tv\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	l = 5; // number of counts
	fwrite(&l, 4, 1, stdout);
	while ((l = kseq_read(seq)) >= 0) {
		int32_t i, z[5];
		uint8_t *s;
		s = calloc(l, 1);
		for (i = 0; i < l; ++i) { // mask lowcase bases and create a copy
			if (islower(seq->seq.s[i])) seq->seq.s[i] = 'N';
			seq->seq.s[i] = s[i] = seq_nt16_table[(int)seq->seq.s[i]];
		}
		for (i = 0; i < l; ++i) { // mask CpG
			int c0 = s[i], c1 = s[i+1];
			if ((c0 == 2 || c0 == 10) && (c1 == 4 || c1 == 5))
				s[i] = s[i+1] = 15;
		}
		i = (l + step - 1) / step;
		fwrite(&i, 4, 1, stdout); // write length
		z[0] = z[1] = z[2] = z[3] = z[4] = 0;
		for (i = 0; i < l; ++i) {
			int b = seq->seq.s[i];
			int c0 = bitcnt_table[b];
			int c1 = bitcnt_table[s[i]];
			if (i && i%step == 0) {
				if (step - z[0] > step * n_ratio) z[0] = z[1] = z[2] = z[3] = z[4] = 0;
				z[0] -= z[2]; z[1] -= z[3];
				if (bin) z[1] = (z[1] >= 1), z[3] = (z[3] >= 1), z[4] = (z[4] >= 1);
				fwrite(z, 4, 5, stdout);
				z[0] = z[1] = z[2] = z[3] = z[4] = 0;
			}
			if (c0 < 3) ++z[0];
			if (c0 == 2) ++z[1];
			if (c1 < 3) ++z[2];
			if (c1 == 2) ++z[3];
			if (c0 == 2 && b != 5 && b != 10) ++z[4]; // transversions
		}
		if (step - z[0] > step * n_ratio) z[0] = z[1] = z[2] = z[3] = z[4] = 0;
		z[0] -= z[2]; z[1] -= z[3];
		if (bin) z[1] = (z[1] >= 1), z[3] = (z[3] >= 1), z[4] = (z[4] >= 1);
		fwrite(z, 4, 5, stdout);
		free(s);
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
