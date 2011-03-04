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

void cnt_het2(const char *seq, const uint8_t *s, int l, int x[4])
{
	int i, cnt[4];
	cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0;
	for (i = 0; i < l; ++i) {
		int c0 = bitcnt_table[seq_nt16_table[(int)seq[i]]];
		int c1 = bitcnt_table[s[i]];
		if (c0 < 3) ++cnt[0];
		if (c0 == 2) ++cnt[1];
		if (c1 < 3) ++cnt[2];
		if (c1 == 2) ++cnt[3];
	}
	x[0] = cnt[0] - cnt[2];
	x[1] = cnt[1] - cnt[3];
	x[2] = cnt[0] - x[0];
	x[3] = cnt[1] - x[1];
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int32_t l, c, step = 100;
	while ((c = getopt(argc, argv, "s")) >= 0) {
		switch (c) {
		case 's': step = atoi(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: cnt_bed_mut [-s 100] <in.fa> [<in.bed>]\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		int32_t i, z[4];
		uint8_t *s;
		s = calloc(l, 1);
		for (i = 0; i < l; ++i) { // create a copy
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
		z[0] = z[1] = z[2] = z[3] = 0;
		for (i = 0; i < l; ++i) {
			int c0 = bitcnt_table[(int)seq->seq.s[i]];
			int c1 = bitcnt_table[s[i]];
			if (i && i%step == 0) {
				z[0] -= z[2]; z[1] -= z[3];
				fwrite(z, 4, 4, stdout);
				z[0] = z[1] = z[2] = z[3] = 0;
			}
			if (c0 < 3) ++z[0];
			if (c0 == 2) ++z[1];
			if (c1 < 3) ++z[2];
			if (c1 == 2) ++z[3];
		}
		z[0] -= z[2]; z[1] -= z[3];
		fwrite(z, 4, 4, stdout);
		free(s);
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
