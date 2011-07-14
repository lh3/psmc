#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

void print_seq(const kseq_t *seq, int start, int end, int id)
{
	int i;
	printf(">%s_%d\n", seq->name.s, id);
	for (i = start; i != end; ++i) {
		int ii = i - start;
		if (ii && ii%60 == 0) putchar('\n');
		putchar(seq->seq.s[i]);
	}
	putchar('\n');
}

void split_psmcfa(int trunk_size, kseq_t *seq)
{
	while (kseq_read(seq) >= 0) {
		int i, k;
		for (i = k = 0; i < seq->seq.l; i += trunk_size) {
			if (seq->seq.l - i < trunk_size * 3 / 2) { // use the full length
				print_seq(seq, i, seq->seq.l, ++k);
				break;
			} else print_seq(seq, i, (i+trunk_size < seq->seq.l)? i+trunk_size : seq->seq.l, ++k);
		}
	}
}

int main(int argc, char *argv[])
{
	int trunk_size = 500000;
	gzFile fp;
	kseq_t *seq;
	if (argc < 2) {
		fprintf(stderr, "Usage: splitfa <in.fa> [trunk_size=%d]\n", trunk_size);
		return 1;
	}
	if (argc >= 3) trunk_size = atoi(argv[2]);
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	split_psmcfa(trunk_size, seq);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
