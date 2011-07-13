#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "psmc.h"

int main(int argc, char *argv[])
{
	int i;
	psmc_par_t *pp;
	psmc_data_t *pd;
	srand48(time(0) ^ getpid());
	pp = psmc_parse_cli(argc, argv);
	pd = psmc_new_data(pp);
	fprintf(pp->fpout, "RD\t0\n");
	psmc_print_data(pp, pd);
	for (i = 0; i != pp->n_iters; ++i) {
		psmc_em(pp, pd);
		fprintf(pp->fpout, "RD\t%d\n", i+1);
		psmc_print_data(pp, pd);
	}
	if ((pp->flag & PSMC_F_DECODE) || pp->fpcnt) psmc_decode(pp, pd);
	if (pp->flag & PSMC_F_SIMU) psmc_simulate(pp, pd);
	psmc_delete_data(pd);
	psmc_delete_par(pp);
	return 0;
}
