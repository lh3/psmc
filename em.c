#include <math.h>
#include <string.h>
#include <assert.h>
#include "psmc.h"
#include "khmm.h"
#include "kmin.h"

typedef struct {
	psmc_par_t *pp;
	psmc_data_t *pd;
	hmm_exp_t *he;
	int cnt;
} em_aux_t;

static double func(int n, double *params, void *data)
{
	em_aux_t *aux = (em_aux_t*)data;
	int i;
	assert(n == aux->pd->n_params);
	++aux->cnt;
	for (i = 0; i != n; ++i)
		aux->pd->params[i] = fabs(params[i]);
	psmc_update_hmm(aux->pp, aux->pd);
	return -hmm_Q(aux->pd->hp, aux->he);
}

double psmc_em(psmc_par_t *pp, psmc_data_t *pd)
{
	hmm_par_t *hp = pd->hp;
	hmm_exp_t *he_sum;
	int i;
	double LL = 0.0;
	he_sum = hmm_new_exp(hp);
	hmm_pre_backward(hp);
	// expectation
	for (i = 0; i != pp->n_seqs; ++i) {
		hmm_exp_t *he;
		hmm_data_t *hd;
		// make the sequence
		psmc_seq_t *s = pp->seqs + i;
		char *seq;
		seq = (char*)calloc(s->L+1, 1);
		memcpy(seq, s->seq, s->L);
		hd = hmm_new_data(s->L, seq, hp);
		// forward-backward, expectation
		hmm_forward(hp, hd);
		hmm_backward(hp, hd);
		LL += hmm_lk(hd);
		he = hmm_expect(hp, hd);
		hmm_add_expect(he, he_sum);
		// free
		hmm_delete_exp(he);
		hmm_delete_data(hd);
		free(seq);
	}
	{ // maximization
		double *params;
		em_aux_t aux;
		hmm_Q0(hp, he_sum);
		pd->lk = LL;
		params = (double*)calloc(pd->n_params, sizeof(double));
		memcpy(params, pd->params, sizeof(double) * pd->n_params);
		aux.pp = pp; aux.pd = pd; aux.he = he_sum; aux.cnt = 0;
		pd->Q0 = hmm_Q(hp, he_sum);
		pd->Q1 = -kmin_hj(func, pd->n_params, params, &aux, KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);
		fprintf(pp->fpout, "IT\t%d\n", aux.cnt);
		free(params);
	}
	{ // update pd->post_sigma
		double sum = 0.0;
		int k;
		for (k = 0; k <= pp->n; ++k) sum += he_sum->E[0][k] + he_sum->E[1][k];
		for (k = 0; k <= pp->n; ++k) pd->post_sigma[k] = (he_sum->E[0][k] + he_sum->E[1][k]) / sum;
	}
	// free
	hmm_delete_exp(he_sum);
	return pd->Q1;
}
