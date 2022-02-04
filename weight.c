#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "weight.h"
#include "mpn.h"
#include "permutations.h"
#include "auxx.h"
#include "weight2.h"

struct amatrix_t *init_amatrix_from_amatrix(struct amatrix_t *amx)
{
	struct amatrix_t *ret=malloc(sizeof(struct amatrix_t));

	assert(ret!=NULL);

	ret->ectx=amx->ectx;
	ret->rng_ctx=amx->rng_ctx;

	ret->nr_occupied=amx->nr_occupied;
	ret->nr_virtual=amx->nr_virtual;

	assert((ret->ectx!=NULL)&&(ret->rng_ctx!=NULL));

	ret->pmxs[0]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);
	ret->pmxs[1]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);

	ret->config=amx->config;

	assert(ret->config!=NULL);

	ret->cached_weight=0.0f;
	ret->cached_weight_is_valid=false;

	return ret;
}

double amatrix_projection_multiplicity(struct amatrix_t *amx)
{
	int nr_occupied_entries,nr_virtual_entries;

	nr_occupied_entries=nr_virtual_entries=0;
	for(int i=0;i<amx->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<amx->pmxs[0]->dimensions;j++)
		{
			if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
			{
				if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
					nr_occupied_entries++;
				else
					nr_virtual_entries++;
			}

			if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
			{
				if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
					nr_occupied_entries++;
				else
					nr_virtual_entries++;
			}
		}
	}

	return pow(amx->nr_virtual,nr_occupied_entries)*pow(amx->nr_occupied,nr_virtual_entries);
}

/*
	Here we calculate the weight associated to a 'amatrix'
*/

double amatrix_weight(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	if(amx->cached_weight_is_valid==true)
	{

#ifndef NDEBUG
		amx->cached_weight_is_valid=false;

		double w1,w2;

		w1=amx->cached_weight;
		w2=amatrix_weight(amx);
		assert(gsl_fcmp(w1,w2,1e-6)==0);

		amx->cached_weight_is_valid=true;
#endif

		return amx->cached_weight;
	}

	/*
		Dimension 1 is a special case that does not need the evaluation
		of the incidence matrix.

		Note that dimension 1 (a bit unexpectedly) corresponds to the Hartree-Fock energy.
	*/

	if(amx->pmxs[0]->dimensions==1)
	{
		int a,b;
		double weight=0.0f,multiplicity;

		a=pmatrix_get_entry(amx->pmxs[0],0,0);
		b=pmatrix_get_entry(amx->pmxs[1],0,0);

		if(a==b)
			weight+=get_hdiag(amx->ectx,a-1);

		weight+=0.5*get_eri(amx->ectx,a-1,b-1,a-1,b-1);
		weight+=get_enuc(amx->ectx)*pow(amx->nr_occupied,-2.0f);

		multiplicity=1.0f;

		amx->cached_weight=weight/multiplicity/amatrix_projection_multiplicity(amx);
		amx->cached_weight_is_valid=true;

		return amx->cached_weight;
	}
	else
	{
		double ret=0;

		struct label_t labels[MAX_LABELS];
		int ilabels=0;

		if(amatrix_check_connectedness(amx)==true)
		{
			gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
			ret=incidence_to_weight(incidence, labels, &ilabels, amx);
			gsl_matrix_int_free(incidence);

			ret/=amatrix_projection_multiplicity(amx);
		}

		amx->cached_weight=ret;
		amx->cached_weight_is_valid=true;

		return amx->cached_weight;
	}

	assert(false);
	return 0.0f;
}
