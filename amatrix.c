#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "pmatrix.h"
#include "loaderis.h"
#include "mpn.h"
#include "multiplicity.h"
#include "cache.h"
#include "auxx.h"

struct amatrix_t *init_amatrix(struct configuration_t *config)
{
	struct amatrix_t *ret=malloc(sizeof(struct amatrix_t));

	assert(ret!=NULL);

	if((config!=NULL)&&(config->erisfile!=NULL))
	{
		FILE *in=fopen(config->erisfile, "r");

		if(!in)
			return NULL;

		ret->ectx=malloc(sizeof(struct energies_ctx_t));
		assert(ret->ectx!=NULL);
		load_energies(in, ret->ectx);

		fclose(in);

		ret->nr_occupied=ret->ectx->nocc;
		ret->nr_virtual=ret->ectx->nvirt;
	}
	else
	{
		ret->ectx=NULL;

		/*
			Here we set some default values. However, if the ERI file is not loaded,
			the number of occupied/virtual orbitals doesn't make a lot of sense.
		*/

		ret->nr_occupied=10;
		ret->nr_virtual=16;
	}

	ret->rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
	assert(ret->rng_ctx!=NULL);

	if((config!=NULL)&&(config->seedrng==true))
		seed_rng(ret->rng_ctx);

	ret->pmxs[0]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);
	ret->pmxs[1]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);

	ret->config=config;

	ret->cached_weight.weight=0.0f;
	ret->cached_weight.l=0;
	ret->cached_weight.h=0;
	ret->cached_weight_is_valid=false;

	return ret;
}

void fini_amatrix(struct amatrix_t *amx)
{
	if(amx)
	{
		if(amx->ectx)
			free(amx->ectx);

		fini_pmatrix(amx->pmxs[0]);
		fini_pmatrix(amx->pmxs[1]);

		free(amx);
	}
}

/*
	Basic matrix operations
*/

int amatrix_get_entry(struct amatrix_t *amx, int i, int j)
{
	assert(amx);
	assert(i>=0);
	assert(i<amx->pmxs[0]->dimensions);
	assert(j>=0);
	assert(j<amx->pmxs[0]->dimensions);

	int cnt=0;

	if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
		cnt++;

	if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
		cnt++;

	return cnt;
}

/*
	We can save the contents of a amatrix in a 'backup' structure using the 'save' operation,
	and then (if needed) we can reload the old contents using the 'restore' operation.
*/

void amatrix_save(struct amatrix_t *amx, struct amatrix_backup_t *backup)
{
	backup->dimensions[0]=amx->pmxs[0]->dimensions;
	backup->dimensions[1]=amx->pmxs[1]->dimensions;

	assert(backup->dimensions[0]==backup->dimensions[1]);
	assert(backup->dimensions[0]<PMATRIX_MAX_DIMENSIONS);

	for(int i=0;i<backup->dimensions[0];i++)
	{
		for(int j=0;j<backup->dimensions[0];j++)
		{
			backup->values[0][i][j]=amx->pmxs[0]->values[i][j];
			backup->values[1][i][j]=amx->pmxs[1]->values[i][j];
		}
	}

	memcpy(&backup->cached_result, &amx->cached_weight, sizeof(struct amatrix_weight_t));
	backup->cached_result_is_valid=amx->cached_weight_is_valid;
}

void amatrix_restore(struct amatrix_t *amx, struct amatrix_backup_t *backup)
{
	amx->pmxs[0]->dimensions=backup->dimensions[0];
	amx->pmxs[1]->dimensions=backup->dimensions[1];

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions<PMATRIX_MAX_DIMENSIONS);

	for(int i=0;i<amx->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<amx->pmxs[0]->dimensions;j++)
		{
			amx->pmxs[0]->values[i][j]=backup->values[0][i][j];
			amx->pmxs[1]->values[i][j]=backup->values[1][i][j];
		}
	}

	memcpy(&amx->cached_weight, &backup->cached_result, sizeof(struct amatrix_weight_t));
	amx->cached_weight_is_valid=backup->cached_result_is_valid;
}

/*
	Printing an 'amatrix'
*/

void amatrix_print(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			int entry=0;

			if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
				entry++;

			if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
				entry++;

			printf("%d ",entry);
		}
		
		printf("\n");
	}
	
	fflush(stdout);
}

/*
	Basic internal consistency check
*/

bool amatrix_check_consistency(struct amatrix_t *amx)
{
	if((amx->pmxs[0]->dimensions<1)||(amx->pmxs[1]->dimensions<1))
		return false;

	if(amx->pmxs[0]->dimensions!=amx->pmxs[1]->dimensions)
		return false;

	return pmatrix_check_consistency(amx->pmxs[0])&&pmatrix_check_consistency(amx->pmxs[1]);
}

/*
	Check if the matrix belongs to the physical sector.

	Note that an amatrix_t object can be unphysical, but still consistent,
	as defined by amatrix_check_consistency().
*/

bool amatrix_is_physical(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	/*
		In dimension 1 the logic is the inverted, the matrices MUST
		have non-zero elements on the diagonal.
	*/

	if(dimensions==1)
	{
		if((pmatrix_get_entry(amx->pmxs[0], 0, 0)!=0)&&
		   (pmatrix_get_entry(amx->pmxs[1], 0, 0)!=0))
			return true;
		else
			return false;
	}

	for(int i=0;i<dimensions;i++)
	{
		if(pmatrix_get_entry(amx->pmxs[0], i, i)!=0)
			return false;

		if(pmatrix_get_entry(amx->pmxs[1], i, i)!=0)
			return false;
	}

	return true;
}

void amatrix_to_python(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	printf("[");

	for(int i=0;i<dimensions;i++)
	{
		printf("[");

		for(int j=0;j<dimensions;j++)
		{
			if(j!=(dimensions-1))
				printf("%d, ", amatrix_get_entry(amx, i, j));
			else
				printf("%d ", amatrix_get_entry(amx, i, j));
		}

		if(i!=(dimensions-1))
			printf("], ");
		else
			printf("]");
	}

	printf("]\n");
}

void amatrix_to_wolfram(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	printf("{");

	for(int i=0;i<dimensions;i++)
	{
		printf("{");

		for(int j=0;j<dimensions;j++)
		{
			if(j!=(dimensions-1))
				printf("%d, ", amatrix_get_entry(amx, i, j));
			else
				printf("%d ", amatrix_get_entry(amx, i, j));
		}

		if(i!=(dimensions-1))
			printf("}, ");
		else
			printf("}");
	}

	printf("}\n");
}

bool actual_amatrix_check_connectedness(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(dimensions>=1);
	if(dimensions==1)
		return true;

	bool result=true;

	gsl_matrix_int *adjacency=gsl_matrix_int_alloc(dimensions, dimensions);

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			int value=0;

			if(amx->pmxs[0]->values[i][j]!=0)
				value++;

			if(amx->pmxs[1]->values[i][j]!=0)
				value++;

			gsl_matrix_int_set(adjacency, i, j, value);
		}
	}

	gsl_matrix_int *powers=gsl_matrix_int_alloc(dimensions, dimensions);
	gsl_matrix_int *tmp=gsl_matrix_int_alloc(dimensions, dimensions);
	gsl_matrix_int *C=gsl_matrix_int_alloc(dimensions, dimensions);

	gsl_matrix_int_set_zero(C);
	gsl_matrix_int_set_identity(powers);

	for(int power=0;power<=(dimensions-1);power++)
	{
		/*
			The power-th power is added to to C
		*/

		gsl_matrix_int_add(C,powers);

		/*
			The (power+1)-th power is calculated
		*/

		gsl_matrix_int_mul(powers,adjacency,tmp);
		gsl_matrix_int_memcpy(powers,tmp);
	}

	/*
		If at least one entry is zero, then the graph is not connected.
	*/

	for(int i=0;i<dimensions;i++)
		for(int j=0;j<dimensions;j++)
			if(gsl_matrix_int_get(C,i,j)==0)
				result=false;

	gsl_matrix_int_free(C);
	gsl_matrix_int_free(tmp);
	gsl_matrix_int_free(powers);
	gsl_matrix_int_free(adjacency);

	return result;
}

bool amatrix_check_connectedness(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	if((dimensions>1)&&(dimensions<=amatrix_cache_max_dimensions)&&(amatrix_cache_is_enabled==true))
	{
		assert(cached_amatrix_check_connectedness(amx)==actual_amatrix_check_connectedness(amx));
		return cached_amatrix_check_connectedness(amx);
	}

	return actual_amatrix_check_connectedness(amx);
}

/*
	Here we calculate the weight associated to a 'amatrix'
*/

struct amatrix_weight_t amatrix_weight(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	if(amx->cached_weight_is_valid==true)
	{

#ifndef NDEBUG
		amx->cached_weight_is_valid=false;
		assert(amx->cached_weight.weight==amatrix_weight(amx).weight);
		amx->cached_weight_is_valid=true;
#endif

		return amx->cached_weight;
	}

	if(amatrix_check_connectedness(amx)==false)
	{
		struct amatrix_weight_t ret={0.0,0,0};
		return ret;
	}

	/*
		Dimension 1 is a special case that does not need the evaluation
		of the incidence matrix.

		Note that dimension 1 corresponds to the Hartree-Fock energy.
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

		amx->cached_weight.weight=amx->config->bias+weight/multiplicity;
		amx->cached_weight.l=0;
		amx->cached_weight.h=2;
		amx->cached_weight_is_valid=true;

		return amx->cached_weight;
	}


	struct amatrix_weight_t ret;

	gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, ret.labels, &ret.ilabels);
	ret=incidence_to_weight(incidence, ret.labels, &ret.ilabels, amx);
	gsl_matrix_int_free(incidence);

	ret.weight=amx->config->bias+ret.weight/amatrix_multiplicity(amx);

	amx->cached_weight=ret;
	amx->cached_weight_is_valid=true;

	return ret;
}
