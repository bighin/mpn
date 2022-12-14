#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "pmatrix.h"
#include "loaderis.h"
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

	ret->cached_weight=0.0f;
	ret->cached_weight_is_valid=false;

	return ret;
}

void fini_amatrix(struct amatrix_t *amx,bool free_ectx)
{
	if(amx)
	{
		if((amx->ectx)&&(free_ectx==true))
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

	backup->cached_result=amx->cached_weight;
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

	amx->cached_weight=backup->cached_result;
	amx->cached_weight_is_valid=backup->cached_result_is_valid;
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

/*
	Printing an 'amatrix', in various formats
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

void amatrix_print_detailed(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			int x=pmatrix_get_entry(amx->pmxs[0], i, j);
			int y=pmatrix_get_entry(amx->pmxs[1], i, j);

			if((x==0)&&(y==0))
				printf("0");
			else if((x==0)&&(y!=0))
				printf("1[%d]",y);
			else if((x!=0)&&(y==0))
				printf("1[%d]",x);
			else if((x!=0)&&(y!=0))
				printf("2[%d,%d]",x,y);
		}

		printf("\n");
	}

	fflush(stdout);
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
	fflush(stdout);
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
	fflush(stdout);
}

/*
	Checks if the amatrix represents a connected graph
*/

bool gsl_matrix_int_check_connectedness(gsl_matrix_int *adjacency,int dimensions)
{
	bool result=true;

	if(dimensions==1)
		return true;

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

	return result;
}

bool actual_amatrix_check_connectedness(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(dimensions>=1);
	if(dimensions==1)
		return true;

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

	bool result=gsl_matrix_int_check_connectedness(adjacency,dimensions);

	gsl_matrix_int_free(adjacency);

	return result;
}

bool amatrix_check_connectedness(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	if(dimensions==1)
	{
		return true;
	}

	if((dimensions>1)&&(dimensions<=amatrix_cache_max_dimensions)&&(amatrix_cache_is_enabled==true))
	{
		assert(cached_amatrix_check_connectedness(amx)==actual_amatrix_check_connectedness(amx));
		return cached_amatrix_check_connectedness(amx);
	}

	return actual_amatrix_check_connectedness(amx);
}
