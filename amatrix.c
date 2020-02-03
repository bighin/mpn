#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "pmatrix.h"
#include "auxx.h"
#include "loaderis.h"
#include "mpn.h"
#include "multiplicity.h"

struct amatrix_t *init_amatrix(char *energies_dot_dat)
{
	struct amatrix_t *ret=malloc(sizeof(struct amatrix_t));

	assert(ret!=NULL);

	FILE *in=fopen(energies_dot_dat,"r");

	if(!in)
		return NULL;

	ret->ectx=malloc(sizeof(struct energies_ctx_t));
	assert(ret->ectx!=NULL);
	load_energies(in, ret->ectx);

	fclose(in);

	ret->nr_occupied=ret->ectx->nocc;
	ret->nr_virtual=ret->ectx->nvirt;

	ret->rng_ctx=gsl_rng_alloc(gsl_rng_mt19937);
	assert(ret->rng_ctx!=NULL);
	seed_rng(ret->rng_ctx);

	ret->pmxs[0]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);
	ret->pmxs[1]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);

	ret->bias=0.0f;
	ret->unphysicalpenalty=1.0f;
	ret->maxorder=16;

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
	assert(backup->dimensions[0]<IMATRIX_MAX_DIMENSIONS);

	for(int i=0;i<backup->dimensions[0];i++)
	{
		for(int j=0;j<backup->dimensions[0];j++)
		{
			backup->values[0][i][j]=amx->pmxs[0]->values[i][j];
			backup->values[1][i][j]=amx->pmxs[1]->values[i][j];
		}
	}
}

void amatrix_restore(struct amatrix_t *amx, struct amatrix_backup_t *backup)
{
	amx->pmxs[0]->dimensions=backup->dimensions[0];
	amx->pmxs[1]->dimensions=backup->dimensions[1];

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions<IMATRIX_MAX_DIMENSIONS);

	for(int i=0;i<amx->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<amx->pmxs[0]->dimensions;j++)
		{
			amx->pmxs[0]->values[i][j]=backup->values[0][i][j];
			amx->pmxs[1]->values[i][j]=backup->values[1][i][j];
		}
	}
}

/*
	Printing and debug functions
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

bool amatrix_check_consistency(struct amatrix_t *amx)
{
	if((amx->pmxs[0]->dimensions<1)||(amx->pmxs[1]->dimensions<1))
		return false;

	if(amx->pmxs[0]->dimensions!=amx->pmxs[1]->dimensions)
		return false;

	return pmatrix_check_consistency(amx->pmxs[0])&&pmatrix_check_consistency(amx->pmxs[1]);
}

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

bool amatrix_check_connectedness(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(dimensions>=1);
	if(dimensions==1)
		return true;

	bool result=true;

	gsl_matrix_int *adjacency=gsl_matrix_int_alloc(dimensions, dimensions);
	gsl_matrix_int *powers=gsl_matrix_int_alloc(dimensions, dimensions);
	gsl_matrix_int *tmp=gsl_matrix_int_alloc(dimensions, dimensions);
	gsl_matrix_int *C=gsl_matrix_int_alloc(dimensions, dimensions);

	gsl_matrix_int_set_zero(C);
	gsl_matrix_int_set_identity(powers);

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

/*
	Here we calculate the weight associated to a 'amatrix'
*/

double amatrix_weight(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	if(amatrix_check_connectedness(amx)==false)
		return 0.0f;

	/*
		Dimension 1 is a special case that does not need the evaluation
		of the incidence matrix.
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

#warning Is this correct?

		multiplicity=1.0f;

		return amx->bias+multiplicity*weight;
	}

	struct label_t labels[MAX_LABELS];
	int ilabels=0;

	gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);

	double weight,multiplicity;

	weight=incidence_to_weight(incidence, labels, &ilabels, amx->ectx, amx->unphysicalpenalty, false);
	multiplicity=amatrix_multiplicity(amx);

	return amx->bias+weight/multiplicity;
}
