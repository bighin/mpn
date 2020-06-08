#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "pmatrix.h"
#include "mpn.h"

struct pmatrix_t *init_pmatrix(int nr_occupied,int nr_virtual,gsl_rng *rngctx)
{
	struct pmatrix_t *ret=malloc(sizeof(struct pmatrix_t));

	assert(ret);

	ret->dimensions=2;
	ret->nr_occupied=nr_occupied;
	ret->nr_virtual=nr_virtual;

	for(int i=0;i<PMATRIX_MAX_DIMENSIONS;i++)
		for(int j=0;j<PMATRIX_MAX_DIMENSIONS;j++)
			ret->values[i][j]=0;

	ret->values[0][0]=pmatrix_get_new_value(ret, rngctx, 0, 0);
	ret->values[1][1]=pmatrix_get_new_value(ret, rngctx, 1, 1);

	return ret;
}

void fini_pmatrix(struct pmatrix_t *pmx)
{
	if(pmx)
		free(pmx);
}

/*
	Basic matrix operations. Note that the matrix format is (i,j):

	(0,0) (0,1) (0,2)
	(1,0) (1,1) (1,2)
	(2,0) (2,1) (2,2)

	so i is the row, while j is the column.
*/

int pmatrix_get_entry(struct pmatrix_t *pmx, int i, int j)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));
	assert(pmx->values[i][j]>=0);

	if(pmx->values[i][j]==0)
		return 0;

	if(pmatrix_entry_type(i,j)==QTYPE_VIRTUAL)
		return 1+((pmx->values[i][j]-1)%pmx->nr_virtual);
	else if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
		return 1+((pmx->values[i][j]-1)%pmx->nr_occupied);

	assert(false);
	return 0;
}

void pmatrix_set_entry(struct pmatrix_t *pmx, int i, int j, int value)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));
	assert(pmx->values[i][j]>=0);
	assert(value>=0);

	if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
		assert(value<=pmx->nr_occupied);
	else if(pmatrix_entry_type(i,j)==QTYPE_VIRTUAL)
		assert(value<=pmx->nr_virtual);

	pmx->values[i][j]=value;
}

int pmatrix_get_raw_entry(struct pmatrix_t *pmx, int i, int j)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));
	assert(pmx->values[i][j]>=0);

	return pmx->values[i][j];
}

void pmatrix_set_raw_entry(struct pmatrix_t *pmx, int i, int j, int value)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));
	assert(pmx->values[i][j]>=0);
	assert(value>=0);

	pmx->values[i][j]=value;
}

void pmatrix_print(struct pmatrix_t *pmx)
{
	for(int i=0;i<pmx->dimensions;i++)
	{
		for(int j=0;j<pmx->dimensions;j++)
		{
			printf("%d ",pmatrix_get_entry(pmx, i, j));
		}

		printf("\n");
	}

	fflush(stdout);
}

void pmatrix_extend(struct pmatrix_t *pmx, gsl_rng *rngctx, int *targeti, int *targetj)
{
	int selector=gsl_rng_uniform_int(rngctx, pmx->dimensions+1);

	assert(selector>=0);
	assert(selector<(pmx->dimensions+1));
	assert(pmx->dimensions<PMATRIX_MAX_DIMENSIONS);

	/*
		There are (dims+1) ways of extending a matrix:

		a) There are dims way of taking one non-zero element and projecting it
		on the newly added row and column

	        | 0 0 1 |     | 0 0 1 0 |
	        | 1 0 0 | --> | 0 0 0 1 |
                | 0 1 0 |     | 0 1 0 0 |
			      | 1 0 0 0 |

		b) In addition to this, one can also add a '1' element in the newly-added corner.

	        | 0 0 1 |     | 0 0 1 0 |
	        | 1 0 0 | --> | 1 0 0 0 |
 	        | 0 1 0 |     | 0 1 0 0 |
        		      | 0 0 0 1 |
	*/

	if(selector<pmx->dimensions)
	{
		/*
			First we increase the dimensions of the matrix
		*/

		pmx->dimensions++;

		/*
			Then we need to cleanup the new entries
		*/

		for(int i=0;i<pmx->dimensions;i++)
		{
			pmatrix_set_entry(pmx, pmx->dimensions-1, i, 0);
			pmatrix_set_entry(pmx, i, pmx->dimensions-1, 0);
		}

		/*
			Then we look for a non-zero element on the j-th column (j=selector)
		*/

		int i,j=selector;

		for(i=0;i<(pmx->dimensions-1);i++)
			if(pmatrix_get_entry(pmx, i, j)!=0)
				break;

		/*
			There should be only at most one non-zero element on each column
		*/

		assert(i<(pmx->dimensions-1));

		for(int i2=i+1;i2<(pmx->dimensions-1);i2++)
			assert(pmatrix_get_entry(pmx, i2, j)==0);

		/*
			Finally we set the two new values.
		*/

		if(pmatrix_entry_type(i,pmx->dimensions-1)==pmatrix_entry_type(i,j))
		{
			pmatrix_set_raw_entry(pmx, i, pmx->dimensions-1, pmatrix_get_raw_entry(pmx, i, j));
		}
		else
		{
			pmatrix_set_entry(pmx, i, pmx->dimensions-1, 1);
			*targeti=i;
			*targetj=pmx->dimensions-1;
		}

		if(pmatrix_entry_type(pmx->dimensions-1, j)==pmatrix_entry_type(i,j))
		{
			pmatrix_set_raw_entry(pmx, pmx->dimensions-1, j, pmatrix_get_raw_entry(pmx, i, j));
		}
		else
		{
			pmatrix_set_entry(pmx, pmx->dimensions-1, j, 1);
			*targeti=pmx->dimensions-1;
			*targetj=j;
		}

		assert((pmatrix_entry_type(i,pmx->dimensions-1)==pmatrix_entry_type(i,j))!=
		       (pmatrix_entry_type(pmx->dimensions-1, j)==pmatrix_entry_type(i,j)));

		pmatrix_set_entry(pmx, i, j, 0);

		return;
	}
	else if(selector==pmx->dimensions)
	{
		/*
			First we increase the dimensions of the matrix
		*/

		pmx->dimensions++;

		/*
			Then we need to cleanup the new entries
		*/

		for(int i=0;i<pmx->dimensions;i++)
		{
			pmatrix_set_entry(pmx, pmx->dimensions-1, i, 0);
			pmatrix_set_entry(pmx, i, pmx->dimensions-1, 0);
		}

		/*
			We just add the new element in the corner.

			Note that it is set to 1, and the caller must take care of
			modifying it, if needed;
		*/

		*targeti=pmx->dimensions-1;
		*targetj=pmx->dimensions-1;

		pmatrix_set_entry(pmx, pmx->dimensions-1, pmx->dimensions-1, 1);

		return;
	}

	assert(false);
}

void pmatrix_squeeze(struct pmatrix_t *pmx, gsl_rng *rngctx)
{
	assert(pmx->dimensions>1);

	int i,j;

	/*
		First we look at the non-zero element in the last column
	*/

	for(i=0;i<pmx->dimensions;i++)
		if(pmatrix_get_entry(pmx, i, pmx->dimensions-1)!=0)
			break;

	/*
		There should be only at most one non-zero element on each column
	*/

	assert(i<pmx->dimensions);

	for(int i2=i+1;i2<pmx->dimensions;i2++)
		assert(pmatrix_get_entry(pmx, i2, pmx->dimensions-1)==0);

	/*
		Then we look at the non-zero element in the last row
	*/

	for(j=0;j<pmx->dimensions;j++)
		if(pmatrix_get_entry(pmx, pmx->dimensions-1, j)!=0)
			break;

	/*
		There should be only at most one non-zero element on each row
	*/

	assert(j<pmx->dimensions);

	for(int j2=j+1;j2<pmx->dimensions;j2++)
		assert(pmatrix_get_entry(pmx, pmx->dimensions-1, j2)==0);

	/*
		Finally we can squeeze the matrix.

		Case b): the element to remove is in the bottom right corner.
	*/

	if(i==(pmx->dimensions-1))
	{
		assert(j==(pmx->dimensions-1));
		pmatrix_set_entry(pmx, i, j, 0);

		pmx->dimensions--;

		return;
	}

	/*
		...otherwise we have case a) two different entries on the last row
		and column are being removed, and a new entry is added to the matrix.
	*/

	assert(pmatrix_get_entry(pmx, i, pmx->dimensions-1)!=0);
	assert(pmatrix_get_entry(pmx, pmx->dimensions-1, j)!=0);

	int newvalue=-1;

	if(pmatrix_entry_type(i, pmx->dimensions-1)==pmatrix_entry_type(i,j))
		newvalue=pmatrix_get_raw_entry(pmx, i, pmx->dimensions-1);
	else if(pmatrix_entry_type(pmx->dimensions-1, j)==pmatrix_entry_type(i,j))
		newvalue=pmatrix_get_raw_entry(pmx, pmx->dimensions-1, j);

	assert(newvalue!=-1);

	pmatrix_set_entry(pmx, i, pmx->dimensions-1, 0);
	pmatrix_set_entry(pmx, pmx->dimensions-1, j, 0);

	pmatrix_set_raw_entry(pmx, i, j, newvalue);

	pmx->dimensions--;
}

void pmatrix_swap_rows(struct pmatrix_t *pmx, int i1, int i2, gsl_rng *rngctx)
{
	assert(i1>=0);
	assert(i1<pmx->dimensions);
	assert(i2>=0);
	assert(i2<pmx->dimensions);
	assert(pmatrix_check_consistency(pmx)==true);

	for(int j=0;j<pmx->dimensions;j++)
	{
		int tmp=pmatrix_get_raw_entry(pmx, i1, j);

		pmatrix_set_raw_entry(pmx, i1, j, pmatrix_get_raw_entry(pmx, i2, j));
		pmatrix_set_raw_entry(pmx, i2, j, tmp);
	}

	assert(pmatrix_check_consistency(pmx)==true);
}

void pmatrix_swap_cols(struct pmatrix_t *pmx, int j1, int j2, gsl_rng *rngctx)
{
	assert(j1>=0);
	assert(j1<pmx->dimensions);
	assert(j2>=0);
	assert(j2<pmx->dimensions);
	assert(pmatrix_check_consistency(pmx)==true);

	for(int i=0;i<pmx->dimensions;i++)
	{
		int tmp=pmatrix_get_raw_entry(pmx, i, j1);

		pmatrix_set_raw_entry(pmx, i, j1, pmatrix_get_raw_entry(pmx, i, j2));
		pmatrix_set_raw_entry(pmx, i, j2, tmp);
	}

	assert(pmatrix_check_consistency(pmx)==true);
}

bool pmatrix_check_consistency(struct pmatrix_t *pmx)
{
	for(int i=0;i<pmx->dimensions;i++)
	{
		int nonzero_entries=0;

		for(int j=0;j<pmx->dimensions;j++)
			if(pmatrix_get_entry(pmx, i, j)!=0)
				nonzero_entries++;

		if(nonzero_entries!=1)
			return false;
	}

	for(int j=0;j<pmx->dimensions;j++)
	{
		int nonzero_entries=0;

		for(int i=0;i<pmx->dimensions;i++)
			if(pmatrix_get_entry(pmx, i, j)!=0)
				nonzero_entries++;

		if(nonzero_entries!=1)
			return false;
	}

	for(int i=0;i<pmx->dimensions;i++)
	{
		for(int j=0;j<pmx->dimensions;j++)
		{
			int entry;

			if((entry=pmatrix_get_entry(pmx, i, j))!=0)
			{
				if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
				{
					if(!((entry>0)&&(entry<=pmx->nr_occupied)))
						return false;
				}
				else
				{
					if(!((entry>0)&&(entry<=pmx->nr_virtual)))
						return false;
				}
			}
		}
	}

	return true;
}

int pmatrix_entry_type(int i,int j)
{
	if(i<=j)
		return QTYPE_OCCUPIED;

	return QTYPE_VIRTUAL;
}

int pmatrix_get_new_value(struct pmatrix_t *pmx, gsl_rng *rngctx, int i, int j)
{
	return 1+gsl_rng_uniform_int(rngctx, pmx->nr_occupied*pmx->nr_virtual);
}
