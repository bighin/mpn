#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gsl/gsl_rng.h>

#include "pmatrix.h"
#include "mpn.h"
#include "auxx.h"

struct pmatrix_t *init_pmatrix(int nr_occupied,int nr_virtual,gsl_rng *rngctx)
{
	struct pmatrix_t *ret=malloc(sizeof(struct pmatrix_t));

	assert(ret);

	ret->dimensions=1;
	ret->nr_occupied=nr_occupied;
	ret->nr_virtual=nr_virtual;

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

	return pmx->values[i][j];
}

void pmatrix_set_entry(struct pmatrix_t *pmx, int i, int j, int value)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));
	assert(pmx->values[i][j]>=0);
	assert(value>=0);

	pmx->values[i][j]=value;
}

void pmatrix_inc_entry(struct pmatrix_t *pmx, int i, int j)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));

	pmx->values[i][j]++;

	int value=pmx->values[i][j];
	assert((value==0)||(value==1)||(value==2));
}

void pmatrix_dec_entry(struct pmatrix_t *pmx, int i, int j)
{
	assert((i>=0)&&(i<pmx->dimensions));
	assert((j>=0)&&(j<pmx->dimensions));

	pmx->values[i][j]--;

	int value=pmx->values[i][j];
	assert((value==0)||(value==1)||(value==2));
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

int pmatrix_sum_row(struct pmatrix_t *pmx, int row)
{
	int ret=0;

	for(int c=0;c<pmx->dimensions;c++)
		ret+=pmatrix_get_entry(pmx, row, c);

	return ret;
}

int pmatrix_sum_column(struct pmatrix_t *pmx, int column)
{
	int ret=0;

	for(int c=0;c<pmx->dimensions;c++)
		ret+=pmatrix_get_entry(pmx, c, column);

	return ret;
}

int pmatrix_trace(struct pmatrix_t *pmx)
{
	int ret=0;

	for(int c=0;c<pmx->dimensions;c++)
		ret+=pmatrix_get_entry(pmx, c, c);

	return ret;
}

double pmatrix_extend(struct pmatrix_t *pmx, gsl_rng *rngctx, int *targeti, int *targetj)
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
			pmatrix_set_entry(pmx, i, pmx->dimensions-1, pmatrix_get_entry(pmx, i, j));
		}
		else
		{
			pmatrix_set_entry(pmx, i, pmx->dimensions-1, 1);
			*targeti=i;
			*targetj=pmx->dimensions-1;
		}

		if(pmatrix_entry_type(pmx->dimensions-1, j)==pmatrix_entry_type(i,j))
		{
			pmatrix_set_entry(pmx, pmx->dimensions-1, j, pmatrix_get_entry(pmx, i, j));
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

		return pmx->dimensions;
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

		return pmx->dimensions;
	}

	assert(false);
	return 0.0f;
}

double pmatrix_squeeze(struct pmatrix_t *pmx, gsl_rng *rngctx)
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

		return 1.0f/(pmx->dimensions+1);
	}

	/*
		...otherwise we have case a) two different entries on the last row
		and column are being removed, and a new entry is added to the matrix.
	*/

	assert(pmatrix_get_entry(pmx, i, pmx->dimensions-1)!=0);
	assert(pmatrix_get_entry(pmx, pmx->dimensions-1, j)!=0);

	int newvalue=-1;

	if(pmatrix_entry_type(i, pmx->dimensions-1)==pmatrix_entry_type(i,j))
		newvalue=pmatrix_get_entry(pmx, i, pmx->dimensions-1);
	else if(pmatrix_entry_type(pmx->dimensions-1, j)==pmatrix_entry_type(i,j))
		newvalue=pmatrix_get_entry(pmx, pmx->dimensions-1, j);

	assert(newvalue!=-1);

	pmatrix_set_entry(pmx, i, pmx->dimensions-1, 0);
	pmatrix_set_entry(pmx, pmx->dimensions-1, j, 0);

	pmatrix_set_entry(pmx, i, j, newvalue);

	pmx->dimensions--;

	return 1.0f/(pmx->dimensions+1);
}

void pmatrix_swap_rows(struct pmatrix_t *pmx, int i1, int i2, int update[2], int reverse[2], gsl_rng *rngctx)
{
	assert(i1>=0);
	assert(i1<pmx->dimensions);
	assert(i2>=0);
	assert(i2<pmx->dimensions);
	assert(pmatrix_check_consistency(pmx)==true);

	update[QTYPE_OCCUPIED]=update[QTYPE_VIRTUAL]=reverse[QTYPE_OCCUPIED]=reverse[QTYPE_VIRTUAL]=0;

#ifndef NDEBUG
	int x1=0,x2=0;
#endif

	/*
		How many occupied and virtual entries are there in the i1-th row?
		How will things change when this becomes the i2-th row?
	*/

	int oldvirtual,oldoccupied,newvirtual,newoccupied;

	oldvirtual=oldoccupied=newvirtual=newoccupied=0;
	for(int j=0;j<pmx->dimensions;j++)
	{
		if(pmatrix_get_entry(pmx, i1, j)!=0)
		{
			if(pmatrix_entry_type(i1,j)==QTYPE_OCCUPIED)
				oldoccupied++;
			else
				oldvirtual++;

			if(pmatrix_entry_type(i2,j)==QTYPE_OCCUPIED)
				newoccupied++;
			else
				newvirtual++;

			/*
				Going from occupied to virtual
			*/

			if((pmatrix_entry_type(i1,j)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i2,j)==QTYPE_VIRTUAL))
				pmatrix_set_entry(pmx, i1, j, pmatrix_get_new_value(pmx, rngctx, i2, j));

			/*
				Going from virtual to occupied
			*/

			if((pmatrix_entry_type(i1,j)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i2,j)==QTYPE_OCCUPIED))
				pmatrix_set_entry(pmx, i1, j, pmatrix_get_new_value(pmx, rngctx, i2, j));

#ifndef NDEBUG
			if((pmatrix_entry_type(i1,j)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i2,j)==QTYPE_VIRTUAL)) x2++;
			if((pmatrix_entry_type(i1,j)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i2,j)==QTYPE_OCCUPIED)) x1++;
#endif
		}
	}

	update[QTYPE_OCCUPIED]+=positive_part(newoccupied-oldoccupied);
	update[QTYPE_VIRTUAL]+=positive_part(newvirtual-oldvirtual);
	reverse[QTYPE_OCCUPIED]+=negative_part(newoccupied-oldoccupied);
	reverse[QTYPE_VIRTUAL]+=negative_part(newvirtual-oldvirtual);

	/*
		Same for the i2-th row.
	*/

	oldvirtual=oldoccupied=newvirtual=newoccupied=0;
	for(int j=0;j<pmx->dimensions;j++)
	{
		if(pmatrix_get_entry(pmx, i2, j)!=0)
		{
			if(pmatrix_entry_type(i2,j)==QTYPE_OCCUPIED)
				oldoccupied++;
			else
				oldvirtual++;

			if(pmatrix_entry_type(i1,j)==QTYPE_OCCUPIED)
				newoccupied++;
			else
				newvirtual++;

			/*
				Going from occupied to virtual
			*/

			if((pmatrix_entry_type(i2,j)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i1,j)==QTYPE_VIRTUAL))
				pmatrix_set_entry(pmx, i2, j, pmatrix_get_new_value(pmx, rngctx, i1, j));

			/*
				Going from virtual to occupied
			*/

			if((pmatrix_entry_type(i2,j)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i1,j)==QTYPE_OCCUPIED))
				pmatrix_set_entry(pmx, i2, j, pmatrix_get_new_value(pmx, rngctx, i1, j));

#ifndef NDEBUG
			if((pmatrix_entry_type(i2,j)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i1,j)==QTYPE_VIRTUAL)) x2++;
			if((pmatrix_entry_type(i2,j)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i1,j)==QTYPE_OCCUPIED)) x1++;
#endif
		}
	}

	update[QTYPE_OCCUPIED]+=positive_part(newoccupied-oldoccupied);
	update[QTYPE_VIRTUAL]+=positive_part(newvirtual-oldvirtual);
	reverse[QTYPE_OCCUPIED]+=negative_part(newoccupied-oldoccupied);
	reverse[QTYPE_VIRTUAL]+=negative_part(newvirtual-oldvirtual);

#ifndef NDEBUG
	assert(update[QTYPE_OCCUPIED]==x1);
	assert(update[QTYPE_VIRTUAL]==x2);
#endif

	/*
		Finally we actually swap the rows!
	*/

	for(int j=0;j<pmx->dimensions;j++)
	{
		int tmp=pmatrix_get_entry(pmx, i1, j);

		pmatrix_set_entry(pmx, i1, j, pmatrix_get_entry(pmx, i2, j));
		pmatrix_set_entry(pmx, i2, j, tmp);
	}

	assert(pmatrix_check_consistency(pmx)==true);
}

void pmatrix_swap_cols(struct pmatrix_t *pmx, int j1, int j2, int update[2], int reverse[2], gsl_rng *rngctx)
{
	assert(j1>=0);
	assert(j1<pmx->dimensions);
	assert(j2>=0);
	assert(j2<pmx->dimensions);
	assert(pmatrix_check_consistency(pmx)==true);

#ifndef NDEBUG
	int x1=0,x2=0;
#endif

	update[QTYPE_OCCUPIED]=update[QTYPE_VIRTUAL]=reverse[QTYPE_OCCUPIED]=reverse[QTYPE_VIRTUAL]=0;

	/*
		How many occupied and virtual entries are there in the i1-th row?
		How will things change when this becomes the i2-th row?
	*/

	int oldvirtual,oldoccupied,newvirtual,newoccupied;

	oldvirtual=oldoccupied=newvirtual=newoccupied=0;
	for(int i=0;i<pmx->dimensions;i++)
	{
		if(pmatrix_get_entry(pmx, i, j1)!=0)
		{
			if(pmatrix_entry_type(i,j1)==QTYPE_OCCUPIED)
				oldoccupied++;
			else
				oldvirtual++;

			if(pmatrix_entry_type(i,j2)==QTYPE_OCCUPIED)
				newoccupied++;
			else
				newvirtual++;

			/*
				Going from occupied to virtual
			*/

			if((pmatrix_entry_type(i,j1)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i,j2)==QTYPE_VIRTUAL))
				pmatrix_set_entry(pmx, i, j1, pmatrix_get_new_value(pmx, rngctx, i, j2));

			/*
				Going from virtual to occupied
			*/

			if((pmatrix_entry_type(i,j1)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i,j2)==QTYPE_OCCUPIED))
				pmatrix_set_entry(pmx, i, j1, pmatrix_get_new_value(pmx, rngctx, i, j2));

#ifndef NDEBUG
			if((pmatrix_entry_type(i,j1)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i,j2)==QTYPE_VIRTUAL)) x2++;
			if((pmatrix_entry_type(i,j1)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i,j2)==QTYPE_OCCUPIED)) x1++;
#endif
		}
	}

	update[QTYPE_OCCUPIED]+=positive_part(newoccupied-oldoccupied);
	update[QTYPE_VIRTUAL]+=positive_part(newvirtual-oldvirtual);
	reverse[QTYPE_OCCUPIED]+=negative_part(newoccupied-oldoccupied);
	reverse[QTYPE_VIRTUAL]+=negative_part(newvirtual-oldvirtual);

	oldvirtual=oldoccupied=newvirtual=newoccupied=0;
	for(int i=0;i<pmx->dimensions;i++)
	{
		if(pmatrix_get_entry(pmx, i, j2)!=0)
		{
			if(pmatrix_entry_type(i,j2)==QTYPE_OCCUPIED)
				oldoccupied++;
			else
				oldvirtual++;

			if(pmatrix_entry_type(i,j1)==QTYPE_OCCUPIED)
				newoccupied++;
			else
				newvirtual++;

			/*
				Going from occupied to virtual
			*/

			if((pmatrix_entry_type(i,j2)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i,j1)==QTYPE_VIRTUAL))
				pmatrix_set_entry(pmx, i, j2, pmatrix_get_new_value(pmx, rngctx, i, j1));

			/*
				Going from virtual to occupied
			*/

			if((pmatrix_entry_type(i,j2)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i,j1)==QTYPE_OCCUPIED))
				pmatrix_set_entry(pmx, i, j2, pmatrix_get_new_value(pmx, rngctx, i, j1));

#ifndef NDEBUG
			if((pmatrix_entry_type(i,j2)==QTYPE_OCCUPIED)&&(pmatrix_entry_type(i,j1)==QTYPE_VIRTUAL)) x2++;
			if((pmatrix_entry_type(i,j2)==QTYPE_VIRTUAL)&&(pmatrix_entry_type(i,j1)==QTYPE_OCCUPIED)) x1++;
#endif
		}
	}

	update[QTYPE_OCCUPIED]+=positive_part(newoccupied-oldoccupied);
	update[QTYPE_VIRTUAL]+=positive_part(newvirtual-oldvirtual);
	reverse[QTYPE_OCCUPIED]+=negative_part(newoccupied-oldoccupied);
	reverse[QTYPE_VIRTUAL]+=negative_part(newvirtual-oldvirtual);

#ifndef NDEBUG
	assert(update[QTYPE_OCCUPIED]==x1);
	assert(update[QTYPE_VIRTUAL]==x2);
#endif
	/*
		Finally we actually swap the columns!
	*/

	for(int i=0;i<pmx->dimensions;i++)
	{
		int tmp=pmatrix_get_entry(pmx, i, j1);

		pmatrix_set_entry(pmx, i, j1, pmatrix_get_entry(pmx, i, j2));
		pmatrix_set_entry(pmx, i, j2, tmp);
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

int pmatrix_get_new_occupied_value(struct pmatrix_t *pmx, gsl_rng *rngctx)
{
	return 1+gsl_rng_uniform_int(rngctx, pmx->nr_occupied);
}

int pmatrix_get_new_virtual_value(struct pmatrix_t *pmx, gsl_rng *rngctx)
{
	return 1+gsl_rng_uniform_int(rngctx, pmx->nr_virtual);
}

int pmatrix_get_new_value(struct pmatrix_t *pmx, gsl_rng *rngctx, int i, int j)
{
	if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
		return pmatrix_get_new_occupied_value(pmx, rngctx);

	return pmatrix_get_new_virtual_value(pmx, rngctx);
}
