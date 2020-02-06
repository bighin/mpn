#include <stdio.h>
#include <assert.h>

#include "multiplicity.h"
#include "auxx.h"

/*
	I follow the algorithmic determination of multiplicity by Quoc
*/

bool create_one_by_swapping_cols(gsl_matrix_int *adjacency, size_t indexi, size_t indexj)
{
	for(size_t j=indexj;j<adjacency->size2;j++)
	{
		/*
			We need to find a '1' in this row
		*/

		if(gsl_matrix_int_get(adjacency,indexi,j)==1)
		{
			if(indexj==j)
				return true;

			gsl_matrix_int_swap_columns(adjacency, indexj, j);
			return true;
		}
	}

	return false;
}

bool create_one_by_swapping_rows(gsl_matrix_int *adjacency, size_t indexi, size_t indexj)
{
	for(size_t i=indexi;i<adjacency->size1;i++)
	{
		/*
			We need to find a '1' in this column
		*/

		if(gsl_matrix_int_get(adjacency,i,indexj)==1)
		{
			if(indexi==i)
				return true;

			gsl_matrix_int_swap_rows(adjacency, indexi, i);
			return true;
		}
	}

	return false;
}

size_t build_block(gsl_matrix_int *adjacency)
{
	size_t indexi,indexj;

	indexi=indexj=0;

	if(create_one_by_swapping_rows(adjacency, indexi, indexj)==false)
		return 0;

	indexi++;

	while((indexi<adjacency->size1)&&(indexj<adjacency->size2))
	{
		if(create_one_by_swapping_rows(adjacency, indexi, indexj)==false)
			break;

		indexj++;

		if(create_one_by_swapping_cols(adjacency, indexi, indexj)==false)
			break;

		indexi++;
	}

	return MIN(indexi+1,indexj+1);
}

gsl_matrix_int *downsize_matrix(gsl_matrix_int *adjacency,gsl_matrix_int *targets)
{
	size_t targeti,targetj,newsize1,newsize2;

	targeti=targetj=0;
	newsize1=newsize2=0;

	for(size_t i=0;i<adjacency->size1;i++)
	{
		bool row_has_elements=false;

		for(size_t j=0;j<adjacency->size2;j++)
		{
			if(gsl_matrix_int_get(targets,i,j)==0)
			{
				gsl_matrix_int_set(adjacency,targeti,targetj,gsl_matrix_int_get(adjacency,i,j));

				targetj++;
				row_has_elements=true;
			}
		}

		if(row_has_elements==true)
		{
			targeti++;
			newsize1++;

			assert((newsize2==0)||(newsize2==targetj));

			if(newsize2==0)
				newsize2=targetj;
		}

		targetj=0;
	}

	/*
		The idiomatic way of resizing a GSL matrix is by alloc'ing a new one and copying,
		see: https://lists.gnu.org/archive/html/help-gsl/2004-09/msg00029.html
	*/

	if((newsize1>0)&&(newsize2>0))
	{
		gsl_matrix_int *newmatrix=gsl_matrix_int_alloc(newsize1, newsize2);

		for(size_t i=0;i<newmatrix->size1;i++)
			for(size_t j=0;j<newmatrix->size2;j++)
				gsl_matrix_int_set(newmatrix, i, j, gsl_matrix_int_get(adjacency, i, j));

		gsl_matrix_int_free(adjacency);

		return newmatrix;
	}

	return NULL;
}

int adjacency_matrix_multiplicity(gsl_matrix_int *adjacency)
{
	/*
		At first we eliminate all the rows and columns containing a '2' entry
	*/

	gsl_matrix_int *targets=gsl_matrix_int_alloc(adjacency->size1, adjacency->size2);
	gsl_matrix_int_set_zero(targets);

	for(size_t i=0;i<adjacency->size1;i++)
	{
		for(size_t j=0;j<adjacency->size2;j++)
		{
			if(gsl_matrix_int_get(adjacency,i,j)==2)
			{
				/*
					The whole row and column is marked for deletion
				*/

				for(size_t k=0;k<adjacency->size1;k++)
					gsl_matrix_int_set(targets,k,j,1);

				for(size_t k=0;k<adjacency->size2;k++)
					gsl_matrix_int_set(targets,i,k,1);
			}
		}
	}

	adjacency=downsize_matrix(adjacency,targets);

	/*
		Then we repeatedly create a block, and we delete it, counting the total
		number of blocks as we do so.
	*/

	int nrblocks=0;

	while(adjacency!=NULL)
	{
		size_t blockdimensions=build_block(adjacency);
		nrblocks++;

		gsl_matrix_int_set_zero(targets);

		for(size_t i=0;i<adjacency->size1;i++)
		{
			for(size_t j=0;j<adjacency->size2;j++)
			{
				if((i<blockdimensions)||(j<blockdimensions))
					gsl_matrix_int_set(targets, i, j, 1);

				/*
					A consistency check: outside the diagonal block,
					one must have only zeroes.
				*/

#ifndef NDEBUG
				int cnt=0;

				if(i<blockdimensions)
					cnt++;

				if(j<blockdimensions)
					cnt++;

				if(cnt==1)
					assert(gsl_matrix_int_get(adjacency, i, j)==0);
#endif
			}
		}

		adjacency=downsize_matrix(adjacency,targets);
	}

	gsl_matrix_int_free(targets);

	/*
		Note that 2**nrblocks = 1<<nrblocks
	*/

	return 1<<nrblocks;
}

double amatrix_multiplicity(struct amatrix_t *amx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(dimensions>=1);
	if(dimensions==1)
		return 2;

	/*
		We create and populate the adjacency matrix as a gsl_matrix_int
	*/

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

	/*
		Note that adjacency_matrix_multiplicity() 'consumes' the adjacency
		matrix and there is no need to free it here.
	*/

	return adjacency_matrix_multiplicity(adjacency);
}
