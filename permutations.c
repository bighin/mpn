#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <assert.h>

#include "permutations.h"
#include "plist.h"
#include "pmatrix.h"

/*
	Given a N-permutation of number from 1 to N, returns its index in the lexicographic
	sorting of all permutations. Note that this is the same order of the permutations listed
	in file permutations.h

	The algorithm is adapted from: http://www.geekviewpoint.com/java/numbers/permutation_index
*/

int get_permutation_index(const int *permutation,int length)
{
	int index=0;

	/*
		Position 0 is paired with factor 0 and so is skipped
	*/

	int position=1;
	int factor=1;

	for(int p=length-1; p>=0; p--)
	{
		int successors=0;
		for(int q=p+1; q<length; q++)
		{
			if (permutation[p] > permutation[q])
				successors++;
		}

		index+=(successors*factor);
		factor*=position;
		position++;
	}

	return index;
}

/*
	Given a permutation in the format of an N-dimensional array of int's, with entries
	going from 1 to N, this function returns the corresponding permutation matrix
*/

gsl_matrix_int *permutation_to_matrix(const int *permutation,int dimensions)
{
	gsl_matrix_int *ret=gsl_matrix_int_alloc(dimensions,dimensions);

	for(int i=0;i<dimensions;i++)
		for(int j=0;j<dimensions;j++)
			gsl_matrix_int_set(ret,i,j,((permutation[i]-1)==j)?(1):(0));

	return ret;
}

void matrix_to_permutation(gsl_matrix_int *m,int *permutation)
{
	assert(m->size1==m->size2);
	int dimensions=m->size1;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(gsl_matrix_int_get(m,i,j)!=0)
			{
				permutation[i]=j+1;
				break;
			}
		}
	}
}

void pmatrix_to_permutation(struct pmatrix_t *m,int *permutation)
{
	int dimensions=m->dimensions;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_get_entry(m,i,j)!=0)
			{
				permutation[i]=j+1;
				break;
			}
		}
	}
}

/*
	The Fisher-Yates algorithm generates a random permutation
*/

void fisher_yates(gsl_rng *rng_ctx, int *array, int length)
{
	for(int c=length-1;c>0;c--)
	{
		int d=gsl_rng_uniform_int(rng_ctx, c);

		int tmp=array[d];
		array[d]=array[c];
		array[c]=tmp;
	}
}
