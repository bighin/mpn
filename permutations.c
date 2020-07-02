#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "permutations.h"
#include "pmatrix.h"
#include "auxx.h"

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

int permutations2[2][2];
int permutations3[6][3];
int permutations4[24][4];
int permutations5[120][5];
int permutations6[720][6];
int permutations7[5040][7];
int permutations8[40320][8];
int permutations9[362880][9];
int permutations10[3628800][10];

void init_permutation_table(int dimensions)
{
	int index,numbers[MAX_ORDER];

	for(int c=0;c<MAX_ORDER;c++)
		numbers[c]=1;

	index=0;
	while(numbers[dimensions]==1)
	{
		bool is_permutation=true;

		for(int c=0;c<dimensions;c++)
			for(int d=c+1;d<dimensions;d++)
				if(numbers[c]==numbers[d])
					is_permutation=false;

		if(is_permutation)
		{
			for(int c=0;c<dimensions;c++)
			{
				switch(dimensions)
				{
					case 2:
					permutations2[index][c]=numbers[dimensions-1-c];
					break;

					case 3:
					permutations3[index][c]=numbers[dimensions-1-c];
					break;

					case 4:
					permutations4[index][c]=numbers[dimensions-1-c];
					break;

					case 5:
					permutations5[index][c]=numbers[dimensions-1-c];
					break;

					case 6:
					permutations6[index][c]=numbers[dimensions-1-c];
					break;

					case 7:
					permutations7[index][c]=numbers[dimensions-1-c];
					break;

					case 8:
					permutations8[index][c]=numbers[dimensions-1-c];
					break;

					case 9:
					permutations9[index][c]=numbers[dimensions-1-c];
					break;

					case 10:
					permutations10[index][c]=numbers[dimensions-1-c];
					break;

					default:
					break;
				}
			}

			index++;
		}

		for(int c=0;c<MAX_ORDER;c++)
		{
			numbers[c]++;

			if(numbers[c]==(dimensions+1))
				numbers[c]=1;
			else
				break;
		}
	}

	assert(index==ifactorial(dimensions));
}

void init_permutation_tables(int max_dimensions)
{
	for(int c=2;c<=max_dimensions;c++)
		init_permutation_table(c);
}

int get_permutation(int dimensions,int pindex,int element)
{
	assert((dimensions>=2)&&(dimensions<=10));
	assert((pindex>=0)&&(pindex<ifactorial(dimensions)));
	assert((element>=0)&&(element<dimensions));

	switch(dimensions)
	{
		case 2:
		return permutations2[pindex][element];

		case 3:
		return permutations3[pindex][element];

		case 4:
		return permutations4[pindex][element];

		case 5:
		return permutations5[pindex][element];

		case 6:
		return permutations6[pindex][element];

		case 7:
		return permutations7[pindex][element];

		case 8:
		return permutations8[pindex][element];

		case 9:
		return permutations9[pindex][element];

		case 10:
		return permutations10[pindex][element];

		default:
		break;
	}

	assert(false);
	return 0;
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
