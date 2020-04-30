#include <string.h>
#include <assert.h>

#include "twins.h"

int assign_label(struct amatrix_weight_t *awt, int pmatrix, int i, int j)
{
	for(int c=0;c<awt->ilabels;c++)
		if((awt->labels[c].pmatrix==pmatrix)&&(awt->labels[c].i==i)&&(awt->labels[c].j==j))
			return c;

	assert(false);
	return -1;
}

void swap_label_values(struct amatrix_weight_t *awt, int index1, int index2)
{
	int tmp;

	tmp=awt->labels[index1].value;
	awt->labels[index1].value=awt->labels[index2].value;
	awt->labels[index2].value=tmp;
}

/*
	Convert an integer to the corresponding Gray code
*/

int integer_to_gray_code(int x)
{
	return x^(x>>1);
}

/*
	Returns the least significant set bit in an integer
*/

int LSB(int x)
{
	return __builtin_ctz(x);
}

/*
	Returns the number of set bit in an integer
*/

int BITCOUNT(int x)
{
	return __builtin_popcount(x);
}

struct permutation_collection_t *init_permutation_collection(void)
{
	struct permutation_collection_t *ret=malloc(sizeof(struct permutation_collection_t));

	ret->ipairs=0;
	ret->index=0;
	ret->permutations_have_ended=false;

	return ret;
}

void fini_permutation_collection(struct permutation_collection_t *pct)
{
	if(pct)
		free(pct);
}

void add_pair(struct permutation_collection_t *pct,int label1,int label2)
{
	assert(pct->ipairs<MAX_PAIRS);

	pct->pairs[pct->ipairs][0]=label1;
	pct->pairs[pct->ipairs][1]=label2;
	pct->ipairs++;
}

bool go_to_next_permutation(struct permutation_collection_t *pct,struct amatrix_weight_t *awt)
{
	if(pct->permutations_have_ended==true)
		return false;

	int current_gray_code,previous_gray_code;

	previous_gray_code=integer_to_gray_code(pct->index);
	pct->index++;
	current_gray_code=integer_to_gray_code(pct->index);

	int difference_mask=previous_gray_code^current_gray_code;

	assert(BITCOUNT(difference_mask)==1);

	int pair_to_swap=LSB(difference_mask);

	/*
		Special care must be taken for the last permutation.

		If according to the the Gray codes we should now swap the (n+1)-th pair,
		this means that we have run through all possible permutations.

		Suppose we have only 3 permutations, then the special case we will triggered
		by setting the bit in 1000 in the following sequence:

		0000
		0001
		0011
		0010
		0110
		0111
		0101
		0100
		1100
		X

		What we do is swapping back the n-th permutation, going back to the original,
		non-permutated state that has never been reached through this function.

		In doing so, we also set the permutations_have_ended flag.
	*/

	if(pair_to_swap>=pct->ipairs)
	{
		assert(pair_to_swap==pct->ipairs);

		pair_to_swap--;

		if(pair_to_swap>=0)
			swap_label_values(awt, pct->pairs[pair_to_swap][0], pct->pairs[pair_to_swap][1]);

		pct->permutations_have_ended=true;

		return true;
	}

	swap_label_values(awt, pct->pairs[pair_to_swap][0], pct->pairs[pair_to_swap][1]);

	return true;
}

struct permutation_collection_t *identify_twins(struct amatrix_t *amx,struct amatrix_weight_t *awt,bool *is_representative,int *combinatorial)
{
	int selected[PMATRIX_MAX_DIMENSIONS][PMATRIX_MAX_DIMENSIONS][2]={0};

	struct
	{
		int i,j,pmatrix;
	}
	entries_above[2], entries_below[2];

	int ientries_above, ientries_below;

	struct permutation_collection_t *pct=init_permutation_collection();

	/*
		We look for lines (i.e. entries in the adjacency matrix) satisfying
		the following conditions:

		1) They cannot be part of the same '2' entry
		2) They must be above the diagonal (j>i)
		3) They must be on the same row
	*/

	int dimensions=amx->pmxs[0]->dimensions;

	for(int i=0;i<dimensions;i++)
	{
		ientries_above=0;

		for(int j=i+1;j<dimensions;j++)
		{
			assert(j>i);

			if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
			{
				entries_above[ientries_above].i=i;
				entries_above[ientries_above].j=j;
				entries_above[ientries_above].pmatrix=0;
				ientries_above++;
			}
			else if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
			{
				entries_above[ientries_above].i=i;
				entries_above[ientries_above].j=j;
				entries_above[ientries_above].pmatrix=1;
				ientries_above++;
			}

			if(ientries_above==2)
			{
				int i0,i1,j0,j1,pmatrix0,pmatrix1;

				i0=entries_above[0].i;
				i1=entries_above[1].i;

				j0=entries_above[0].j;
				j1=entries_above[1].j;

				pmatrix0=entries_above[0].pmatrix;
				pmatrix1=entries_above[1].pmatrix;

				if((selected[i0][j0][pmatrix0]==0)&&(selected[i1][j1][pmatrix1]==0))
				{
					selected[i0][j0][pmatrix0]=1;
					selected[i1][j1][pmatrix1]=1;

					int label0=assign_label(awt,pmatrix0,i0,j0);
					int label1=assign_label(awt,pmatrix1,i1,j1);

					add_pair(pct,label0,label1);
				}

				break;
			}
		}
	}

	for(int j=0;j<dimensions;j++)
	{
		ientries_above=0;

		for(int i=0;i<j;i++)
		{
			assert(j>i);

			if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
			{
				entries_above[ientries_above].i=i;
				entries_above[ientries_above].j=j;
				entries_above[ientries_above].pmatrix=0;
				ientries_above++;
			}
			else if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
			{
				entries_above[ientries_above].i=i;
				entries_above[ientries_above].j=j;
				entries_above[ientries_above].pmatrix=1;
				ientries_above++;
			}

			if(ientries_above==2)
			{
				int i0,i1,j0,j1,pmatrix0,pmatrix1;

				i0=entries_above[0].i;
				i1=entries_above[1].i;

				j0=entries_above[0].j;
				j1=entries_above[1].j;

				pmatrix0=entries_above[0].pmatrix;
				pmatrix1=entries_above[1].pmatrix;

				if((selected[i0][j0][pmatrix0]==0)&&(selected[i1][j1][pmatrix1]==0))
				{
					selected[i0][j0][pmatrix0]=1;
					selected[i1][j1][pmatrix1]=1;

					int label0=assign_label(awt,pmatrix0,i0,j0);
					int label1=assign_label(awt,pmatrix1,i1,j1);

					add_pair(pct,label0,label1);
				}

				break;
			}
		}
	}

	for(int i=0;i<dimensions;i++)
	{
		ientries_below=0;

		for(int j=0;j<i;j++)
		{
			assert(i>j);

			if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
			{
				entries_below[ientries_below].i=i;
				entries_below[ientries_below].j=j;
				entries_below[ientries_below].pmatrix=0;
				ientries_below++;
			}
			else if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
			{
				entries_below[ientries_below].i=i;
				entries_below[ientries_below].j=j;
				entries_below[ientries_below].pmatrix=1;
				ientries_below++;
			}

			if(ientries_below==2)
			{
				int i0,i1,j0,j1,pmatrix0,pmatrix1;

				i0=entries_below[0].i;
				i1=entries_below[1].i;

				j0=entries_below[0].j;
				j1=entries_below[1].j;

				pmatrix0=entries_below[0].pmatrix;
				pmatrix1=entries_below[1].pmatrix;

				if((selected[i0][j0][pmatrix0]==0)&&(selected[i1][j1][pmatrix1]==0))
				{
					selected[i0][j0][pmatrix0]=1;
					selected[i1][j1][pmatrix1]=1;

					int label0=assign_label(awt,pmatrix0,i0,j0);
					int label1=assign_label(awt,pmatrix1,i1,j1);

					add_pair(pct,label0,label1);
				}

				break;
			}
		}
	}

	for(int j=0;j<dimensions;j++)
	{
		ientries_below=0;

		for(int i=j+1;i<dimensions;i++)
		{
			assert(i>j);

			if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
			{
				entries_below[ientries_below].i=i;
				entries_below[ientries_below].j=j;
				entries_below[ientries_below].pmatrix=0;
				ientries_below++;
			}
			else if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
			{
				entries_below[ientries_below].i=i;
				entries_below[ientries_below].j=j;
				entries_below[ientries_below].pmatrix=1;
				ientries_below++;
			}

			if(ientries_below==2)
			{
				int i0,i1,j0,j1,pmatrix0,pmatrix1;

				i0=entries_below[0].i;
				i1=entries_below[1].i;

				j0=entries_below[0].j;
				j1=entries_below[1].j;

				pmatrix0=entries_below[0].pmatrix;
				pmatrix1=entries_below[1].pmatrix;

				if((selected[i0][j0][pmatrix0]==0)&&(selected[i1][j1][pmatrix1]==0))
				{
					selected[i0][j0][pmatrix0]=1;
					selected[i1][j1][pmatrix1]=1;

					int label0=assign_label(awt,pmatrix0,i0,j0);
					int label1=assign_label(awt,pmatrix1,i1,j1);

					add_pair(pct,label0,label1);
				}

				break;
			}
		}
	}

	*is_representative=true;
	*combinatorial=1;

	for(int c=0;c<pct->ipairs;c++)
	{
		int label0=pct->pairs[c][0];
		int label1=pct->pairs[c][1];

		if(awt->labels[label0].value>awt->labels[label1].value)
			*is_representative=false;

		if(awt->labels[label0].value==awt->labels[label1].value)
			*combinatorial=*combinatorial*2;
	}

	return pct;
}
