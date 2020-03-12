#include <stdbool.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "cache.h"
#include "auxx.h"
#include "multiplicity.h"
#include "permutations.h"

/*
	The global variables where the actual cache content is kept
*/

uint8_t *amatrix_cache[32];
int amatrix_cache_max_dimensions=-1;
bool amatrix_cache_is_enabled=true;

/*
	The function amatrix_to_index(), given a amatrix_t struct, returns a unique index,
	representative of the arrangement of the zero/non-zero entries.

	Keep in mind that there's nothing about the quantum numbers, since it would be
	so expensive to cache that information. Because of this we cannot cache the weight,
	that depends on the quantum numbers.
*/

int amatrix_to_index(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	int dimensions=amx->pmxs[0]->dimensions;

	int pa[16],pb[16];

	pmatrix_to_permutation(amx->pmxs[0],pa);
	pmatrix_to_permutation(amx->pmxs[1],pb);

	return ifactorial(dimensions)*get_permutation_index(pb,dimensions)+get_permutation_index(pa,dimensions);
}

/*
	Returns the value of the largest index you can get from amatrix_to_index()
*/

int cache_largest_index(int dimensions)
{
	return ifactorial(dimensions)*ifactorial(dimensions);
}

/*
	The multiplicity can take only values that are powers of two, and working
	with matrices of dimension at most 6, the maximum value is 8.

	Here we convert the double value return by amatrix_multiplicity to an int.
*/

int multiplicity_to_int(double multiplicity)
{
	if(gsl_fcmp(multiplicity,1.0,1e-6)==0)
		return 1;

	if(gsl_fcmp(multiplicity,2.0,1e-6)==0)
		return 2;

	if(gsl_fcmp(multiplicity,4.0,1e-6)==0)
		return 4;

	if(gsl_fcmp(multiplicity,8.0,1e-6)==0)
		return 8;

	assert(false);
	return 0;
}

void fill_cache(int dimensions,int expected_connected,int expected_not_connected)
{
	struct amatrix_t *amx=init_amatrix(NULL);

	int connected,not_connected,nr_permutations;

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=dimensions;

	connected=not_connected=0;
	nr_permutations=ifactorial(dimensions);

	for(int i=0;i<nr_permutations;i++)
	{
		for(int j=0;j<nr_permutations;j++)
		{
			gsl_matrix_int *a,*b;

			switch(dimensions)
			{
				case 2:
				a=permutation_to_matrix(permutations2[i], dimensions);
				b=permutation_to_matrix(permutations2[j], dimensions);
				break;

				case 3:
				a=permutation_to_matrix(permutations3[i], dimensions);
				b=permutation_to_matrix(permutations3[j], dimensions);
				break;

				case 4:
				a=permutation_to_matrix(permutations4[i], dimensions);
				b=permutation_to_matrix(permutations4[j], dimensions);
				break;

				case 5:
				a=permutation_to_matrix(permutations5[i], dimensions);
				b=permutation_to_matrix(permutations5[j], dimensions);
				break;

				case 6:
				a=permutation_to_matrix(permutations6[i], dimensions);
				b=permutation_to_matrix(permutations6[j], dimensions);
				break;

				default:
				a=b=NULL;
				assert(false);
			}

			for(int k=0;k<dimensions;k++)
			{
				for(int l=0;l<dimensions;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			bool is_connected=actual_amatrix_check_connectedness(amx);
			double multiplicity=actual_amatrix_multiplicity(amx);

			if(is_connected==true)
				connected++;
			else
				not_connected++;

			cache_set_entry(amatrix_to_index(amx), dimensions, multiplicity_to_int(multiplicity), is_connected);

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	fini_amatrix(amx);

	assert((connected==expected_connected)&&(not_connected==expected_not_connected));
}

bool init_cache(int max_dimensions)
{
	assert((max_dimensions>1)&&(max_dimensions<=6));
	int total_alloced=0;

	for(int dimensions=0;dimensions<32;dimensions++)
		amatrix_cache[dimensions]=NULL;

	for(int dimensions=2;dimensions<=max_dimensions;dimensions++)
	{
		int size_of_current_allocation=cache_largest_index(dimensions);

		amatrix_cache[dimensions]=malloc(size_of_current_allocation);
		total_alloced+=size_of_current_allocation;

		for(int c=0;c<size_of_current_allocation;c++)
			amatrix_cache[dimensions][c]=0;
	}

	printf("Cache size: ");
	print_file_size(stdout,total_alloced);
	printf("\n");

	amatrix_cache_max_dimensions=max_dimensions;

	/*
		Now we just fill the caches.

		TODO: this could be extended at least to order 8.
		However, in that case the multiplicity will no longer fit in three bits!
	*/

	if(max_dimensions>=2) fill_cache(2,3,1);
	if(max_dimensions>=3) fill_cache(3,26,10);
	if(max_dimensions>=4) fill_cache(4,426,150);
	if(max_dimensions>=5) fill_cache(5,11064,3336);
	if(max_dimensions>=6) fill_cache(6,413640,104760);
}

void free_cache(void)
{
	for(int dimensions=0;dimensions<32;dimensions++)
		if(amatrix_cache[dimensions]!=NULL)
			free(amatrix_cache[dimensions]);
}

/*
	Low level functions to get/set a cache entry
*/

uint8_t cache_get_entry(int index, int dimensions)
{
	assert((dimensions>1)&&(dimensions<=amatrix_cache_max_dimensions));
	assert(amatrix_cache[dimensions]!=NULL);
	assert(index<cache_largest_index(dimensions));

	/*
		TODO: if RAM is a problem, each entry could be squeezed in half a byte.

		We only have to modify cache_get_entry(), cache_set_entry() and the
		malloc() call in init_cache().
	*/

	return amatrix_cache[dimensions][index];
}

void cache_set_entry(int index, int dimensions, int multiplicity, bool isconnected)
{
	assert((dimensions>1)&&(dimensions<=amatrix_cache_max_dimensions));
	assert(amatrix_cache[dimensions]!=NULL);
	assert(index<cache_largest_index(dimensions));

	/*
		With diagrams up to order 6, the multiplicity can take only the values
		1, 2, 4, 8. We subtract one so that it can fit in three bits.
	*/

	uint8_t byte=(multiplicity-1)&0x7;

	if(isconnected==true)
		byte|=0x8;

	amatrix_cache[dimensions][index]=byte;
}

/*
	High-level functions for dealing with cache entries given an 'amatrix'
*/

int cached_amatrix_multiplicity(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	int dimensions=amx->pmxs[0]->dimensions;

	uint8_t result=cache_get_entry(amatrix_to_index(amx), dimensions);

	/*
		The lowest three bits contain multiplicity-1
	*/

	return 1+(result&0x7);
}

bool cached_amatrix_check_connectedness(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	int dimensions=amx->pmxs[0]->dimensions;

	uint8_t result=cache_get_entry(amatrix_to_index(amx), dimensions);

	/*
		The connectedness information is in the fourth bit.
	*/

	return (result&0x8)?(true):(false);
}
