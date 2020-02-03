#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>

#include "amatrix.h"
#include "mpn.h"
#include "mc.h"
#include "multiplicity.h"
#include "config.h"

void run_debug_tests(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	amatrix_print(amx);
	printf("\n");

	assert(amatrix_check_consistency(amx)==true);

	for(int c=0;c<5;c++)
	{
		update_extend(amx, true);

		amatrix_print(amx);
		printf("Weight: %f\n\n",amatrix_weight(amx));

		assert(amatrix_check_consistency(amx)==true);
	}

	for(int c=0;c<1000000;c++)
	{
		update_shuffle(amx, true);
		printf("%f\n",amatrix_weight(amx));
		assert(amatrix_check_consistency(amx)==true);
	}

	for(int c=0;c<1000000;c++)
	{
		update_shuffle(amx, true);

		if(amatrix_is_physical(amx)==true)
		{
			printf("Found a physical matrix!\n");

			amatrix_print(amx);
			printf("\n");

			struct label_t labels[MAX_LABELS];
			int ilabels;

			amatrix_to_python(amx);
			amatrix_to_wolfram(amx);
			gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
			incidence_to_weight(incidence, labels, &ilabels, amx->ectx, amx->unphysicalpenalty, true);
			gsl_matrix_int_free(incidence);

			exit(0);
		}

		assert(amatrix_check_consistency(amx)==true);
	}

	for(int c=0;c<5;c++)
	{
		update_squeeze(amx, true);
		amatrix_print(amx);
		printf("\n");

		assert(amatrix_check_consistency(amx)==true);
	}

	fini_amatrix(amx);
}

void run_debug_tests2(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	amx->pmxs[0]->dimensions=6;

	pmatrix_set_entry(amx->pmxs[0], 0, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 3, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 0, 3));
	pmatrix_set_entry(amx->pmxs[0], 0, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 1, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 2, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 1, 2));
	pmatrix_set_entry(amx->pmxs[0], 1, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 2, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 5, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 2, 5));

	pmatrix_set_entry(amx->pmxs[0], 3, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 4, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 3, 4));
	pmatrix_set_entry(amx->pmxs[0], 3, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 4, 0, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 4, 0));
	pmatrix_set_entry(amx->pmxs[0], 4, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 5, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 1, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 5, 1));
	pmatrix_set_entry(amx->pmxs[0], 5, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 5, 0);

	amx->pmxs[1]->dimensions=6;

	pmatrix_set_entry(amx->pmxs[1], 0, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 1, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 0, 1));
	pmatrix_set_entry(amx->pmxs[1], 0, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 1, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 3, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 1, 3));
	pmatrix_set_entry(amx->pmxs[1], 1, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 2, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 4, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 2, 4));
	pmatrix_set_entry(amx->pmxs[1], 2, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 3, 0, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 3, 0));
	pmatrix_set_entry(amx->pmxs[1], 3, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 4, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 5, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 4, 5));

	pmatrix_set_entry(amx->pmxs[1], 5, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 2, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 5, 2));
	pmatrix_set_entry(amx->pmxs[1], 5, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 5, 0);

	assert(amatrix_check_consistency(amx)==true);

	amatrix_print(amx);
	printf("\n");

	struct label_t labels[MAX_LABELS];
	int ilabels;

	amatrix_to_python(amx);
	amatrix_to_wolfram(amx);
	gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
	incidence_to_weight(incidence, labels, &ilabels, amx->ectx, amx->unphysicalpenalty, true);
	gsl_matrix_int_free(incidence);
}

/*
	Here we debug the connectedness test
*/

gsl_matrix_int *permutation_to_matrix(int *permutation,int dimensions)
{
	gsl_matrix_int *ret=gsl_matrix_int_alloc(dimensions,dimensions);

	for(int i=0;i<dimensions;i++)
		for(int j=0;j<dimensions;j++)
			gsl_matrix_int_set(ret,i,j,((permutation[i]-1)==j)?(1):(0));

	return ret;
}

#include "permutations.h"

/*
	Here we check the number of connected diagrams at order 3,4,5,6 and compare
	them with values obtained from Wolfram Mathematica, please see the notebook
	Connectedness.nb
*/

void run_debug_tests3(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	int connected,not_connected;

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=3;
	connected=not_connected=0;

	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations3[i],3);
			b=permutation_to_matrix(permutations3[j],3);

			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==26)&&(not_connected==10));

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=4;
	connected=not_connected=0;

	for(int i=0;i<24;i++)
	{
		for(int j=0;j<24;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations4[i],4);
			b=permutation_to_matrix(permutations4[j],4);

			for(int k=0;k<4;k++)
			{
				for(int l=0;l<4;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==426)&&(not_connected==150));

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=5;
	connected=not_connected=0;

	for(int i=0;i<120;i++)
	{
		for(int j=0;j<120;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations5[i],5);
			b=permutation_to_matrix(permutations5[j],5);

			for(int k=0;k<5;k++)
			{
				for(int l=0;l<5;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==11064)&&(not_connected==3336));

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=6;
	connected=not_connected=0;

	for(int i=0;i<720;i++)
	{
		for(int j=0;j<720;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations6[i],6);
			b=permutation_to_matrix(permutations6[j],6);

			for(int k=0;k<6;k++)
			{
				for(int l=0;l<6;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==413640)&&(not_connected==104760));
}

void run_debug_tests3b(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	int connected,not_connected;

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=3;
	connected=not_connected=0;

	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations3[i],3);
			b=permutation_to_matrix(permutations3[j],3);

			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=4;
	connected=not_connected=0;

	for(int i=0;i<24;i++)
	{
		for(int j=0;j<24;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations4[i],4);
			b=permutation_to_matrix(permutations4[j],4);

			for(int k=0;k<4;k++)
			{
				for(int l=0;l<4;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=5;
	connected=not_connected=0;

	for(int i=0;i<120;i++)
	{
		for(int j=0;j<120;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations5[i],5);
			b=permutation_to_matrix(permutations5[j],5);

			for(int k=0;k<5;k++)
			{
				for(int l=0;l<5;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=6;
	connected=not_connected=0;

	for(int i=0;i<720;i++)
	{
		for(int j=0;j<720;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations6[i],6);
			b=permutation_to_matrix(permutations6[j],6);

			for(int k=0;k<6;k++)
			{
				for(int l=0;l<6;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
}

void run_debug_tests4(void)
{
	/*
		'Manual' MP2 calculation
	*/

	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	double result=0.0;

	pmatrix_set_entry(amx->pmxs[0],0,0,0);
	pmatrix_set_entry(amx->pmxs[0],0,1,0);
	pmatrix_set_entry(amx->pmxs[0],1,0,0);
	pmatrix_set_entry(amx->pmxs[0],1,1,0);
	pmatrix_set_entry(amx->pmxs[1],0,0,0);
	pmatrix_set_entry(amx->pmxs[1],0,1,0);
	pmatrix_set_entry(amx->pmxs[1],1,0,0);
	pmatrix_set_entry(amx->pmxs[1],1,1,0);

	for(int i=1;i<=amx->nr_virtual;i++)
	{
		for(int j=1;j<=amx->nr_virtual;j++)
		{
			for(int a=1;a<=amx->nr_occupied;a++)
			{
				for(int b=1;b<=amx->nr_occupied;b++)
				{
					pmatrix_set_entry(amx->pmxs[0],0,1,a);
					pmatrix_set_entry(amx->pmxs[0],1,0,i);
					pmatrix_set_entry(amx->pmxs[1],0,1,b);
					pmatrix_set_entry(amx->pmxs[1],1,0,j);

					result+=amatrix_weight(amx);

					bool verbose=false;

					if(verbose==true)
					{

						printf("%d %d %d %d || %f\n", pmatrix_get_entry(amx->pmxs[0], 0, 1),
						       pmatrix_get_entry(amx->pmxs[0], 1, 0),
						       pmatrix_get_entry(amx->pmxs[1], 0, 1),
						       pmatrix_get_entry(amx->pmxs[1], 1, 0),
						       amatrix_weight(amx));
					}
				}
			}
		}
	}

	printf("%f\n",result);
}

void run_debug_tests5(void)
{
	/*
		'Manual' MP2 calculation, with stochastic sampling of quantum numbers
	*/

	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	pmatrix_set_entry(amx->pmxs[0],0,0,0);
	pmatrix_set_entry(amx->pmxs[0],0,1,0);
	pmatrix_set_entry(amx->pmxs[0],1,0,0);
	pmatrix_set_entry(amx->pmxs[0],1,1,0);
	pmatrix_set_entry(amx->pmxs[1],0,0,0);
	pmatrix_set_entry(amx->pmxs[1],0,1,0);
	pmatrix_set_entry(amx->pmxs[1],1,0,0);
	pmatrix_set_entry(amx->pmxs[1],1,1,0);

	double result=0.0;
	int iterations=50000000;

	for(int c=0;c<iterations;c++)
	{
		int i=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_virtual);
		int j=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_virtual);
		int a=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_occupied);
		int b=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_occupied);

		pmatrix_set_entry(amx->pmxs[0],0,1,a);
		pmatrix_set_entry(amx->pmxs[0],1,0,i);
		pmatrix_set_entry(amx->pmxs[1],0,1,b);
		pmatrix_set_entry(amx->pmxs[1],1,0,j);

		result+=amatrix_weight(amx);
	}

	printf("%f\n",result/iterations*pow(amx->nr_virtual*amx->nr_occupied,2.0f));
}

void run_debug_tests6(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	double result=0.0;

	for(int a=0;a<amx->nr_occupied;a++)
	{
		for(int b=0;b<amx->nr_occupied;b++)
		{
			if(a==b)
			{
				printf("[%f]\n",get_occupied_energy(amx->ectx, a));
				result+=get_occupied_energy(amx->ectx, a);
			}

			//result+=0.5*get_eri(amx->ectx,a,b,a,b);
		}
	}

	printf("%f\n",result);
}

void run_debug_tests7(void)
{
	/*
		'Manual' Hartree-Fock calculation
	*/

	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	double result=0.0;

	pmatrix_set_entry(amx->pmxs[0],0,0,0);
	pmatrix_set_entry(amx->pmxs[1],0,0,0);

	for(int a=1;a<=amx->nr_occupied;a++)
	{
		for(int b=1;b<=amx->nr_occupied;b++)
		{
			pmatrix_set_entry(amx->pmxs[0], 0, 0, a);
			pmatrix_set_entry(amx->pmxs[1], 0, 0, b);

			result+=amatrix_weight(amx);

			bool verbose=false;

			if(verbose==true)
			{

				printf("%d %d %d %d || %f\n", pmatrix_get_entry(amx->pmxs[0], 0, 1),
				       pmatrix_get_entry(amx->pmxs[0], 1, 0), pmatrix_get_entry(amx->pmxs[1], 0, 1),
				       pmatrix_get_entry(amx->pmxs[1], 1, 0), amatrix_weight(amx));
			}
		}
	}

	printf("%f\n",result);
}

void run_debug_tests3c(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	double connected,not_connected;

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=3;
	connected=not_connected=0;

	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations3[i],3);
			b=permutation_to_matrix(permutations3[j],3);

			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected+=1.0f/amatrix_multiplicity(amx);
				else
					not_connected+=1.0f/amatrix_multiplicity(amx);
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%f %f\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=4;
	connected=not_connected=0;

	for(int i=0;i<24;i++)
	{
		for(int j=0;j<24;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations4[i],4);
			b=permutation_to_matrix(permutations4[j],4);

			for(int k=0;k<4;k++)
			{
				for(int l=0;l<4;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected+=1.0f/amatrix_multiplicity(amx);
				else
					not_connected+=1.0f/amatrix_multiplicity(amx);
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%f %f\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=5;
	connected=not_connected=0;

	for(int i=0;i<120;i++)
	{
		for(int j=0;j<120;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations5[i],5);
			b=permutation_to_matrix(permutations5[j],5);

			for(int k=0;k<5;k++)
			{
				for(int l=0;l<5;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected+=1.0f/amatrix_multiplicity(amx);
				else
					not_connected+=1.0f/amatrix_multiplicity(amx);
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%f %f\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=6;
	connected=not_connected=0;

	for(int i=0;i<720;i++)
	{
		for(int j=0;j<720;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations6[i],6);
			b=permutation_to_matrix(permutations6[j],6);

			for(int k=0;k<6;k++)
			{
				for(int l=0;l<6;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected+=1.0f/amatrix_multiplicity(amx);
				else
					not_connected+=1.0f/amatrix_multiplicity(amx);
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%f %f\n",connected,not_connected);
}

int old_main(void)
{
	//run_debug_tests();
	//run_debug_tests3();
	//run_debug_tests3b();
	//run_debug_tests4();
	//run_debug_tests5();
	//run_debug_tests6();
	//run_debug_tests7();
	//run_debug_tests3c();
	//do_diagmc("/Users/zakk/Desktop/MPn/psi4/H2O.dat","/Users/zakk/Desktop/MPn/test.dat",1000000,100.0f,false);

	return 0;
}

void usage(char *argv0)
{
	printf("Usage: %s <inifile> [<otherinifiles> ...]\n",argv0);

	exit(0);
}

int main(int argc,char *argv[])
{
	if(argc<2)
		usage(argv[0]);

	bool first=true;

	for(int c=1;c<argc;c++)
	{
		struct configuration_t config;

		if(first==true)
			fprintf(stderr,"Diagrammatic Monte Carlo for Møller-Plesset theory.\n");

		load_config_defaults(&config);

		if(load_configuration(argv[c],&config)==false)
			continue;

		do_diagmc(&config);
		first=false;

		if((c+1)!=argc)
			printf("\n");
	}
}
