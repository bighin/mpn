#include <math.h>
#include <curses.h>
#include <assert.h>
#include <string.h>

#include "auxx.h"

/*
        Routine to seed GSL's random number generator
*/

void seed_rng(gsl_rng *rng)
{
        char *devname="/dev/urandom";
        FILE *dev;
        unsigned long seed;

        if((dev=fopen(devname,"r"))!=NULL)
        {
                fread(&seed,sizeof(unsigned long),1,dev);
                fclose(dev);

                gsl_rng_set(rng,seed);
        }
        else
        {
                printf("Warning: couldn't read from %s to seed the RNG.\n",devname);
        }
}

/*
	Factorial of an integer, using the builtin gamma function
*/

double factorial(int i)
{
	return tgamma(i+1);
}

/*
	Factorial of an integer, using only integer arithmetic
*/

int ifactorial(int n)
{
	int result = 1;

	for (int i = 1; i <= n; ++i)
		result *= i;

	return result;
}

/*
	Power operation using integer arithmetic
*/

int ipow(int base,int exp)
{
	int result=1;

	for (;;)
	{
		if (exp&1)
			result *= base;
		exp>>=1;

		if (!exp)
			break;

		base*=base;
	}

	return result;
}

/*
	Simply prints a matrix, with a newline after each row.
	It seems like this simple function is missing in GSL.
*/

void gsl_matrix_int_print(gsl_matrix_int *m)
{
	for(size_t i=0;i<m->size1;i++)
	{
		for(size_t j=0;j<m->size2;j++)
			printf("%d ", gsl_matrix_int_get(m, i, j));

		printf("\n");
	}
}

/*
	Matrix multiplication: C = AB
*/

void gsl_matrix_int_mul(gsl_matrix_int *A,gsl_matrix_int *B,gsl_matrix_int *C)
{
	assert(A->size2==B->size1);
	assert(A->size1==C->size1);
	assert(B->size2==C->size2);

	for(size_t i=0;i<C->size1;i++)
	{
		for(size_t j=0;j<C->size2;j++)
		{
			int result=0;

			for(size_t k=0;k<A->size2;k++)
				result+=gsl_matrix_int_get(A,i,k)*gsl_matrix_int_get(B,k,j);

			gsl_matrix_int_set(C,i,j,result);
		}
	}
}

/*
        Matrix exponentiation: N = M^i
*/

void gsl_matrix_int_power(gsl_matrix_int *M,gsl_matrix_int *N,int i)
{
	assert(M->size1==M->size2);
	assert(N->size1==N->size2);
	assert(M->size1==N->size1);
	assert(i>=0);

	/*
		Special cases
	*/

	if(i==0)
	{
		gsl_matrix_int_set_identity(N);
		return;
	}
	else if(i==1)
	{
		gsl_matrix_int_memcpy(N,M);
		return;
	}

	/*
		General case
	*/

	gsl_matrix_int *tmp=gsl_matrix_int_alloc(M->size1,M->size2);

	gsl_matrix_int_set_identity(N);

	while(i-->0)
	{
		gsl_matrix_int_mul(N,M,tmp);
		gsl_matrix_int_memcpy(N,tmp);
	}

	gsl_matrix_int_free(tmp);
}

/*
	Returns true if the col1-th and col2-th rows coincide.
*/

bool columns_are_identical(gsl_matrix_int *m, size_t col1, size_t col2)
{
	for(size_t row=0;row<m->size1;row++)
	{
		int x1,x2;

		x1=gsl_matrix_int_get(m, row, col1);
		x2=gsl_matrix_int_get(m, row, col2);

		if(x1!=x2)
			return false;
	}

	return true;
}

/*
	Positive and negative part of an integer number
*/

int positive_part(int x)
{
	return (x>0)?(x):(0);
}

int negative_part(int x)
{
	return (x<0)?(-x):(0);
}

/*
	Returns the n-th character of the string s if it is long enough, otherwise an underscore.
*/

char get_nth_character(char *s,size_t n)
{
	if(n<strlen(s))
		return s[n];

	return '_';
}

/*
	Print a file size using the appropriate units.
*/

void print_file_size(FILE *out,int size)
{
	char *symbols[]={"B","KiB","MiB","GiB","TiB","PiB"};

	for(int c=0;c<6;c++)
	{
		if(size<1024)
		{
			fprintf(out,"%d %s", size, symbols[c]);
			return;
		}

		size/=1024;
	}

	fprintf(out,"%d %s", size, "EiB");
}
