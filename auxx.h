#ifndef __AUXX_H__
#define __AUXX_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

void seed_rng(gsl_rng *rng);

double factorial(int i);
int ifactorial(int n);
int ipow(int base,int exp);

void gsl_matrix_int_print(gsl_matrix_int *m);
void gsl_matrix_int_mul(gsl_matrix_int *A,gsl_matrix_int *B,gsl_matrix_int *C);
void gsl_matrix_int_power(gsl_matrix_int *M,gsl_matrix_int *N,int i);
bool columns_are_identical(gsl_matrix_int *m, size_t col1, size_t col2);

int positive_part(int x);
int negative_part(int x);

#define MIN(x,y)	(((x)<(y))?(x):(y))
#define MAX(x,y)	(((x)>(y))?(x):(y))

char get_nth_character(char *s,size_t n);

void print_file_size(FILE *out,int size);

char *find_and_replace(const char *src,const char *from,const char *to);

void normalize_distribution(double *dists, int nrstates);
void to_cumulative_distribution(const double *dists, double *cdists, int nrstates);
int cdist_search(const double *cdists, int lo, int hi, double selector);
int cdist_linear_search(const double *cdists, int lo, int hi, double selector);

#endif //__AUXX_H__
