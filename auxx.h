#ifndef __AUXX_H__
#define __AUXX_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

void seed_rng(gsl_rng *rng);

double factorial(int i);

char get_nth_character(char *s,size_t n);

void gsl_matrix_int_print(gsl_matrix_int *m);
void gsl_matrix_int_mul(gsl_matrix_int *A,gsl_matrix_int *B,gsl_matrix_int *C);
void gsl_matrix_int_power(gsl_matrix_int *M,gsl_matrix_int *N,int i);
bool columns_are_identical(gsl_matrix_int *m, size_t col1, size_t col2);

#endif //__AUXX_H__
