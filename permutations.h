#ifndef __PERMUTATIONS_H__
#define __PERMUTATIONS_H__

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

#include "pmatrix.h"

int get_permutation_index(const int *permutation,int length);

gsl_matrix_int *permutation_to_matrix(const int *permutation,int dimensions);
void matrix_to_permutation(gsl_matrix_int *m,int *permutation);

void pmatrix_to_permutation(struct pmatrix_t *m,int *permutation);

extern int permutations2[2][2];
extern int permutations3[6][3];
extern int permutations4[24][4];
extern int permutations5[120][5];
extern int permutations6[720][6];
extern int permutations7[5040][7];
extern int permutations8[40320][8];
extern int permutations9[362880][9];
extern int permutations10[3628800][10];

void init_permutation_tables(int max_dimensions);
int get_permutation(int dimensions,int pindex,int element);

void fisher_yates(gsl_rng *rng_ctx, int *array, int length);

#endif //__PERMUTATIONS_H__
