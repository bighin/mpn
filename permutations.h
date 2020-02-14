#ifndef __PERMUTATIONS_H__
#define __PERMUTATIONS_H__

#include <gsl/gsl_matrix.h>

int get_permutation_index(int *permutation,int length);
gsl_matrix_int *permutation_to_matrix(const int *permutation,int dimensions);
void matrix_to_permutation(gsl_matrix_int *m,int *permutation);

extern int permutations2[2][2];
extern int permutations3[6][3];
extern int permutations4[24][4];
extern int permutations5[120][5];
extern int permutations6[720][6];

#endif //__PERMUTATIONS_H__
