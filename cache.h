#ifndef __CACHE_H__
#define __CACHE_H__

#include <stdbool.h>
#include <stdint.h>

#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"

extern int amatrix_cache_max_dimensions;
extern bool amatrix_cache_is_enabled;

int amatrix_to_index(struct amatrix_t *amx);
int cache_largest_index(int dimensions);

gsl_matrix_int *permutation_to_matrix(const int *permutation,int dimensions);
void matrix_to_permutation(gsl_matrix_int *m,int *permutation);

bool init_cache(int max_dimensions);

uint8_t cache_get_entry(int index, int dimensions);
void cache_set_entry(int index, int dimensions, int multiplicity, bool isconnected);

int cached_amatrix_multiplicity(struct amatrix_t *amx);
bool cached_amatrix_check_connectedness(struct amatrix_t *amx);

#endif //__CACHE_H__
