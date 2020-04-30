#ifndef __TWINS_H__
#define __TWINS_H__

#include <stdbool.h>

#include "amatrix.h"

#define MAX_PAIRS	(32)

struct permutation_collection_t
{
	int pairs[MAX_PAIRS][2];
	int ipairs;

	int index;
	bool permutations_have_ended;
};

struct permutation_collection_t *init_permutation_collection(void);
void fini_permutation_collection(struct permutation_collection_t *pct);

bool go_to_next_permutation(struct permutation_collection_t *pct,struct amatrix_weight_t *awt);

struct permutation_collection_t *identify_twins(struct amatrix_t *amx,struct amatrix_weight_t *awt,bool *is_representative,int *combinatorial);

#endif //__TWINS_H__
