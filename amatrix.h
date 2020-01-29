#ifndef __CMATRIX_H__
#define __CMATRIX_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "pmatrix.h"
#include "reader.h"

/*
	The 'amatrix' struct: a matrix following a certain set of rules,
	in which every non-zero entry is associated to some quantum numbers.
*/

struct amatrix_t
{
	int nr_occupied,nr_virtual;

	/*
		The two permutation matrices composing the current matrix
	*/

	struct pmatrix_t *pmxs[2];

	/*
		A GSL RNG context
	*/

	gsl_rng *rng_ctx;

	/*
		A context defined in reader.h, containing the orbital energies
		and the H tensor.
	*/

	struct energies_ctx_t *ectx;

	/*
		Various parameters
	*/

	double bias,unphysical_penalty;
	int max_order;
};

struct amatrix_t *init_amatrix(char *energies_dot_dat);
void fini_amatrix(struct amatrix_t *amx);

int amatrix_get_entry(struct amatrix_t *amx, int i, int j);

void amatrix_print(struct amatrix_t *amx);

struct amatrix_stack_t
{
	int dimensions[2];
	int values[2][IMATRIX_MAX_DIMENSIONS][IMATRIX_MAX_DIMENSIONS];
};

void amatrix_push(struct amatrix_t *amx, struct amatrix_stack_t *stack);
void amatrix_pop(struct amatrix_t *amx, struct amatrix_stack_t *stack);

bool amatrix_check_consistency(struct amatrix_t *amx);
bool amatrix_is_physical(struct amatrix_t *amx);

void amatrix_to_python(struct amatrix_t *amx);
void amatrix_to_wolfram(struct amatrix_t *amx);

bool amatrix_check_connectedness(struct amatrix_t *amx);
double amatrix_weight(struct amatrix_t *amx);

#endif //__CMATRIX_H__
