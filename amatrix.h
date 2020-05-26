#ifndef __AMATRIX_H__
#define __AMATRIX_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "pmatrix.h"
#include "config.h"
#include "limits.h"

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
		Various parameters that are read from the configuration file are store here
	*/

	struct configuration_t *config;

	/*
		The cached weight
	*/

	double cached_weight;
	bool cached_weight_is_valid;
};

struct amatrix_t *init_amatrix(struct configuration_t *config);
void fini_amatrix(struct amatrix_t *amx,bool free_ectx);

int amatrix_get_entry(struct amatrix_t *amx, int i, int j);

struct amatrix_backup_t
{
	int dimensions[2];
	int values[2][PMATRIX_MAX_DIMENSIONS][PMATRIX_MAX_DIMENSIONS];

	double cached_result;
	bool cached_result_is_valid;
};

void amatrix_save(struct amatrix_t *amx, struct amatrix_backup_t *backup);
void amatrix_restore(struct amatrix_t *amx, struct amatrix_backup_t *backup);

bool amatrix_check_consistency(struct amatrix_t *amx);
bool amatrix_is_physical(struct amatrix_t *amx);

void amatrix_print(struct amatrix_t *amx);
void amatrix_to_python(struct amatrix_t *amx);
void amatrix_to_wolfram(struct amatrix_t *amx);

bool actual_amatrix_check_connectedness(struct amatrix_t *amx);
bool amatrix_check_connectedness(struct amatrix_t *amx);

#endif //__AMATRIX_H__
