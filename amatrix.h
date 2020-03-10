#ifndef __AMATRIX_H__
#define __AMATRIX_H__

#include <stdbool.h>
#include <gsl/gsl_rng.h>

#include "pmatrix.h"
#include "loaderis.h"
#include "config.h"
#include "mpn.h"

/*
	This structure saves the result of the evaluation of the weight of a matrix,
	as well as a couple of other quantities (number of hole states and number of loops).
*/

struct amatrix_weight_t
{
	/*
		The actual weight.
	*/

	double weight;

	/*
		This information can be used to reconstruct the weight, given the labels.
	*/

	int l,h;

	int nr_denominators;
	struct
	{
		int labels[32];
		int qtypes[32];
		int ilabels;
	}
	denominators[16];

	int nr_numerators;
	struct
	{
		int labels[4];
	}
	numerators[16];

	double inversefactor,unphysical_penalty;

#warning Creare un file limits, qui ci va MAXLABELS

	struct label_t labels[32];
	int ilabels;
};

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

	struct amatrix_weight_t cached_weight;
	bool cached_weight_is_valid;
};

struct amatrix_t *init_amatrix(struct configuration_t *config);
void fini_amatrix(struct amatrix_t *amx);

int amatrix_get_entry(struct amatrix_t *amx, int i, int j);

void amatrix_print(struct amatrix_t *amx);

struct amatrix_backup_t
{
	int dimensions[2];
	int values[2][PMATRIX_MAX_DIMENSIONS][PMATRIX_MAX_DIMENSIONS];

	struct amatrix_weight_t cached_result;
	bool cached_result_is_valid;
};

void amatrix_save(struct amatrix_t *amx, struct amatrix_backup_t *backup);
void amatrix_restore(struct amatrix_t *amx, struct amatrix_backup_t *backup);

bool amatrix_check_consistency(struct amatrix_t *amx);
bool amatrix_is_physical(struct amatrix_t *amx);

void amatrix_to_python(struct amatrix_t *amx);
void amatrix_to_wolfram(struct amatrix_t *amx);

bool actual_amatrix_check_connectedness(struct amatrix_t *amx);
bool amatrix_check_connectedness(struct amatrix_t *amx);
struct amatrix_weight_t amatrix_weight(struct amatrix_t *amx);

#endif //__AMATRIX_H__
