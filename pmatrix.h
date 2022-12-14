#ifndef __PMATRIX_H__
#define __PMATRIX_H__

#include <gsl/gsl_rng.h>

#include "loaderis.h"
#include "limits.h"

struct pmatrix_t
{
	/*
		The dimensions, maximum value and the actual values in the matrix
	*/

	int dimensions,nr_occupied,nr_virtual;
	int values[PMATRIX_MAX_DIMENSIONS][PMATRIX_MAX_DIMENSIONS];
};

struct pmatrix_t *init_pmatrix(int nr_occupied,int nr_virtual,gsl_rng *rngctx);
void fini_pmatrix(struct pmatrix_t *pmx);

int pmatrix_get_entry(struct pmatrix_t *pmx, int i, int j);
void pmatrix_set_entry(struct pmatrix_t *pmx, int i, int j, int value);

int pmatrix_get_raw_entry(struct pmatrix_t *pmx, int i, int j);
void pmatrix_set_raw_entry(struct pmatrix_t *pmx, int i, int j, int value);

void pmatrix_print(struct pmatrix_t *pmx);

void pmatrix_extend(struct pmatrix_t *pmx, gsl_rng *rngctx, int *targeti, int *targetj);
void pmatrix_squeeze(struct pmatrix_t *pmx, gsl_rng *rngctx);

void pmatrix_swap_rows(struct pmatrix_t *pmx, int i1, int i2, gsl_rng *rngctx);
void pmatrix_swap_cols(struct pmatrix_t *pmx, int i1, int i2, gsl_rng *rngctx);

bool pmatrix_check_consistency(struct pmatrix_t *pmx);

int pmatrix_entry_type(int i,int j);
int pmatrix_get_new_value(struct pmatrix_t *pmx, gsl_rng *rngctx, int i, int j);

#endif //__PMATRIX_H__
