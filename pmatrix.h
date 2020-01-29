#ifndef __IMATRIX_H__
#define __IMATRIX_H__

#include <gsl/gsl_rng.h>

#define IMATRIX_MAX_DIMENSIONS	(512)

struct pmatrix_t
{
	/*
		The dimensions, maximum value and the actual values in the matrix
	*/

	int dimensions,nr_occupied,nr_virtual;
	int values[IMATRIX_MAX_DIMENSIONS][IMATRIX_MAX_DIMENSIONS];
};

struct pmatrix_t *init_pmatrix(int nr_occupied,int nr_virtual,gsl_rng *rngctx);
void fini_pmatrix(struct pmatrix_t *pmx);

int pmatrix_get_entry(struct pmatrix_t *pmx, int i, int j);
void pmatrix_set_entry(struct pmatrix_t *pmx, int i, int j, int value);
void pmatrix_inc_entry(struct pmatrix_t *pmx, int i, int j);
void pmatrix_dec_entry(struct pmatrix_t *pmx, int i, int j);
void pmatrix_print(struct pmatrix_t *pmx);
int pmatrix_sum_row(struct pmatrix_t *pmx, int row);
int pmatrix_sum_column(struct pmatrix_t *pmx, int column);
int pmatrix_trace(struct pmatrix_t *pmx);

int pmatrix_extend(struct pmatrix_t *pmx, gsl_rng *rngctx, int *targeti, int *targetj);
int pmatrix_squeeze(struct pmatrix_t *pmx, gsl_rng *rngctx, int *targeti, int *targetj);

void pmatrix_swap_rows(struct pmatrix_t *pmx, int i1, int i2, int update[2], int reverse[2], gsl_rng *rngctx);
void pmatrix_swap_cols(struct pmatrix_t *pmx, int i1, int i2, int update[2], int reverse[2], gsl_rng *rngctx);

bool pmatrix_check_consistency(struct pmatrix_t *pmx);

int pmatrix_entry_type(int i,int j);

int pmatrix_get_new_occupied_value(struct pmatrix_t *pmx, gsl_rng *rngctx);
int pmatrix_get_new_virtual_value(struct pmatrix_t *pmx, gsl_rng *rngctx);
int pmatrix_get_new_value(struct pmatrix_t *pmx, gsl_rng *rngctx, int i, int j);

#endif //__IMATRIX_H__
