#ifndef __WEIGHT2_H__
#define __WEIGHT2_H__

#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "mpn.h"
#include "cache.h"
#include "limits.h"

/*
	This structure saves most of the intermediate results that can be used to recalculate the weight after
	having modified a quantum number
*/

struct weight_info_t
{
	/*
		This information can be used to reconstruct the weight, given the labels.
	*/

	int l,h;

	struct
	{
		int labels[MAX_LABELS];
		int qtypes[MAX_LABELS];
		int ilabels;
	}
	denominators[MAX_DENOMINATORS];
	int nr_denominators;

	struct
	{
		int labels[4];
	}
	numerators[MAX_NUMERATORS];
	int nr_numerators;

	double inversefactor,unphysical_penalty,weight;

	struct label_t labels[MAX_LABELS];
	int ilabels;
};

struct weight_info_t incidence_to_weight_info(gsl_matrix_int *B, struct label_t *labels, int *ilabels, struct amatrix_t *amx);

double reconstruct_weight(struct amatrix_t *amx, struct weight_info_t *awt);

int coordinate_to_label_index(struct label_t *labels,int ilabels,int i,int j,int pmatrix);

#endif //__WEIGHT2_H__
