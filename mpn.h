#ifndef __MPN_H__
#define __MPN_H__

#include <stdbool.h>

#include <gsl/gsl_matrix.h>

#include "loaderis.h"
#include "amatrix.h"

struct label_t
{
	int id;
	int i,j;
	char mnemonic;
	int value;

#define QTYPE_OCCUPIED	(0)
#define QTYPE_VIRTUAL	(1)

	int qtype;
	bool visited;
	bool selfloop;
};

/*
	These two limits could be estimated in a better way!
*/

#define MAX_LABELS		(512)
#define MAX_MATRIX_ELEMENTS	(512)

struct amatrix_weight_t incidence_to_weight(gsl_matrix_int *B, struct label_t *labels, int *ilabels, struct amatrix_t *amx);
gsl_matrix_int *amatrix_calculate_incidence(struct amatrix_t *amx, struct label_t labels[MAX_LABELS], int *ilabels);

#endif //__MPN_H__
