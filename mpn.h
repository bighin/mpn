#ifndef __MPN_H__
#define __MPN_H__

#include <stdbool.h>

#include <gsl/gsl_matrix.h>

#include "reader.h"
#include "amatrix.h"

struct label_t
{
	int id;
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

double incidence_to_weight(gsl_matrix_int *B, struct label_t *labels, int *ilabels, struct energies_ctx_t *ectx,double unphysical_penalty,bool verbose);
gsl_matrix_int *amatrix_calculate_incidence(struct amatrix_t *amx, struct label_t labels[MAX_LABELS], int *ilabels);

int do_diagmc(char *energies_dot_dat,char *output,int iterations,double timelimit,bool show_progressbar);

#endif //__MPN_H__
