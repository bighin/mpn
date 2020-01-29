#ifndef __READER_H__
#define __READER_H__

#include <stdio.h>
#include <stdbool.h>

struct energies_ctx_t
{
	int nso,nocc,nvirt;
	double *eocc,*evirt;
	double *htensor;
};

bool load_energies(FILE *in, struct energies_ctx_t *ctx);

double get_occupied_energy(struct energies_ctx_t *ctx,int n);
double get_virtual_energy(struct energies_ctx_t *ctx,int n);
double get_htensor(struct energies_ctx_t *ctx,int i,int j,int a,int b);

#endif //__READER_H__
