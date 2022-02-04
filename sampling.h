#ifndef __SAMPLING_H__
#define __SAMPLING_H__

#include <stdio.h>
#include <stdbool.h>

#include "amatrix.h"
#include "config.h"

struct sampling_ctx_t;

struct sampling_ctx_t *init_sampling_ctx(int maxdimensions);
void fini_sampling_ctx(struct sampling_ctx_t *sctx);

void sampling_ctx_measure(struct sampling_ctx_t *sctx,struct amatrix_t *amx,struct configuration_t *config,long int counter);
double sampling_ctx_get_physical_pct(struct sampling_ctx_t *sctx);
void sampling_ctx_print_report(struct sampling_ctx_t *sctx,struct amatrix_t *amx,FILE *out,bool finalize);

void rfactors_sample_sign(struct amatrix_t *amx, int sign);
void rfactors_output_summary(const char *filename);

#endif //__SAMPLING_H__
