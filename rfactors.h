#ifndef __RFACTORS_H__
#define __RFACTORS_H__

#include "amatrix.h"

void rfactors_init(void);
void rfactors_sample_sign(struct amatrix_t *amx, int sign);
void rfactors_output_summary(const char *filename);

#endif //__RFACTORS_H__
