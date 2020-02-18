#ifndef __REGULARIZATION_H__
#define __REGULARIZATION_H__

#include "amatrix.h"

double regularize(struct amatrix_t *amx,double denominator);

char *regularization_type_description(int regularization);
char *resummation_type_description(int resummation);

#endif //__REGULARIZATION_H__
