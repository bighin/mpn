#ifndef __JIT_H__
#define __JIT_H__

#include <jit/jit.h>

#include "amatrix.h"

jit_function_t weight_to_jit(struct amatrix_t *amx,struct amatrix_weight_t *awt,int label1,int label2,jit_context_t context);

#endif //__JIT_H__
