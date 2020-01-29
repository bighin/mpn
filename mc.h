#ifndef __MC_H__
#define __MC_H__

#include "amatrix.h"

#define UPDATE_UNPHYSICAL       (0)
#define UPDATE_REJECTED         (1)
#define UPDATE_ACCEPTED         (2)
#define UPDATE_ERROR            (3)

int update_shuffle(struct amatrix_t *amx, bool always_accept);
int update_extend(struct amatrix_t *amx, bool always_accept);
int update_squeeze(struct amatrix_t *amx, bool always_accept);

#endif //__MC_H__
