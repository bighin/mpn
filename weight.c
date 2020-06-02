#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "weight.h"
#include "mpn.h"
#include "permutations.h"
#include "auxx.h"

struct amatrix_collection_t
{
	struct amatrix_t **amatrices;
	bool *mustfree;

	int iamatrices,nralloced;
};

struct amatrix_collection_t *init_collection(void)
{
	struct amatrix_collection_t *ret;

	if(!(ret=malloc(sizeof(struct amatrix_collection_t))))
		return NULL;

	ret->amatrices=NULL;
	ret->mustfree=NULL;
	ret->iamatrices=0;
	ret->nralloced=0;

	return ret;
}

void fini_collection(struct amatrix_collection_t *acl)
{
	if(acl)
	{
		for(int c=0;c<acl->iamatrices;c++)
			if(acl->mustfree[c]==true)
				fini_amatrix(acl->amatrices[c],false);

		if(acl->amatrices)
			free(acl->amatrices);

		if(acl->mustfree)
			free(acl->mustfree);

		free(acl);
	}
}

void collection_add(struct amatrix_collection_t *acl,struct amatrix_t *amx,bool mustfree)
{
	if(acl->iamatrices>=acl->nralloced)
	{
		assert(acl->iamatrices==acl->nralloced);

		acl->nralloced+=16;
		acl->amatrices=realloc(acl->amatrices,sizeof(struct amatrix_t *)*acl->nralloced);
		acl->mustfree=realloc(acl->mustfree,sizeof(bool)*acl->nralloced);
	}

	acl->amatrices[acl->iamatrices]=amx;
	acl->mustfree[acl->iamatrices]=mustfree;
	acl->iamatrices++;
}

struct amatrix_t *init_amatrix_from_amatrix(struct amatrix_t *amx)
{
	struct amatrix_t *ret=malloc(sizeof(struct amatrix_t));

	assert(ret!=NULL);

	ret->ectx=amx->ectx;
	ret->rng_ctx=amx->rng_ctx;

	ret->nr_occupied=amx->nr_occupied;
	ret->nr_virtual=amx->nr_virtual;

	assert((ret->ectx!=NULL)&&(ret->rng_ctx!=NULL));

	ret->pmxs[0]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);
	ret->pmxs[1]=init_pmatrix(ret->nr_occupied, ret->nr_virtual, ret->rng_ctx);

	ret->config=amx->config;

	assert(ret->config!=NULL);

	ret->cached_weight=0.0f;
	ret->cached_weight_is_valid=false;

	return ret;
}

#define ORDERING_GT	(1)
#define ORDERING_LT	(-1)
#define ORDERING_EQ	(0)

int amatrix_ordering(struct amatrix_t *first,struct amatrix_t *second)
{
	assert(first->pmxs[0]->dimensions==second->pmxs[0]->dimensions);
	assert(first->pmxs[1]->dimensions==second->pmxs[1]->dimensions);

	for(int i=0;i<first->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<first->pmxs[0]->dimensions;j++)
		{
			if(first->pmxs[0]->values[i][j]>second->pmxs[0]->values[i][j])
				return ORDERING_GT;
			else if(first->pmxs[0]->values[i][j]<second->pmxs[0]->values[i][j])
				return ORDERING_LT;

			if(first->pmxs[1]->values[i][j]>second->pmxs[1]->values[i][j])
				return ORDERING_GT;
			else if(first->pmxs[1]->values[i][j]<second->pmxs[1]->values[i][j])
				return ORDERING_LT;
		}
	}

	return ORDERING_EQ;
}

void cdet_generate_grouping(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl,amx,false);

	/*
		The transposed adjacency matrix
	*/

	struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

	twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
	twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

	assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

	for(int i=0;i<twin->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<twin->pmxs[0]->dimensions;j++)
		{
			twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[j][i];
			twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[j][i];
		}
	}

	collection_add(acl,twin,true);

	/*
		The anti-transposed adjacency matrix
	*/

	twin=init_amatrix_from_amatrix(amx);

	twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
	twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

	assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

	for(int i=0;i<twin->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<twin->pmxs[0]->dimensions;j++)
		{
			int dimensions=twin->pmxs[0]->dimensions;

			twin->pmxs[0]->values[dimensions-j-1][dimensions-i-1]=amx->pmxs[0]->values[i][j];
			twin->pmxs[1]->values[dimensions-j-1][dimensions-i-1]=amx->pmxs[1]->values[i][j];
		}
	}

	collection_add(acl,twin,true);

	/*
		Finally we apply both operations at once
	*/

	twin=init_amatrix_from_amatrix(amx);

	twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
	twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

	assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

	for(int i=0;i<twin->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<twin->pmxs[0]->dimensions;j++)
		{
			int dimensions=twin->pmxs[0]->dimensions;

			twin->pmxs[0]->values[dimensions-i-1][dimensions-j-1]=amx->pmxs[0]->values[i][j];
			twin->pmxs[1]->values[dimensions-i-1][dimensions-j-1]=amx->pmxs[1]->values[i][j];
		}
	}

	collection_add(acl,twin,true);
}

void cdet_generate_grouping_alt(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The grouping can have any form, but it has to satisfy the property
		that given a diagram D, and the grouping F(D)

		D' \in f(D) ==> D \in f(D')

		so that the groupings form a partition of the set of all diagram, and
		that the is-in-the-same-grouping relationship is commutative.
	*/

	int flags=0x00;

	/*
		The original adjacency matrix
	*/

	collection_add(acl,amx,false);

	/*
		The transposed adjacency matrix
	*/

	if(flags&0x01)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[j][i];
				twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[j][i];
			}
		}

		collection_add(acl, twin,true);
	}

	/*
		The anti-transposed adjacency matrix
	*/

	if(flags&0x02)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				int dimensions=twin->pmxs[0]->dimensions;

				twin->pmxs[0]->values[dimensions-j-1][dimensions-i-1]=amx->pmxs[0]->values[i][j];
				twin->pmxs[1]->values[dimensions-j-1][dimensions-i-1]=amx->pmxs[1]->values[i][j];
			}
		}

		collection_add(acl, twin, true);
	}

	/*
		Finally we apply both operations at once
	*/

	if(flags&0x04)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				int dimensions=twin->pmxs[0]->dimensions;

				twin->pmxs[0]->values[dimensions-i-1][dimensions-j-1]=amx->pmxs[0]->values[i][j];
				twin->pmxs[1]->values[dimensions-i-1][dimensions-j-1]=amx->pmxs[1]->values[i][j];
			}
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_alt2(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The grouping can have any form, but it has to satisfy the property
		that given a diagram D, and the grouping F(D)

		D' \in f(D) ==> D \in f(D')

		so that the groupings form a partition of the set of all diagram, and
		that the is-in-the-same-grouping relationship is commutative.
	*/

	int flags=0x02;

	/*
		The original adjacency matrix
	*/

	if(flags&0x01)
	{
		collection_add(acl, amx, false);
	}

	/*
		The transposed adjacency matrix
	*/

	if(flags&0x02)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[j][i];
				twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[j][i];
			}
		}

		collection_add(acl, twin,true);
	}

	/*
		The anti-transposed adjacency matrix
	*/

	if(flags&0x04)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				int dimensions=twin->pmxs[0]->dimensions;

				twin->pmxs[0]->values[dimensions-j-1][dimensions-i-1]=amx->pmxs[0]->values[i][j];
				twin->pmxs[1]->values[dimensions-j-1][dimensions-i-1]=amx->pmxs[1]->values[i][j];
			}
		}

		collection_add(acl, twin, true);
	}

	/*
		Finally we apply both operations at once
	*/

	if(flags&0x08)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				int dimensions=twin->pmxs[0]->dimensions;

				twin->pmxs[0]->values[dimensions-i-1][dimensions-j-1]=amx->pmxs[0]->values[i][j];
				twin->pmxs[1]->values[dimensions-i-1][dimensions-j-1]=amx->pmxs[1]->values[i][j];
			}
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altS2(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	/*
		The twin
	*/

	if(amx->pmxs[0]->dimensions==4)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[j][i];
				twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[j][i];
			}
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altS3(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	/*
		The twin
	*/

	if(amx->pmxs[0]->dimensions==4)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

#define FIRST_HALF(x)	(((x)==0)||((x)==1))
#define SECOND_HALF(x)	(((x)==2)||((x)==3))

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				if(FIRST_HALF(i)&&FIRST_HALF(j))
				{
					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][j];
				}

				if(FIRST_HALF(i)&&SECOND_HALF(j))
				{
					int jprime=(j==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][jprime];
				}

				if(SECOND_HALF(i)&&FIRST_HALF(j))
				{
					int iprime=(i==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][j];
				}

				if(SECOND_HALF(i)&&SECOND_HALF(j))
				{
					int iprime=(i==3)?(2):(3);
					int jprime=(j==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][jprime];
				}
			}
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altS4(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	/*
		The twin
	*/

	if(amx->pmxs[0]->dimensions==4)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

#define FIRST_HALF(x)	(((x)==0)||((x)==1))
#define SECOND_HALF(x)	(((x)==2)||((x)==3))

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				if(FIRST_HALF(i)&&FIRST_HALF(j))
				{
					twin->pmxs[0]->values[j][i]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[j][i]=amx->pmxs[1]->values[i][j];
				}

				if(FIRST_HALF(i)&&SECOND_HALF(j))
				{
					int iprime=(i==0)?(1):(0);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][j];
				}

				if(SECOND_HALF(i)&&FIRST_HALF(j))
				{
					int jprime=(j==0)?(1):(0);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][jprime];
				}

				if(SECOND_HALF(i)&&SECOND_HALF(j))
				{
					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][j];
				}
			}
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altSx(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	/*
		The second twin
	*/

	if(amx->pmxs[0]->dimensions==4)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[j][i];
				twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[j][i];
			}
		}

		collection_add(acl, twin, true);
	}

	/*
		The third twin
	*/

	if(amx->pmxs[0]->dimensions==4)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

#define FIRST_HALF(x)	(((x)==0)||((x)==1))
#define SECOND_HALF(x)	(((x)==2)||((x)==3))

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				if(FIRST_HALF(i)&&FIRST_HALF(j))
				{
					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][j];
				}

				if(FIRST_HALF(i)&&SECOND_HALF(j))
				{
					int jprime=(j==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][jprime];
				}

				if(SECOND_HALF(i)&&FIRST_HALF(j))
				{
					int iprime=(i==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][j];
				}

				if(SECOND_HALF(i)&&SECOND_HALF(j))
				{
					int iprime=(i==3)?(2):(3);
					int jprime=(j==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][jprime];
				}
			}
		}

		collection_add(acl, twin, true);
	}

	/*
		The fourth twin
	*/

	if(amx->pmxs[0]->dimensions==4)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		assert(twin->pmxs[0]->dimensions==twin->pmxs[1]->dimensions);

#define FIRST_HALF(x)	(((x)==0)||((x)==1))
#define SECOND_HALF(x)	(((x)==2)||((x)==3))

		for(int i=0;i<twin->pmxs[0]->dimensions;i++)
		{
			for(int j=0;j<twin->pmxs[0]->dimensions;j++)
			{
				if(FIRST_HALF(i)&&FIRST_HALF(j))
				{
					twin->pmxs[0]->values[j][i]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[j][i]=amx->pmxs[1]->values[i][j];
				}

				if(FIRST_HALF(i)&&SECOND_HALF(j))
				{
					int iprime=(i==0)?(1):(0);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][j];
				}

				if(SECOND_HALF(i)&&FIRST_HALF(j))
				{
					int jprime=(j==0)?(1):(0);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][jprime];
				}

				if(SECOND_HALF(i)&&SECOND_HALF(j))
				{
					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][j];
				}
			}
		}

		collection_add(acl, twin, true);
	}
}

int get_permutation(int order,int index,int element)
{
	switch(order)
	{
		case 2:
		return permutations2[index][element];

		case 3:
		return permutations3[index][element];

		case 4:
		return permutations4[index][element];

		case 5:
		return permutations5[index][element];

		case 6:
		return permutations6[index][element];
	}

	return 0;
}

void cdet_generate_grouping_altR(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	/*
		We now generate the twins
	*/

	int dimensions=amx->pmxs[0]->dimensions;

	struct
	{
		int i, j, pmatrix;
	} candidates[32];

	int icandidates=0;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
			{
				candidates[icandidates].i=i;
				candidates[icandidates].j=j;
				candidates[icandidates].pmatrix=0;
				icandidates++;
			}

			if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
			{
				candidates[icandidates].i=i;
				candidates[icandidates].j=j;
				candidates[icandidates].pmatrix=1;
				icandidates++;
			}
		}
	}

	/*
		We need to choose a candidate label according to a deterministic procedure:
		we use the simple Fowler–Noll–Vo 1a hash function...
	*/

	uint64_t hash=0xcbf29ce484222325;

	const uint8_t bytes[2][32][32]={{{171, 212, 172, 236, 140, 73,  207, 209, 67,  218, 70,  175, 66,  59, 53,  231, 161, 89,  141, 255, 207, 192, 191, 73, 210, 244, 11,  30, 93,  19,  154, 71}, {210, 17, 66,  82,  184, 115, 206, 139, 33,  233, 191, 195, 84, 153, 66,  202, 112, 196, 127, 137, 139, 47, 109, 129, 103, 56,  173, 136, 212, 137, 0,  154}, {89, 233, 69,  92, 105, 165, 36, 234, 204, 198, 162, 71, 157, 235, 33,  247, 169, 18, 209, 224, 153, 250, 99, 196, 40,  109, 223, 196, 114, 102, 159, 174}, {251, 116, 18,  124, 99, 31, 222, 19,  255, 4,   247, 210, 86, 29,  45,  241, 101, 155, 197, 25,  224, 144, 155, 35,  221, 6,  136, 125, 197, 208, 243, 239}, {107, 20,  46,  56, 249, 164, 246, 42, 129, 254, 121, 211, 33,  128, 251, 180, 89, 220, 39,  45, 211, 49,  21,  164, 222, 194, 142, 91,  234, 214, 161, 146}, {105, 87,  48, 92,  127, 85, 232, 159, 165, 210, 109, 99, 176, 3,   215, 73,  235, 84,  177, 105, 14,  19, 215, 217, 47,  98, 169, 96,  53, 178, 111, 5},   {46,  10,  145, 180, 232, 51,  21, 145, 225, 84, 218, 132, 47, 117, 128, 95, 77, 191, 145, 103, 51, 158, 19, 122, 114, 158, 179, 198, 25, 56, 51, 255}, {93, 209, 129, 17, 90, 39, 219, 242, 136, 236, 27, 85, 58, 39, 108, 166, 70, 189, 100, 128, 243, 87, 217, 179, 77, 169, 142, 145, 1, 117, 253, 190}, {132, 224, 148, 239, 11, 12, 160, 153, 144, 128, 38, 49, 189, 231, 126, 77, 246, 114, 70, 76, 118, 72, 103, 70, 200, 66, 244, 103, 246, 192, 61, 117}, {85, 1, 9, 196, 235, 130, 13, 43, 59, 167, 88, 112, 232, 216, 209, 133, 150, 86, 237, 231, 76, 80, 208, 207, 37, 129, 246, 5, 0, 104, 214, 178}, {194, 51, 199, 165, 202, 237, 157, 17, 81, 64, 236, 135, 100, 229, 151, 247, 182, 159, 162, 194, 115, 186, 104, 173, 245, 162, 214, 229, 204, 92, 227, 111}, {146, 46, 77, 225, 255, 181, 139, 159, 224, 232, 65, 49, 195, 110, 15, 100, 176, 230, 10, 135, 11, 46, 79, 227, 246, 188, 206, 57, 37, 134, 27, 67}, {1, 122, 217, 122, 32, 31, 77, 170, 238, 234, 89, 67, 112, 203, 71, 83, 156, 119, 4, 178, 82, 238, 98, 219, 129, 166, 69, 128, 253, 163, 221, 52}, {149, 77, 99, 118, 133, 240, 122, 52, 104, 114, 157, 84, 217, 250, 17, 131, 29, 228, 126, 245, 227, 78, 193, 216, 250, 220, 10, 42, 95, 107, 69, 109}, {240, 16, 48, 93, 208, 219, 240, 141, 217, 251, 130, 93, 52, 221, 28, 210, 237, 61, 21, 234, 114, 90, 181, 37, 179, 54, 8, 129, 121, 87, 187, 53}, {106, 236, 70, 233, 226, 224, 187, 181, 226, 12, 101, 23, 213, 58, 83, 242, 68, 89, 101, 135, 97, 173, 230, 82, 215, 215, 72, 56, 182, 13, 193, 152}, {120, 1, 152, 99, 145, 19, 197, 185, 100, 133, 60, 125, 8, 113, 188, 165, 180, 176, 132, 7, 100, 52, 118, 197, 71, 57, 213, 10, 29, 4, 202, 230}, {218, 2, 200, 67, 71, 211, 65, 212, 234, 29, 37, 103, 191, 12, 179, 134, 253, 52, 156, 79, 225, 17, 188, 231, 238, 76, 242, 202, 219, 254, 6, 133}, {100, 4, 17, 127, 77, 192, 150, 37, 72, 240, 128, 108, 207, 144, 63, 47, 12, 200, 143, 157, 105, 170, 206, 28, 234, 51, 133, 125, 195, 120, 181, 49}, {224, 166, 191, 152, 214, 83, 193, 177, 18, 218, 115, 163, 254, 249, 177, 93, 64, 15, 84, 6, 93, 228, 228, 212, 150, 38, 228, 165, 164, 82, 134, 105}, {202, 70, 239, 186, 105, 110, 156, 80, 234, 234, 6, 50, 18, 78, 77, 254, 154, 79, 235, 29, 3, 67, 80, 176, 181, 74, 41, 108, 8, 122, 108, 235}, {113, 201, 202, 169, 97, 118, 143, 223, 84, 95, 61, 135, 133, 146, 252, 102, 126, 242, 95, 253, 181, 101, 255, 197, 147, 223, 213, 137, 228, 101, 65, 147}, {185, 36, 221, 167, 34, 147, 201, 8, 2, 59, 68, 72, 69, 65, 210, 207, 104, 197, 243, 210, 25, 31, 217, 174, 97, 123, 107, 27, 178, 171, 21, 203}, {209, 87, 77, 70, 236, 113, 65, 25, 92, 52, 184, 110, 227, 78, 101, 44, 90, 239, 35, 143, 20, 239, 214, 248, 224, 247, 214, 44, 206, 48, 32, 96}, {136, 99, 80, 39, 158, 193, 108, 216, 45, 122, 237, 173, 55, 64, 148, 14, 205, 21, 28, 55, 235, 151, 3, 248, 219, 60, 150, 52, 224, 221, 203, 99}, {119, 152, 76, 113, 140, 176, 227, 71, 125, 192, 235, 105, 186, 82, 52, 47, 195, 158, 32, 199, 70, 187, 200, 204, 242, 29, 152, 162, 70, 210, 45, 69}, {29, 166, 227, 57, 217, 21, 75, 20, 149, 18, 172, 132, 161, 181, 108, 148, 86, 192, 14, 223, 16, 79, 122, 112, 125, 85, 176, 21, 217, 192, 168, 81}, {27, 95, 4, 95, 226, 143, 174, 134, 21, 178, 130, 122, 177, 121, 39, 98, 148, 31, 101, 224, 19, 119, 163, 91, 113, 110, 191, 64, 173, 251, 160, 226}, {106, 107, 74, 141, 70, 34, 220, 159, 26, 113, 65, 18, 89, 213, 189, 86, 34, 110, 180, 44, 118, 21, 167, 169, 132, 105, 248, 160, 245, 35, 178, 230}, {97, 166, 157, 88, 88, 230, 93, 7, 22, 72, 196, 201, 74, 219, 120, 161, 161, 228, 24, 115, 184, 196, 149, 46, 145, 48, 114, 251, 16, 106, 47, 40}, {3, 84, 157, 225, 113, 214, 148, 10, 183, 25, 156, 116, 73, 129, 130, 174, 24, 232, 232, 121, 32, 111, 14, 172, 202, 67, 95, 255, 239, 220, 165, 147}, {83, 93, 3, 192, 14, 4, 125, 194, 21, 113, 74, 46, 126, 151, 52, 253, 221, 126, 205, 241, 28, 25, 73, 122, 164, 162, 101, 83, 19, 14, 234, 229}},
				    {{243, 150, 20,  241, 234, 201, 39,  58,  153, 56,  184, 9,   121, 64, 203, 191, 82,  125, 5,   124, 119, 233, 193, 42, 162, 143, 243, 38, 195, 137, 57,  61}, {200, 5,  206, 198, 126, 153, 240, 84,  181, 81,  246, 152, 15, 209, 229, 42,  59,  39,  20,  180, 142, 9,  165, 221, 125, 130, 109, 154, 113, 43,  54, 103}, {53, 165, 221, 97, 55,  148, 35, 117, 39,  154, 249, 76, 139, 178, 191, 70,  45,  65, 132, 33,  250, 67,  24, 161, 121, 151, 245, 212, 41,  165, 248, 246}, {35,  88,  112, 89,  62, 40, 108, 157, 73,  212, 129, 155, 5,  236, 183, 247, 77,  153, 247, 155, 41,  196, 139, 146, 187, 31, 191, 181, 206, 131, 32,  203}, {253, 109, 230, 96, 51,  164, 249, 35, 242, 218, 135, 216, 114, 121, 82,  127, 10, 220, 223, 98, 45,  106, 250, 0,   211, 255, 116, 137, 130, 141, 85,  9},   {86,  135, 41, 251, 20,  10, 82,  121, 209, 214, 55,  12, 183, 249, 66,  119, 186, 227, 126, 202, 158, 38, 225, 147, 234, 86, 157, 186, 10, 34,  73,  135}, {247, 209, 152, 99,  210, 141, 45, 136, 81, 30, 87, 64, 130, 32, 129, 134, 231, 202, 239, 82, 34, 102, 164, 89, 133, 207, 104, 147, 219, 109, 157, 16}, {251, 215, 8, 166, 95, 27, 72, 144, 25, 165, 18, 36, 23, 190, 101, 122, 217, 192, 23, 134, 51, 51, 173, 54, 117, 36, 56, 215, 159, 86, 41, 60}, {82, 206, 177, 168, 165, 31, 108, 65, 36, 180, 62, 115, 23, 49, 128, 86, 86, 161, 241, 165, 149, 0, 17, 183, 58, 18, 186, 51, 194, 203, 65, 233}, {200, 194, 169, 149, 242, 241, 73, 245, 117, 112, 41, 143, 123, 0, 36, 149, 83, 83, 143, 102, 97, 214, 129, 247, 41, 186, 23, 16, 30, 91, 18, 98}, {48, 115, 44, 194, 22, 112, 11, 191, 102, 116, 221, 154, 40, 34, 84, 42, 187, 224, 46, 81, 11, 118, 2, 81, 37, 146, 26, 39, 253, 244, 90, 189}, {173, 195, 128, 116, 23, 131, 210, 70, 229, 153, 172, 204, 117, 217, 174, 125, 149, 79, 165, 222, 9, 87, 58, 224, 14, 238, 60, 7, 212, 196, 45, 64}, {27, 128, 3, 55, 183, 153, 176, 122, 212, 98, 214, 26, 115, 156, 92, 181, 224, 14, 10, 57, 111, 40, 194, 173, 9, 170, 31, 150, 136, 191, 84, 215}, {253, 68, 139, 23, 54, 70, 178, 231, 28, 212, 229, 29, 158, 186, 103, 33, 8, 102, 41, 251, 111, 202, 235, 120, 200, 103, 255, 187, 10, 184, 156, 31}, {83, 239, 213, 245, 244, 201, 215, 203, 201, 38, 248, 166, 190, 46, 252, 137, 107, 102, 152, 153, 211, 106, 143, 30, 188, 176, 202, 68, 198, 203, 95, 171}, {76, 22, 63, 172, 92, 247, 200, 98, 215, 109, 113, 235, 172, 163, 242, 254, 238, 219, 62, 67, 75, 4, 216, 138, 191, 98, 185, 78, 18, 14, 100, 52}, {49, 169, 138, 32, 104, 50, 187, 247, 225, 140, 64, 119, 132, 113, 28, 150, 3, 43, 56, 59, 69, 176, 22, 250, 246, 254, 7, 159, 212, 197, 240, 244}, {211, 8, 218, 192, 38, 60, 111, 56, 226, 37, 46, 207, 229, 209, 11, 158, 3, 142, 183, 99, 90, 65, 206, 5, 4, 2, 133, 20, 173, 117, 107, 134}, {72, 177, 224, 19, 166, 240, 157, 47, 115, 33, 181, 194, 2, 48, 79, 192, 83, 204, 5, 143, 227, 91, 102, 211, 186, 180, 61, 215, 251, 127, 191, 159}, {232, 100, 222, 17, 227, 148, 20, 101, 252, 71, 23, 118, 190, 215, 105, 61, 91, 175, 223, 217, 245, 6, 32, 73, 51, 5, 153, 238, 207, 152, 155, 89}, {123, 122, 170, 169, 178, 28, 182, 42, 22, 141, 19, 136, 10, 143, 13, 108, 229, 28, 67, 236, 15, 22, 111, 66, 44, 200, 124, 106, 9, 23, 138, 21}, {254, 131, 85, 194, 222, 0, 151, 199, 200, 115, 12, 237, 169, 33, 97, 235, 63, 20, 37, 231, 188, 109, 97, 233, 104, 213, 220, 123, 33, 80, 86, 86}, {216, 172, 164, 202, 8, 23, 50, 56, 108, 107, 234, 184, 139, 80, 13, 154, 37, 79, 117, 117, 47, 118, 37, 136, 120, 115, 30, 160, 243, 16, 229, 48}, {78, 219, 254, 237, 36, 158, 161, 45, 246, 104, 191, 197, 117, 27, 124, 1, 45, 134, 152, 15, 196, 145, 192, 131, 114, 86, 185, 24, 67, 37, 113, 222}, {91, 215, 42, 253, 102, 10, 185, 21, 121, 199, 64, 212, 1, 30, 40, 45, 241, 91, 147, 20, 65, 248, 49, 134, 215, 105, 16, 236, 215, 111, 117, 221}, {163, 184, 137, 234, 221, 6, 25, 158, 188, 199, 150, 170, 114, 124, 204, 140, 126, 79, 7, 217, 230, 187, 136, 233, 159, 79, 199, 19, 255, 212, 95, 161}, {5, 92, 168, 45, 124, 74, 155, 59, 212, 60, 6, 204, 28, 47, 210, 97, 252, 239, 236, 38, 87, 3, 121, 29, 64, 226, 153, 249, 46, 136, 243, 199}, {50, 102, 110, 162, 35, 133, 196, 139, 237, 45, 206, 231, 83, 126, 5, 232, 11, 145, 119, 246, 188, 6, 72, 34, 76, 48, 11, 105, 234, 230, 76, 11}, {243, 153, 74, 138, 62, 99, 128, 138, 142, 123, 163, 97, 131, 95, 61, 82, 247, 156, 218, 79, 89, 90, 185, 94, 151, 149, 104, 156, 186, 148, 106, 10}, {137, 250, 75, 205, 222, 43, 166, 10, 151, 114, 108, 149, 22, 136, 141, 25, 67, 10, 175, 20, 159, 55, 207, 94, 166, 218, 53, 18, 252, 213, 128, 73}, {90, 205, 75, 93, 83, 84, 239, 148, 89, 246, 72, 14, 47, 234, 85, 97, 229, 109, 47, 106, 229, 167, 29, 39, 23, 38, 153, 195, 72, 204, 61, 146}, {247, 218, 227, 122, 232, 204, 31, 22, 30, 164, 199, 235, 207, 151, 202, 89, 241, 197, 22, 43, 127, 40, 182, 125, 103, 50, 169, 219, 146, 32, 59, 101}}};

	for(int c=0;c<icandidates;c++)
	{
		int i,j,pmatrix;

		i=candidates[c].i;
		j=candidates[c].j;
		pmatrix=candidates[c].pmatrix;

		hash^=bytes[pmatrix][i][j];
		hash*=0x100000001b3;
	}

	int chosen_candidate=hash%icandidates;
	int qtype=pmatrix_entry_type(candidates[chosen_candidate].i,candidates[chosen_candidate].j);

	/*
		...and then we sum over all possible values of the chosen candidate label.
	*/

	for(int value=1;value<((qtype==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));value++)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		for(int i=0;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[0], i, j);
				twin->pmxs[1]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[1], i, j);
			}
		}

		pmatrix_set_entry(twin->pmxs[candidates[chosen_candidate].pmatrix],
				  candidates[chosen_candidate].i,
				  candidates[chosen_candidate].j,
				  value);

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altRplus(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	/*
		We now generate the twins
	*/

	int dimensions=amx->pmxs[0]->dimensions;

	struct
	{
		int i, j, pmatrix;
	} candidates[32];

	int icandidates=0;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
			{
				candidates[icandidates].i=i;
				candidates[icandidates].j=j;
				candidates[icandidates].pmatrix=0;
				icandidates++;
			}

			if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
			{
				candidates[icandidates].i=i;
				candidates[icandidates].j=j;
				candidates[icandidates].pmatrix=1;
				icandidates++;
			}
		}
	}

	/*
		We need to choose a candidate label according to a deterministic procedure:
		we use the simple Fowler–Noll–Vo 1a hash function...
	*/

	uint64_t hash=0xcbf29ce484222325;

	const uint8_t bytes[2][32][32]={{{171, 212, 172, 236, 140, 73,  207, 209, 67,  218, 70,  175, 66,  59, 53,  231, 161, 89,  141, 255, 207, 192, 191, 73, 210, 244, 11,  30, 93,  19,  154, 71}, {210, 17, 66,  82,  184, 115, 206, 139, 33,  233, 191, 195, 84, 153, 66,  202, 112, 196, 127, 137, 139, 47, 109, 129, 103, 56,  173, 136, 212, 137, 0,  154}, {89, 233, 69,  92, 105, 165, 36, 234, 204, 198, 162, 71, 157, 235, 33,  247, 169, 18, 209, 224, 153, 250, 99, 196, 40,  109, 223, 196, 114, 102, 159, 174}, {251, 116, 18,  124, 99, 31, 222, 19,  255, 4,   247, 210, 86, 29,  45,  241, 101, 155, 197, 25,  224, 144, 155, 35,  221, 6,  136, 125, 197, 208, 243, 239}, {107, 20,  46,  56, 249, 164, 246, 42, 129, 254, 121, 211, 33,  128, 251, 180, 89, 220, 39,  45, 211, 49,  21,  164, 222, 194, 142, 91,  234, 214, 161, 146}, {105, 87,  48, 92,  127, 85, 232, 159, 165, 210, 109, 99, 176, 3,   215, 73,  235, 84,  177, 105, 14,  19, 215, 217, 47,  98, 169, 96,  53, 178, 111, 5},   {46,  10,  145, 180, 232, 51,  21, 145, 225, 84, 218, 132, 47, 117, 128, 95, 77, 191, 145, 103, 51, 158, 19, 122, 114, 158, 179, 198, 25, 56, 51, 255}, {93, 209, 129, 17, 90, 39, 219, 242, 136, 236, 27, 85, 58, 39, 108, 166, 70, 189, 100, 128, 243, 87, 217, 179, 77, 169, 142, 145, 1, 117, 253, 190}, {132, 224, 148, 239, 11, 12, 160, 153, 144, 128, 38, 49, 189, 231, 126, 77, 246, 114, 70, 76, 118, 72, 103, 70, 200, 66, 244, 103, 246, 192, 61, 117}, {85, 1, 9, 196, 235, 130, 13, 43, 59, 167, 88, 112, 232, 216, 209, 133, 150, 86, 237, 231, 76, 80, 208, 207, 37, 129, 246, 5, 0, 104, 214, 178}, {194, 51, 199, 165, 202, 237, 157, 17, 81, 64, 236, 135, 100, 229, 151, 247, 182, 159, 162, 194, 115, 186, 104, 173, 245, 162, 214, 229, 204, 92, 227, 111}, {146, 46, 77, 225, 255, 181, 139, 159, 224, 232, 65, 49, 195, 110, 15, 100, 176, 230, 10, 135, 11, 46, 79, 227, 246, 188, 206, 57, 37, 134, 27, 67}, {1, 122, 217, 122, 32, 31, 77, 170, 238, 234, 89, 67, 112, 203, 71, 83, 156, 119, 4, 178, 82, 238, 98, 219, 129, 166, 69, 128, 253, 163, 221, 52}, {149, 77, 99, 118, 133, 240, 122, 52, 104, 114, 157, 84, 217, 250, 17, 131, 29, 228, 126, 245, 227, 78, 193, 216, 250, 220, 10, 42, 95, 107, 69, 109}, {240, 16, 48, 93, 208, 219, 240, 141, 217, 251, 130, 93, 52, 221, 28, 210, 237, 61, 21, 234, 114, 90, 181, 37, 179, 54, 8, 129, 121, 87, 187, 53}, {106, 236, 70, 233, 226, 224, 187, 181, 226, 12, 101, 23, 213, 58, 83, 242, 68, 89, 101, 135, 97, 173, 230, 82, 215, 215, 72, 56, 182, 13, 193, 152}, {120, 1, 152, 99, 145, 19, 197, 185, 100, 133, 60, 125, 8, 113, 188, 165, 180, 176, 132, 7, 100, 52, 118, 197, 71, 57, 213, 10, 29, 4, 202, 230}, {218, 2, 200, 67, 71, 211, 65, 212, 234, 29, 37, 103, 191, 12, 179, 134, 253, 52, 156, 79, 225, 17, 188, 231, 238, 76, 242, 202, 219, 254, 6, 133}, {100, 4, 17, 127, 77, 192, 150, 37, 72, 240, 128, 108, 207, 144, 63, 47, 12, 200, 143, 157, 105, 170, 206, 28, 234, 51, 133, 125, 195, 120, 181, 49}, {224, 166, 191, 152, 214, 83, 193, 177, 18, 218, 115, 163, 254, 249, 177, 93, 64, 15, 84, 6, 93, 228, 228, 212, 150, 38, 228, 165, 164, 82, 134, 105}, {202, 70, 239, 186, 105, 110, 156, 80, 234, 234, 6, 50, 18, 78, 77, 254, 154, 79, 235, 29, 3, 67, 80, 176, 181, 74, 41, 108, 8, 122, 108, 235}, {113, 201, 202, 169, 97, 118, 143, 223, 84, 95, 61, 135, 133, 146, 252, 102, 126, 242, 95, 253, 181, 101, 255, 197, 147, 223, 213, 137, 228, 101, 65, 147}, {185, 36, 221, 167, 34, 147, 201, 8, 2, 59, 68, 72, 69, 65, 210, 207, 104, 197, 243, 210, 25, 31, 217, 174, 97, 123, 107, 27, 178, 171, 21, 203}, {209, 87, 77, 70, 236, 113, 65, 25, 92, 52, 184, 110, 227, 78, 101, 44, 90, 239, 35, 143, 20, 239, 214, 248, 224, 247, 214, 44, 206, 48, 32, 96}, {136, 99, 80, 39, 158, 193, 108, 216, 45, 122, 237, 173, 55, 64, 148, 14, 205, 21, 28, 55, 235, 151, 3, 248, 219, 60, 150, 52, 224, 221, 203, 99}, {119, 152, 76, 113, 140, 176, 227, 71, 125, 192, 235, 105, 186, 82, 52, 47, 195, 158, 32, 199, 70, 187, 200, 204, 242, 29, 152, 162, 70, 210, 45, 69}, {29, 166, 227, 57, 217, 21, 75, 20, 149, 18, 172, 132, 161, 181, 108, 148, 86, 192, 14, 223, 16, 79, 122, 112, 125, 85, 176, 21, 217, 192, 168, 81}, {27, 95, 4, 95, 226, 143, 174, 134, 21, 178, 130, 122, 177, 121, 39, 98, 148, 31, 101, 224, 19, 119, 163, 91, 113, 110, 191, 64, 173, 251, 160, 226}, {106, 107, 74, 141, 70, 34, 220, 159, 26, 113, 65, 18, 89, 213, 189, 86, 34, 110, 180, 44, 118, 21, 167, 169, 132, 105, 248, 160, 245, 35, 178, 230}, {97, 166, 157, 88, 88, 230, 93, 7, 22, 72, 196, 201, 74, 219, 120, 161, 161, 228, 24, 115, 184, 196, 149, 46, 145, 48, 114, 251, 16, 106, 47, 40}, {3, 84, 157, 225, 113, 214, 148, 10, 183, 25, 156, 116, 73, 129, 130, 174, 24, 232, 232, 121, 32, 111, 14, 172, 202, 67, 95, 255, 239, 220, 165, 147}, {83, 93, 3, 192, 14, 4, 125, 194, 21, 113, 74, 46, 126, 151, 52, 253, 221, 126, 205, 241, 28, 25, 73, 122, 164, 162, 101, 83, 19, 14, 234, 229}},
					{{243, 150, 20,  241, 234, 201, 39,  58,  153, 56,  184, 9,   121, 64, 203, 191, 82,  125, 5,   124, 119, 233, 193, 42, 162, 143, 243, 38, 195, 137, 57,  61}, {200, 5,  206, 198, 126, 153, 240, 84,  181, 81,  246, 152, 15, 209, 229, 42,  59,  39,  20,  180, 142, 9,  165, 221, 125, 130, 109, 154, 113, 43,  54, 103}, {53, 165, 221, 97, 55,  148, 35, 117, 39,  154, 249, 76, 139, 178, 191, 70,  45,  65, 132, 33,  250, 67,  24, 161, 121, 151, 245, 212, 41,  165, 248, 246}, {35,  88,  112, 89,  62, 40, 108, 157, 73,  212, 129, 155, 5,  236, 183, 247, 77,  153, 247, 155, 41,  196, 139, 146, 187, 31, 191, 181, 206, 131, 32,  203}, {253, 109, 230, 96, 51,  164, 249, 35, 242, 218, 135, 216, 114, 121, 82,  127, 10, 220, 223, 98, 45,  106, 250, 0,   211, 255, 116, 137, 130, 141, 85,  9},   {86,  135, 41, 251, 20,  10, 82,  121, 209, 214, 55,  12, 183, 249, 66,  119, 186, 227, 126, 202, 158, 38, 225, 147, 234, 86, 157, 186, 10, 34,  73,  135}, {247, 209, 152, 99,  210, 141, 45, 136, 81, 30, 87, 64, 130, 32, 129, 134, 231, 202, 239, 82, 34, 102, 164, 89, 133, 207, 104, 147, 219, 109, 157, 16}, {251, 215, 8, 166, 95, 27, 72, 144, 25, 165, 18, 36, 23, 190, 101, 122, 217, 192, 23, 134, 51, 51, 173, 54, 117, 36, 56, 215, 159, 86, 41, 60}, {82, 206, 177, 168, 165, 31, 108, 65, 36, 180, 62, 115, 23, 49, 128, 86, 86, 161, 241, 165, 149, 0, 17, 183, 58, 18, 186, 51, 194, 203, 65, 233}, {200, 194, 169, 149, 242, 241, 73, 245, 117, 112, 41, 143, 123, 0, 36, 149, 83, 83, 143, 102, 97, 214, 129, 247, 41, 186, 23, 16, 30, 91, 18, 98}, {48, 115, 44, 194, 22, 112, 11, 191, 102, 116, 221, 154, 40, 34, 84, 42, 187, 224, 46, 81, 11, 118, 2, 81, 37, 146, 26, 39, 253, 244, 90, 189}, {173, 195, 128, 116, 23, 131, 210, 70, 229, 153, 172, 204, 117, 217, 174, 125, 149, 79, 165, 222, 9, 87, 58, 224, 14, 238, 60, 7, 212, 196, 45, 64}, {27, 128, 3, 55, 183, 153, 176, 122, 212, 98, 214, 26, 115, 156, 92, 181, 224, 14, 10, 57, 111, 40, 194, 173, 9, 170, 31, 150, 136, 191, 84, 215}, {253, 68, 139, 23, 54, 70, 178, 231, 28, 212, 229, 29, 158, 186, 103, 33, 8, 102, 41, 251, 111, 202, 235, 120, 200, 103, 255, 187, 10, 184, 156, 31}, {83, 239, 213, 245, 244, 201, 215, 203, 201, 38, 248, 166, 190, 46, 252, 137, 107, 102, 152, 153, 211, 106, 143, 30, 188, 176, 202, 68, 198, 203, 95, 171}, {76, 22, 63, 172, 92, 247, 200, 98, 215, 109, 113, 235, 172, 163, 242, 254, 238, 219, 62, 67, 75, 4, 216, 138, 191, 98, 185, 78, 18, 14, 100, 52}, {49, 169, 138, 32, 104, 50, 187, 247, 225, 140, 64, 119, 132, 113, 28, 150, 3, 43, 56, 59, 69, 176, 22, 250, 246, 254, 7, 159, 212, 197, 240, 244}, {211, 8, 218, 192, 38, 60, 111, 56, 226, 37, 46, 207, 229, 209, 11, 158, 3, 142, 183, 99, 90, 65, 206, 5, 4, 2, 133, 20, 173, 117, 107, 134}, {72, 177, 224, 19, 166, 240, 157, 47, 115, 33, 181, 194, 2, 48, 79, 192, 83, 204, 5, 143, 227, 91, 102, 211, 186, 180, 61, 215, 251, 127, 191, 159}, {232, 100, 222, 17, 227, 148, 20, 101, 252, 71, 23, 118, 190, 215, 105, 61, 91, 175, 223, 217, 245, 6, 32, 73, 51, 5, 153, 238, 207, 152, 155, 89}, {123, 122, 170, 169, 178, 28, 182, 42, 22, 141, 19, 136, 10, 143, 13, 108, 229, 28, 67, 236, 15, 22, 111, 66, 44, 200, 124, 106, 9, 23, 138, 21}, {254, 131, 85, 194, 222, 0, 151, 199, 200, 115, 12, 237, 169, 33, 97, 235, 63, 20, 37, 231, 188, 109, 97, 233, 104, 213, 220, 123, 33, 80, 86, 86}, {216, 172, 164, 202, 8, 23, 50, 56, 108, 107, 234, 184, 139, 80, 13, 154, 37, 79, 117, 117, 47, 118, 37, 136, 120, 115, 30, 160, 243, 16, 229, 48}, {78, 219, 254, 237, 36, 158, 161, 45, 246, 104, 191, 197, 117, 27, 124, 1, 45, 134, 152, 15, 196, 145, 192, 131, 114, 86, 185, 24, 67, 37, 113, 222}, {91, 215, 42, 253, 102, 10, 185, 21, 121, 199, 64, 212, 1, 30, 40, 45, 241, 91, 147, 20, 65, 248, 49, 134, 215, 105, 16, 236, 215, 111, 117, 221}, {163, 184, 137, 234, 221, 6, 25, 158, 188, 199, 150, 170, 114, 124, 204, 140, 126, 79, 7, 217, 230, 187, 136, 233, 159, 79, 199, 19, 255, 212, 95, 161}, {5, 92, 168, 45, 124, 74, 155, 59, 212, 60, 6, 204, 28, 47, 210, 97, 252, 239, 236, 38, 87, 3, 121, 29, 64, 226, 153, 249, 46, 136, 243, 199}, {50, 102, 110, 162, 35, 133, 196, 139, 237, 45, 206, 231, 83, 126, 5, 232, 11, 145, 119, 246, 188, 6, 72, 34, 76, 48, 11, 105, 234, 230, 76, 11}, {243, 153, 74, 138, 62, 99, 128, 138, 142, 123, 163, 97, 131, 95, 61, 82, 247, 156, 218, 79, 89, 90, 185, 94, 151, 149, 104, 156, 186, 148, 106, 10}, {137, 250, 75, 205, 222, 43, 166, 10, 151, 114, 108, 149, 22, 136, 141, 25, 67, 10, 175, 20, 159, 55, 207, 94, 166, 218, 53, 18, 252, 213, 128, 73}, {90, 205, 75, 93, 83, 84, 239, 148, 89, 246, 72, 14, 47, 234, 85, 97, 229, 109, 47, 106, 229, 167, 29, 39, 23, 38, 153, 195, 72, 204, 61, 146}, {247, 218, 227, 122, 232, 204, 31, 22, 30, 164, 199, 235, 207, 151, 202, 89, 241, 197, 22, 43, 127, 40, 182, 125, 103, 50, 169, 219, 146, 32, 59, 101}}};

	for(int c=0;c<icandidates;c++)
	{
		int i,j,pmatrix;

		i=candidates[c].i;
		j=candidates[c].j;
		pmatrix=candidates[c].pmatrix;

		hash^=bytes[pmatrix][i][j];
		hash*=0x100000001b3;
	}

	int chosen_candidate=hash%icandidates;
	int qtype=pmatrix_entry_type(candidates[chosen_candidate].i,candidates[chosen_candidate].j);

	/*
		...and then we sum over all possible values of the chosen candidate label.
	*/

	for(int value=1;value<((qtype==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));value++)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		for(int i=0;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[0], i, j);
				twin->pmxs[1]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[1], i, j);
			}
		}

		pmatrix_set_entry(twin->pmxs[candidates[chosen_candidate].pmatrix],
				  candidates[chosen_candidate].i,
				  candidates[chosen_candidate].j,
				  value);

		collection_add(acl, twin, true);
	}

	/*
		S1 <-> S3 symmetry
	*/

	struct amatrix_backup_t backup;
	amatrix_save(amx,&backup);

	bool candidate_found=false;

#define FIRST_HALF(x)	(((x)==0)||((x)==1))
#define SECOND_HALF(x)	(((x)==2)||((x)==3))

	for(int i=0;i<amx->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<amx->pmxs[0]->dimensions;j++)
		{
			if(FIRST_HALF(i)&&FIRST_HALF(j))
			{
				amx->pmxs[0]->values[i][j]=backup.values[0][i][j];
				amx->pmxs[1]->values[i][j]=backup.values[1][i][j];

				if((candidate_found==false)&&(i==candidates[chosen_candidate].i)&&(j==candidates[chosen_candidate].j))
				{
					candidates[chosen_candidate].i=i;
					candidates[chosen_candidate].j=j;
					candidate_found=true;
				}
			}

			if(FIRST_HALF(i)&&SECOND_HALF(j))
			{
				int jprime=(j==3)?(2):(3);

				amx->pmxs[0]->values[i][j]=backup.values[0][i][jprime];
				amx->pmxs[1]->values[i][j]=backup.values[1][i][jprime];

				if((i==candidates[chosen_candidate].i)&&(jprime==candidates[chosen_candidate].j))
				{
					candidates[chosen_candidate].i=i;
					candidates[chosen_candidate].j=j;
					candidate_found=true;
				}
			}

			if(SECOND_HALF(i)&&FIRST_HALF(j))
			{
				int iprime=(i==3)?(2):(3);

				amx->pmxs[0]->values[i][j]=backup.values[0][iprime][j];
				amx->pmxs[1]->values[i][j]=backup.values[1][iprime][j];

				if((iprime==candidates[chosen_candidate].i)&&(j==candidates[chosen_candidate].j))
				{
					candidates[chosen_candidate].i=i;
					candidates[chosen_candidate].j=j;
					candidate_found=true;
				}
			}

			if(SECOND_HALF(i)&&SECOND_HALF(j))
			{
				int iprime=(i==3)?(2):(3);
				int jprime=(j==3)?(2):(3);

				amx->pmxs[0]->values[i][j]=backup.values[0][iprime][jprime];
				amx->pmxs[1]->values[i][j]=backup.values[1][iprime][jprime];

				if((iprime==candidates[chosen_candidate].i)&&(jprime==candidates[chosen_candidate].j))
				{
					candidates[chosen_candidate].i=i;
					candidates[chosen_candidate].j=j;
					candidate_found=true;
				}
			}
		}
	}

	assert(candidate_found==true);

	/*
		Remember that 'qtype' could have changed applying the symmetry.
	*/

	qtype=pmatrix_entry_type(candidates[chosen_candidate].i,candidates[chosen_candidate].j);

	for(int value=1;value<((qtype==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));value++)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		for(int i=0;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[0], i, j);
				twin->pmxs[1]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[1], i, j);
			}
		}

		pmatrix_set_entry(twin->pmxs[candidates[chosen_candidate].pmatrix],
				  candidates[chosen_candidate].i,
				  candidates[chosen_candidate].j,
				  value);

		collection_add(acl, twin, true);
	}

	amatrix_restore(amx,&backup);
}

void cdet_generate_grouping_altQ(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	int dimensions=amx->pmxs[0]->dimensions;
	int target_type=QTYPE_VIRTUAL;

	struct
	{
		int i, j, pmatrix;
	} candidates[32];

	int icandidates=0;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_entry_type(i, j)==target_type)
			{
				if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
				{
					candidates[icandidates].i=i;
					candidates[icandidates].j=j;
					candidates[icandidates].pmatrix=0;
					icandidates++;
				}

				if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
				{
					candidates[icandidates].i=i;
					candidates[icandidates].j=j;
					candidates[icandidates].pmatrix=1;
					icandidates++;
				}
			}
		}
	}

	if((icandidates<2)||(icandidates>6))
		return;

	assert(icandidates<32);

	/*
		The twins
	*/

	int values[32];

	for(int c=0;c<icandidates;c++)
		values[c]=pmatrix_get_raw_entry(amx->pmxs[candidates[c].pmatrix],candidates[c].i,candidates[c].j);

	for(int c=1;c<ifactorial(icandidates);c++)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		for(int i=0;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[0], i, j);
				twin->pmxs[1]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[1], i, j);
			}
		}

		for(int d=0;d<icandidates;d++)
		{
			int dprime=get_permutation(icandidates,c,d)-1;

			struct pmatrix_t *pprime=twin->pmxs[candidates[d].pmatrix];
			int iprime=candidates[d].i;
			int jprime=candidates[d].j;
			int value=values[dprime];

			pmatrix_set_raw_entry(pprime, iprime, jprime, value);
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altQplus(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	/*
		The original adjacency matrix
	*/

	collection_add(acl, amx, false);

	int dimensions=amx->pmxs[0]->dimensions;
	int target_type=QTYPE_VIRTUAL;

	struct
	{
		int i, j, pmatrix;
	} candidates[32];

	int icandidates=0;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_entry_type(i, j)==target_type)
			{
				if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
				{
					candidates[icandidates].i=i;
					candidates[icandidates].j=j;
					candidates[icandidates].pmatrix=0;
					icandidates++;
				}

				if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
				{
					candidates[icandidates].i=i;
					candidates[icandidates].j=j;
					candidates[icandidates].pmatrix=1;
					icandidates++;
				}
			}
		}
	}

	if((icandidates<2)||(icandidates>6))
		return;

	assert(icandidates<32);

	/*
		The twins, part 1
	*/

	int values[32];

	for(int c=0;c<icandidates;c++)
		values[c]=pmatrix_get_raw_entry(amx->pmxs[candidates[c].pmatrix],candidates[c].i,candidates[c].j);

	for(int c=1;c<ifactorial(icandidates);c++)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

		for(int i=0;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				twin->pmxs[0]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[0], i, j);
				twin->pmxs[1]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[1], i, j);
			}
		}

		for(int d=0;d<icandidates;d++)
		{
			int dprime=get_permutation(icandidates,c,d)-1;

			struct pmatrix_t *pprime=twin->pmxs[candidates[d].pmatrix];
			int iprime=candidates[d].i;
			int jprime=candidates[d].j;
			int value=values[dprime];

			pmatrix_set_raw_entry(pprime, iprime, jprime, value);
		}

		collection_add(acl, twin, true);
	}

	/*
		The twins, part 2 (S1 <-> S3 symmetry)
	*/

	for(int c=1;c<ifactorial(icandidates);c++)
	{
		struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

		twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
		twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

#define FIRST_HALF(x)	(((x)==0)||((x)==1))
#define SECOND_HALF(x)	(((x)==2)||((x)==3))

		for(int i=0;i<dimensions;i++)
		{
			for(int j=0;j<dimensions;j++)
			{
				if(FIRST_HALF(i)&&FIRST_HALF(j))
				{
					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][j];
				}

				if(FIRST_HALF(i)&&SECOND_HALF(j))
				{
					int jprime=(j==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][jprime];
				}

				if(SECOND_HALF(i)&&FIRST_HALF(j))
				{
					int iprime=(i==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][j];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][j];
				}

				if(SECOND_HALF(i)&&SECOND_HALF(j))
				{
					int iprime=(i==3)?(2):(3);
					int jprime=(j==3)?(2):(3);

					twin->pmxs[0]->values[i][j]=amx->pmxs[0]->values[iprime][jprime];
					twin->pmxs[1]->values[i][j]=amx->pmxs[1]->values[iprime][jprime];
				}
			}
		}

		for(int d=0;d<icandidates;d++)
		{
			int dprime=get_permutation(icandidates,c,d)-1;

			struct pmatrix_t *pprime=twin->pmxs[candidates[d].pmatrix];
			int iprime=candidates[d].i;
			int jprime=candidates[d].j;
			int value=values[dprime];

			pmatrix_set_raw_entry(pprime, iprime, jprime, value);
		}

		collection_add(acl, twin, true);
	}
}

void cdet_generate_grouping_altL(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	int dimensions=amx->pmxs[0]->dimensions;

	struct
	{
		int i, j, pmatrix;
	} candidates[32];

	int icandidates=0;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
			{
				candidates[icandidates].i=i;
				candidates[icandidates].j=j;
				candidates[icandidates].pmatrix=0;
				icandidates++;
			}

			if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
			{
				candidates[icandidates].i=i;
				candidates[icandidates].j=j;
				candidates[icandidates].pmatrix=1;
				icandidates++;
			}
		}
	}

	if(icandidates<2)
	{
		collection_add(acl, amx, false);
		return;
	}

	int qtypes[2];

	qtypes[0]=pmatrix_entry_type(candidates[0].i,candidates[0].j);
	qtypes[1]=pmatrix_entry_type(candidates[1].i,candidates[1].j);

	/*
		...and then we sum over all possible values of the chosen candidate labels.
	*/

	int values[2];

	for(values[0]=0;values[0]<((qtypes[0]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[0]++)
	{
		for(values[1]=0;values[1]<((qtypes[1]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[1]++)
		{
			if((values[0]==0)&&(values[1]==0))
				continue;

			struct amatrix_t *twin=init_amatrix_from_amatrix(amx);

			twin->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
			twin->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

			for(int i=0;i<dimensions;i++)
			{
				for(int j=0;j<dimensions;j++)
				{
					twin->pmxs[0]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[0], i, j);
					twin->pmxs[1]->values[i][j]=pmatrix_get_raw_entry(amx->pmxs[1], i, j);
				}
			}

			pmatrix_set_entry(twin->pmxs[candidates[0].pmatrix],
					  candidates[0].i,
					  candidates[0].j,
					  values[0]);

			pmatrix_set_entry(twin->pmxs[candidates[1].pmatrix],
					  candidates[1].i,
					  candidates[1].j,
					  values[1]);

			collection_add(acl, twin, true);
		}
	}
}

double amatrix_projection_multiplicity(struct amatrix_t *amx)
{
	int nr_occupied_entries,nr_virtual_entries;

	nr_occupied_entries=nr_virtual_entries=0;
	for(int i=0;i<amx->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<amx->pmxs[0]->dimensions;j++)
		{
			if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
			{
				if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
					nr_occupied_entries++;
				else
					nr_virtual_entries++;
			}

			if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
			{
				if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
					nr_occupied_entries++;
				else
					nr_virtual_entries++;
			}
		}
	}

	return pow(amx->nr_virtual,nr_occupied_entries)*pow(amx->nr_occupied,nr_virtual_entries);
}

/*
	Here we calculate the weight associated to a 'amatrix'
*/

double amatrix_weight(struct amatrix_t *amx)
{
	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	if(amx->cached_weight_is_valid==true)
	{

#ifndef NDEBUG
		amx->cached_weight_is_valid=false;

		double w1,w2;

		w1=amx->cached_weight;
		w2=amatrix_weight(amx);
		assert(gsl_fcmp(w1,w2,1e-6)==0);

		amx->cached_weight_is_valid=true;
#endif

		return amx->cached_weight;
	}

	/*
		Dimension 1 is a special case that does not need the evaluation
		of the incidence matrix.

		Note that dimension 1 (a bit unexpectedly) corresponds to the Hartree-Fock energy.
	*/

	if(amx->pmxs[0]->dimensions==1)
	{
		int a,b;
		double weight=0.0f,multiplicity;

		a=pmatrix_get_entry(amx->pmxs[0],0,0);
		b=pmatrix_get_entry(amx->pmxs[1],0,0);

		if(a==b)
			weight+=get_hdiag(amx->ectx,a-1);

		weight+=0.5*get_eri(amx->ectx,a-1,b-1,a-1,b-1);
		weight+=get_enuc(amx->ectx)*pow(amx->nr_occupied,-2.0f);

		multiplicity=1.0f;

		amx->cached_weight=weight/multiplicity/amatrix_projection_multiplicity(amx);
		amx->cached_weight_is_valid=true;

		return amx->cached_weight;
	}
	else
	{
		struct amatrix_collection_t *acl=init_collection();

		double ret,combinatorial,Rdenominator;
		int nr_weights;

		cdet_generate_grouping_altL(amx, acl);

		combinatorial=1.0f;
		ret=Rdenominator=0.0f;
		nr_weights=0;
		for(int c=0;c<acl->iamatrices;c++)
		{
			assert(amatrix_check_consistency(acl->amatrices[c])==true);

			struct label_t labels[MAX_LABELS];
			int ilabels=0;

			double thisweight;

			if(amatrix_check_connectedness(acl->amatrices[c])==true)
			{
				gsl_matrix_int *incidence=amatrix_calculate_incidence(acl->amatrices[c], labels, &ilabels);
				thisweight=incidence_to_weight(incidence, labels, &ilabels, acl->amatrices[c]);
				gsl_matrix_int_free(incidence);

				thisweight/=amatrix_projection_multiplicity(amx);
			}
			else
			{
				thisweight=0.0f;
			}

			ret+=thisweight;
			Rdenominator+=fabs(thisweight);
			nr_weights++;
		}

		fini_collection(acl);

		if(false)
		{
			double R=fabs(ret)/Rdenominator;

			{
				//char *banner=(amatrix_is_physical(amx)==true) ? ("physical") : ("unphysical");
				//printf("R=%e, sum=%e [%s, order=%d]\n", R, ret, banner, amx->pmxs[0]->dimensions);
			}
		}

		for(int c=0;c<acl->iamatrices;c++)
			for(int d=0;c<d;d++)
				if(amatrix_ordering(acl->amatrices[c],acl->amatrices[d])==ORDERING_EQ)
					combinatorial*=2.0f;

		if(nr_weights==0)
		{
			amx->cached_weight=0.0f;
			amx->cached_weight_is_valid=true;
		}
		else
		{
			amx->cached_weight=ret/combinatorial/((double)(nr_weights));
			amx->cached_weight_is_valid=true;
		}

		return amx->cached_weight;
	}

	assert(false);
	return 0.0f;
}

void cdet_debug(struct configuration_t *cfg)
{
	struct amatrix_t *amx=init_amatrix(cfg);

	amx->pmxs[0]->dimensions=4;
	amx->pmxs[1]->dimensions=4;

	for(int i=0;i<amx->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<amx->pmxs[0]->dimensions;j++)
		{
			amx->pmxs[0]->values[i][j]=0;
			amx->pmxs[1]->values[i][j]=0;

#define INDEX_MATCH(a,b)	(((a)==i)&&((b)==j))

			if(INDEX_MATCH(1,0)||INDEX_MATCH(0,1)||INDEX_MATCH(3,2)||INDEX_MATCH(2,3))
				amx->pmxs[0]->values[i][j]=1;

			if(INDEX_MATCH(1,0)||INDEX_MATCH(0,2)||INDEX_MATCH(3,1)||INDEX_MATCH(2,3))
				amx->pmxs[1]->values[i][j]=1;
		}
	}

	if(amatrix_check_consistency(amx)==false)
	{
		printf("Consistency check failed!\n");
		goto cleanup;
	}

	struct amatrix_collection_t *acl=init_collection();

	cdet_generate_grouping_altSx(amx, acl);

	for(int c=0;c<acl->iamatrices;c++)
	{
		if(amatrix_check_consistency(acl->amatrices[c])==false)
		{
			printf("Consistency check failed!\n");
			goto cleanup;
		}

		amatrix_print(acl->amatrices[c]);
		printf("\n");
	}

	fini_collection(acl);

	cleanup:

	fini_amatrix(amx,true);
}
