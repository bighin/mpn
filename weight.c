#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "weight.h"
#include "mpn.h"

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

	int flags=0x07;

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

#undef NO_CDET_AT_ALL
#ifdef NO_CDET_AT_ALL
	{
		struct label_t labels[MAX_LABELS];
		int ilabels=0;

		if(amatrix_check_connectedness(amx)==false)
			return 0.0f;

		gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
		double thisweight=incidence_to_weight(incidence, labels, &ilabels, amx);
		gsl_matrix_int_free(incidence);

		amx->cached_weight=thisweight/get_mul2(amx);
		amx->cached_weight_is_valid=true;

		return amx->cached_weight;
	}
#else
	{
		struct amatrix_collection_t *acl=init_collection();

		double ret,combinatorial,Rdenominator;
		int nr_weights;

		cdet_generate_grouping_alt(amx, acl);

		//printf("[%d] ",amx->pmxs[0]->dimensions);

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

			//printf("%f ",thisweight);

			ret+=thisweight;
			Rdenominator+=fabs(thisweight);
			nr_weights++;
		}

		//printf("\n");

		fini_collection(acl);

		if((fabs(ret)>1e-6)&&(Rdenominator>1e-6)&&(amx->pmxs[0]->dimensions>3)&&(nr_weights>=2))
		{
			double R=fabs(ret)/Rdenominator;

			//char *banner=(amatrix_is_physical(amx)==true)?("physical"):("unphysical");
			//printf("%f [%s, order=%d]\n", R, banner, amx->pmxs[0]->dimensions);
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
#endif

	assert(false);
	return 0.0f;
}
