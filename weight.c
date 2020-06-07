#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix_int.h>

#include "amatrix.h"
#include "weight.h"
#include "mpn.h"
#include "permutations.h"
#include "auxx.h"
#include "weight2.h"

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

struct amatrix_collection_t
{
	struct amatrix_t **amatrices;
	struct weight_info_t **winfos;

	int iamatrices,nralloced;
};

struct amatrix_collection_t *init_collection(void)
{
	struct amatrix_collection_t *ret;

	if(!(ret=malloc(sizeof(struct amatrix_collection_t))))
		return NULL;

	ret->amatrices=NULL;
	ret->winfos=NULL;
	ret->iamatrices=0;
	ret->nralloced=0;

	return ret;
}

void fini_collection(struct amatrix_collection_t *acl)
{
	if(acl)
	{
		for(int c=0;c<acl->iamatrices;c++)
		{
			if(acl->amatrices[c])
				fini_amatrix(acl->amatrices[c], false);

			if(acl->winfos[c])
				free(acl->winfos[c]);
		}

		if(acl->amatrices)
			free(acl->amatrices);

		if(acl->winfos)
			free(acl->winfos);

		free(acl);
	}
}

/*
	A 'amatrix_t' struct is added to the amatrix collection and it will be freed
	when the collection will be destroyed.
*/

void collection_add(struct amatrix_collection_t *acl,struct amatrix_t *amx, struct weight_info_t *winfo)
{
	if(acl->iamatrices>=acl->nralloced)
	{
		assert(acl->iamatrices==acl->nralloced);

		acl->nralloced+=16;
		acl->amatrices=realloc(acl->amatrices,sizeof(struct amatrix_t *)*acl->nralloced);
		acl->winfos=realloc(acl->winfos,sizeof(struct weight_info_t *)*acl->nralloced);
	}

	acl->amatrices[acl->iamatrices]=amx;
	acl->winfos[acl->iamatrices]=NULL;

	if(winfo!=NULL)
	{
		acl->winfos[acl->iamatrices]=malloc(sizeof(struct weight_info_t));
		memcpy(acl->winfos[acl->iamatrices],winfo,sizeof(struct weight_info_t));
	}

	acl->iamatrices++;
}

/*
	The copy of an 'amatrix_t' struct is added to the amatrix collection; this means
	that the original amatrix will NOT be freed when the collection will be destroyed.
*/

void collection_add_copy(struct amatrix_collection_t *acl, struct amatrix_t *amx,struct weight_info_t *winfo)
{
	struct amatrix_t *copy=init_amatrix_from_amatrix(amx);

	copy->pmxs[0]->dimensions=amx->pmxs[0]->dimensions;
	copy->pmxs[1]->dimensions=amx->pmxs[1]->dimensions;

	assert(copy->pmxs[0]->dimensions==copy->pmxs[1]->dimensions);

	for(int i=0;i<copy->pmxs[0]->dimensions;i++)
	{
		for(int j=0;j<copy->pmxs[0]->dimensions;j++)
		{
			copy->pmxs[0]->values[i][j]=amx->pmxs[0]->values[i][j];
			copy->pmxs[1]->values[i][j]=amx->pmxs[1]->values[i][j];
		}
	}

	collection_add(acl, copy, winfo);
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

void cdet_generate_grouping_altL0(struct amatrix_t *amx, struct amatrix_collection_t *acl)
{
	collection_add_copy(acl, amx, NULL);
}

void cdet_generate_grouping_altL1(struct amatrix_t *amx, struct amatrix_collection_t *acl)
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
		collection_add_copy(acl, amx, NULL);
		return;
	}

	struct weight_info_t winfo;

	{
		struct label_t labels[MAX_LABELS];
		int ilabels=0;

		gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
		winfo=incidence_to_weight_info(incidence, labels, &ilabels, amx);
		gsl_matrix_int_free(incidence);
	}

	int qtypes[1];
	int ilabels[1];

	qtypes[0]=pmatrix_entry_type(candidates[0].i,candidates[0].j);
	ilabels[0]=coordinate_to_label_index(winfo.labels,winfo.ilabels,candidates[0].i,candidates[0].j,candidates[0].pmatrix);

	/*
		...and then we sum over all possible values of the chosen candidate labels.
	*/

	int values[0];

	for(values[0]=1;values[0]<=((qtypes[0]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[0]++)
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

		pmatrix_set_entry(twin->pmxs[candidates[0].pmatrix],
				  candidates[0].i,
				  candidates[0].j,
				  values[0]);

		winfo.labels[ilabels[0]].value=values[0];

		collection_add(acl, twin, &winfo);
	}
}

void cdet_generate_grouping_altL2(struct amatrix_t *amx, struct amatrix_collection_t *acl)
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
		collection_add_copy(acl, amx, NULL);
		return;
	}

	struct weight_info_t winfo;

	{
		struct label_t labels[MAX_LABELS];
		int ilabels=0;

		gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
		winfo=incidence_to_weight_info(incidence, labels, &ilabels, amx);
		gsl_matrix_int_free(incidence);
	}

	int qtypes[2];
	int ilabels[2];

	qtypes[0]=pmatrix_entry_type(candidates[0].i,candidates[0].j);
	qtypes[1]=pmatrix_entry_type(candidates[1].i,candidates[1].j);

	ilabels[0]=coordinate_to_label_index(winfo.labels,winfo.ilabels,candidates[0].i,candidates[0].j,candidates[0].pmatrix);
	ilabels[1]=coordinate_to_label_index(winfo.labels,winfo.ilabels,candidates[1].i,candidates[1].j,candidates[1].pmatrix);

	/*
		...and then we sum over all possible values of the chosen candidate labels.
	*/

	int values[2];

	for(values[0]=1;values[0]<=((qtypes[0]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[0]++)
	{
		for(values[1]=1;values[1]<=((qtypes[1]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[1]++)
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

			pmatrix_set_entry(twin->pmxs[candidates[0].pmatrix],
					  candidates[0].i,
					  candidates[0].j,
					  values[0]);

			pmatrix_set_entry(twin->pmxs[candidates[1].pmatrix],
					  candidates[1].i,
					  candidates[1].j,
					  values[1]);

			winfo.labels[ilabels[0]].value=values[0];
			winfo.labels[ilabels[1]].value=values[1];

			collection_add(acl, twin, &winfo);
		}
	}
}

void cdet_generate_grouping_altL2x(struct amatrix_t *amx, struct amatrix_collection_t *acl)
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

	int first_candidate=-1,second_candidate=-1;

	for(int c=0;c<icandidates;c++)
	{
		if(pmatrix_entry_type(candidates[c].i,candidates[c].j)==QTYPE_OCCUPIED)
		{
			first_candidate=c;
			break;
		}
	}

	for(int c=0;c<icandidates;c++)
	{
		if(pmatrix_entry_type(candidates[c].i,candidates[c].j)==QTYPE_VIRTUAL)
		{
			second_candidate=c;
			break;
		}
	}

	if((first_candidate==-1)||(second_candidate==-1))
	{
		collection_add_copy(acl, amx, NULL);
		return;
	}

	struct weight_info_t winfo;

	{
		struct label_t labels[MAX_LABELS];
		int ilabels=0;

		gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
		winfo=incidence_to_weight_info(incidence, labels, &ilabels, amx);
		gsl_matrix_int_free(incidence);
	}

	int qtypes[2];
	int ilabels[2];

	qtypes[0]=pmatrix_entry_type(candidates[first_candidate].i,candidates[first_candidate].j);
	qtypes[1]=pmatrix_entry_type(candidates[second_candidate].i,candidates[second_candidate].j);

	ilabels[0]=coordinate_to_label_index(winfo.labels,winfo.ilabels,candidates[first_candidate].i,candidates[first_candidate].j,candidates[first_candidate].pmatrix);
	ilabels[1]=coordinate_to_label_index(winfo.labels,winfo.ilabels,candidates[second_candidate].i,candidates[second_candidate].j,candidates[second_candidate].pmatrix);

	/*
		...and then we sum over all possible values of the chosen candidate labels.
	*/

	int values[2];

	for(values[0]=1;values[0]<=((qtypes[0]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[0]++)
	{
		for(values[1]=1;values[1]<=((qtypes[1]==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied));values[1]++)
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

			pmatrix_set_entry(twin->pmxs[candidates[first_candidate].pmatrix],
					  candidates[first_candidate].i,
					  candidates[first_candidate].j,
					  values[0]);

			pmatrix_set_entry(twin->pmxs[candidates[second_candidate].pmatrix],
					  candidates[second_candidate].i,
					  candidates[second_candidate].j,
					  values[1]);

			winfo.labels[ilabels[0]].value=values[0];
			winfo.labels[ilabels[1]].value=values[1];

			collection_add(acl, twin, &winfo);
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

	bool use_cdet=true;

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
	else if(use_cdet)
	{
		struct amatrix_collection_t *acl=init_collection();

		double ret,combinatorial,Rdenominator;
		int nr_weights;

		cdet_generate_grouping_altL1(amx, acl);

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
				if(acl->winfos[c]==NULL)
				{
					gsl_matrix_int *incidence=amatrix_calculate_incidence(acl->amatrices[c], labels, &ilabels);
					thisweight=incidence_to_weight(incidence, labels, &ilabels, acl->amatrices[c]);
					gsl_matrix_int_free(incidence);
				}
				else
				{
					thisweight=reconstruct_weight(acl->amatrices[c],acl->winfos[c]);
				}

				thisweight/=amatrix_projection_multiplicity(acl->amatrices[c]);
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
	else
	{
		double ret=0;

		struct label_t labels[MAX_LABELS];
		int ilabels=0;

		if(amatrix_check_connectedness(amx)==true)
		{
			gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
			ret=incidence_to_weight(incidence, labels, &ilabels, amx);
			gsl_matrix_int_free(incidence);

			ret/=amatrix_projection_multiplicity(amx);
		}

		amx->cached_weight=ret;
		amx->cached_weight_is_valid=true;

		return amx->cached_weight;
	}

	assert(false);
	return 0.0f;
}
