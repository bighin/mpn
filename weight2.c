#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "weight2.h"
#include "mpn.h"
#include "multiplicity.h"
#include "auxx.h"

void add_denominator_entry(struct weight_info_t *awt, int label, int qtype)
{
	int index=awt->denominators[awt->nr_denominators].ilabels;

	awt->denominators[awt->nr_denominators].labels[index]=label;
	awt->denominators[awt->nr_denominators].qtypes[index]=qtype;
	awt->denominators[awt->nr_denominators].ilabels++;
}

void add_numerator(struct weight_info_t *awt, int l1, int l2, int l3, int l4)
{
	awt->numerators[awt->nr_numerators].labels[0]=l1;
	awt->numerators[awt->nr_numerators].labels[1]=l2;
	awt->numerators[awt->nr_numerators].labels[2]=l3;
	awt->numerators[awt->nr_numerators].labels[3]=l4;
	awt->nr_numerators++;
}

void add_unphysical_penalty(struct weight_info_t *awt, double penalty)
{
	awt->unphysical_penalty*=penalty;
}

/*
	This function calculates a diagram's weight given the incidence matrix
*/

struct weight_info_t incidence_to_weight_info(gsl_matrix_int *B, struct label_t *labels, int *ilabels, struct amatrix_t *amx)
{
	int mels[MAX_MATRIX_ELEMENTS][4];
	assert(B->size1<=MAX_MATRIX_ELEMENTS);

	double inversefactor=1.0f;

	/*
		Rule 2: each row in the incidence matrix corresponds to a matrix element.
	*/

	for(size_t i=0;i<B->size1;i++)
	{
		bool has_selfloop=false;

		for(int l=0;l<4;l++)
			mels[i][l]=-1;

		for(size_t j=0;j<B->size2;j++)
		{
			if(gsl_matrix_int_get(B, i, j)==1)
			{
				if(mels[i][0]==-1)
					mels[i][0]=j;
				else if(mels[i][1]==-1)
					mels[i][1]=j;
				else
					assert(false);
			}

			if(gsl_matrix_int_get(B, i, j)==-1)
			{
				if(mels[i][2]==-1)
					mels[i][2]=j;
				else if(mels[i][3]==-1)
					mels[i][3]=j;
				else
					assert(false);
			}

			if(gsl_matrix_int_get(B, i, j)==2)
				has_selfloop=true;
		}

		if(has_selfloop==false)
			for(int l=0;l<4;l++)
				assert(mels[i][l]!=-1);
	}

	/*
		Rule 3: identical columns give rise to a factorial factor
	*/

	for(size_t i=0;i<B->size2;i++)
	{
		int cnt=0;

		for(size_t j=i+1;j<B->size2;j++)
			if(columns_are_identical(B, i, j)==true)
				cnt++;

		if(cnt!=0)
			inversefactor*=factorial(cnt+1);
	}

	/*
		This structure will contain the final calculated weight, along with plenty of
		additional information, allowing one to reconstruct the weight from scratch.
	*/

	struct weight_info_t ret;

	ret.nr_denominators=0;
	ret.nr_numerators=0;
	ret.unphysical_penalty=1.0f;

	/*
		Rule 4: denominators
	*/

	double denominators=1.0f;
	int excitation_level=0;

	for(size_t i=0;i<(B->size1-1);i++)
	{
		double denominator=0.0f;
		int energies_in_denominator=0;

		ret.denominators[ret.nr_denominators].ilabels=0;

		for(size_t j=0;j<B->size2;j++)
		{
			int sum=0;

			for(size_t k=0;k<=i;k++)
				sum+=gsl_matrix_int_get(B,k,j);

			if(sum!=0)
			{
				/*
					The j-th label contributes to the i-th denominator
				*/

				add_denominator_entry(&ret, j, labels[j].qtype);

				switch(labels[j].qtype)
				{
					case QTYPE_OCCUPIED:
					denominator+=get_occupied_energy(amx->ectx,labels[j].value-1);
					energies_in_denominator++;
					break;

					case QTYPE_VIRTUAL:
					denominator-=get_virtual_energy(amx->ectx,labels[j].value-1);
					energies_in_denominator++;
					break;

					default:
					assert(false);
				}
			}
		}

		if(energies_in_denominator>0)
		{
			denominators*=denominator;
			ret.nr_denominators++;

			excitation_level=MAX(excitation_level,energies_in_denominator/2);
		}
	}

	ret.excitation_level=excitation_level;

	/*
		Additional rule: phase factor
	*/

	int h,l;

	l=count_loops(labels, ilabels, mels, B->size1);

	h=0;
	for(int i=0;i<*ilabels;i++)
		if(labels[i].qtype==QTYPE_OCCUPIED)
			h++;

	inversefactor*=pow(-1.0f,l+h);

	/*
		Finally we build up the total weight.
	*/

	double numerators=1.0f;

	for(size_t i=0;i<B->size1;i++)
	{
		/*
			If these four entries are not set, it means that the i-th line of the adjacency matrix
			contains a '2' entry, coming from the unphysical sector.

			Therefore, it cannot give rise to a proper matrix element.
		*/

		if((mels[i][0]==-1)||(mels[i][1]==-1)||(mels[i][2]==-1)||(mels[i][3]==-1))
		{
			/*
				We assign an unphysical penalty for 'selfloops', i.e. edges on the
				graph connecting a vertex with itself, appearing only in the unphysical sector.
			*/

			numerators*=amx->config->unphysicalpenalty;
			add_unphysical_penalty(&ret,amx->config->unphysicalpenalty);

			continue;
		}

		/*
			We have to convert the quantum numbers into the format used by get_eri()
		*/

		int i1,i2,i3,i4;

		i1=labels[mels[i][0]].value-1+((labels[mels[i][0]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i2=labels[mels[i][1]].value-1+((labels[mels[i][1]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i3=labels[mels[i][2]].value-1+((labels[mels[i][2]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i4=labels[mels[i][3]].value-1+((labels[mels[i][3]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));

		numerators*=get_eri(amx->ectx, i1, i2, i3, i4);

		add_numerator(&ret,mels[i][0],mels[i][1],mels[i][2],mels[i][3]);
	}

	/*
		Finally we put all the results into an 'weight_info_t' struct and
		we return it.
	*/

	memcpy(ret.labels,labels,sizeof(struct label_t)*MAX_LABELS);
	ret.ilabels=*ilabels;

	ret.l=l;
	ret.h=h;
	ret.inversefactor=inversefactor;
	ret.weight=pow(inversefactor,-1.0f)*numerators/denominators/amatrix_multiplicity(amx);

#ifndef NDEBUG
	{
		double w1=reconstruct_weight(amx,&ret);

		/*
			We have to reset this, since count_loops() will be invoked again in incidence_to_weight()
		*/

		for(int c=0;c<*ilabels;c++)
			labels[c].visited=false;

		double w2=incidence_to_weight(B,labels,ilabels,amx);

		assert(gsl_fcmp(w1,w2,1e-6)==0);
	}
#endif

	return ret;
}

double reconstruct_weight(struct amatrix_t *amx, struct weight_info_t *awt)
{
	/*
		Keep in mind that here we do not check for connectedness.
	*/

	double denominators=1.0f;

	for(int c=0;c<awt->nr_denominators;c++)
	{
		double denominator=0.0f;

		for(int d=0;d<awt->denominators[c].ilabels;d++)
		{
			int label=awt->denominators[c].labels[d];

			switch(awt->denominators[c].qtypes[d])
			{
				case QTYPE_OCCUPIED:
				denominator+=get_occupied_energy(amx->ectx, awt->labels[label].value-1);
				break;

				case QTYPE_VIRTUAL:
				denominator-=get_virtual_energy(amx->ectx, awt->labels[label].value-1);
				break;
			}
		}

		denominators*=denominator;
	}

	double numerators=1.0f;

	for(int c=0;c<awt->nr_numerators;c++)
	{
		int l1,l2,l3,l4;

		l1=awt->numerators[c].labels[0];
		l2=awt->numerators[c].labels[1];
		l3=awt->numerators[c].labels[2];
		l4=awt->numerators[c].labels[3];

		int i1,i2,i3,i4;

		i1=awt->labels[l1].value-1+((awt->labels[l1].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i2=awt->labels[l2].value-1+((awt->labels[l2].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i3=awt->labels[l3].value-1+((awt->labels[l3].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i4=awt->labels[l4].value-1+((awt->labels[l4].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));

		numerators*=get_eri(amx->ectx, i1, i2, i3, i4);
	}

	numerators*=awt->unphysical_penalty;

	double weight=pow(awt->inversefactor,-1.0f)*numerators/denominators/amatrix_multiplicity(amx);

	/*
		Lindelöf resummation should happen here, if needed.
	*/

	return weight;
}

int coordinate_to_label_index(struct label_t *labels,int ilabels,int i,int j,int pmatrix)
{
	for(int c=0;c<ilabels;c++)
		if((labels[c].i==i)&&(labels[c].j==j)&&(labels[c].pmatrix==pmatrix))
			return c;

	assert(false);
	return 0;
}

int get_excitation_level(struct amatrix_t *amx)
{
	struct label_t labels[MAX_LABELS];
	int ilabels=0;

#warning This result could be cached with the weight

	int dimensions=amx->pmxs[0]->dimensions;

	if(dimensions==1)
		return 0;
	else if(dimensions==2)
		return 2;

	gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
	struct weight_info_t w=incidence_to_weight_info(incidence, labels, &ilabels, amx);
	gsl_matrix_int_free(incidence);

	return w.excitation_level;
}
