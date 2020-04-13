#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "mpn.h"
#include "amatrix.h"
#include "pmatrix.h"
#include "auxx.h"
#include "loaderis.h"
#include "multiplicity.h"

/*
	For this function and the next one, see Szabo-Ostlund, page 360.

	This function could be optimized by writing down a lookup table when
	we write the matrix elements... Is this a bottleneck?
*/

int find_successor(int node,int mels[MAX_MATRIX_ELEMENTS][4],int imels)
{
	for(int i=0;i<imels;i++)
	{
		if(mels[i][0]==node)
			return mels[i][2];

		if(mels[i][1]==node)
			return mels[i][3];
	}

	assert(false);
	return 0;
}

int count_loops(struct label_t *labels, int *ilabels, int mels[MAX_MATRIX_ELEMENTS][4], int nrmels)
{
	int loops=0;

	for(int c=0;c<*ilabels;c++)
	{
		/*
			We don't consider 'selfloops', i.e. edges on the graph connecting a vertex
			with itself. These appear only in the unphysical sector.
		*/

		if((labels[c].visited==false)&&(labels[c].selfloop==false))
		{
			labels[c].visited=true;

			/*
				We just go through the graph following the bra/kets until we're back
				to where we started, as explained in Szabo-Ostlund, page 360.
			*/

			for(int node=find_successor(c,mels,nrmels);node!=c;node=find_successor(node,mels,nrmels))
				labels[node].visited=true;

			loops++;
		}
	}

#ifndef NDEBUG

	for(int c=0;c<*ilabels;c++)
	{
		if(labels[c].selfloop==false)
			assert(labels[c].visited==true);
		else
			assert(labels[c].visited==false);
	}

#endif

	return loops;
}

void add_denominator_entry(struct amatrix_weight_t *awt, int label, int qtype)
{
	int index=awt->denominators[awt->nr_denominators].ilabels;

	awt->denominators[awt->nr_denominators].labels[index]=label;
	awt->denominators[awt->nr_denominators].qtypes[index]=qtype;
	awt->denominators[awt->nr_denominators].ilabels++;
}

void add_numerator(struct amatrix_weight_t *awt,int l1,int l2,int l3,int l4)
{
	awt->numerators[awt->nr_numerators].labels[0]=l1;
	awt->numerators[awt->nr_numerators].labels[1]=l2;
	awt->numerators[awt->nr_numerators].labels[2]=l3;
	awt->numerators[awt->nr_numerators].labels[3]=l4;
	awt->nr_numerators++;
}

void add_unphysical_penalty(struct amatrix_weight_t *awt,double penalty)
{
	awt->unphysical_penalty*=penalty;
}

double reconstruct_single_weight(struct amatrix_t *amx, struct amatrix_weight_t *awt)
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

	return weight;
}

double reconstruct_weight(struct amatrix_t *amx, struct amatrix_weight_t *awt)
{
	return reconstruct_single_weight(amx,awt);
}

/*
	This function calculates a diagram's weight given the incidence matrix
*/

struct amatrix_weight_t incidence_to_weight(gsl_matrix_int *B, struct label_t *labels, int *ilabels, struct amatrix_t *amx)
{
	bool verbose=false;

	int mels[MAX_MATRIX_ELEMENTS][4];
	assert(B->size1<=MAX_MATRIX_ELEMENTS);

	double inversefactor=1.0f;

	if(verbose==true)
	{
		printf("Labels: ");

		for(int i=0;i<*ilabels;i++)
			printf("%c (%d) ", labels[i].mnemonic, labels[i].value);

		printf("\n");

		gsl_matrix_int_print(B);
	}

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

	struct amatrix_weight_t ret;

	ret.nr_denominators=0;
	ret.nr_numerators=0;
	ret.unphysical_penalty=1.0f;

	/*
		Rule 4: denominators
	*/

	double denominators=1.0f;

	for(size_t i=0;i<(B->size1-1);i++)
	{
		if(verbose==true)
			printf("{");

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

				if(verbose==true)
					printf("%c",labels[j].mnemonic);

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
		}

		if(verbose==true)
			printf("}\n");
	}

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

	if(verbose==true)
		printf("H: %d, L: %d\n",h,l);

	/*
		Finally print out all of this, in a computer-readable form, while also
		calculating the total weight.
	*/

	if(verbose==true)
		printf("1/%f\n",inversefactor);

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

		if(verbose==true)
		{
			printf("<%c", labels[mels[i][0]].mnemonic);
			printf("%c|H|", labels[mels[i][1]].mnemonic);
			printf("%c", labels[mels[i][2]].mnemonic);
			printf("%c>\n", labels[mels[i][3]].mnemonic);
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

	if(verbose==true)
		printf("Final weight: %f\n",pow(inversefactor,-1.0f)*numerators/denominators);

	/*
		Finally we put all the results into an 'amatrix_weight_t' struct and
		we return it.
	*/

	ret.weight=pow(inversefactor,-1.0f)*numerators/denominators/amatrix_multiplicity(amx);
	ret.l=l;
	ret.h=h;
	ret.inversefactor=inversefactor;

	memcpy(ret.labels,labels,sizeof(struct label_t)*MAX_LABELS);
	ret.ilabels=*ilabels;

	if(fabs(ret.weight)>1e-4)
		assert(gsl_fcmp(ret.weight, reconstruct_single_weight(amx, &ret), 1e-8)==0);

	return ret;
}

gsl_matrix_int *amatrix_calculate_incidence(struct amatrix_t *amx, struct label_t labels[MAX_LABELS], int *ilabels)
{
	int dimensions=amx->pmxs[0]->dimensions;
	assert(dimensions>0);
	assert(dimensions<=PMATRIX_MAX_DIMENSIONS);

	gsl_matrix_int *incidence=gsl_matrix_int_alloc(dimensions,2*dimensions);

	char *mnemonics[2]={"nopqrstuvwxyz","abcdefghijklm"};
	int imnemonics[2]={0, 0};

	*ilabels=0;

	/*
		Rule 1: label assignment
	*/

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(i==j)
			{
				if(amatrix_is_physical(amx)==true)
					assert(pmatrix_get_entry(amx->pmxs[0], i, j)==0);

				/*
					If we are here, we must be in the unphysical sector
				*/

				if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
				{
					assert(*ilabels<(2*dimensions));

					for(int k=0;k<dimensions;k++)
						gsl_matrix_int_set(incidence, k, *ilabels, 0);

					for(int k=0;k<dimensions;k++)
						if(k==i)
							gsl_matrix_int_set(incidence, k, *ilabels, 2);

					int qtype=pmatrix_entry_type(i,j);

					labels[*ilabels].id=*ilabels;
					labels[*ilabels].i=i;
					labels[*ilabels].j=j;
					labels[*ilabels].pmatrix=0;
					labels[*ilabels].mnemonic=get_nth_character(mnemonics[qtype],imnemonics[qtype]++);
					labels[*ilabels].qtype=qtype;
					labels[*ilabels].value=pmatrix_get_entry(amx->pmxs[0], i, j);
					labels[*ilabels].visited=false;
					labels[*ilabels].selfloop=true;
					(*ilabels)++;

					if(qtype==QTYPE_OCCUPIED)
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_occupied));
					else
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_virtual));

					assert(*ilabels<=(2*dimensions));
					assert(*ilabels<MAX_LABELS);
				}

				if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
				{
					assert(*ilabels<(2*dimensions));

					for(int k=0;k<dimensions;k++)
						gsl_matrix_int_set(incidence, k, *ilabels, 0);

					for(int k=0;k<dimensions;k++)
						if(k==i)
							gsl_matrix_int_set(incidence, k, *ilabels, 2);

					int qtype=pmatrix_entry_type(i,j);

					labels[*ilabels].id=*ilabels;
					labels[*ilabels].i=i;
					labels[*ilabels].j=j;
					labels[*ilabels].pmatrix=1;
					labels[*ilabels].mnemonic=get_nth_character(mnemonics[qtype],imnemonics[qtype]++);
					labels[*ilabels].qtype=qtype;
					labels[*ilabels].value=pmatrix_get_entry(amx->pmxs[1], i, j);
					labels[*ilabels].visited=false;
					labels[*ilabels].selfloop=true;
					(*ilabels)++;

					if(qtype==QTYPE_OCCUPIED)
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_occupied));
					else
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_virtual));

					assert(*ilabels<=(2*dimensions));
					assert(*ilabels<MAX_LABELS);
				}
			}
			else
			{
				if(pmatrix_get_entry(amx->pmxs[0], i, j)!=0)
				{
					assert(*ilabels<(2*dimensions));

					for(int k=0;k<dimensions;k++)
						gsl_matrix_int_set(incidence, k, *ilabels, 0);

					for(int k=0;k<dimensions;k++)
					{
						if(k==i)
							gsl_matrix_int_set(incidence, k, *ilabels, 1);
						else if(k==j)
							gsl_matrix_int_set(incidence, k, *ilabels, -1);
					}

					int qtype=pmatrix_entry_type(i,j);

					labels[*ilabels].id=*ilabels;
					labels[*ilabels].i=i;
					labels[*ilabels].j=j;
					labels[*ilabels].pmatrix=0;
					labels[*ilabels].mnemonic=get_nth_character(mnemonics[qtype], imnemonics[qtype]++);
					labels[*ilabels].qtype=qtype;
					labels[*ilabels].value=pmatrix_get_entry(amx->pmxs[0], i, j);
					labels[*ilabels].visited=false;
					labels[*ilabels].selfloop=false;
					(*ilabels)++;

					if(qtype==QTYPE_OCCUPIED)
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_occupied));
					else
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_virtual));

					assert(*ilabels<=(2*dimensions));
					assert(*ilabels<MAX_LABELS);
				}

				if(pmatrix_get_entry(amx->pmxs[1], i, j)!=0)
				{
					assert(*ilabels<(2*dimensions));

					for(int k=0;k<dimensions;k++)
						gsl_matrix_int_set(incidence, k, *ilabels, 0);

					for(int k=0;k<dimensions;k++)
					{
						if(k==i)
							gsl_matrix_int_set(incidence, k, *ilabels, 1);
						else if(k==j)
							gsl_matrix_int_set(incidence, k, *ilabels, -1);
					}

					int qtype=pmatrix_entry_type(i,j);

					labels[*ilabels].id=*ilabels;
					labels[*ilabels].i=i;
					labels[*ilabels].j=j;
					labels[*ilabels].pmatrix=1;
					labels[*ilabels].mnemonic=get_nth_character(mnemonics[qtype],imnemonics[qtype]++);
					labels[*ilabels].qtype=qtype;
					labels[*ilabels].value=pmatrix_get_entry(amx->pmxs[1], i, j);
					labels[*ilabels].visited=false;
					labels[*ilabels].selfloop=false;
					(*ilabels)++;

					if(qtype==QTYPE_OCCUPIED)
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_occupied));
					else
						assert((labels[*ilabels-1].value>0)&&(labels[*ilabels-1].value<=amx->nr_virtual));

					assert(*ilabels<=(2*dimensions));
					assert(*ilabels<MAX_LABELS);
				}
			}
		}
	}

	/*
		Now that we have created the incidence matrix, we can verify that it satisfies some
		properties. These checks are not performed when compiling in 'Release' mode.
	*/

#ifndef NDEBUG

	for(int i=0;i<dimensions;i++)
	{
		int sumrow=0,abssumrow=0;

		for(int j=0;j<2*dimensions;j++)
		{
			if(gsl_matrix_int_get(incidence, i, j)!=2)
				sumrow+=gsl_matrix_int_get(incidence, i, j);

			abssumrow+=abs(gsl_matrix_int_get(incidence, i, j));
		}

		assert(sumrow==0);
		assert(abssumrow==4);
	}

	for(int j=0;j<2*dimensions;j++)
	{
		int sumcolumn=0,abssumcolumn=0;

		for(int i=0;i<dimensions;i++)
		{
			if(gsl_matrix_int_get(incidence, i, j)!=2)
				sumcolumn+=gsl_matrix_int_get(incidence, i, j);

			abssumcolumn+=abs(gsl_matrix_int_get(incidence, i, j));
		}

		assert(sumcolumn==0);
		assert(abssumcolumn==2);
	}

#endif

	return incidence;
}
