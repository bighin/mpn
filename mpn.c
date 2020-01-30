#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "mpn.h"
#include "amatrix.h"
#include "pmatrix.h"
#include "auxx.h"
#include "reader.h"
#include "mc.h"

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

/*
	This function calculates a diagram's weight given the incidence matrix
*/

double incidence_to_weight(gsl_matrix_int *B, struct label_t *labels, int *ilabels, struct energies_ctx_t *ectx,double unphysical_penalty,bool verbose)
{
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
		Rule 4: denominators
	*/

	double denominators=1.0f;

	for(size_t i=0;i<(B->size1-1);i++)
	{
		if(verbose==true)
			printf("{");

		double denominator=0.0f;
		int energies_in_denominator=0;

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

				switch(labels[j].qtype)
				{
					case QTYPE_OCCUPIED:
					denominator+=get_occupied_energy(ectx,labels[j].value-1);
					energies_in_denominator++;
					break;

					case QTYPE_VIRTUAL:
					denominator-=get_virtual_energy(ectx,labels[j].value-1);
					energies_in_denominator++;
					break;

					default:
					assert(false);
				}
			}
		}

		if(energies_in_denominator>0)
			denominators*=denominator;

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
				We assing an unphysical penalty for 'selfloops', i.e. edges on the
				graph connecting a vertex with itself, appearing only in the unphysical sector.
			*/

			numerators*=unphysical_penalty;

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
			We have to convert the quantum numbers into the format used by get_htensor()
		*/

		int i1,i2,i3,i4;

		i1=labels[mels[i][0]].value-1+((labels[mels[i][0]].qtype==QTYPE_VIRTUAL)?(ectx->nocc):(0));
		i2=labels[mels[i][1]].value-1+((labels[mels[i][1]].qtype==QTYPE_VIRTUAL)?(ectx->nocc):(0));
		i3=labels[mels[i][2]].value-1+((labels[mels[i][2]].qtype==QTYPE_VIRTUAL)?(ectx->nocc):(0));
		i4=labels[mels[i][3]].value-1+((labels[mels[i][3]].qtype==QTYPE_VIRTUAL)?(ectx->nocc):(0));

		numerators*=get_htensor(ectx,i1,i2,i3,i4);
	}

	if(verbose==true)
		printf("Final weight: %f\n",pow(inversefactor,-1.0f)*numerators/denominators);

	return pow(inversefactor,-1.0f)*numerators/denominators;
}

char get_nth_character(char *s,size_t n)
{
	if(n<strlen(s))
		return s[n];

	return '_';
}

gsl_matrix_int *amatrix_calculate_incidence(struct amatrix_t *amx, struct label_t labels[MAX_LABELS], int *ilabels)
{
	int dimensions=amx->pmxs[0]->dimensions;
	assert(dimensions>0);
	assert(dimensions<=IMATRIX_MAX_DIMENSIONS);

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
		properties. This checks are not performed when compiling in 'Release' mode.
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

void run_debug_tests(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	amatrix_print(amx);
	printf("\n");

	assert(amatrix_check_consistency(amx)==true);

	for(int c=0;c<5;c++)
	{
		update_extend(amx, true);

		amatrix_print(amx);
		printf("Weight: %f\n\n",amatrix_weight(amx));

		assert(amatrix_check_consistency(amx)==true);
	}

	for(int c=0;c<1000000;c++)
	{
		update_shuffle(amx, true);
		printf("%f\n",amatrix_weight(amx));
		assert(amatrix_check_consistency(amx)==true);
	}

	for(int c=0;c<1000000;c++)
	{
		update_shuffle(amx, true);

		if(amatrix_is_physical(amx)==true)
		{
			printf("Found a physical matrix!\n");

			amatrix_print(amx);
			printf("\n");

			struct label_t labels[MAX_LABELS];
			int ilabels;

			amatrix_to_python(amx);
			amatrix_to_wolfram(amx);
			gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
			incidence_to_weight(incidence, labels, &ilabels, amx->ectx, amx->unphysical_penalty, true);
			gsl_matrix_int_free(incidence);

			exit(0);
		}

		assert(amatrix_check_consistency(amx)==true);
	}

	for(int c=0;c<5;c++)
	{
		update_squeeze(amx, true);
		amatrix_print(amx);
		printf("\n");

		assert(amatrix_check_consistency(amx)==true);
	}

	fini_amatrix(amx);
}

void run_debug_tests2(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	amx->pmxs[0]->dimensions=6;

	pmatrix_set_entry(amx->pmxs[0], 0, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 3, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 0, 3));
	pmatrix_set_entry(amx->pmxs[0], 0, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 0, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 1, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 2, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 1, 2));
	pmatrix_set_entry(amx->pmxs[0], 1, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 1, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 2, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 2, 5, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 2, 5));

	pmatrix_set_entry(amx->pmxs[0], 3, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 3, 4, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 3, 4));
	pmatrix_set_entry(amx->pmxs[0], 3, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 4, 0, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 4, 0));
	pmatrix_set_entry(amx->pmxs[0], 4, 1, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 4, 5, 0);

	pmatrix_set_entry(amx->pmxs[0], 5, 0, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 1, pmatrix_get_new_value(amx->pmxs[0], amx->rng_ctx, 5, 1));
	pmatrix_set_entry(amx->pmxs[0], 5, 2, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 3, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 4, 0);
	pmatrix_set_entry(amx->pmxs[0], 5, 5, 0);

	amx->pmxs[1]->dimensions=6;

	pmatrix_set_entry(amx->pmxs[1], 0, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 1, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 0, 1));
	pmatrix_set_entry(amx->pmxs[1], 0, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 0, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 1, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 3, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 1, 3));
	pmatrix_set_entry(amx->pmxs[1], 1, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 1, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 2, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 2, 4, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 2, 4));
	pmatrix_set_entry(amx->pmxs[1], 2, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 3, 0, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 3, 0));
	pmatrix_set_entry(amx->pmxs[1], 3, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 3, 5, 0);

	pmatrix_set_entry(amx->pmxs[1], 4, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 2, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 4, 5, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 4, 5));

	pmatrix_set_entry(amx->pmxs[1], 5, 0, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 1, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 2, pmatrix_get_new_value(amx->pmxs[1], amx->rng_ctx, 5, 2));
	pmatrix_set_entry(amx->pmxs[1], 5, 3, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 4, 0);
	pmatrix_set_entry(amx->pmxs[1], 5, 5, 0);

	assert(amatrix_check_consistency(amx)==true);

	amatrix_print(amx);
	printf("\n");

	struct label_t labels[MAX_LABELS];
	int ilabels;

	amatrix_to_python(amx);
	amatrix_to_wolfram(amx);
	gsl_matrix_int *incidence=amatrix_calculate_incidence(amx, labels, &ilabels);
	incidence_to_weight(incidence, labels, &ilabels, amx->ectx, amx->unphysical_penalty, true);
	gsl_matrix_int_free(incidence);
}

/*
	Here we debug the connectedness test
*/

gsl_matrix_int *permutation_to_matrix(int *permutation,int dimensions)
{
	gsl_matrix_int *ret=gsl_matrix_int_alloc(dimensions,dimensions);

	for(int i=0;i<dimensions;i++)
		for(int j=0;j<dimensions;j++)
			gsl_matrix_int_set(ret,i,j,((permutation[i]-1)==j)?(1):(0));

	return ret;
}

#include "permutations.h"

/*
	Here we check the number of connected diagrams at order 3,4,5,6 and compare
	them with values obtained from Wolfram Mathematica, please see the notebook
	Connectedness.nb
*/

void run_debug_tests3(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	int connected,not_connected;

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=3;
	connected=not_connected=0;

	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations3[i],3);
			b=permutation_to_matrix(permutations3[j],3);

			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==26)&&(not_connected==10));

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=4;
	connected=not_connected=0;

	for(int i=0;i<24;i++)
	{
		for(int j=0;j<24;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations4[i],4);
			b=permutation_to_matrix(permutations4[j],4);

			for(int k=0;k<4;k++)
			{
				for(int l=0;l<4;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==426)&&(not_connected==150));

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=5;
	connected=not_connected=0;

	for(int i=0;i<120;i++)
	{
		for(int j=0;j<120;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations5[i],5);
			b=permutation_to_matrix(permutations5[j],5);

			for(int k=0;k<5;k++)
			{
				for(int l=0;l<5;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==11064)&&(not_connected==3336));

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=6;
	connected=not_connected=0;

	for(int i=0;i<720;i++)
	{
		for(int j=0;j<720;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations6[i],6);
			b=permutation_to_matrix(permutations6[j],6);

			for(int k=0;k<6;k++)
			{
				for(int l=0;l<6;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_check_connectedness(amx)==true)
				connected++;
			else
				not_connected++;

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
	assert((connected==413640)&&(not_connected==104760));
}

void run_debug_tests3b(void)
{
	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	int connected,not_connected;

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=3;
	connected=not_connected=0;

	for(int i=0;i<6;i++)
	{
		for(int j=0;j<6;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations3[i],3);
			b=permutation_to_matrix(permutations3[j],3);

			for(int k=0;k<3;k++)
			{
				for(int l=0;l<3;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=4;
	connected=not_connected=0;

	for(int i=0;i<24;i++)
	{
		for(int j=0;j<24;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations4[i],4);
			b=permutation_to_matrix(permutations4[j],4);

			for(int k=0;k<4;k++)
			{
				for(int l=0;l<4;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=5;
	connected=not_connected=0;

	for(int i=0;i<120;i++)
	{
		for(int j=0;j<120;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations5[i],5);
			b=permutation_to_matrix(permutations5[j],5);

			for(int k=0;k<5;k++)
			{
				for(int l=0;l<5;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);

	amx->pmxs[0]->dimensions=amx->pmxs[1]->dimensions=6;
	connected=not_connected=0;

	for(int i=0;i<720;i++)
	{
		for(int j=0;j<720;j++)
		{
			gsl_matrix_int *a,*b;

			a=permutation_to_matrix(permutations6[i],6);
			b=permutation_to_matrix(permutations6[j],6);

			for(int k=0;k<6;k++)
			{
				for(int l=0;l<6;l++)
				{
					amx->pmxs[0]->values[k][l]=gsl_matrix_int_get(a, k, l);
					amx->pmxs[1]->values[k][l]=gsl_matrix_int_get(b, k, l);
				}
			}

			if(amatrix_is_physical(amx)==true)
			{
				if(amatrix_check_connectedness(amx)==true)
					connected++;
				else
					not_connected++;
			}

			gsl_matrix_int_free(a);
			gsl_matrix_int_free(b);
		}
	}

	printf("%d %d\n",connected,not_connected);
}

void run_debug_tests4(void)
{
	/*
		'Manual' MP2 calculation
	*/

	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	double result=0.0;

	pmatrix_set_entry(amx->pmxs[0],0,0,0);
	pmatrix_set_entry(amx->pmxs[0],0,1,0);
	pmatrix_set_entry(amx->pmxs[0],1,0,0);
	pmatrix_set_entry(amx->pmxs[0],1,1,0);
	pmatrix_set_entry(amx->pmxs[1],0,0,0);
	pmatrix_set_entry(amx->pmxs[1],0,1,0);
	pmatrix_set_entry(amx->pmxs[1],1,0,0);
	pmatrix_set_entry(amx->pmxs[1],1,1,0);

	for(int i=1;i<=amx->nr_virtual;i++)
	{
		for(int j=1;j<=amx->nr_virtual;j++)
		{
			for(int a=1;a<=amx->nr_occupied;a++)
			{
				for(int b=1;b<=amx->nr_occupied;b++)
				{
					pmatrix_set_entry(amx->pmxs[0],0,1,a);
					pmatrix_set_entry(amx->pmxs[0],1,0,i);
					pmatrix_set_entry(amx->pmxs[1],0,1,b);
					pmatrix_set_entry(amx->pmxs[1],1,0,j);

					result+=amatrix_weight(amx);

					bool verbose=false;

					if(verbose==true)
					{

						printf("%d %d %d %d || %f\n", pmatrix_get_entry(amx->pmxs[0], 0, 1),
							                      pmatrix_get_entry(amx->pmxs[0], 1, 0),
						                              pmatrix_get_entry(amx->pmxs[1], 0, 1),
						                              pmatrix_get_entry(amx->pmxs[1], 1, 0),
						                              amatrix_weight(amx));
					}
				}
			}
		}
	}

	printf("%f\n",result);
}


void run_debug_tests5(void)
{
	/*
		'Manual' MP2 calculation, with stochastic sampling of quantum numbers
	*/

	struct amatrix_t *amx=init_amatrix("/Users/zakk/Desktop/MPn/psi4/H2O.dat");

	pmatrix_set_entry(amx->pmxs[0],0,0,0);
	pmatrix_set_entry(amx->pmxs[0],0,1,0);
	pmatrix_set_entry(amx->pmxs[0],1,0,0);
	pmatrix_set_entry(amx->pmxs[0],1,1,0);
	pmatrix_set_entry(amx->pmxs[1],0,0,0);
	pmatrix_set_entry(amx->pmxs[1],0,1,0);
	pmatrix_set_entry(amx->pmxs[1],1,0,0);
	pmatrix_set_entry(amx->pmxs[1],1,1,0);

	double result=0.0;
	int iterations=50000000;

	for(int c=0;c<iterations;c++)
	{
		int i=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_virtual);
		int j=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_virtual);
		int a=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_occupied);
		int b=1+gsl_rng_uniform_int(amx->rng_ctx,amx->nr_occupied);

		pmatrix_set_entry(amx->pmxs[0],0,1,a);
		pmatrix_set_entry(amx->pmxs[0],1,0,i);
		pmatrix_set_entry(amx->pmxs[1],0,1,b);
		pmatrix_set_entry(amx->pmxs[1],1,0,j);

		result+=amatrix_weight(amx);
	}

	printf("%f\n",result/iterations*pow(amx->nr_virtual*amx->nr_occupied,2.0f));
}

int main(void)
{
	//run_debug_tests();
	run_debug_tests4();
	run_debug_tests5();
	//do_diagmc("/Users/zakk/Desktop/MPn/psi4/H2O.dat","/Users/zakk/Desktop/MPn/test.dat",10000000,100.0f,false);

	return 0;
}
