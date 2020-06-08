/*
	Note that this file is almost entirely C, we compile it as C++ since
	we want to use the ALEA library from ALPSCore.
*/

#include <alps/alea.hpp>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "amatrix.h"
#include "weight.h"
#include "permutations.h"
#include "cache.h"
#include "sampling.h"
#include "auxx.h"

#ifdef __cplusplus
}
#endif

int gsl_matrix_int_compare(gsl_matrix_int *a,const gsl_matrix_int *b)
{
	if(a->size1>b->size1)
		return 1;
	else if(a->size1<b->size1)
		return -1;

	if(a->size2>b->size2)
		return 1;
	else if(a->size2<b->size2)
		return -1;

	for(size_t i=0;i<a->size1;i++)
	{
		for(size_t j=0;j<a->size2;j++)
		{
			if(gsl_matrix_int_get(a,i,j)>gsl_matrix_int_get(b,i,j))
				return 1;
			else if(gsl_matrix_int_get(a,i,j)<gsl_matrix_int_get(b,i,j))
				return -1;
		}
	}

	return 0;
}

int *get_permutation_array(int order,int c)
{
	switch(order)
	{
		case 2:
		return permutations2[c];

		case 3:
		return permutations3[c];

		case 4:
		return permutations4[c];

		case 5:
		return permutations5[c];

		case 6:
		return permutations6[c];

		default:
		break;
	}

	assert(false);
	return NULL;
}

bool matrices_in_same_bin(gsl_matrix_int *a,gsl_matrix_int *b)
{
	bool result=false;

	gsl_matrix_int *aT=gsl_matrix_int_alloc(a->size2,a->size1);

	for(size_t i=0;i<aT->size1;i++)
		for(size_t j=0;j<aT->size2;j++)
			gsl_matrix_int_set(aT,i,j,gsl_matrix_int_get(a,j,i));

	if(gsl_matrix_int_compare(a,b)==0)
		result=true;
	else if(gsl_matrix_int_compare(aT,b)==0)
		result=true;

	gsl_matrix_int_free(aT);

	return result;
}

struct sampling_ctx_t
{
	alps::alea::autocorr_acc<double> *autocorrelation;

	alps::alea::batch_acc<double> *overall_sign, *signs[MAX_ORDER];
	alps::alea::batch_acc<double> *plus[MAX_ORDER], *minus[MAX_ORDER];
	alps::alea::batch_acc<double> *orders[MAX_ORDER], **topologies[MAX_ORDER];

	long int nr_samples, nr_physical_samples, nr_samples_by_order[MAX_ORDER], nr_positive_samples[MAX_ORDER], nr_negative_samples[MAX_ORDER];

	int *bin_mappings[MAX_ORDER];
	int nr_of_bins[MAX_ORDER];

	int maxdimensions;
};

/*
	Each diagram topology can be identified by an identifier, as returned by amatrix_to_index().

	However, these identifiers are not unique: two different id's can correspond to
	the same adjacency matrix, due to multiplicity considerations. Also the identifiers
	are sequential when taking into account adjacency matrices corresponding to unconnected
	graphs; once these are removed the identifiers are sparse.

	Here we establish a lookup table -- order by order -- between the index returned
	by amatrix_to_index() and a different index which has these properties:

	1) It only labels connected, physical matrices, and it is sequential.
	Matrices representing unconnected graphs are not assigned an index.

	2) Adjacency matrices connected by simple transformations (conjugation, transposition)
	are assigned the same indices.

	This new index will then used to label bins -- again, order by order -- for sampling.
*/

void init_topology_bins(struct sampling_ctx_t *sctx,int dimensions)
{
	gsl_matrix_int **mxs=(gsl_matrix_int **)(malloc(sizeof(gsl_matrix_int *)*ifactorial(dimensions)*ifactorial(dimensions)));

	int *bin_mapping=(int *)(malloc(sizeof(int)*ifactorial(dimensions)*ifactorial(dimensions)));
	bool *valid=(bool *)(malloc(sizeof(bool)*ifactorial(dimensions)*ifactorial(dimensions)));

	for(int i=0;i<ifactorial(dimensions)*ifactorial(dimensions);i++)
		mxs[i]=gsl_matrix_int_alloc(dimensions,dimensions);

	for(int d=0;d<ifactorial(dimensions);d++)
	{
		for(int c=0;c<ifactorial(dimensions);c++)
		{
			int index=ifactorial(dimensions)*d+c;

			valid[index]=true;

			gsl_matrix_int *m=permutation_to_matrix(get_permutation_array(dimensions,c),dimensions);
			gsl_matrix_int *n=permutation_to_matrix(get_permutation_array(dimensions,d),dimensions);

			gsl_matrix_int_add(m,n);
			gsl_matrix_int_memcpy(mxs[index],m);

			bool is_unphysical=false;

			for(int l=0;l<dimensions;l++)
				if(gsl_matrix_int_get(mxs[index],l,l)!=0)
					is_unphysical=true;

			if((is_unphysical)||(gsl_matrix_int_check_connectedness(mxs[index],dimensions)==false))
			{
				bin_mapping[index]=-1;
				valid[index]=false;
			}

			gsl_matrix_int_free(m);
			gsl_matrix_int_free(n);
		}
	}

	int nextid=0;

	for(int i=0;i<ifactorial(dimensions)*ifactorial(dimensions);i++)
	{
		if(valid[i]==false)
			continue;

		bin_mapping[i]=nextid;
		valid[i]=false;

		for(int j=i+1;j<ifactorial(dimensions)*ifactorial(dimensions);j++)
		{
			if((valid[j]==true)&&(matrices_in_same_bin(mxs[i],mxs[j])==true))
			{
				bin_mapping[j]=nextid;
				valid[j]=false;
			}
		}

		nextid++;
	}

	for(int i=0;i<ifactorial(dimensions)*ifactorial(dimensions);i++)
		gsl_matrix_int_free(mxs[i]);

	sctx->nr_of_bins[dimensions]=nextid;

	/*
		Finally we initialize the bin accumulators
	*/

	sctx->topologies[dimensions]=(alps::alea::batch_acc<double> **)(malloc(sizeof(alps::alea::batch_acc<double> *)*sctx->nr_of_bins[dimensions]));

	for(int c=0;c<sctx->nr_of_bins[dimensions];c++)
		sctx->topologies[dimensions][c]=new alps::alea::batch_acc<double>;

	sctx->bin_mappings[dimensions]=bin_mapping;

	if(valid)
		free(valid);
}

char *get_bin_description(struct sampling_ctx_t *sctx,int order,int id,char buf[256])
{
	bool first=true;

	buf[0]='\0';

	for(int i=0;i<ifactorial(order)*ifactorial(order);i++)
	{
		if(sctx->bin_mappings[order][i]==id)
		{
			char tmp[256];

			if(first)
				snprintf(tmp, 256, "%d",i);
			else
				snprintf(tmp, 256, "%s,%d",buf,i);

			strncpy(buf,tmp,256);
			first=false;
		}
	}

	return buf;
}

int amatrix_to_topology_index(struct amatrix_t *amx,struct sampling_ctx_t *sctx)
{
	int dimensions=amx->pmxs[0]->dimensions;

	if(dimensions>=3)
		return sctx->bin_mappings[dimensions][amatrix_to_index(amx)];

	return 0;
}

struct sampling_ctx_t *init_sampling_ctx(int maxdimensions)
{
	struct sampling_ctx_t *ret=(struct sampling_ctx_t *)(malloc(sizeof(struct sampling_ctx_t)));

	assert(ret!=NULL);

	ret->maxdimensions=maxdimensions;

	for(int dimensions=0;dimensions<MAX_ORDER;dimensions++)
	{
		ret->nr_of_bins[dimensions]=0;

		if((dimensions>=3)&&(dimensions<=maxdimensions))
			init_topology_bins(ret, dimensions);
	}

	/*
		We initialize the accumulators and variables we will use for performing measurements.
	*/

	ret->autocorrelation=new alps::alea::autocorr_acc<double>(1);
	ret->overall_sign=new alps::alea::batch_acc<double>;

	for(int c=0;c<MAX_ORDER;c++)
	{
		ret->signs[c]=new alps::alea::batch_acc<double>;
		ret->plus[c]=new alps::alea::batch_acc<double>;
		ret->minus[c]=new alps::alea::batch_acc<double>;
		ret->orders[c]=new alps::alea::batch_acc<double>;
	}

	ret->nr_samples=ret->nr_physical_samples=0;

	for(int c=0;c<MAX_ORDER;c++)
		ret->nr_samples_by_order[c]=ret->nr_positive_samples[c]=ret->nr_negative_samples[c]=0;

	return ret;
}

void fini_sampling_ctx(struct sampling_ctx_t *sctx)
{
	if(sctx)
	{
		delete sctx->autocorrelation;
		delete sctx->overall_sign;

		for(int c=0;c<MAX_ORDER;c++)
		{
			delete sctx->signs[c];
			delete sctx->plus[c];
			delete sctx->minus[c];
			delete sctx->orders[c];
		}

		for(int dimensions=3;dimensions<=sctx->maxdimensions;dimensions++)
		{
			if(sctx->bin_mappings[dimensions])
				free(sctx->bin_mappings[dimensions]);

			if(sctx->topologies[dimensions])
			{
				for(int c=0;c<sctx->nr_of_bins[dimensions];c++)
					delete sctx->topologies[dimensions][c];

				free(sctx->topologies[dimensions]);
			}
		}

		free(sctx);
	}
}

void sampling_ctx_measure(struct sampling_ctx_t *sctx,struct amatrix_t *amx,struct configuration_t *config,long int counter)
{
	sctx->nr_samples++;

	if(amatrix_is_physical(amx))
	{
		sctx->nr_physical_samples++;

		if((counter>config->thermalization)&&((sctx->nr_physical_samples%config->decorrelation)==0))
		{
			double weight=amatrix_weight(amx);
			double sign=(weight>0.0f) ? (1.0f) : (-1.0f);

			int order=amx->pmxs[0]->dimensions;
			int topology=amatrix_to_topology_index(amx,sctx);

			(*sctx->autocorrelation) << weight;
			(*sctx->overall_sign) << sign;
			(*sctx->signs[order]) << sign;

			for(int c=0;c<=sctx->maxdimensions;c++)
			{
				(*sctx->orders[c]) << ((c!=order) ? (0.0) : (sign));
				(*sctx->plus[c]) << ((c!=order) ? (0.0) : (fpositive_part(sign)));
				(*sctx->minus[c]) << ((c!=order) ? (0.0) : (fnegative_part(sign)));

				if(c>=3)
				{
					for(int d=0;d<sctx->nr_of_bins[c];d++)
						(*sctx->topologies[c][d]) << (((c==order)&&(d==topology)) ? (sign) : (0.0));
				}
			}

			sctx->nr_samples_by_order[order]++;

			if(weight>=0.0)
			{
				sctx->nr_positive_samples[order]++;
			}
			else
			{
				sctx->nr_negative_samples[order]++;
			}
		}
	}
}

double sampling_ctx_get_physical_pct(struct sampling_ctx_t *sctx)
{
	return 100.0f*((double)(sctx->nr_physical_samples))/((double)(sctx->nr_samples));
}

void order_description(char *buf,int length,int order)
{
	if(order==1)
		snprintf(buf,length,"HF");
	else
		snprintf(buf,length,"MP%d",order);

	buf[length-1]='\0';
}

void sampling_ctx_print_report(struct sampling_ctx_t *sctx,struct amatrix_t *amx,FILE *out,bool finalize)
{
	fprintf(out,"# <Order> <Positive physical samples> <Negative physical samples> <Percentage> <Sign> <Positive fraction> <Negative fraction> <Sign (from ALEA)>\n");

	alps::alea::autocorr_result<double> result_autocorrelation;

	if(finalize==true)
		result_autocorrelation=sctx->autocorrelation->finalize();
	else
		result_autocorrelation=sctx->autocorrelation->result();

	alps::alea::batch_result<double> result_overall_sign, result_signs[MAX_ORDER];
	alps::alea::batch_result<double> result_plus[MAX_ORDER], result_minus[MAX_ORDER];
	alps::alea::batch_result<double> result_orders[MAX_ORDER], *result_topologies[MAX_ORDER];

	for(int c=0;c<MAX_ORDER;c++)
	{
		result_topologies[c]=(alps::alea::batch_result<double> *)(malloc(sizeof(alps::alea::batch_result<double>)*sctx->nr_of_bins[c]));

		for(int d=0;d<sctx->nr_of_bins[c];d++)
			new (result_topologies[c]+d) alps::alea::batch_result<double>;
	}

	if(finalize==true)
	{
		result_overall_sign=sctx->overall_sign->finalize();

		for(int c=0;c<MAX_ORDER;c++)
			result_signs[c]=sctx->signs[c]->finalize();

		for(int c=0;c<MAX_ORDER;c++)
			result_plus[c]=sctx->plus[c]->finalize();

		for(int c=0;c<MAX_ORDER;c++)
			result_minus[c]=sctx->minus[c]->finalize();

		for(int c=0;c<MAX_ORDER;c++)
			result_orders[c]=sctx->orders[c]->finalize();

		for(int c=3;c<=sctx->maxdimensions;c++)
			for(int d=0;d<sctx->nr_of_bins[c];d++)
				result_topologies[c][d]=sctx->topologies[c][d]->finalize();
	}
	else
	{
		result_overall_sign=sctx->overall_sign->result();

		for(int c=0;c<MAX_ORDER;c++)
			result_signs[c]=sctx->signs[c]->result();

		for(int c=0;c<MAX_ORDER;c++)
			result_plus[c]=sctx->plus[c]->result();

		for(int c=0;c<MAX_ORDER;c++)
			result_minus[c]=sctx->minus[c]->result();

		for(int c=0;c<MAX_ORDER;c++)
			result_orders[c]=sctx->orders[c]->result();

		for(int c=3;c<=sctx->maxdimensions;c++)
			for(int d=0;d<sctx->nr_of_bins[c];d++)
				result_topologies[c][d]=sctx->topologies[c][d]->result();
	}

	long int total_positive,total_negative;

	total_positive=total_negative=0;
	for(int order=amx->config->minorder;order<=amx->config->maxorder;order++)
	{
		total_positive+=sctx->nr_positive_samples[order];
		total_negative+=sctx->nr_negative_samples[order];
	}

	for(int order=amx->config->minorder;order<=amx->config->maxorder;order++)
	{
		double pct,sign;

		if((total_positive-total_negative)!=0)
			pct=100.0f*((double)(sctx->nr_positive_samples[order]-sctx->nr_negative_samples[order]))/(total_positive-total_negative);
		else
			pct=NAN;

		if((sctx->nr_positive_samples[order]+sctx->nr_negative_samples[order])!=0)
			sign=((double)(sctx->nr_positive_samples[order]-sctx->nr_negative_samples[order]))/((double)(sctx->nr_positive_samples[order]+sctx->nr_negative_samples[order]));
		else
			sign=NAN;

		fprintf(out, "%d %ld %ld %f %f ", order, sctx->nr_positive_samples[order], sctx->nr_negative_samples[order], pct, sign);

		fprintf(out,"%f+-%f %f+-%f ",result_plus[order].mean()(0),result_plus[order].stderror()(0),
			result_minus[order].mean()(0),result_minus[order].stderror()(0));

		/*
			Remember that result_signs[order].mean() has type alps::alea::column
			which is a shorthand for Eigen::column
		*/

		fprintf(out, "%f+-%f\n",result_signs[order].mean()(0),result_signs[order].stderror()(0));
	}

	fprintf(out,"# Overall sign: %f+-%f\n",result_overall_sign.mean()(0),result_overall_sign.stderror()(0));
	fprintf(out,"# Order-by-order ratios:\n");

	for(int order1=amx->config->minorder;order1<=amx->config->maxorder;order1++)
	{
		for(int order2=amx->config->minorder;order2<=amx->config->maxorder;order2++)
		{
			if(order1==order2)
				continue;

			/*
				We calculate the contributions at each order, their standard error,
				and then we calculation the ratio and its error using the usual
				error propagation formula.
			*/

			double phi1,phi2,sigmaphi1,sigmaphi2;

			phi1=result_orders[order1].mean()(0);
			phi2=result_orders[order2].mean()(0);

			sigmaphi1=result_orders[order1].stderror()(0);
			sigmaphi2=result_orders[order2].stderror()(0);

			char desc1[128],desc2[128];

			order_description(desc1,128,order1);
			order_description(desc2,128,order2);

			double ratio, sigmaratio;

			ratio=phi1/phi2;
			sigmaratio=ratio*sqrt(pow(sigmaphi1/phi1,2.0f)+pow(sigmaphi2/phi2,2.0f));

			fprintf(out,"%s/%s %f +- %f (%f%%)\n", desc1, desc2, ratio, sigmaratio, 100.0f*sigmaratio/ratio);
		}
	}

	fprintf(out,"# Topology-by-topology ratios:\n");

	double ratios=0.0f,sigmaratios=0.0f;

	for(int d=0;d<sctx->nr_of_bins[4];d++)
	{
		double phi1,phi2,sigmaphi1,sigmaphi2;

		phi1=result_topologies[4][d].mean()(0);
		phi2=result_orders[2].mean()(0);

		sigmaphi1=result_topologies[4][d].stderror()(0);
		sigmaphi2=result_orders[2].stderror()(0);

		double ratio, sigmaratio;

		ratio=phi1/phi2;
		sigmaratio=ratio*sqrt(pow(sigmaphi1/phi1,2.0f)+pow(sigmaphi2/phi2,2.0f));

		ratios+=ratio;
		sigmaratios+=sigmaratio*sigmaratio;

		char buf[256];

		fprintf(out,"MP4(%s)/MP2 %f +- %f (%f%%)\n",get_bin_description(sctx,4,d,buf),ratio,sigmaratio,100.0f*sigmaratio/ratio);
	}

	sigmaratios=sqrt(sigmaratios);

	fprintf(out,"MP4(all)/MP2 %f +- %f (%f%%)\n",ratios,sigmaratios,100.0f*sigmaratios/ratios);

	for(int c=0;c<MAX_ORDER;c++)
	{
		if(result_topologies[c])
		{
			for(int d=0;d<sctx->nr_of_bins[c];d++)
				(result_topologies[c]+d)->~batch_result();

			free(result_topologies[c]);
		}
	}

	fprintf(out,"# Measured autocorrelation time = %f\n",result_autocorrelation.tau()(0));
	fflush(out);
}
