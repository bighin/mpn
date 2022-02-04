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
#include <gsl/gsl_math.h>

#include "amatrix.h"
#include "weight.h"
#include "weight2.h"
#include "permutations.h"
#include "cache.h"
#include "sampling.h"
#include "auxx.h"
#include "limits.h"

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

struct sampling_ctx_t
{
	alps::alea::autocorr_acc<double> *autocorrelation;

	alps::alea::batch_acc<double> *overall_sign, *signs[MAX_ORDER];
	alps::alea::batch_acc<double> *plus[MAX_ORDER], *minus[MAX_ORDER];
	alps::alea::batch_acc<double> *orders[MAX_ORDER];

	long int nr_samples, nr_physical_samples, nr_samples_by_order[MAX_ORDER], nr_positive_samples[MAX_ORDER], nr_negative_samples[MAX_ORDER];

	int maxdimensions;
};

struct sampling_ctx_t *init_sampling_ctx(int maxdimensions)
{
	struct sampling_ctx_t *ret=(struct sampling_ctx_t *)(malloc(sizeof(struct sampling_ctx_t)));

	assert(ret!=NULL);

	ret->maxdimensions=maxdimensions;

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

			(*sctx->autocorrelation) << weight;
			(*sctx->overall_sign) << sign;
			(*sctx->signs[order]) << sign;

			for(int c=0;c<=sctx->maxdimensions;c++)
			{
				(*sctx->orders[c]) << ((c!=order) ? (0.0) : (sign));
				(*sctx->plus[c]) << ((c!=order) ? (0.0) : (fpositive_part(sign)));
				(*sctx->minus[c]) << ((c!=order) ? (0.0) : (fnegative_part(sign)));
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
	alps::alea::batch_result<double> result_orders[MAX_ORDER];

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
	}

	long int total_positive, total_negative;

	total_positive=total_negative=0;
	for(int order=amx->config->minorder;order<=amx->config->maxorder;order++)
	{
		total_positive+=sctx->nr_positive_samples[order];
		total_negative+=sctx->nr_negative_samples[order];
	}

	for(int order=amx->config->minorder;order<=amx->config->maxorder;order++)
	{
		double pct, sign;

		if((total_positive-total_negative)!=0)
			pct=100.0f*((double)(sctx->nr_positive_samples[order]-sctx->nr_negative_samples[order]))/
				(total_positive-total_negative);
		else
			pct=NAN;

		if((sctx->nr_positive_samples[order]+sctx->nr_negative_samples[order])!=0)
			sign=((double)(sctx->nr_positive_samples[order]-sctx->nr_negative_samples[order]))/
				((double)(sctx->nr_positive_samples[order]+sctx->nr_negative_samples[order]));
		else
			sign=NAN;

		fprintf(out, "%d %ld %ld %f %f ", order, sctx->nr_positive_samples[order],
			sctx->nr_negative_samples[order], pct, sign);

		fprintf(out, "%f+-%f %f+-%f ", result_plus[order].mean()(0), result_plus[order].stderror()(0),
			result_minus[order].mean()(0), result_minus[order].stderror()(0));

		/*
			Remember that result_signs[order].mean() has type alps::alea::column
			which is a shorthand for Eigen::column
		*/

		fprintf(out, "%f+-%f\n", result_signs[order].mean()(0), result_signs[order].stderror()(0));
	}

	fprintf(out, "# Overall sign: %f +- %f\n", result_overall_sign.mean()(0), result_overall_sign.stderror()(0));
	fprintf(out, "# Order-by-order ratios:\n");

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

			double phi1, phi2, sigmaphi1, sigmaphi2;

			phi1=result_orders[order1].mean()(0);
			phi2=result_orders[order2].mean()(0);

			sigmaphi1=result_orders[order1].stderror()(0);
			sigmaphi2=result_orders[order2].stderror()(0);

			char desc1[128], desc2[128];

			order_description(desc1, 128, order1);
			order_description(desc2, 128, order2);

			double ratio, sigmaratio;

			ratio=phi1/phi2;
			sigmaratio=ratio*sqrt(pow(sigmaphi1/phi1, 2.0f)+pow(sigmaphi2/phi2, 2.0f));

			fprintf(out, "%s/%s %f +- %f (%f%%)\n", desc1, desc2, ratio, sigmaratio,
				100.0f*sigmaratio/ratio);
		}
	}

	fprintf(out,"# Measured autocorrelation time = %f\n",result_autocorrelation.tau()(0));
	fflush(out);
}
