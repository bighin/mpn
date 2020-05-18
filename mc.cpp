/*
	Note that this file is almost entirely C, we compile it as C++ since
	we want to use the ALEA library from ALPSCore.
*/

#include <alps/alea.hpp>

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <signal.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "mc.h"
#include "amatrix.h"
#include "auxx.h"
#include "mpn.h"
#include "config.h"
#include "permutations.h"

#include "libprogressbar/progressbar.h"

/*
	The updates
*/

int coordinate_to_label_index(struct label_t *labels,int ilabels,int i,int j,int pmatrix)
{
	for(int c=0;c<ilabels;c++)
		if((labels[c].i==i)&&(labels[c].j==j)&&(labels[c].pmatrix==pmatrix))
			return c;

	assert(false);
	return 0;
}

double extend_pdf(struct amatrix_t *amx,struct amatrix_weight_t *awt)
{
	double weight=reconstruct_weight(amx, awt);

	return fabs(weight);
}

void pmatrix_set_entry_with_distribution(struct pmatrix_t *pmx, int i, int j, int value, gsl_rng *rngctx)
{
	if(pmatrix_entry_type(i,j)==QTYPE_VIRTUAL)
	{
		int actualvalue=value+pmx->nr_virtual*gsl_rng_uniform_int(rngctx, pmx->nr_occupied);

		pmatrix_set_entry(pmx,i,j,actualvalue);

		assert(pmatrix_get_entry(pmx,i,j)==value);
		assert((actualvalue>=1)&&(actualvalue<=(pmx->nr_virtual*pmx->nr_occupied)));
		return;
	}
	else if(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)
	{
		int actualvalue=value+pmx->nr_occupied*gsl_rng_uniform_int(rngctx, pmx->nr_virtual);

		pmatrix_set_entry(pmx,i,j,actualvalue);

		assert(pmatrix_get_entry(pmx,i,j)==value);
		assert((actualvalue>=1)&&(actualvalue<=(pmx->nr_virtual*pmx->nr_occupied)));
		return;
	}

	assert(false);
}

int update_extend(struct amatrix_t *amx, bool always_accept)
{
	double weightratio=1.0f/fabs(amatrix_weight(amx).weight);
	double probability=1.0f;

	if(amx->pmxs[0]->dimensions>=amx->config->maxorder)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	int i1,j1,i2,j2;

	probability*=pmatrix_extend(amx->pmxs[0], amx->rng_ctx, &i1, &j1);
	probability*=pmatrix_extend(amx->pmxs[1], amx->rng_ctx, &i2, &j2);

	if(amatrix_check_connectedness(amx)==false)
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	amx->cached_weight_is_valid=false;

	struct amatrix_weight_t w=amatrix_weight(amx);

	/*
		Now we have to select the two missing values from a suitable distribution, that
		depends, therefore, on two variables.

		We map the (discrete) joint probability distribution function to a unidimensional
		one, using p(index)=p(c,d) with:

		index=c+d*nr_states1

		the obvious inverses being:

		c=index%nr_states1
		d=index/nr_states1
	*/

	int nr_states1=(pmatrix_entry_type(i1,j1)==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied);
	int nr_states2=(pmatrix_entry_type(i2,j2)==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied);

	int label1=coordinate_to_label_index(w.labels,w.ilabels,i1,j1,0);
	int label2=coordinate_to_label_index(w.labels,w.ilabels,i2,j2,1);

	double *dists=(double *)malloc(sizeof(double)*nr_states1*nr_states2);
	double *cdists=(double *)malloc(sizeof(double)*nr_states1*nr_states2);

	for(int c=0;c<nr_states1;c++)
	{
		for(int d=0;d<nr_states2;d++)
		{
			int index=c+d*nr_states1;

			assert(index==((index%nr_states1)+(index/nr_states1)*nr_states1));

			w.labels[label1].value=1+c;
			w.labels[label2].value=1+d;

			dists[index]=extend_pdf(amx, &w);
		}
	}

	normalize_distribution(dists,nr_states1*nr_states2);
	to_cumulative_distribution(dists,cdists,nr_states1*nr_states2);

	/*
		Now that we have the probability distribution, we sample one random value from it.
	*/

	double selector=gsl_rng_uniform(amx->rng_ctx);

	int c,d,index;

	index=cdist_search(cdists, 0, nr_states1*nr_states2, selector);
	c=index%nr_states1;
	d=index/nr_states1;

	assert(index==(c+d*nr_states1));
	assert(cdist_search(cdists, 0, nr_states1*nr_states2, selector)==cdist_linear_search(cdists, 0, nr_states1*nr_states2, selector));

	pmatrix_set_entry_with_distribution(amx->pmxs[0], i1, j1, 1+c, amx->rng_ctx);
	pmatrix_set_entry_with_distribution(amx->pmxs[1], i2, j2, 1+d, amx->rng_ctx);
	amx->cached_weight_is_valid=false;

	w.labels[label1].value=1+c;
	w.labels[label2].value=1+d;

	probability*=1.0f/dists[index];

	free(cdists);
	free(dists);

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	double acceptance_ratio;
	double currentweight=reconstruct_weight(amx, &w);

	amx->cached_weight=w;
	amx->cached_weight.weight=currentweight;
	amx->cached_weight_is_valid=true;

	weightratio*=fabs(currentweight);
	acceptance_ratio=weightratio*probability;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio) ? (true) : (false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

void find_squeeze_target(struct pmatrix_t *pmx,int *iprime,int *jprime)
{
#ifndef NDEBUG
	*iprime=*jprime=-1;
#endif

	int dimensions=pmx->dimensions;

	for(int i=0;i<dimensions;i++)
	{
		if(pmatrix_get_entry(pmx,i,dimensions-1)!=0)
		{
			*iprime=i;
			break;
		}
	}

	for(int j=0;j<dimensions;j++)
	{
		if(pmatrix_get_entry(pmx,dimensions-1,j)!=0)
		{
			*jprime=j;
			break;
		}
	}

	assert((*iprime!=-1)&&(*jprime!=-1));

	/*
		If the Squeeze update will be targeting the corner entry in
		the bottom-right corner, we can return straight away.
	*/

	if((*iprime==(dimensions-1))&&(*jprime==(dimensions-1)))
		return;

	/*
		Otherwise the update will be targeting an entry either in the rightmost
		column or in the bottom row. Let's find the coordinates and return.
	*/

	assert(pmatrix_entry_type(*iprime,dimensions-1)!=pmatrix_entry_type(dimensions-1,*jprime));

	if(pmatrix_entry_type(*iprime,*jprime)!=pmatrix_entry_type(*iprime,dimensions-1))
	{
		*jprime=dimensions-1;
		return;
	}
	else if(pmatrix_entry_type(*iprime,*jprime)!=pmatrix_entry_type(dimensions-1,*jprime))
	{
		*iprime=dimensions-1;
		return;
	}

	assert(false);
}

int update_squeeze(struct amatrix_t *amx, bool always_accept)
{
	struct amatrix_weight_t w=amatrix_weight(amx);

	double weightratio=1.0f/fabs(w.weight);
	double probability=1.0f;

	if(amx->pmxs[0]->dimensions<=amx->config->minorder)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	int i1,j1,i2,j2;

	find_squeeze_target(amx->pmxs[0],&i1,&j1);
	find_squeeze_target(amx->pmxs[1],&i2,&j2);

	int label1=coordinate_to_label_index(w.labels,w.ilabels,i1,j1,0);
	int label2=coordinate_to_label_index(w.labels,w.ilabels,i2,j2,1);

	int value1=pmatrix_get_entry(amx->pmxs[0],i1,j1);
	int value2=pmatrix_get_entry(amx->pmxs[1],i2,j2);

	int nr_states1=(pmatrix_entry_type(i1,j1)==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied);
	int nr_states2=(pmatrix_entry_type(i2,j2)==QTYPE_VIRTUAL)?(amx->nr_virtual):(amx->nr_occupied);

	double *dists=(double *)malloc(sizeof(double)*nr_states1*nr_states2);

	for(int c=0;c<nr_states1;c++)
	{
		for(int d=0;d<nr_states2;d++)
		{
			int index=c+d*nr_states1;

			w.labels[label1].value=1+c;
			w.labels[label2].value=1+d;

			dists[index]=extend_pdf(amx,&w);
		}
	}

	normalize_distribution(dists,nr_states1*nr_states2);

	int index=(value1-1)+(value2-1)*nr_states1;
	probability*=dists[index];

	free(dists);

	probability*=pmatrix_squeeze(amx->pmxs[0], amx->rng_ctx);
	probability*=pmatrix_squeeze(amx->pmxs[1], amx->rng_ctx);
	amx->cached_weight_is_valid=false;

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx).weight);
	acceptance_ratio=weightratio*probability;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio)?(true):(false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_shuffle(struct amatrix_t *amx, bool always_accept)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	/*
		If we are in dimension 1, there's nothing to shuffle.
	*/

	if(dimensions==1)
		return UPDATE_UNPHYSICAL;

	double weightratio=1.0f/fabs(amatrix_weight(amx).weight);

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	/*
		We select which one of the permutation matrices we want to play with
	*/

	struct pmatrix_t *target=amx->pmxs[gsl_rng_uniform_int(amx->rng_ctx, 2)];

	/*
		We select the rows/columns we want to swap.

		Note that we are choosing the first one among dims possible choices,
		and the second one among (dims-1) possible choices
	*/

	int i=gsl_rng_uniform_int(amx->rng_ctx, dimensions);
	int j=gsl_rng_uniform_int(amx->rng_ctx, dimensions-1);

	if(j>=i)
		j++;

	assert(i!=j);

	/*
		Are we swapping rows or columns?
	*/

	int update[2],reverse[2];

	switch(gsl_rng_uniform_int(amx->rng_ctx, 2))
	{
		case 0:
		pmatrix_swap_rows(target, i, j, update, reverse, amx->rng_ctx);
		break;

		case 1:
		pmatrix_swap_cols(target, i, j, update, reverse, amx->rng_ctx);
		break;
	}

	amx->cached_weight_is_valid=false;

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx).weight);
	acceptance_ratio=weightratio;
	acceptance_ratio*=pow(amx->nr_occupied, update[QTYPE_OCCUPIED]);
	acceptance_ratio*=pow(amx->nr_virtual, update[QTYPE_VIRTUAL]);
	acceptance_ratio/=pow(amx->nr_occupied, reverse[QTYPE_OCCUPIED]);
	acceptance_ratio/=pow(amx->nr_virtual, reverse[QTYPE_VIRTUAL]);

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio)?(true):(false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_modify(struct amatrix_t *amx, bool always_accept)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	double weightratio=1.0f/fabs(amatrix_weight(amx).weight);

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	/*
		We select which one of the permutation matrices we want to play with
	*/

	struct pmatrix_t *target=amx->pmxs[gsl_rng_uniform_int(amx->rng_ctx, 2)];

	/*
		We select the row the element we want to modify lies in. Since there's one
		and only one element per row, we do not need to select a column.
	*/

	for(int repeats=0;repeats<1;repeats++)
	{
		int i=gsl_rng_uniform_int(amx->rng_ctx, dimensions);

		for(int j=0;j<dimensions;j++)
			if(pmatrix_get_entry(target, i, j)!=0)
				pmatrix_set_entry(target, i, j, pmatrix_get_new_value(target, amx->rng_ctx, i, j));
	}

	amx->cached_weight_is_valid=false;

	/*
		The update is balanced with itself, the acceptance ratio is simply given
		by the (modulus of the) weights ratio.
	*/

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx).weight);
	acceptance_ratio=weightratio;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio)?(true):(false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_swap(struct amatrix_t *amx, bool always_accept)
{
	int dimensions=amx->pmxs[0]->dimensions;

	assert(amx->pmxs[0]->dimensions==amx->pmxs[1]->dimensions);
	assert(amx->pmxs[0]->dimensions>=1);

	double weightratio=1.0f/fabs(amatrix_weight(amx).weight);

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	/*
		Do we play with virtual or with occupied states?
	*/

	int target_type=(gsl_rng_uniform_int(amx->rng_ctx, 2)==0) ? (QTYPE_OCCUPIED) : (QTYPE_VIRTUAL);

	struct
	{
		int i,j,pmatrix;
	}
	candidates[32];

	int icandidates=0;

	for(int i=0;i<dimensions;i++)
	{
		for(int j=0;j<dimensions;j++)
		{
			if(pmatrix_entry_type(i,j)==target_type)
			{
				if(pmatrix_get_entry(amx->pmxs[0],i,j)!=0)
				{
					candidates[icandidates].i=i;
					candidates[icandidates].j=j;
					candidates[icandidates].pmatrix=0;
					icandidates++;
				}

				if(pmatrix_get_entry(amx->pmxs[1],i,j)!=0)
				{
					candidates[icandidates].i=i;
					candidates[icandidates].j=j;
					candidates[icandidates].pmatrix=1;
					icandidates++;
				}
			}
		}
	}

	if(icandidates<2)
		return UPDATE_UNPHYSICAL;

	int values[32];

	for(int c=0;c<icandidates;c++)
		values[c]=pmatrix_get_entry(amx->pmxs[candidates[c].pmatrix],candidates[c].i,candidates[c].j);

	/*
		TODO: Probably it would be better to use perfect probability distributions.
	*/

	fisher_yates(amx->rng_ctx,values,icandidates);

	for(int c=0;c<icandidates;c++)
		pmatrix_set_entry(amx->pmxs[candidates[c].pmatrix],candidates[c].i,candidates[c].j,values[c]);

	amx->cached_weight_is_valid=false;

	/*
		The update is balanced with itself, the acceptance ratio is simply given
		by the (modulus of the) weights ratio.
	*/

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx).weight);
	acceptance_ratio=weightratio;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio)?(true):(false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

/*
	Auxiliary functions
*/

void show_update_statistics(FILE *out,long int proposed,long int accepted,long int rejected)
{
	double accepted_pct,rejected_pct;

	if(proposed>0)
	{
		accepted_pct=100.0f*((double)(accepted))/((double)(proposed));
		rejected_pct=100.0f*((double)(rejected))/((double)(proposed));
	}
	else
	{
		accepted_pct=rejected_pct=0.0f;
	}

	fprintf(out,"proposed %ld, accepted %ld (%f%%), rejected %ld (%f%%).\n",proposed,accepted,accepted_pct,rejected,rejected_pct);
}

void order_description(char *buf,int length,int order)
{
	if(order==1)
		snprintf(buf,length,"HF");
	else
		snprintf(buf,length,"MP%d",order);

	buf[length-1]='\0';
}

static volatile int keep_running=1;

void interrupt_handler(int dummy __attribute__((unused)))
{
	keep_running=0;
}

#ifdef __cplusplus
}
#endif

/*
	The actual DiagMC routine.
*/

int do_diagmc(struct configuration_t *config)
{

#define DIAGRAM_NR_UPDATES        (5)

	int (*updates[DIAGRAM_NR_UPDATES])(struct amatrix_t *amx, bool always_accept);
	const char *update_names[DIAGRAM_NR_UPDATES];

	long int proposed[DIAGRAM_NR_UPDATES], accepted[DIAGRAM_NR_UPDATES], rejected[DIAGRAM_NR_UPDATES];
	int update_probability[DIAGRAM_NR_UPDATES],cumulative_probability[DIAGRAM_NR_UPDATES];

	/*
		We set up the updates we will be using
	*/

	updates[0]=update_extend;
	updates[1]=update_squeeze;
	updates[2]=update_shuffle;
	updates[3]=update_modify;
	updates[4]=update_swap;

	update_names[0]="Extend";
	update_names[1]="Squeeze";
	update_names[2]="Shuffle";
	update_names[3]="Modify";
	update_names[4]="Swap";

	/*
		Update probabilities: note that they must be the same for symmetric updates,
		see the assert() below.
	*/

	update_probability[0]=1;
	update_probability[1]=1;
	update_probability[2]=1;
	update_probability[3]=1;
	update_probability[4]=1;

	assert(update_probability[0]==update_probability[1]);

	/*
		Here we calculate the cumulative probabilities from the update probabilities.
	*/

	cumulative_probability[0]=update_probability[0];

	for(int c=1;c<DIAGRAM_NR_UPDATES;c++)
		cumulative_probability[c]=update_probability[c]+cumulative_probability[c-1];

	/*
		We reset the update statistics
	*/

	for(int d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;

	/*
		We initialize the accumulators and variables we will use for performing measurements.
	*/

	assert(config->maxorder>config->minorder);
	assert(config->maxorder<MAX_ORDER);

	alps::alea::autocorr_acc<double> autocorrelation(1);
	alps::alea::batch_acc<double> overall_sign, signs[MAX_ORDER], plus[MAX_ORDER], minus[MAX_ORDER], orders[MAX_ORDER];

#define	MAX_EXCITATIONS (16)

	alps::alea::batch_acc<double> excitations[MAX_EXCITATIONS];

	long int nr_samples, nr_physical_samples, nr_samples_by_order[MAX_ORDER], nr_positive_samples[MAX_ORDER], nr_negative_samples[MAX_ORDER];

	nr_samples=nr_physical_samples=0;

	for(int c=0;c<MAX_ORDER;c++)
		nr_samples_by_order[c]=nr_positive_samples[c]=nr_negative_samples[c]=0;

	/*
		We print some informative message, and then we open the log file
	*/

	printf("Performing %ld iterations\n",config->iterations);

	FILE *out;
	char output[1024];

	snprintf(output,1024,"%s.dat",config->prefix);
	output[1023]='\0';

	if(!(out=fopen(output,"w+")))
	{
		fprintf(stderr,"Error: couldn't open %s for writing\n",output);
		return 0;
	}

	printf("Writing results to '%s'\n",output);

	/*
		The diagram parameters are loaded from the configuration, as a new 'amatrix' is created
	*/

	struct amatrix_t *amx=init_amatrix(config);

	if(!amx)
	{
		fprintf(stderr,"Error: couldn't load the ERIs file (%s).\n",config->erisfile);
		return 0;
	}

#warning Do we have a better way of doing this?

	while(amx->pmxs[0]->dimensions<config->minorder)
		update_extend(amx,true);

	/*
		We setup an interrupt handler to gracefully handle a CTRL-C.
	*/

	keep_running=1;
	signal(SIGINT,interrupt_handler);

	progressbar *progress;

	if(config->progressbar==true)
		progress=progressbar_new("Progress",config->iterations/262144);
	else
		progress=NULL;

	/*
		We save the start time
	*/

	struct timeval starttime;
	gettimeofday(&starttime,NULL);

	/*
		This is the main DiagMC loop
	*/

	long int counter;
	for(counter=0;(counter<config->iterations)&&(keep_running==1);counter++)
	{
		int update_type,status,selector;

		selector=gsl_rng_uniform_int(amx->rng_ctx, cumulative_probability[DIAGRAM_NR_UPDATES-1]);
		update_type=-1;

		for(int c=0;c<DIAGRAM_NR_UPDATES;c++)
		{
			if(cumulative_probability[c]>selector)
			{
				update_type=c;
				break;
			}
		}

		assert(update_type!=-1);

		status=updates[update_type](amx, false);
		proposed[update_type]++;

		switch(status)
		{
			case UPDATE_ACCEPTED:
			accepted[update_type]++;
			break;

			case UPDATE_UNPHYSICAL:
			case UPDATE_REJECTED:
			rejected[update_type]++;
			break;

			case UPDATE_ERROR:
			default:
			assert(false);
		}

		nr_samples++;

		if(false)
		{
			printf("================================\n");
			amatrix_print(amx);
			printf("--------------------------------\n");
			pmatrix_print(amx->pmxs[0]);
			printf("--------------------------------\n");
			pmatrix_print(amx->pmxs[1]);
			printf("--------------------------------\n");

			if(amatrix_is_physical(amx)==true)
				printf("PHYSICAL\n");
			else
				printf("UNPHYSICAL\n");

			printf("================================\n");

			char nil[1024];
			gets(nil);
		}

		if(amatrix_is_physical(amx))
		{
			nr_physical_samples++;

			if((counter>config->thermalization)&&((nr_physical_samples%config->decorrelation)==0))
			{
				struct amatrix_weight_t ws=amatrix_weight(amx);
				double sign=(ws.weight>0.0f)?(1.0f):(-1.0f);
				int order=amx->pmxs[0]->dimensions;

				autocorrelation << ws.weight;
				overall_sign << sign;
				signs[order] << sign;

				for(int c=0;c<MAX_ORDER;c++)
				{
					orders[c] << ((c!=order) ? (0.0) : (sign));
					plus[c] << ((c!=order) ? (0.0) : (positive_part(sign)));
					minus[c] << ((c!=order) ? (0.0) : (negative_part(sign)));
				}

				{
					int excitation_level=(order<=2)?(0):(ws.excitation_level);

					excitations[excitation_level] << sign;
				}

				nr_samples_by_order[order]++;

				if(ws.weight>=0.0)
				{
					nr_positive_samples[order]++;
				}
				else
				{
					nr_negative_samples[order]++;
				}
			}
		}

		if((counter%262144)==0)
		{
			struct timeval now;
			double elapsedtime;

			if(config->progressbar==true)
				progressbar_inc(progress);

			gettimeofday(&now,NULL);

			elapsedtime=(now.tv_sec-starttime.tv_sec)*1000.0;
			elapsedtime+=(now.tv_usec-starttime.tv_usec)/1000.0;
			elapsedtime/=1000;

			if((config->timelimit>0.0f)&&(elapsedtime>config->timelimit))
				keep_running=0;
		}
	}

	if(keep_running==0)
	{
		printf("Caught SIGINT or time limit exceeded, exiting earlier.\n");
	}

	if(config->progressbar)
		progressbar_finish(progress);

	/*
		Now we print the statistics we collected to the output file in a nice way.
	*/

	fprintf(out,"# Diagrammatic Monte Carlo for Møller-Plesset theory\n");
	fprintf(out,"#\n");
	fprintf(out,"# Electron repulsion integrals loaded from '%s'\n",config->erisfile);
	fprintf(out,"# Output file is '%s'\n",output);
	fprintf(out,"# Binary compiled from git commit %s\n",GITCOMMIT);
	fprintf(out,"#\n");
	fprintf(out,"# Unphysical penalty: %f\n",config->unphysicalpenalty);
	fprintf(out,"# Minimum order: %d\n",config->minorder);
	fprintf(out,"# Maximum order: %d\n",config->maxorder);
	fprintf(out,"# Epsilon (for Lindelöf resummation): %f\n",config->epsilon);
	fprintf(out,"#\n");

	fprintf(out,"# Iterations (done/planned): %ld/%ld\n",counter,config->iterations);
	fprintf(out,"# Thermalization: %ld\n",config->thermalization);
	fprintf(out,"# Decorrelation: %d\n",config->decorrelation);
	fprintf(out,"# Iterations in the physical sector: %f%%\n",100.0f*((double)(nr_physical_samples))/((double)(nr_samples)));
	fprintf(out,"#\n");

	/*
		Here we calculate the elapsed time.
	*/

	struct timeval now;
	double elapsedtime;

	gettimeofday(&now,NULL);

	elapsedtime=(now.tv_sec-starttime.tv_sec)*1000.0;
	elapsedtime+=(now.tv_usec-starttime.tv_usec)/1000.0;
	elapsedtime/=1000;

	fprintf(out,"# Total time: %f seconds\n",elapsedtime);
	fprintf(out,"#\n");

	/*
		Now we print some update statistics
	*/

	long int total_proposed,total_accepted,total_rejected;
	total_proposed=total_accepted=total_rejected=0;

	fprintf(out,"# Update statistics:\n");

	for(int d=0;d<DIAGRAM_NR_UPDATES;d++)
	{
		fprintf(out,"# Update #%d (%s): ",d,update_names[d]);
		show_update_statistics(out,proposed[d],accepted[d],rejected[d]);

		total_proposed+=proposed[d];
		total_accepted+=accepted[d];
		total_rejected+=rejected[d];
	}

	fprintf(out,"# Total: ");
	show_update_statistics(out,total_proposed,total_accepted,total_rejected);
	fprintf(out,"#\n");

	/*
		Finally, we output the actual results.
	*/

	fprintf(out,"# <Order> <Positive physical samples> <Negative physical samples> <Percentage> <Sign> <Positive fraction> <Negative fraction> <Sign (from ALEA)>\n");

	alps::alea::autocorr_result<double> result_autocorrelation=autocorrelation.finalize();
	alps::alea::batch_result<double> result_overall_sign, result_signs[MAX_ORDER], result_plus[MAX_ORDER],result_minus[MAX_ORDER], result_orders[MAX_ORDER];
	alps::alea::batch_result<double> result_excitations[MAX_EXCITATIONS];

	result_overall_sign=overall_sign.finalize();

	for(int c=0;c<MAX_ORDER;c++)
		result_signs[c]=signs[c].finalize();

	for(int c=0;c<MAX_ORDER;c++)
		result_plus[c]=plus[c].finalize();

	for(int c=0;c<MAX_ORDER;c++)
		result_minus[c]=minus[c].finalize();

	for(int c=0;c<MAX_ORDER;c++)
		result_orders[c]=orders[c].finalize();

	for(int c=0;c<MAX_EXCITATIONS;c++)
		result_excitations[c]=excitations[c].finalize();

	long int total_positive,total_negative;

	total_positive=total_negative=0;
	for(int order=amx->config->minorder;order<=amx->config->maxorder;order++)
	{
		total_positive+=nr_positive_samples[order];
		total_negative+=nr_negative_samples[order];
	}

	for(int order=amx->config->minorder;order<=amx->config->maxorder;order++)
	{
		double pct,sign;

		if((total_positive-total_negative)!=0)
			pct=100.0f*((double)(nr_positive_samples[order]-nr_negative_samples[order]))/(total_positive-total_negative);
		else
			pct=NAN;

		if((nr_positive_samples[order]+nr_negative_samples[order])!=0)
			sign=((double)(nr_positive_samples[order]-nr_negative_samples[order]))/((double)(nr_positive_samples[order]+nr_negative_samples[order]));
		else
			sign=NAN;

		fprintf(out, "%d %ld %ld %f %f ", order, nr_positive_samples[order], nr_negative_samples[order], pct, sign);

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

	for(int c=0;c<MAX_EXCITATIONS;c++)
	{
		fprintf(out,"Excitation level(%d): %f +- %f\n",c,result_excitations[c].mean()(0),result_excitations[c].stderror()(0));
	}

	fprintf(out,"# Measured autocorrelation time = %f\n",result_autocorrelation.tau()(0));

	//auto divide = [] (double x,double y) -> double { return x/y; };
	//auto joined_data=alps::alea::join(result_signs[1],result_signs[2]);
	//auto transformer=alps::alea::make_transformer(std::function<double(double,double)>(divide));
	//auto ratio=alps::alea::transform(alps::alea::jackknife_prop(),transformer,joined_data);

	/*
		From the tutorial (https://github.com/ALPSCore/ALPSCore/wiki/Tutorial:-ALEA):

		Functions of multiple random variables (X,Y) can be realized by grouping the arguments
		together using alps::alea::join and then applying the transform on the combined result.
	*/

	//fprintf(out,"\nRatio is: %f\n",ratio.mean()(0));
	//std::cout << ratio.mean() << std::endl;

	fini_amatrix(amx);

	if(out)
		fclose(out);

	return 0;
}
