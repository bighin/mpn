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
#include "weight.h"
#include "sampling.h"

#include "libprogressbar/progressbar.h"

/*
	The updates
*/

int update_extend(struct amatrix_t *amx, bool always_accept)
{
	double weightratio=1.0f/fabs(amatrix_weight(amx));
	double extend_probability,squeeze_probability;

	if(amx->pmxs[0]->dimensions>=amx->config->maxorder)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	extend_probability=1.0f/pow(amx->pmxs[0]->dimensions+1, 2.0f)/pow(amx->nr_occupied*amx->nr_virtual,2.0f);
	squeeze_probability=1.0f;

	int i1, j1, i2, j2;

	pmatrix_extend(amx->pmxs[0], amx->rng_ctx, &i1, &j1);
	pmatrix_extend(amx->pmxs[1], amx->rng_ctx, &i2, &j2);

	pmatrix_set_raw_entry(amx->pmxs[0],i1,j1,pmatrix_get_new_value(amx->pmxs[0],amx->rng_ctx,i1,j1));
	pmatrix_set_raw_entry(amx->pmxs[1],i2,j2,pmatrix_get_new_value(amx->pmxs[1],amx->rng_ctx,i2,j2));

	amx->cached_weight_is_valid=false;

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	double acceptance_ratio;
	double currentweight=amatrix_weight(amx);

	weightratio*=fabs(currentweight);
	acceptance_ratio=weightratio/extend_probability*squeeze_probability;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio)?(true):(false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_squeeze(struct amatrix_t *amx, bool always_accept)
{
	double weightratio=fabs(amatrix_weight(amx));
	double extend_probability,squeeze_probability;

	if(amx->pmxs[0]->dimensions<=amx->config->minorder)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	extend_probability=1.0f/pow(amx->pmxs[0]->dimensions, 2.0f)/pow(amx->nr_occupied*amx->nr_virtual,2.0f);
	squeeze_probability=1.0f;

	pmatrix_squeeze(amx->pmxs[0], amx->rng_ctx);
	pmatrix_squeeze(amx->pmxs[1], amx->rng_ctx);

	amx->cached_weight_is_valid=false;

	double acceptance_ratio;

	weightratio/=fabs(amatrix_weight(amx));
	acceptance_ratio=weightratio/extend_probability*squeeze_probability;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<(1.0f/acceptance_ratio))?(true):(false);

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

	double weightratio=1.0f/fabs(amatrix_weight(amx));

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

	switch(gsl_rng_uniform_int(amx->rng_ctx, 2))
	{
		case 0:
		pmatrix_swap_rows(target, i, j, amx->rng_ctx);
		break;

		case 1:
		pmatrix_swap_cols(target, i, j, amx->rng_ctx);
		break;
	}

	amx->cached_weight_is_valid=false;

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx));
	acceptance_ratio=weightratio;

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

	double weightratio=1.0f/fabs(amatrix_weight(amx));

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
				pmatrix_set_raw_entry(target, i, j, pmatrix_get_new_value(target, amx->rng_ctx, i, j));
	}

	amx->cached_weight_is_valid=false;

	/*
		The update is balanced with itself, the acceptance ratio is simply given
		by the (modulus of the) weights ratio.
	*/

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx));
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

	double weightratio=1.0f/fabs(amatrix_weight(amx));

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	/*
		Do we play with virtual or with occupied states?
	*/

	int target_type=(gsl_rng_uniform_int(amx->rng_ctx, 2)==0)?(QTYPE_OCCUPIED):(QTYPE_VIRTUAL);

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

	assert(icandidates<32);

	int values[32];

	for(int c=0;c<icandidates;c++)
		values[c]=pmatrix_get_raw_entry(amx->pmxs[candidates[c].pmatrix],candidates[c].i,candidates[c].j);

	fisher_yates(amx->rng_ctx,values,icandidates);

	for(int c=0;c<icandidates;c++)
		pmatrix_set_raw_entry(amx->pmxs[candidates[c].pmatrix],candidates[c].i,candidates[c].j,values[c]);

	amx->cached_weight_is_valid=false;

	/*
		The update is balanced with itself, the acceptance ratio is simply given
		by the (modulus of the) weights ratio.
	*/

	double acceptance_ratio;

	weightratio*=fabs(amatrix_weight(amx));
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

static volatile int keep_running=1;
static volatile int print_summary=0;

static void signal_handler(int signo)
{
	switch(signo)
	{

		case SIGINT:
		keep_running=0;
		break;

		case SIGUSR1:
		print_summary=1;
		break;

		default:
		break;
	}
}

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
		Update probabilities: note that they must be the same for complementary update pairs,
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
		We reset the update statistics and prepare a sampling context for the measurements
	*/

	for(int d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;

	struct sampling_ctx_t *sctx=init_sampling_ctx(config->maxorder);

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

	/*
		Let's extend the matrix until we hit the minimum allowed dimensions.
	*/

	assert(config->maxorder>config->minorder);
	assert(config->maxorder<MAX_ORDER);

	while(amx->pmxs[0]->dimensions<config->minorder)
		update_extend(amx,true);

	/*
		We setup a signal handler to gracefully handle a CTRL-C (i.e. SIGINT),
		and to print a short summary on SIGUSR1.
	*/

	keep_running=1;
	print_summary=0;

	signal(SIGINT,signal_handler);
	signal(SIGUSR1,signal_handler);

	/*
		We initialize the progress bar
	*/

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

		sampling_ctx_measure(sctx,amx,config,counter);

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

			if(print_summary==1)
			{
				sampling_ctx_print_report(sctx,amx,stdout,false);
				print_summary=0;
			}
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
	fprintf(out,"# Iterations in the physical sector: %f%%\n",sampling_ctx_get_physical_pct(sctx));
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
		Finally, we output the actual results...
	*/

	sampling_ctx_print_report(sctx,amx,out,true);

	/*
		...and we perform some final cleanups!
	*/

	fini_amatrix(amx,true);
	fini_sampling_ctx(sctx);

	if(out)
		fclose(out);

	return 0;
}
