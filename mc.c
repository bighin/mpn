#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

#include "mc.h"
#include "amatrix.h"
#include "auxx.h"
#include "mpn.h"

#include "libprogressbar/progressbar.h"

/*
	The updates
*/

int update_extend(struct amatrix_t *amx, bool always_accept)
{
	double weightratio=1.0f/amatrix_weight(amx);
	double probability=1.0f;

	if((amx->pmxs[0]->dimensions+1)>=amx->max_order)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	int dimensions=amx->pmxs[0]->dimensions;
	int i,j;

	switch(pmatrix_extend(amx->pmxs[0], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability/=((dimensions+1)*amx->nr_virtual*amx->nr_occupied);
		probability*=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)?(amx->nr_occupied):(amx->nr_virtual);
		break;

		case 2:
		assert(pmatrix_entry_type(i,j)==QTYPE_VIRTUAL);
		probability/=((dimensions+1)*amx->nr_virtual);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	switch(pmatrix_extend(amx->pmxs[1], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability/=((dimensions+1)*amx->nr_virtual*amx->nr_occupied);
		probability*=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED) ? (amx->nr_occupied) : (amx->nr_virtual);
		break;

		case 2:
		assert(pmatrix_entry_type(i,j)==QTYPE_VIRTUAL);
		probability/=((dimensions+1)*amx->nr_virtual);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	double acceptance_ratio;

	weightratio*=amatrix_weight(amx);
	acceptance_ratio=fabs(weightratio)*probability;

	bool is_accepted=(gsl_rng_uniform(amx->rng_ctx)<acceptance_ratio) ? (true) : (false);

	if((is_accepted==false)&&(always_accept==false))
	{
		amatrix_restore(amx, &backup);
		return UPDATE_REJECTED;
	}

	return UPDATE_ACCEPTED;
}

int update_squeeze(struct amatrix_t *amx, bool always_accept)
{
	double weightratio=1.0f/amatrix_weight(amx);
	double probability=1.0f;

	if(amx->pmxs[0]->dimensions<=2)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	int dimensions=amx->pmxs[0]->dimensions;
	int i,j;

	switch(pmatrix_squeeze(amx->pmxs[0], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability*=(dimensions*amx->nr_virtual*amx->nr_occupied);
		probability/=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)?(amx->nr_occupied):(amx->nr_virtual);
		break;

		case 2:
		probability*=(dimensions*amx->nr_virtual);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	switch(pmatrix_squeeze(amx->pmxs[1], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability*=(dimensions*amx->nr_virtual*amx->nr_occupied);
		probability/=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED) ? (amx->nr_occupied) : (amx->nr_virtual);
		break;

		case 2:
		probability*=(dimensions*amx->nr_virtual);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	double acceptance_ratio;

	weightratio*=amatrix_weight(amx);
	acceptance_ratio=fabs(weightratio)*probability;

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

	double weightratio=1.0f/amatrix_weight(amx);

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

	double acceptance_ratio;

	weightratio*=amatrix_weight(amx);
	acceptance_ratio=fabs(weightratio);
	acceptance_ratio/=pow(amx->nr_occupied, update[QTYPE_OCCUPIED]);
	acceptance_ratio/=pow(amx->nr_virtual, update[QTYPE_VIRTUAL]);
	acceptance_ratio*=pow(amx->nr_occupied, reverse[QTYPE_OCCUPIED]);
	acceptance_ratio*=pow(amx->nr_virtual, reverse[QTYPE_VIRTUAL]);

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

	double weightratio=1.0f/amatrix_weight(amx);

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

	/*
		The update is balanced with itself, the acceptance ratio is simply given
		by the (modulus of the) weights ratio.
	*/

	double acceptance_ratio;

	weightratio*=amatrix_weight(amx);
	acceptance_ratio=fabs(weightratio);

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

void show_update_statistics(FILE *out,int proposed,int accepted,int rejected)
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

	fprintf(out,"proposed %d, accepted %d (%f%%), rejected %d (%f%%).\n",proposed,accepted,accepted_pct,rejected,rejected_pct);
}

static volatile int keep_running=1;

void interrupt_handler(int dummy __attribute__((unused)))
{
	keep_running=0;
}

/*
	The actual DiagMC routine.
*/

int do_diagmc(char *energies_dot_dat,char *output,int iterations,double timelimit,bool show_progressbar)
{

#define DIAGRAM_NR_UPDATES	(4)

	int (*updates[DIAGRAM_NR_UPDATES])(struct amatrix_t *amx,bool always_accept);
	char *update_names[DIAGRAM_NR_UPDATES];

	int proposed[DIAGRAM_NR_UPDATES],accepted[DIAGRAM_NR_UPDATES],rejected[DIAGRAM_NR_UPDATES];
	int total_proposed,total_accepted,total_rejected;

	/*
		We set up the updates we will be using
	*/

	updates[0]=update_extend;
	updates[1]=update_squeeze;
	updates[2]=update_shuffle;
	updates[3]=update_modify;

	update_names[0]="Extend";
	update_names[1]="Squeeze";
	update_names[2]="Shuffle";
	update_names[3]="Modify";

	/*
		We reset the update statistics
	*/

	for(int d=0;d<DIAGRAM_NR_UPDATES;d++)
		proposed[d]=accepted[d]=rejected[d]=0;

	/*
		We initialize the variables we will use for sampling.
	*/

	double average_weight,average_squared_weight;
	int nr_positive,nr_negative,nr_samples,nr_samples_including_unphysical;

	average_weight=average_squared_weight=0.0f;
	nr_positive=nr_negative=nr_samples=nr_samples_including_unphysical=0;

	/*
		We print some informative message, and then we open the log file
	*/

	fprintf(stderr,"Performing %d iterations\n",iterations);

	FILE *out;
	if(!(out=fopen(output,"w+")))
	{
		fprintf(stderr,"Error: couldn't open %s for writing\n",output);
		return 0;
	}

	fprintf(stderr,"Writing results to '%s'\n",output);

	/*
		The diagram parameters are loaded from the configuration, and then a new 'amatrix' is created
	*/

	struct amatrix_t *amx=init_amatrix(energies_dot_dat);

	amx->bias=0.0f;
	amx->unphysical_penalty=0.01f;
	amx->max_order=2;

	/*
		We setup an interrupt handler to gracefully handle a CTRL-C, and initialize a structure needed
		by the ncurses library to return info about the current terminal.
	*/

	keep_running=1;
	signal(SIGINT,interrupt_handler);

	progressbar *progress;

	if(show_progressbar==true)
		progress=progressbar_new("Progress",iterations/16384);
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

	for(int c=0;(c<iterations)&&(keep_running==1);c++)
	{
		int update_type,status;

		update_type=gsl_rng_uniform_int(amx->rng_ctx, DIAGRAM_NR_UPDATES);
		status=updates[update_type](amx, false);
		proposed[update_type]++;

		switch(status)
		{
			case UPDATE_ACCEPTED:
			accepted[update_type]++;
			break;

			/*
				FIXME:

				Here unphysical and rejected updates are treated in exactly the same way.
				Is this OK?
			*/

			case UPDATE_UNPHYSICAL:
			case UPDATE_REJECTED:
			rejected[update_type]++;
			break;

			case UPDATE_ERROR:
			default:
			assert(false);
		}

		nr_samples_including_unphysical++;

		if(amatrix_is_physical(amx))
		{
			/*
				TODO: one can get the weight for free since it has already been calculated!
				We can just cache the result!
			*/

			double weight=amatrix_weight(amx);

			//printf("%d %d %d %d || %f\n",pmatrix_get_entry(amx->pmxs[0],0,1),
			//                             pmatrix_get_entry(amx->pmxs[0],1,0),
			//                             pmatrix_get_entry(amx->pmxs[1],0,1),
			//                             pmatrix_get_entry(amx->pmxs[1],1,0),weight);

			average_weight+=weight;
			average_squared_weight+=weight*weight;
			nr_samples++;

			if(weight>0.0f)
				nr_positive++;
			else
				nr_negative++;
		}

		if((c%16384)==0)
		{
			struct timeval now;
			double elapsedtime;

			if(show_progressbar==true)
				progressbar_inc(progress);

			gettimeofday(&now,NULL);

			elapsedtime=(now.tv_sec-starttime.tv_sec)*1000.0;
			elapsedtime+=(now.tv_usec-starttime.tv_usec)/1000.0;
			elapsedtime/=1000;

			if((timelimit>0.0f)&&(elapsedtime>timelimit))
				keep_running=0;
		}
	}

	if(keep_running==0)
	{
		printf("Caught SIGINT or time limit exceeded, exiting earlier.\n");
	}

	if(show_progressbar)
		progressbar_finish(progress);

	fini_amatrix(amx);

	/*
		Now we print the statistics we collected to the output file in a nice way.
	*/

	fprintf(out,"# Diagrammatic Monte Carlo for Møller-Plesset theory\n");
	fprintf(out,"#\n");
	fprintf(out,"# Electron repulsion integrals loaded from '%s'\n",energies_dot_dat);
	fprintf(out,"# Output file is '%s'\n",output);
	fprintf(out,"#\n");
	fprintf(out,"#\n");

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
	fprintf(out,"# Sampled quantity: Correlation energy\n");
	fprintf(out,"#\n");
	fprintf(out,"# Iterations: %d\n",iterations);
	fprintf(out,"# Samples in the physical sector: %f%%\n",100.0f*((double)(nr_samples))/((double)(nr_samples_including_unphysical)));

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
		Finally, we output the actual results.
	*/

	average_weight/=nr_samples;
	average_squared_weight/=nr_samples;

#warning Check the following two formulas!

	double error_on_average=(average_squared_weight-average_weight*average_weight)/sqrt(nr_samples);
	double sign=((double)(nr_positive-nr_negative))/((double)(nr_samples));

	fprintf(out,"# <Correlation energy> <Error on correlation energy> <Sign>\n");
	fprintf(out,"%f %f %f\n",average_weight,error_on_average,sign);

	if(out)
		fclose(out);

	return 0;
}
