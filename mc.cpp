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
#include <gsl/gsl_rng.h>

#include "mc.h"
#include "amatrix.h"
#include "auxx.h"
#include "mpn.h"
#include "config.h"

#include "libprogressbar/progressbar.h"

/*
	The updates
*/

int update_extend(struct amatrix_t *amx, bool always_accept)
{
	double weightratio=1.0f/amatrix_weight(amx).weight;
	double probability=1.0f;

	if(amx->pmxs[0]->dimensions>=amx->config->maxorder)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	int dimensions=amx->pmxs[0]->dimensions;
	int i,j;

	switch(pmatrix_extend(amx->pmxs[0], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability*=((dimensions+1)*amx->nr_virtual*amx->nr_occupied);
		probability/=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)?(amx->nr_occupied):(amx->nr_virtual);
		break;

		case 2:
		assert(pmatrix_entry_type(dimensions,dimensions)==QTYPE_OCCUPIED);
		probability*=((dimensions+1)*amx->nr_occupied);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	switch(pmatrix_extend(amx->pmxs[1], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability*=((dimensions+1)*amx->nr_virtual*amx->nr_occupied);
		probability/=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED) ? (amx->nr_occupied) : (amx->nr_virtual);
		break;

		case 2:
		assert(pmatrix_entry_type(dimensions,dimensions)==QTYPE_OCCUPIED);
		probability*=((dimensions+1)*amx->nr_occupied);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	amx->cached_weight_is_valid=false;

	/*
		Finally we calculate the acceptance ratio for the update.
	*/

	double acceptance_ratio;

	weightratio*=amatrix_weight(amx).weight;
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
	double weightratio=1.0f/amatrix_weight(amx).weight;
	double probability=1.0f;

	if(amx->pmxs[0]->dimensions<=amx->config->minorder)
		return UPDATE_UNPHYSICAL;

	struct amatrix_backup_t backup;
	amatrix_save(amx, &backup);

	int dimensions=amx->pmxs[0]->dimensions;
	int i,j;

	switch(pmatrix_squeeze(amx->pmxs[0], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability/=(dimensions*amx->nr_virtual*amx->nr_occupied);
		probability*=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED)?(amx->nr_occupied):(amx->nr_virtual);
		break;

		case 2:
		assert(pmatrix_entry_type(dimensions,dimensions)==QTYPE_OCCUPIED);
		probability/=(dimensions*amx->nr_occupied);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	switch(pmatrix_squeeze(amx->pmxs[1], amx->rng_ctx, &i, &j))
	{
		case 1:
		probability/=(dimensions*amx->nr_virtual*amx->nr_occupied);
		probability*=(pmatrix_entry_type(i,j)==QTYPE_OCCUPIED) ? (amx->nr_occupied) : (amx->nr_virtual);
		break;

		case 2:
		assert(pmatrix_entry_type(dimensions,dimensions)==QTYPE_OCCUPIED);
		probability/=(dimensions*amx->nr_occupied);
		break;

		case -1:
		default:
		assert(false);
		return UPDATE_ERROR;
	}

	amx->cached_weight_is_valid=false;

	double acceptance_ratio;

	weightratio*=amatrix_weight(amx).weight;
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

	/*
		If we are in dimension 1, there's nothing to shuffle.
	*/

	if(dimensions==1)
		return UPDATE_UNPHYSICAL;

	double weightratio=1.0f/amatrix_weight(amx).weight;

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

	weightratio*=amatrix_weight(amx).weight;
	acceptance_ratio=fabs(weightratio);
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

	double weightratio=1.0f/amatrix_weight(amx).weight;

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

	weightratio*=amatrix_weight(amx).weight;
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

#define DIAGRAM_NR_UPDATES        (4)

	int (*updates[DIAGRAM_NR_UPDATES])(struct amatrix_t *amx, bool always_accept);
	const char *update_names[DIAGRAM_NR_UPDATES];

	long int proposed[DIAGRAM_NR_UPDATES], accepted[DIAGRAM_NR_UPDATES], rejected[DIAGRAM_NR_UPDATES];

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
		We initialize the accumulators and variables we will use for performing measurements.
	*/

	assert(config->maxorder>config->minorder);
	assert(config->maxorder<256);

	alps::alea::autocorr_acc<double> autocorrelation(1);
	alps::alea::batch_acc<double> signs[256], lsigns[16][16], hsigns[16][16], lhsigns[16][16];

	long int nr_samples, nr_physical_samples, nr_samples_by_order[256], nr_positive_samples[256], nr_negative_samples[256];

	long int nr_positive_samples_l[256][256], nr_negative_samples_l[256][256];
	long int nr_positive_samples_h[256][256], nr_negative_samples_h[256][256];
	long int nr_positive_samples_lh[256][256], nr_negative_samples_lh[256][256];

	nr_samples=nr_physical_samples=0;

	for(int c=0;c<256;c++)
		nr_samples_by_order[c]=nr_positive_samples[c]=nr_negative_samples[c]=0;

	for(int c=0;c<256;c++)
	{
		for(int d=0;d<256;d++)
		{
			nr_positive_samples_l[c][d]=nr_negative_samples_l[c][d]=0;
			nr_positive_samples_h[c][d]=nr_negative_samples_h[c][d]=0;
			nr_positive_samples_lh[c][d]=nr_negative_samples_lh[c][d]=0;
		}
	}

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

		nr_samples++;

		if(amatrix_is_physical(amx))
		{
			nr_physical_samples++;

			if((counter>config->thermalization)&&((nr_physical_samples%config->decorrelation)==0))
			{
				struct amatrix_weight_t ws=amatrix_weight(amx);
				double sign=(ws.weight>0.0f)?(1.0f):(-1.0f);
				int order=amx->pmxs[0]->dimensions;
				int l=ws.l;
				int h=ws.h;

				autocorrelation << ws.weight;
				signs[order] << sign;
				lsigns[order][l] << sign;
				hsigns[order][h] << sign;
				lhsigns[order][l+h] << sign;

				nr_samples_by_order[order]++;

				if(ws.weight>=0.0)
				{
					nr_positive_samples[order]++;
					nr_positive_samples_l[order][l]++;
					nr_positive_samples_h[order][h]++;
					nr_positive_samples_lh[order][l+h]++;
				}
				else
				{
					nr_negative_samples[order]++;
					nr_negative_samples_l[order][l]++;
					nr_negative_samples_h[order][h]++;
					nr_negative_samples_lh[order][l+h]++;
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
	fprintf(out,"#\n");
	fprintf(out,"# Bias: %f\n",config->bias);
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

	fprintf(out,"# <Order> <Positive physical samples> <Negative physical samples> <Percentage> <Sign> <Sign (from ALEA)>\n");

	alps::alea::autocorr_result<double> result_autocorrelation=autocorrelation.finalize();
	alps::alea::batch_result<double> result_signs[256], result_lsigns[16][16],result_hsigns[16][16], result_lhsigns[16][16];

	for(int c=0;c<256;c++)
		result_signs[c]=signs[c].finalize();

	for(int c=0;c<16;c++)
	{
		for(int d=0;d<16;d++)
		{
			result_lsigns[c][d]=lsigns[c][d].finalize();
			result_hsigns[c][d]=hsigns[c][d].finalize();
			result_lhsigns[c][d]=lhsigns[c][d].finalize();
		}
	}

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

		/*
			Remember that result_signs[order].mean() has type alps::alea::column
			which is a shorthand for Eigen::column
		*/

		fprintf(out, "%f+-%f\n",result_signs[order].mean()(0),result_signs[order].stderror()(0));
	}

	fprintf(out,"# Order-by-order ratios:\n");

	for(int order1=amx->config->minorder;order1<=amx->config->maxorder;order1++)
	{
		for(int order2=amx->config->minorder;order2<=amx->config->maxorder;order2++)
		{
			if(order1==order2)
				continue;

			double pct1,pct2;

			pct1=((double)(nr_positive_samples[order1]-nr_negative_samples[order1]))/(total_positive-total_negative);
			pct2=((double)(nr_positive_samples[order2]-nr_negative_samples[order2]))/(total_positive-total_negative);

			char desc1[128],desc2[128];

			order_description(desc1,128,order1);
			order_description(desc2,128,order2);

			fprintf(out,"%s/%s %f\n", desc1, desc2, pct1/pct2);
		}
	}

	fprintf(out,"# Measured autocorrelation time = %f\n",result_autocorrelation.tau()(0));

	fprintf(out,"L:\n");
	for(int c=0;c<16;c++)
	{
		double sign=((double)(nr_positive_samples[c]-nr_negative_samples[c]))/((double)(nr_positive_samples[c]+nr_negative_samples[c]));
		fprintf(out, "(%f[%f]) ", result_signs[c].stderror()(0)/result_signs[c].mean()(0), sign);

		for(int d=0;d<16;d++)
		{
			sign=((double)(nr_positive_samples_l[c][d]-nr_negative_samples_l[c][d]))/((double)(nr_positive_samples_l[c][d]+nr_negative_samples_l[c][d]));
			fprintf(out, "%f[%f] ", result_lsigns[c][d].stderror()(0)/result_lsigns[c][d].mean()(0), sign);
		}

		fprintf(out,"\n");
	}

	fprintf(out,"H:\n");
	for(int c=0;c<16;c++)
	{
		double sign=((double)(nr_positive_samples[c]-nr_negative_samples[c]))/((double)(nr_positive_samples[c]+nr_negative_samples[c]));
		fprintf(out, "(%f[%f]) ", result_signs[c].stderror()(0)/result_signs[c].mean()(0), sign);

		for(int d=0;d<16;d++)
		{
			sign=((double)(nr_positive_samples_h[c][d]-nr_negative_samples_h[c][d]))/((double)(nr_positive_samples_h[c][d]+nr_negative_samples_h[c][d]));
			fprintf(out, "%f[%f] ", result_hsigns[c][d].stderror()(0)/result_hsigns[c][d].mean()(0), sign);
		}

		fprintf(out,"\n");
	}

	fprintf(out,"LH:\n");
	for(int c=0;c<16;c++)
	{
		double sign=((double)(nr_positive_samples[c]-nr_negative_samples[c]))/((double)(nr_positive_samples[c]+nr_negative_samples[c]));
		fprintf(out, "(%f[%f]) ", result_signs[c].stderror()(0)/result_signs[c].mean()(0), sign);

		for(int d=0;d<16;d++)
		{
			sign=((double)(nr_positive_samples_lh[c][d]-nr_negative_samples_lh[c][d]))/((double)(nr_positive_samples_lh[c][d]+nr_negative_samples_lh[c][d]));
			fprintf(out, "%f[%f] ", result_lhsigns[c][d].stderror()(0)/result_lhsigns[c][d].mean()(0), sign);
		}

		fprintf(out,"\n");
	}

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
