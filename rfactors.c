#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "amatrix.h"
#include "cache.h"
#include "rfactors.h"

/*
        This is not thread safe!
*/

long int rfactor4global[2];
long int rfactors4[576][2];

void rfactors_init(void)
{
}

void rfactors_sample_sign(struct amatrix_t *amx, int sign)
{
	assert((sign==1)||(sign==-1));

	if(amx->pmxs[0]->dimensions!=4)
		return;

	int index=amatrix_to_index(amx);

	if(sign==1)
	{
		rfactor4global[0]++;
		rfactors4[index][0]++;
	}
	else
	{
		rfactor4global[1]++;
		rfactors4[index][1]++;
	}
}

void rfactors_output_summary(const char *filename)
{
	FILE *out=fopen(filename,"w+");

	{
		double num,den;

		num=labs(rfactor4global[0]-rfactor4global[1]);
		den=rfactor4global[0]+rfactor4global[1];

		fprintf(out,"# GLOBAL: %f\n",num/den);
	}

	for(int c=0;c<576;c++)
	{
		if((rfactors4[c][0]+rfactors4[c][1])==0)
			continue;

		double num,den;

		num=labs(rfactors4[c][0]-rfactors4[c][1]);
		den=rfactors4[c][0]+rfactors4[c][1];

		fprintf(out,"%d %f # %s\n",c,num/den,(rfactors4[c][0]>rfactors4[c][1])?("positive"):("negative"));
	}

	if(out)
		fclose(out);
}
