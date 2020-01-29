#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct energies_ctx_t
{
	int nso,nocc,nvirt;
	double *eocc,*evirt;
	double *htensor;
};

int htensor_index(int i,int j,int a,int b,int nocc,int nvirt)
{
	int ntot=nocc+nvirt;

	assert((i>=0)&&(i<ntot));
	assert((j>=0)&&(j<ntot));
	assert((a>=0)&&(a<ntot));
	assert((b>=0)&&(b<ntot));

	return ((i*ntot+j)*ntot+a)*ntot+b;
}

int htensor_size(int nocc,int nvirt)
{
	int ntot=nocc+nvirt;

	return 1+htensor_index(ntot-1,ntot-1,ntot-1,ntot-1,nocc,nvirt);
}

#define MAX_TOKENS		(1024)
#define TOKEN_MAX_LENGTH	(1024)

bool parse_tokens(char tokens[MAX_TOKENS][TOKEN_MAX_LENGTH],int nrtokens,struct energies_ctx_t *ctx)
{
	if(nrtokens<1)
		return false;

	if(strcmp(tokens[0],"nso")==0)
	{
		if(nrtokens!=2)
			return false;
		
		ctx->nso=atoi(tokens[1]);
	}
	else if(strcmp(tokens[0],"nocc")==0)
	{
		if(nrtokens!=2)
			return false;
		
		ctx->nocc=atoi(tokens[1]);
	}
	else if(strcmp(tokens[0],"nvirt")==0)
	{
		if(nrtokens!=2)
			return false;

		ctx->nvirt=atoi(tokens[1]);
	}
	else if(strcmp(tokens[0],"eocc")==0)
	{
		if(ctx->nocc==-1)
			return false;
		
		if(nrtokens!=(ctx->nocc+1))
			return false;

		if(ctx->eocc!=NULL)
			return false;
		
		ctx->eocc=malloc(sizeof(double)*ctx->nocc);

		for(int c=1;c<=ctx->nocc;c++)
			ctx->eocc[c-1]=atof(tokens[c]);
	}
	else if(strcmp(tokens[0],"evirt")==0)
	{
		if(ctx->nvirt==-1)
			return false;
		
		if(nrtokens!=(ctx->nvirt+1))
			return false;

		if(ctx->evirt!=NULL)
			return false;
		
		ctx->evirt=malloc(sizeof(double)*ctx->nvirt);

		for(int c=1;c<=ctx->nvirt;c++)
			ctx->evirt[c-1]=atof(tokens[c]);
	}
	else if(strcmp(tokens[0],"htensor")==0)
	{
		if(nrtokens!=6)
			return false;

		if((ctx->nocc==-1)||(ctx->nvirt==-1))
			return false;

		if(ctx->htensor==NULL)
		{
			/*
				First time, we allocate the tensor
			*/

			printf("Tensor size: %ld MiB\n",sizeof(double)*htensor_size(ctx->nocc,ctx->nvirt)/1024/1024);

			ctx->htensor=malloc(sizeof(double)*htensor_size(ctx->nocc,ctx->nvirt));
		}

		int i,j,a,b;

		i=atoi(tokens[1]);
		j=atoi(tokens[2]);
		a=atoi(tokens[3]);
		b=atoi(tokens[4]);

		ctx->htensor[htensor_index(i,j,a,b,ctx->nocc,ctx->nvirt)]=atof(tokens[5]);
	}
	else
	{
		return false;
	}

	return true;
}

bool load_energies(FILE *in, struct energies_ctx_t *ctx)
{
	int nrlines=0;
	
	ctx->nso=-1;
	ctx->nocc=-1;
	ctx->nvirt=-1;

	ctx->eocc=NULL;
	ctx->evirt=NULL;
	
	ctx->htensor=NULL;

	while((!feof(in))&&(!ferror(in)))
	{
		char line[1024];

		fgets(line,1024,in);
		line[1023]='\0';
		nrlines++;

		if(line[0]=='#')
			continue;

		if(strlen(line)==0)
			continue;

		/*
			We tokenize the string
		*/

		char *string,*tofree,*token;
		char tokens[MAX_TOKENS][TOKEN_MAX_LENGTH];
		int nrtokens=0;
		
		tofree=string=strdup(line);

		while((token=strsep(&string," "))!=NULL)
		{
			strncpy(tokens[nrtokens++],token,TOKEN_MAX_LENGTH);
		
			if(nrtokens>=MAX_TOKENS)
			{
				printf("Error parsing line %d, maximum number of tokens (%d) exceeded!",nrlines,MAX_TOKENS);
				break;
			}
		}

		free(tofree);

		if(parse_tokens(tokens,nrtokens,ctx)!=true)
		{
			printf("Error parsing line %d, skipping it!",nrlines);	
		}
	}
	
	/*
		On error, the caller is responsible for freeing non-null pointers.
	*/
	
	if((ctx->nso==-1)||(ctx->nocc==-1)||(ctx->nvirt==-1))
		return false;
	
	if(ctx->nso!=(ctx->nocc+ctx->nvirt))
		return false;
	
	if((ctx->eocc==NULL)||(ctx->evirt==NULL))
		return false;

	if(ctx->htensor==NULL)
		return false;

	return true;
}

/*
	Remember that in this context the indices can take the following values:

	i 	[0,nocc-1]
	j 	[0,nocc-1]
	a 	[nocc,nocc+nvirt-1]
	b 	[nocc,nocc+virt-1]

	This is different from the rest of the code, where all the indices are shifted by 1,
	or even arranged in different ways.
*/

double get_occupied_energy(struct energies_ctx_t *ctx,int n)
{
	assert((n>=0)&&(n<ctx->nocc));

	return ctx->eocc[n];
}

double get_virtual_energy(struct energies_ctx_t *ctx,int n)
{
	assert((n>=0)&&(n<ctx->nvirt));

	return ctx->evirt[n];
}

double get_htensor(struct energies_ctx_t *ctx,int i,int j,int a,int b)
{
	assert((i>=0)&&(i<(ctx->nocc+ctx->nvirt)));
	assert((j>=0)&&(j<(ctx->nocc+ctx->nvirt)));
	assert((a>=0)&&(a<(ctx->nocc+ctx->nvirt)));
	assert((b>=0)&&(b<(ctx->nocc+ctx->nvirt)));

	return ctx->htensor[htensor_index(i,j,a,b,ctx->nocc,ctx->nvirt)];
}
