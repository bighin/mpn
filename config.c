#include <string.h>
#include <stdlib.h>

#include "config.h"
#include "auxx.h"
#include "inih/ini.h"

int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);

#define MATCH(s,n) ((strcmp(section,s)==0)&&(strcmp(name,n)==0))

	if(MATCH("general","prefix"))
	{
		if(strcmp(value,"auto")==0)
		{
			if(!strstr(pconfig->inipath,".ini"))
			{
				printf("Error: using automatic prefix, but the configuration file path does not contain '.ini'\n");
				exit(0);
			}

			pconfig->prefix=find_and_replace(pconfig->inipath,".ini","");
		}
		else
		{
			pconfig->prefix=strdup(value);
		}
	}
	else if(MATCH("general","progressbar"))
	{
		if(!strcmp(value,"true"))
			pconfig->progressbar=true;
		else if(!strcmp(value,"false"))
			pconfig->progressbar=false;
		else
			return 0;
	}
	else if(MATCH("general","erisfile"))
	{
		pconfig->erisfile=strdup(value);
	}
	else if(MATCH("general","seedrng"))
	{
		if(!strcmp(value,"true"))
			pconfig->seedrng=true;
		else if(!strcmp(value,"false"))
			pconfig->seedrng=false;
		else
			return 0;
	}
	else if(MATCH("parameters","unphysicalpenalty"))
	{
		pconfig->unphysicalpenalty=atof(value);
	}
	else if(MATCH("parameters","minorder"))
	{
		pconfig->minorder=atoi(value);
	}
	else if(MATCH("parameters","maxorder"))
	{
		pconfig->maxorder=atoi(value);
	}
	else if(MATCH("parameters","epsilon"))
	{
		pconfig->epsilon=atof(value);
	}
	else if(MATCH("parameters","maxtau"))
	{
		pconfig->maxtau=atof(value);
	}
	else if(MATCH("parameters","chempot"))
	{
		pconfig->chempot=atof(value);
	}
	else if(MATCH("sampling","iterations"))
	{
		pconfig->iterations=(long int)(strtol(value,(char **)NULL,10));
	}
	else if(MATCH("sampling","thermalization"))
	{
		pconfig->thermalization=(long int)(strtol(value,(char **)NULL,10));
	}
	else if(MATCH("sampling","timelimit"))
	{
		pconfig->timelimit=atof(value);
	}
	else if(MATCH("sampling","decorrelation"))
	{
		pconfig->decorrelation=atoi(value);
	}
	else if(MATCH("sampling","nrbins"))
	{
		pconfig->nrbins=atoi(value);
	}
	else if(MATCH("sampling","binwidth"))
	{
		pconfig->binwidth=atof(value);
	}
	else
	{
		/* Unknown section/name, error */
		return 0;
	}

	return 1;
}

void load_config_defaults(struct configuration_t *config)
{
	config->prefix=strdup("default");
	config->progressbar=false;
	config->erisfile=NULL;
	config->seedrng=true;

	config->unphysicalpenalty=0.01f;
	config->minorder=1;
	config->minorder=8;
	config->epsilon=0.0f;
	config->maxtau=100.0f;
	config->chempot=0.0f;

	config->iterations=10000000;
	config->thermalization=config->iterations/100;
	config->timelimit=0.0f;
	config->decorrelation=10;
	config->nrbins=100;
	config->binwidth=0.1;

	config->inipath=NULL;
}

bool load_configuration(char *configfile,struct configuration_t *config)
{
	load_config_defaults(config);
	config->inipath=strdup(configfile);

	if(ini_parse(configfile,configuration_handler,config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		return false;
	}

	printf("Loaded '%s'\n",configfile);
	return true;
}
