#include <string.h>
#include <stdlib.h>

#include "config.h"
#include "inih/ini.h"

int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);

#define MATCH(s,n) ((strcmp(section,s)==0)&&(strcmp(name,n)==0))

	if(MATCH("general","prefix"))
	{
		pconfig->prefix=strdup(value);
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
	else if(MATCH("parameters","bias"))
	{
		pconfig->bias=atof(value);
	}
	else if(MATCH("parameters","unphysicalpenalty"))
	{
		pconfig->unphysicalpenalty=atof(value);
	}
	else if(MATCH("parameters","minorder"))
	{
		pconfig->minorder=atof(value);
	}
	else if(MATCH("parameters","maxorder"))
	{
		pconfig->maxorder=atof(value);
	}
	else if(MATCH("sampling","iterations"))
	{
		pconfig->iterations=atoi(value);
	}
	else if(MATCH("sampling","thermalization"))
	{
		pconfig->thermalization=atoi(value);
	}
	else if(MATCH("sampling","timelimit"))
	{
		pconfig->timelimit=atof(value);
	}
	else if(MATCH("sampling","decorrelation"))
	{
		pconfig->decorrelation=atoi(value);
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
	config->erisfile=strdup("H2O.dat");

	config->bias=0.0f;
	config->unphysicalpenalty=0.01f;
	config->minorder=1;
	config->minorder=16;

	config->iterations=10000000;
	config->thermalization=config->iterations/100;
	config->timelimit=0.0f;
	config->decorrelation=10;

	config->ininame=NULL;
}

bool load_configuration(char *configfile,struct configuration_t *config)
{
	load_config_defaults(config);

	if(ini_parse(configfile,configuration_handler,config)<0)
	{
		fprintf(stderr,"Couldn't read or parse '%s'\n",configfile);
		return false;
	}

	config->ininame=strdup(configfile);

	fprintf(stderr,"Loaded '%s'\n",configfile);
	return true;
}
