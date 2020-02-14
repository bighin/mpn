#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdbool.h>

struct configuration_t
{
	/* "general" section */

	char *prefix;
	bool progressbar;
	char *erisfile;

	/* "parameters" section */

	double bias;
	double unphysicalpenalty;
	int minorder,maxorder;

	/* "sampling" section */

	long int iterations;
	int thermalization;
	double timelimit;
	int decorrelation;

	/* The name of the file the configuration has been loaded from */

	char *ininame;
};

int configuration_handler(void *user,const char *section,const char *name,const char *value);
void load_config_defaults(struct configuration_t *config);
bool load_configuration(char *configfile,struct configuration_t *config);

#endif //__CONFIG_H__
