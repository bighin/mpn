#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdbool.h>

struct configuration_t
{
	/* "general" section */

	char *prefix;
	bool progressbar;
	char *erisfile;
	bool seedrng;

	/* "parameters" section */

	double unphysicalpenalty;
	int minorder,maxorder;
	double epsilon;

	/* "sampling" section */

	long int iterations;
	long int thermalization;
	double timelimit;
	int decorrelation;

	/* The name of the file the configuration has been loaded from */

	char *inipath;
};

int configuration_handler(void *user,const char *section,const char *name,const char *value);
void load_config_defaults(struct configuration_t *config);
bool load_configuration(char *configfile,struct configuration_t *config);

#endif //__CONFIG_H__
