#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "amatrix.h"
#include "mc.h"
#include "config.h"
#include "cache.h"

void usage(char *argv0)
{
	printf("Usage: %s <inifile> [<otherinifiles> ...]\n",argv0);

	exit(0);
}

int main(int argc,char *argv[])
{
	if(argc<2)
		usage(argv[0]);

	bool first=true;

	for(int c=1;c<argc;c++)
	{
		struct configuration_t config;

		if(first==true)
		{
			printf("Diagrammatic Monte Carlo for Møller-Plesset theory.\n");

			amatrix_cache_is_enabled=true;
			init_cache(6);
		}

		load_config_defaults(&config);

		if(load_configuration(argv[c],&config)==false)
			continue;

		do_diagmc(&config);
		first=false;

		if((c+1)!=argc)
			printf("\n");
	}
}
