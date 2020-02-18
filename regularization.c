#include <math.h>
#include <assert.h>

#include "regularization.h"

double regularize(struct amatrix_t *amx,double denominator)
{
	switch(amx->regularization)
	{
		case REGULARIZATION_TYPE_NONE:
		return denominator;

		/*
			alpha^p regularization
		*/

		case REGULARIZATION_TYPE_ALPHA:
		return denominator+pow(amx->alpha+denominator,-amx->p);

		/*
			sigma^p regularization
		*/

		case REGULARIZATION_TYPE_SIGMA:
		return denominator/(1.0f-exp(-amx->sigma*pow(denominator,amx->p)));

		default:
		assert(false);
	}

	return denominator;
}

char *regularization_type_description(int regularization)
{
	switch(regularization)
	{
		case REGULARIZATION_TYPE_NONE:
		return "none";

		case REGULARIZATION_TYPE_ALPHA:
		return "alpha^p";

		case REGULARIZATION_TYPE_SIGMA:
		return "sigma^p";

		default:
		assert(false);
	}

	return "unknown (error)";
}

char *resummation_type_description(int resummation)
{
	switch(resummation)
	{
		case RESUMMATION_TYPE_NONE:
		return "none";

		case RESUMMATION_TYPE_LINDELOEF:
		return "Lindelöf";

		default:
		assert(false);
	}

	return "unknown (error)";
}
