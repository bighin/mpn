#include <math.h>
#include <stdbool.h>
#include <jit/jit.h>

#include "amatrix.h"
#include "multiplicity.h"

#define CONST_FLOAT64(v) 	(jit_value_create_float64_constant(F, jit_type_float64, v))
#define CONST_INT(v) 		(jit_value_create_nint_constant(F, jit_type_int, v))

/*
	This is essentially the function reconstruct_weight() from mpn.c, in JIT form.
*/

jit_function_t weight_to_jit(struct amatrix_t *amx, struct amatrix_weight_t *awt,int labelindex1, int labelindex2, jit_context_t context)
{
	jit_context_build_start(context);

	/*
		Create function signature and object. Signature is: float64 (*)(void *, int, int)
	*/

	jit_type_t params[3] = {jit_type_void_ptr, jit_type_int, jit_type_int};
	jit_type_t signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, params, 3, 1);
	jit_function_t F = jit_function_create(context, signature);

	/*
		The variables that we are going to use
	*/

	jit_value_t ctxpointer, value1, value2, denominators, denominator, numerators;
	ctxpointer = jit_value_get_param(F, 0);
	value1 = jit_value_get_param(F, 1);
	value2 = jit_value_get_param(F, 2);
	denominators = jit_value_create(F, jit_type_float64);
	denominator = jit_value_create(F, jit_type_float64);
	numerators = jit_value_create(F, jit_type_float64);

	/*
		The signatures of the external functions we are going to call
	*/

	jit_type_t get_occupied_energy_params[] = {jit_type_void_ptr, jit_type_int};
	jit_type_t get_occupied_energy_signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, get_occupied_energy_params, 2, 1);

	jit_type_t get_virtual_energy_params[] = {jit_type_void_ptr, jit_type_int};
	jit_type_t get_virtual_energy_signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, get_virtual_energy_params, 2, 1);

	jit_type_t get_eri_params[] = {jit_type_void_ptr, jit_type_int, jit_type_int, jit_type_int, jit_type_int};
	jit_type_t get_eri_signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, get_eri_params, 5, 1);

	/*
		The actual calculations
	*/

	// denominators = 1.0f
	jit_insn_store(F, denominators, CONST_FLOAT64(1.0));

	for(int c=0;c<awt->nr_denominators;c++)
	{
		// denominator = 0.0f
		jit_insn_store(F, denominator, CONST_FLOAT64(0.0));

		for(int d=0;d<awt->denominators[c].ilabels;d++)
		{
			int label=awt->denominators[c].labels[d];
			int qtype=awt->denominators[c].qtypes[d];

			switch(qtype)
			{
				case QTYPE_OCCUPIED:

				if((label==labelindex1)||(label==labelindex2))
				{
					jit_value_t get_occupied_energy_args[2], tmp;

					get_occupied_energy_args[0]=ctxpointer;
					get_occupied_energy_args[1]=jit_insn_sub(F, (label==labelindex1) ? (value1) : (value2), CONST_INT(1));

					tmp = jit_insn_call_native(F, "get_occupied_energy", get_occupied_energy, get_occupied_energy_signature,
							           get_occupied_energy_args, 2, JIT_CALL_NOTHROW);

					// denominator += tmp
					denominator = jit_insn_add(F, denominator, tmp);
				}
				else
				{
					// denominator += get_occupied_energy(amx->ectx, awt->labels[label].value-1)
					denominator = jit_insn_add(F, denominator, CONST_FLOAT64(get_occupied_energy(amx->ectx, awt->labels[label].value-1)));
				}

				break;

				case QTYPE_VIRTUAL:

				if((label==labelindex1)||(label==labelindex2))
				{
					jit_value_t get_virtual_energy_args[2], tmp;

					get_virtual_energy_args[0]=ctxpointer;
					get_virtual_energy_args[1]=jit_insn_sub(F, (label==labelindex1) ? (value1) : (value2), CONST_INT(1));

					tmp = jit_insn_call_native(F, "get_virtual_energy", get_virtual_energy, get_virtual_energy_signature,
						                   get_virtual_energy_args, 2, JIT_CALL_NOTHROW);

					// denominator -= tmp
					denominator = jit_insn_sub(F, denominator, tmp);
				}
				else
				{
					// denominator -= get_virtual_energy(amx->ectx, awt->labels[label].value-1)
					denominator = jit_insn_sub(F, denominator, CONST_FLOAT64(get_virtual_energy(amx->ectx, awt->labels[label].value-1)));
				}

				break;
			}
		}


		// denominators *= denominator
		denominators = jit_insn_mul(F, denominators, denominator);
	}

	// numerators = 1.0f
	jit_insn_store(F, numerators, CONST_FLOAT64(1.0));

	for(int c=0;c<awt->nr_numerators;c++)
	{
		/*
			First we fetch the four label indices that enter into this ERI
		*/

		int ls[4];

		ls[0]=awt->numerators[c].labels[0];
		ls[1]=awt->numerators[c].labels[1];
		ls[2]=awt->numerators[c].labels[2];
		ls[3]=awt->numerators[c].labels[3];

		/*
			The we convert them to the format used by get_eri()
		*/

		int is[4];

		is[0]=awt->labels[ls[0]].value-1+((awt->labels[ls[0]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		is[1]=awt->labels[ls[1]].value-1+((awt->labels[ls[1]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		is[2]=awt->labels[ls[2]].value-1+((awt->labels[ls[2]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		is[3]=awt->labels[ls[3]].value-1+((awt->labels[ls[3]].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));

		bool is_numeric=true;

		for(int i=0;i<=3;i++)
			if((ls[i]==labelindex1)||(ls[i]==labelindex2))
				is_numeric=false;

		if(is_numeric==true)
		{
			double numerator=get_eri(amx->ectx, is[0], is[1], is[2], is[3]);

			// numerators *= numerator
			numerators = jit_insn_mul(F, numerators, CONST_FLOAT64(numerator));
		}
		else
		{
			jit_value_t vs[4];

			for(int i=0;i<=3;i++)
				vs[i]=jit_value_create(F, jit_type_int);

			for(int i=0;i<=3;i++)
			{
				if((ls[i]!=labelindex1)&&(ls[i]!=labelindex2))
				{
					jit_insn_store(F, vs[i], CONST_INT(is[i]));
				}
				else
				{
					if(ls[i]==labelindex1)
						jit_insn_store(F, vs[i], value1);
					else
						jit_insn_store(F, vs[i], value2);

					//vs[i] = vs[i] - 1
					vs[i]=jit_insn_sub(F, vs[i], CONST_INT(1));

					// vs[i] += (qtype==QTYPE_VIRTUAL)?(nr_occupied):(0)
					if(awt->labels[ls[i]].qtype==QTYPE_VIRTUAL)
						vs[i]=jit_insn_add(F, vs[i], CONST_INT(amx->ectx->nocc));
				}
			}

			jit_value_t get_eri_args[5], tmp;

			get_eri_args[0]=ctxpointer;
			get_eri_args[1]=vs[0];
			get_eri_args[2]=vs[1];
			get_eri_args[3]=vs[2];
			get_eri_args[4]=vs[3];

			tmp = jit_insn_call_native(F, "get_eri", get_eri, get_eri_signature,
				                   get_eri_args, 5, JIT_CALL_NOTHROW);

			// numerators *= tmp
			numerators = jit_insn_mul(F, numerators, tmp);
		}
	}

	// numerators *= awt->unphysical_penalty
	numerators = jit_insn_mul(F, numerators, CONST_FLOAT64(awt->unphysical_penalty));

	double lindeloef_factor=exp(amx->config->epsilon*awt->h*log(awt->h));
	double extrafactors=pow(awt->inversefactor,-1.0f)*lindeloef_factor;

	// weight = numerators/denominators*pow(awt->inversefactor,-1.0f)*lindeloef_factor
	jit_value_t weight = jit_insn_div(F, numerators, denominators);
	weight = jit_insn_mul(F, weight, CONST_FLOAT64(extrafactors));

	// weight = amx->config->bias + weight/amatrix_multiplicity(amx)
	weight = jit_insn_div(F, weight, CONST_FLOAT64(amatrix_multiplicity(amx)));
	weight = jit_insn_add(F, weight, CONST_FLOAT64(amx->config->bias));

	/*
		Final cleanup
	*/

	jit_insn_return(F, weight);
	jit_context_build_end(context);
	return F;
}
