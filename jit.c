#include <math.h>
#include <jit/jit.h>

#include "amatrix.h"
#include "multiplicity.h"

#define CONST_FLOAT64(v) (jit_value_create_float64_constant(F, jit_type_float64, v))
#define CONST_INT(v) (jit_value_create_float64_constant(F, jit_type_int, v))
#define CONST_VOID_PTR(v) (jit_value_create_void_ptr_constant(F, jit_type_void_ptr, v))

jit_function_t weight_to_jit(struct amatrix_t *amx,struct amatrix_weight_t *awt,int label1,int label2,jit_context_t context)
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

	jit_value_t pointer, value1, value2, denominators, denominator, numerators, f64tmp, itmp, result;
	pointer = jit_value_get_param(F, 0);
	value1 = jit_value_get_param(F, 1);
	value2 = jit_value_get_param(F, 2);
	denominators = jit_value_create(F, jit_type_float64);
	denominator = jit_value_create(F, jit_type_float64);
	numerators = jit_value_create(F, jit_type_float64);
	f64tmp = jit_value_create(F, jit_type_float64);
	itmp = jit_value_create(F, jit_type_int);
	// Can I remove the next line?
	result = jit_value_create(F, jit_type_float64);

	/*
		The signatures of the external functions we are going to call
	*/

	jit_type_t get_occupied_energy_params[] = {jit_type_void_ptr, jit_type_int};
	jit_type_t get_occupied_energy_signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, get_occupied_energy_params, 2, 1);

	jit_type_t get_virtual_energy_params[] = {jit_type_void_ptr, jit_type_int};
	jit_type_t get_virtual_energy_signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, get_virtual_energy_params, 2, 1);

	jit_type_t get_eri_params[] = {jit_type_void_ptr, jit_type_int, jit_type_int, jit_type_int, jit_type_int};
	jit_type_t get_eri_signature = jit_type_create_signature(jit_abi_cdecl, jit_type_float64, get_eri_params, 2, 1);

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

			switch(awt->denominators[c].qtypes[d])
			{
				case QTYPE_OCCUPIED:

				if((label==label1)||(label==label2))
				{
					jit_value_t get_occupied_energy_args[2];

					itmp = jit_insn_sub(F, (label==label1)?(value1):(value2), CONST_INT(1));

					get_occupied_energy_args[0]=pointer;
					get_occupied_energy_args[1]=itmp;

					f64tmp = jit_insn_call_native(F, "get_occupied_energy", get_occupied_energy, get_occupied_energy_signature,
								      get_occupied_energy_args, 2, JIT_CALL_NOTHROW);

					// denominator += f64tmp
					denominator = jit_insn_add(F, denominator, f64tmp);
				}
				else
				{
					// denominator += get_occupied_energy(amx->ectx, awt->labels[label].value-1)
					denominator = jit_insn_add(F, denominator, CONST_FLOAT64(get_occupied_energy(amx->ectx, awt->labels[label].value-1)));
				}

				break;

				case QTYPE_VIRTUAL:

				if((label==label1)||(label==label2))
				{
					jit_value_t get_virtual_energy_args[2];

					itmp = jit_insn_sub(F, (label==label1)?(value1):(value2), CONST_INT(1));

					get_virtual_energy_args[0]=pointer;
					get_virtual_energy_args[1]=itmp;

					f64tmp = jit_insn_call_native(F, "get_virtual_energy", get_virtual_energy, get_virtual_energy_signature,
								      get_virtual_energy_args, 2, JIT_CALL_NOTHROW);

					// denominator -= f64tmp
					denominator = jit_insn_sub(F, denominator, f64tmp);
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
		int l1,l2,l3,l4;

		l1=awt->numerators[c].labels[0];
		l2=awt->numerators[c].labels[1];
		l3=awt->numerators[c].labels[2];
		l4=awt->numerators[c].labels[3];

		int i1,i2,i3,i4;

		i1=awt->labels[l1].value-1+((awt->labels[l1].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i2=awt->labels[l2].value-1+((awt->labels[l2].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i3=awt->labels[l3].value-1+((awt->labels[l3].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));
		i4=awt->labels[l4].value-1+((awt->labels[l4].qtype==QTYPE_VIRTUAL)?(amx->ectx->nocc):(0));

		if((l1!=label1)&&(l1!=label2)&&(l2!=label1)&&(l2!=label2)&&
		   (l3!=label1)&&(l3!=label2)&&(l4!=label1)&&(l4!=label2))
		{
			double numerator=get_eri(amx->ectx, i1, i2, i3, i4);

			// numerators *= numerator
			numerators = jit_insn_mul(F, numerators, CONST_FLOAT64(numerator));
		}
		else
		{
			jit_value_t get_eri_args[5];
#if 0
			get_eri_args[0]=pointer;
			get_eri_args[1]=label_to_jit_value(F, l1, label1, value1, label2, value2);
			get_eri_args[2]=label_to_jit_value(F, l2, label1, value1, label2, value2);
			get_eri_args[3]=label_to_jit_value(F, l3, label1, value1, label2, value2);
			get_eri_args[4]=label_to_jit_value(F, l4, label1, value1, label2, value2);
#endif

			f64tmp = jit_insn_call_native(F, "get_eri", get_eri, get_eri_signature,
						      get_eri_args, 5, JIT_CALL_NOTHROW);

			// numerators *= f64mp
			numerators = jit_insn_mul(F, numerators, f64tmp);
		}
	}

	double lindeloef_factor=exp(amx->config->epsilon*awt->h*log(awt->h));
	double extrafactors=awt->unphysical_penalty*lindeloef_factor*pow(awt->inversefactor,-1.0f);

	result = jit_insn_div(F, numerators, denominators);
	result = jit_insn_mul(F, result, CONST_FLOAT64(extrafactors));

	// result = amx->config->bias + result/amatrix_multiplicity(amx)
	result = jit_insn_div(F, result, CONST_FLOAT64(amatrix_multiplicity(amx)));
	result = jit_insn_add(F, result, CONST_FLOAT64(amx->config->bias));

	/*
		Final cleanup
	*/

	jit_insn_return(F, result);

	jit_context_build_end(context);

	return F;
}
