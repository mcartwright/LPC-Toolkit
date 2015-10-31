/*----------------------------------------------------------------------------------
Filename:		mbc.areaScale~.c
Project:		LPC Toolkit
Author:			Mark Cartwright
Created:		7/23/07
Updated:		7/23/07
Description:	Scale area coefficients by an array of weights
-------------------------------------------------------------------------------------*/

#include "ext.h"
#include "z_dsp.h"

#define MAX_ORDER 200

void *areaScale_class;

typedef struct _areaScale
{
    t_pxobject	a_obj;
    t_float*	a_mults;
	int			a_order;
} t_areaScale;

t_int *areaScale_perform(t_int *w);
void areaScale_list(t_areaScale *x, t_symbol *msg, short argc, t_atom *argv);
void areaScale_order(t_areaScale *x, long order);
void areaScale_init(t_areaScale *x);
void areaScale_free(t_areaScale *x);
void areaScale_dsp(t_areaScale *x, t_signal **sp, short *count);
void areaScale_assist(t_areaScale *x, void *b, long m, long a, char *s);
void *areaScale_new(long order);

void main(void)
{
	setup((t_messlist **)&areaScale_class, (method)areaScale_new, (method)areaScale_free, (short)sizeof(t_areaScale), 0L, A_DEFLONG, 0);
	dsp_initclass();
	addmess((method)areaScale_order, "order", A_LONG, 0);
	addmess((method)areaScale_dsp, "dsp", A_CANT, 0);
	addmess((method)areaScale_list, "list", A_GIMME, 0);
	addmess((method)areaScale_assist,"assist",A_CANT,0);
}

t_int *areaScale_perform(t_int *w)
{
    t_float *areaIn = (t_float *)(w[1]);
    t_float *areaIdxIn = (t_float *)(w[2]);
	t_areaScale *x = (t_areaScale *)(w[3]);
	t_float *areaOut = (t_float *)(w[4]);
	t_float *areaIdxOut = (t_float *)(w[5]);
	int n = (int)(w[6]);
	t_float* mults = x->a_mults;
	t_float in, out, mult;
	
	if (x->a_obj.z_disabled)
		goto out;
	
    while (n--) {
		if ((int)(*areaIdxIn) > 0) {
			in = *areaIn;
			mult = mults[(int)(*areaIdxIn)];
			out = in * mult;
			*areaOut  = *areaIn * mults[(int)(*areaIdxIn)];
			*areaIdxOut = *areaIdxIn;
		} else {
			*areaOut = 0.0;
			*areaIdxOut = 0.0;
		}
		
		areaIn++;
		areaIdxIn++;
		areaOut++;
		areaIdxOut++;
	}
	
out:
    return (w+7);
}

void areaScale_list(t_areaScale *x, t_symbol *msg, short argc, t_atom *argv) {

	int i;
	if (argc == x->a_order) {
		int min = (argc > x->a_order) ? x->a_order : argc; //use the minimum of the two to avoid any pointer erors
		t_float* mults = x->a_mults;
		
		for (i = 0; i < min; i++) {
			*mults = atom_getfloatarg(i,argc,argv);
			mults++;
		}
	} else {
		error("mbc.areaScale~: number of items in the input list must equal the specified order");
	}
	
	int min = (argc > x->a_order) ? x->a_order : argc; //use the minimum of the two to avoid any pointer erors
	/*for (i=0; i < min; i++) {
		cpost("\nmults(%d) = %f",i+1,x->a_mults[i]);
	}*/
	
}

void areaScale_order(t_areaScale *x, long order) {
	x->a_order = order;
}

void areaScale_dsp(t_areaScale *x, t_signal **sp, short *count)
{
	dsp_add(areaScale_perform, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

void areaScale_assist(t_areaScale *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		switch (a) {
			case 0: sprintf(s,"(signal) Areas"); break;
			case 1: sprintf(s,"(signal) Coeff Index"); break;
		}
	}
	else {
		switch (a) {
			case 0: sprintf(s,"(signal) Scaled Areas"); break;
			case 1: sprintf(s,"(signal) Coeff Index"); break;
		}
	}
}

void areaScale_init(t_areaScale *x) {
	int i;
	
	x->a_mults = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	
	for (i=0; i < x->a_order; i++) {
		x->a_mults[i] = 1.0;
	}
}

void areaScale_free(t_areaScale *x) {
	dsp_free((t_pxobject *) x);
	freebytes(x->a_mults, MAX_ORDER * sizeof(t_float));
}

void *areaScale_new(long order)
{
    t_areaScale *x = (t_areaScale *)newobject(areaScale_class);
    dsp_setup((t_pxobject *)x,2);
    outlet_new((t_pxobject *)x, "signal");
	outlet_new((t_pxobject *)x, "signal");
	x->a_obj.z_misc |= Z_NO_INPLACE;
	
    x->a_order = order;
	
	areaScale_init(x);
	
    return (x);
}



