/*----------------------------------------------------------------------------------
Filename:		mbc.areaSwap~.c
Project:		LPC Toolkit
Author:			Mark Cartwright
Created:		7/23/07
Updated:		7/23/07
Description:	Swap area coefficients by an array of weights
-------------------------------------------------------------------------------------*/

#include "ext.h"
#include "z_dsp.h"

#define MAX_ORDER 200

void *areaSwap_class;

typedef struct _areaSwap
{
    t_pxobject	a_obj;
    int*		a_pos;
	t_float*	a_areaInBuff;
	t_float*	a_areaOutBuff;
	int			a_order;
	int			a_out_idx;
} t_areaSwap;

t_int *areaSwap_perform(t_int *w);
void areaSwap_list(t_areaSwap *x, t_symbol *msg, short argc, t_atom *argv);
void areaSwap_order(t_areaSwap *x, long order);
void areaSwap_init(t_areaSwap *x);
void areaSwap_free(t_areaSwap *x);
void areaSwap_dsp(t_areaSwap *x, t_signal **sp, short *count);
void areaSwap_assist(t_areaSwap *x, void *b, long m, long a, char *s);
void *areaSwap_new(long order);

void main(void)
{
	setup((t_messlist **)&areaSwap_class, (method)areaSwap_new, (method)areaSwap_free, (short)sizeof(t_areaSwap), 0L, A_DEFLONG, 0);
	dsp_initclass();
	addmess((method)areaSwap_order, "order", A_LONG, 0);
	addmess((method)areaSwap_dsp, "dsp", A_CANT, 0);
	addmess((method)areaSwap_list, "list", A_GIMME, 0);
	addmess((method)areaSwap_assist,"assist",A_CANT,0);
}

t_int *areaSwap_perform(t_int *w)
{
    t_float *areaIn = (t_float *)(w[1]);
    t_float *areaIdxIn = (t_float *)(w[2]);
	t_areaSwap *x = (t_areaSwap *)(w[3]);
	t_float *areaOut = (t_float *)(w[4]);
	t_float *areaIdxOut = (t_float *)(w[5]);
	int n = (int)(w[6]);
	
	int order = x->a_order;
	int out_idx = 0; //x->a_out_idx;
	int in_idx = 0;
	int i;
	int* pos = x->a_pos;
	
	if (x->a_obj.z_disabled)
		goto out;
	
    while (n--) {
		if (x->a_out_idx < order) {
			areaOut[out_idx] = x->a_areaOutBuff[x->a_out_idx];
			areaIdxOut[out_idx] = x->a_out_idx + 1;
			cpost("\nareaOut[%d] = %f",out_idx,areaOut[out_idx]);
			x->a_out_idx++;
		} else {
			areaOut[out_idx] = 0.0;
			areaIdxOut[out_idx] = 0.0;
		}
		if ((int)(areaIdxIn[in_idx]) > 0) {
			x->a_areaInBuff[(int)(areaIdxIn[in_idx])-1] = areaIn[in_idx];
			if ((int)(areaIdxIn[in_idx]) == order) {
				
				for (i=0; i < order; i++) {
					x->a_areaOutBuff[pos[i] - 1] = x->a_areaInBuff[i];
					//cpost("\npos[%d] = %d",i,pos[i]);
				}
												
				x->a_out_idx = 0;
				n++; //to make up for last decrement
				while (n--) {
					if (x->a_out_idx < order) {
						areaOut[out_idx] = x->a_areaOutBuff[x->a_out_idx];
						areaIdxOut[out_idx] = x->a_out_idx + 1;
						cpost("\nareaOut[%d] = %f",out_idx,areaOut[out_idx]);
						x->a_out_idx++;
					} else {
						areaOut[out_idx] = 0.0;
						areaIdxOut[out_idx] = 0.0;
					}
					if ((int)(areaIdxIn[in_idx]) > 0) {
						x->a_areaInBuff[(int)(areaIdxIn[in_idx])-1] = areaIn[in_idx];
					}
					in_idx++;
					out_idx++;
				}
				break;
			}
		}
		in_idx++;
		out_idx++;
	}
	
out:
    return (w+7);
}

void areaSwap_list(t_areaSwap *x, t_symbol *msg, short argc, t_atom *argv) {

	int i, pos_in;
	if (argc == x->a_order) {
		int min = (argc > x->a_order) ? x->a_order : argc; //use the minimum of the two to avoid any pointer erors
		int* pos = x->a_pos;
		int order = x->a_order;
		
		for (i = 0; i < min; i++) {
			pos_in = atom_getintarg(i,argc,argv);
			if (pos_in > 0 && pos_in <= order) {
				*pos = atom_getintarg(i,argc,argv);
			} else {
				error("mbc.areaSwap~: position input is out of bounds for the given order.");
				*pos = i+1;
			}
			pos++;
		}
	} else {
		error("mbc.areaSwap~: number of items in the input list must equal the specified order");
	}
	
	int min = (argc > x->a_order) ? x->a_order : argc; //use the minimum of the two to avoid any pointer erors
	for (i=0; i < min; i++) {
		cpost("\npos(%d) = %d",i,x->a_pos[i]);
	}
	
}

void areaSwap_order(t_areaSwap *x, long order) {
	x->a_order = order;
}

void areaSwap_dsp(t_areaSwap *x, t_signal **sp, short *count)
{
	dsp_add(areaSwap_perform, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}

void areaSwap_assist(t_areaSwap *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		switch (a) {
			case 0: sprintf(s,"(signal) Areas"); break;
			case 1: sprintf(s,"(signal) Coeff Index"); break;
		}
	}
	else {
		switch (a) {
			case 0: sprintf(s,"(signal) Swapd Areas"); break;
			case 1: sprintf(s,"(signal) Coeff Index"); break;
		}
	}
}

void areaSwap_init(t_areaSwap *x) {
	int i;
	
	x->a_pos = (int *) getbytes( MAX_ORDER * sizeof(int));
	x->a_areaInBuff = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->a_areaOutBuff = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	
	for (i=0; i < x->a_order; i++) {
		x->a_pos[i] = i+1;
		x->a_areaInBuff[i] = 0.0;
		x->a_areaOutBuff[i] = 0.0;
	}
}

void areaSwap_free(t_areaSwap *x) {
	dsp_free((t_pxobject *) x);
	freebytes(x->a_pos, MAX_ORDER * sizeof(int));
	freebytes(x->a_areaInBuff, MAX_ORDER * sizeof(t_float));
	freebytes(x->a_areaOutBuff, MAX_ORDER * sizeof(t_float));
}

void *areaSwap_new(long order)
{
    t_areaSwap *x = (t_areaSwap *)newobject(areaSwap_class);
    dsp_setup((t_pxobject *)x,2);
    outlet_new((t_pxobject *)x, "signal");
	outlet_new((t_pxobject *)x, "signal");
	x->a_obj.z_misc |= Z_NO_INPLACE;
	
    x->a_order = order;
	
	areaSwap_init(x);
	
    return (x);
}



