/*----------------------------------------------------------------------------------
Filename:		mbc.areaGen~.c
Project:		LPC Toolkit
Author:			Mark Cartwright
Created:		5/30/07
Updated:		5/30/07
Description:	Generation of area coefficients for an allpole vocal tract model from
				an input list.
-------------------------------------------------------------------------------------*/

#include "ext.h"
#include "z_dsp.h"

void *areaGen_class;

typedef struct _areaGen
{
    t_pxobject	a_obj;			//object
	t_float*	a_areas;		//area of area coefficients
	t_float*	a_areaBuff;		//input buffer for area coefficients
    int			a_order;		//filter order
	int			a_outIdx;		//out index
	int			a_buffFull;		//buffFull flag
} t_areaGen;

t_int *areaGen_perform(t_int *w);
void areaGen_dsp(t_areaGen *x, t_signal **sp, short *count);
void areaGen_list(t_areaGen *x, t_symbol *msg, short argc, t_atom *argv);
void areaGen_assist(t_areaGen *x, void *b, long m, long a, char *s);
void areaGen_init(t_areaGen *x);
void areaGen_free(t_areaGen *x);
void *areaGen_new(long order);

void main(void)
{
	setup((t_messlist **)&areaGen_class, (method)areaGen_new, (method)areaGen_free, (short)sizeof(t_areaGen), 0L, A_DEFLONG, 0);
	dsp_initclass();
	addmess((method)areaGen_dsp, "dsp", A_CANT, 0);
	addmess((method)areaGen_list, "list", A_GIMME, 0);
	addmess((method)areaGen_assist,"assist",A_CANT,0);
}

t_int *areaGen_perform(t_int *w)
{
	t_areaGen *x = (t_areaGen *)(w[1]);
	t_float *areaOut = (t_float *)(w[2]);
	t_float *coeffIdxOut = (t_float *)(w[3]);
	int n = (int)(w[4]);	
	int outidx = x->a_outIdx;
	int order = x->a_order;
	int buffFull = x->a_buffFull;
	int i;
	
	while(n--) {
		if (outidx < order) {
			*areaOut = x->a_areas[outidx];
			*coeffIdxOut = (t_float)(outidx+1);
			outidx++;
			areaOut++; 
			coeffIdxOut++;
		} else if(buffFull) {
			for(i=0;i<order;i++) {
				x->a_areas[i] = x->a_areaBuff[i];
			}
			outidx = 0;
			buffFull = 0;
			
			*areaOut = 0.0;
			*coeffIdxOut = 0.0;
			areaOut++;
			coeffIdxOut++;
		} else {
			*areaOut = 0.0;
			*coeffIdxOut = 0.0;
			areaOut++;
			coeffIdxOut++;
		}
		
	}
	
	x->a_buffFull = buffFull;
	x->a_outIdx = outidx;
	
	return (w+5);
}

void areaGen_list(t_areaGen *x, t_symbol *msg, short argc, t_atom *argv) {

	int i;
	if (argc == x->a_order) {
		int min = (argc > x->a_order) ? x->a_order : argc; //use the minimum of the two to avoid any pointer erors
		t_float* areaBuff = x->a_areaBuff;
		
		for (i = 0; i < min; i++) {
			*areaBuff = atom_getfloatarg(i,argc,argv);
			areaBuff++;
		}
		x->a_buffFull = 1;
	} else {
		error("mbc.areaGen~: number of items in the input list must equal to the specified order");
	}
	
}

void areaGen_dsp(t_areaGen *x, t_signal **sp, short *count)
{
	x->a_buffFull = 1;
	dsp_add(areaGen_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

void areaGen_assist(t_areaGen *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		sprintf(s,"(list) Area Coefficients");
	}
	else {
		switch(a) {
			case 0: sprintf(s,"(signal) Area Coefficients"); break;
			case 1: sprintf(s,"(signal) Coefficient Index"); break;
		}
	}
}

void areaGen_init(t_areaGen *x) 
{
	x->a_outIdx = x->a_order + 1;
	x->a_areas = (t_float *) sysmem_newptrclear( x->a_order * sizeof(t_float));
	x->a_areaBuff = (t_float *) sysmem_newptrclear( x->a_order * sizeof(t_float));
}

void areaGen_free(t_areaGen *x)
{
	dsp_free((t_pxobject *) x);
	
	sysmem_freeptr(x->a_areas);
	sysmem_freeptr(x->a_areaBuff);
}

void *areaGen_new(long order)
{
    t_areaGen *x = (t_areaGen *)newobject(areaGen_class);
    dsp_setup((t_pxobject *)x,0); //0 signal objects
    outlet_new((t_pxobject *)x, "signal");
	outlet_new((t_pxobject *)x, "signal");
	
	x->a_order = (int)order;
	areaGen_init(x);
	
	return (x);
}



