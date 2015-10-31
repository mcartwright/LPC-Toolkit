#include "ext.h"
#include "z_dsp.h"

#define MAX_ORDER 200

void *coeffConv_class;

typedef struct _coeffConv
{
    t_pxobject c_obj;
	t_float* c_a;			//filer coefficients
	t_float* c_inBuff;		//filter coeff input buffer
	t_float* c_Ar;			//tube areas
	double **c_A;			//for parcor to coeff conversion
	double *c_K;			//parcor coefficients
	long c_order;
	int c_sCoeffType;		//type of source coefficients: CT_PARCOR,CT_FILTER,CT_AREA
	int c_tCoeffType;		//type of destiation coefficients: CT_PARCOR,CT_FILTER,CT_AREA
	int c_out_idx;			//current outframe buffer index
} t_coeffConv;

enum coeffTypes {
	CT_PARCOR,
	CT_FILTER,
	CT_AREA
};

t_int *coeffConv_perf_f2a(t_int *w);
t_int *coeffConv_perf_f2p(t_int *w);
t_int *coeffConv_perf_p2a(t_int *w);
t_int *coeffConv_perf_p2f(t_int *w);
t_int *coeffConv_perf_a2f(t_int *w);
t_int *coeffConv_perf_a2p(t_int *w);
void coeffConv_order(t_coeffConv *x, int order);
void coeffConv_init(t_coeffConv *x);
void coeffConv_free(t_coeffConv *x);
void coeffConv_clear(t_coeffConv *x);
void coeffConv_dsp(t_coeffConv *x, t_signal **sp, short *count);
void coeffConv_assist(t_coeffConv *x, void *b, long m, long a, char *s);
void *coeffConv_new(t_symbol *s, int argc, t_atom *argv);

void main(void)
{
	setup((t_messlist **)&coeffConv_class, (method)coeffConv_new, (method)coeffConv_free, (short)sizeof(t_coeffConv), 0L, A_GIMME, 0); 
	dsp_initclass();
	addmess((method)coeffConv_dsp, "dsp", A_CANT, 0);
	addmess((method)coeffConv_order,"order",A_DEFLONG,0);
	addmess((method)coeffConv_assist,"assist",A_CANT,0);
}

t_int *coeffConv_perf_f2a(t_int *w) {
}

t_int *coeffConv_perf_f2p(t_int *w) {
}

t_int *coeffConv_perf_p2a(t_int *w)
{
	t_float *coeffIn = (t_float *)(w[1]);
	t_float *coeffIdxIn = (t_float *)(w[2]);
	t_coeffConv *x = (t_coeffConv *)(w[3]);
	t_float *coeffOut = (t_float *)(w[4]);
	t_float *coeffIdxOut = (t_float *)(w[5]);
	int n = (int)(w[6]);
	int order = x->c_order;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	
	while(n--) {
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[6]);
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if (x->c_out_idx < order) {
			coeffOut[out_idx] = x->c_Ar[x->c_out_idx];
			coeffIdxOut[out_idx] = x->c_out_idx + 1;
			x->c_out_idx++;
		} else {
			coeffOut[out_idx] = 0.0;
			coeffIdxOut[out_idx] = 0.0;
		}
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->c_K[i+1] = x->c_inBuff[i];
				}
								
				x->c_Ar[0] = 1.0;
				for (i = 1; i < order; i++) { // droping a coefficient or not?  this is the way that vtoa does it...
					x->c_Ar[i] = x->c_Ar[i-1]*(1 - x->c_K[i])/(1 + x->c_K[i]);
				}
				
				x->c_out_idx = 0;
				n++; //to make up for last decrement
				while (n--) {
					if (x->c_out_idx < order) {
						coeffOut[out_idx] = x->c_Ar[x->c_out_idx];
						coeffIdxOut[out_idx] = x->c_out_idx + 1;
						x->c_out_idx++;
					} else {
						coeffOut[out_idx] = 0.0;
						coeffIdxOut[out_idx] = 0.0;
					}
					if ((int)(coeffIdxIn[in_idx]) > 0) {
						x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
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
	
	return (w+7);
}


t_int *coeffConv_perf_p2f(t_int *w)
{

	t_float *coeffIn = (t_float *)(w[1]);
	t_float *coeffIdxIn = (t_float *)(w[2]);
	t_coeffConv *x = (t_coeffConv *)(w[3]);
	t_float *coeffOut = (t_float *)(w[4]);
	t_float *coeffIdxOut = (t_float *)(w[5]);
	int n = (int)(w[6]);
	int order = x->c_order;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	
	while(n--) {
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[6]);
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if (x->c_out_idx < order) {
			coeffOut[out_idx] = x->c_a[x->c_out_idx];
			coeffIdxOut[out_idx] = x->c_out_idx + 1;
			x->c_out_idx++;
		} else {
			coeffOut[out_idx] = 0; 
			coeffIdxOut[out_idx] = 0;
		}
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->c_K[i+1] = x->c_inBuff[i];
				}
				
				//solve for filter coefficients
				for (i = 1; i <= order; i++) {
					x->c_A[i][i] = x->c_K[i];
					i1 = i - 1;
					if (i1 >= 1) {
						for (j = 1; j <=i1; j++) {
							ji = i - j;
							x->c_A[j][i] = x->c_A[j][i1] - x->c_K[i] * x->c_A[ji][i1];
						}				
					}
				}
				
				for (j = 1; j <=order; j++) {
					x->c_a[j-1] = x->c_A[j][order];
				}
				x->c_out_idx = 0;
				n++; //to make up for last decrement
				while (n--) {
					if (x->c_out_idx < order) {
						coeffOut[out_idx] = x->c_a[x->c_out_idx];
						coeffIdxOut[out_idx] = x->c_out_idx + 1;
						x->c_out_idx++;
					} else {
						coeffOut[out_idx] = 0.0;
						coeffIdxOut[out_idx] = 0.0;
					}
					if ((int)(coeffIdxIn[in_idx]) > 0) {
						x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
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
	
	return (w+7);
}

t_int *coeffConv_perf_a2f(t_int *w)
{
	t_float *coeffIn = (t_float *)(w[1]);
	t_float *coeffIdxIn = (t_float *)(w[2]);
	t_coeffConv *x = (t_coeffConv *)(w[3]);
	t_float *coeffOut = (t_float *)(w[4]);
	t_float *coeffIdxOut = (t_float *)(w[5]);
	int n = (int)(w[6]);
	int order = x->c_order;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	
	while(n--) {
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[6]);
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if (x->c_out_idx < order) {
			coeffOut[out_idx] = x->c_a[x->c_out_idx];
			coeffIdxOut[out_idx] = x->c_out_idx + 1;
			x->c_out_idx++;
		} else {
			coeffOut[out_idx] = 0.0;
			coeffIdxOut[out_idx] = 0.0;
		}
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->c_Ar[i+1] = x->c_inBuff[i];
				}
				
				//solve for filter coefficients from area coefficients
				for (i = 1; i <= order; i++) {
					//convert from area coefficients to reflection coefficients (remember reflection coefficients are -PARCOR);
					if (i != order)
						x->c_K[i] = -(x->c_Ar[i+1] - x->c_Ar[i])/(x->c_Ar[i+1] + x->c_Ar[i]);
					else
						x->c_K[i] = -0.7; //no losses at lips NOTE: need to make this changeable
					
					x->c_A[i][i] = x->c_K[i];
								
					i1 = i - 1;
					if (i1 >= 1) {
						for (j = 1; j <=i1; j++) {
							ji = i - j;
							x->c_A[j][i] = x->c_A[j][i1] - x->c_K[i] * x->c_A[ji][i1];
						}				
					}
				}
				
				for (j = 1; j <=order; j++) {
					x->c_a[j-1] = x->c_A[j][order];
				}
				x->c_out_idx = 0;
				n++; //to make up for last decrement
				while (n--) {
					if (x->c_out_idx < order) {
						coeffOut[out_idx] = x->c_a[x->c_out_idx];
						coeffIdxOut[out_idx] = x->c_out_idx + 1;
						x->c_out_idx++;
					} else {
						coeffOut[out_idx] = 0.0;
						coeffIdxOut[out_idx] = 0.0;
					}
					if ((int)(coeffIdxIn[in_idx]) > 0) {
						x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
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
	
	return (w+7);
}

t_int *coeffConv_perf_a2p(t_int *w)
{
	t_float *coeffIn = (t_float *)(w[1]);
	t_float *coeffIdxIn = (t_float *)(w[2]);
	t_coeffConv *x = (t_coeffConv *)(w[3]);
	t_float *coeffOut = (t_float *)(w[4]);
	t_float *coeffIdxOut = (t_float *)(w[5]);
	int n = (int)(w[6]);
	int order = x->c_order;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	
	while(n--) {
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[6]);
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if (x->c_out_idx < order) {
			coeffOut[out_idx] = x->c_K[x->c_out_idx+1];
			coeffIdxOut[out_idx] = x->c_out_idx + 1;
			x->c_out_idx++;
		} else {
			coeffOut[out_idx] = 0.0;
			coeffIdxOut[out_idx] = 0.0;
		}
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->c_Ar[i+1] = x->c_inBuff[i];
				}
				
				//solve for filter coefficients from area coefficients
				for (i = 1; i <= order; i++) {
					//convert from area coefficients to reflection coefficients (remember reflection coefficients are -PARCOR);
					if (i != order)
						x->c_K[i] = -(x->c_Ar[i+1] - x->c_Ar[i])/(x->c_Ar[i+1] + x->c_Ar[i]);
					else
						x->c_K[i] = -0.7; //no losses at lips NOTE: need to make this changeable
				}
			
				x->c_out_idx = 0;
				n++; //to make up for last decrement
				while (n--) {
					if (x->c_out_idx < order) {
						coeffOut[out_idx] = x->c_K[x->c_out_idx + 1];
						coeffIdxOut[out_idx] = x->c_out_idx + 1;
						x->c_out_idx++;
					} else {
						coeffOut[out_idx] = 0.0;
						coeffIdxOut[out_idx] = 0.0;
					}
					if ((int)(coeffIdxIn[in_idx]) > 0) {
						x->c_inBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
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
	
	return (w+7);
}


void coeffConv_init(t_coeffConv *x) {
	int i;
	
	x->c_a = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->c_inBuff = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->c_Ar = (t_float *) getbytes( (MAX_ORDER + 1) * sizeof(t_float));
	x->c_K = (double *) getbytes( (MAX_ORDER + 1) * sizeof(double));
	x->c_A = (double **) getbytes( (MAX_ORDER + 1) * sizeof(double*));
	for(i=0; i<MAX_ORDER; i++) {
		x->c_A[i] = (double *)getbytes( (MAX_ORDER + 1) * sizeof(double));
	}
	coeffConv_clear(x);
}

void coeffConv_free(t_coeffConv *x) {
	int i;
	
	dsp_free((t_pxobject *) x);
	freebytes(x->c_a, MAX_ORDER * sizeof(t_float));
	freebytes(x->c_inBuff, MAX_ORDER * sizeof(t_float));
	freebytes(x->c_Ar, (MAX_ORDER + 1) * sizeof(t_float));
	freebytes(x->c_K, (MAX_ORDER + 1) * sizeof(double));
	for(i=0; i<MAX_ORDER; i++) {
		freebytes(x->c_A[i], (MAX_ORDER + 1) * sizeof(double));
	}
	freebytes(x->c_A, MAX_ORDER * sizeof(double*));

}

void coeffConv_clear(t_coeffConv *x) {
	int i, j;
	for(i = 0; i < MAX_ORDER; i++) {
		x->c_a[i] = 0.0;
		x->c_inBuff[i] = 0.0;
		for (j = 0; j < (MAX_ORDER+1); j++) {
			x->c_A[i][j] = 0.0;
		}
	}
	
	for(i = 0; i < (MAX_ORDER+1); i++) {
		x->c_Ar[i] = 0.0;
		x->c_K[i] = 0.0;
	}
	x->c_out_idx = 0;
}

void coeffConv_dsp(t_coeffConv *x, t_signal **sp, short *count)
{
	//NOTE: need to specify parcor, filter, or area perform string !!!!!!
	switch(x->c_sCoeffType) {
		case CT_FILTER:
			switch (x->c_tCoeffType) {
				case CT_PARCOR:
					dsp_add(coeffConv_perf_f2p, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
					break;
					
				case CT_AREA:
					dsp_add(coeffConv_perf_f2a, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
					break;
				
				default:
					error("mbc.coeffConv~: no target coefficient type selected");
					break;
			}
			break;

		case CT_PARCOR:
			switch (x->c_tCoeffType) {
				case CT_FILTER:
					dsp_add(coeffConv_perf_p2f, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
					break;
				
				case CT_AREA:
					dsp_add(coeffConv_perf_p2a, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
					break;
				
				default:
					error("mbc.coeffConv~: no target coefficient type selected");
					break;
			}
			break;
		
		case CT_AREA:
			switch (x->c_tCoeffType) {
				case CT_FILTER:
					dsp_add(coeffConv_perf_a2f, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
					break;
				
				case CT_PARCOR:
					dsp_add(coeffConv_perf_a2p, 6, sp[0]->s_vec, sp[1]->s_vec, x, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
					break;
				
				default:
					error("mbc.coeffConv~: no target coefficient type selected");
					break;
			}
			break;
		
		default:
			error("mbc.coeffConv~: no coefficient type selected");
			break;
	}	
	coeffConv_clear(x);
}

void coeffConv_order(t_coeffConv *x, int order) {
	if (order < 1) {
		error("mbc.coeffConv~: mbc.coeffConv~ needs a positive integer value for order");
		order = 1;
	} else if (order > 199) {
		error("mbc.coeffConv~: max order is 199");
		order = 199;
	}
	
	int i = x->c_order;
	x->c_order = order;
	for ( i = i+1; i < order; i++) {
		x->c_a[i] = 0.0;
		x->c_inBuff[i] = 0.0;
		x->c_Ar[i] = 0.0;
	}
}

void coeffConv_assist(t_coeffConv *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		switch (a) {
			case 0: sprintf(s,"(signal) In Coefficients"); break;
			case 1: sprintf(s,"(signal) In Coeff Index "); break;
		}
	}
	else {
		switch (a) {
			case 0: sprintf(s,"(signal) Out Coefficients"); break;
			case 1: sprintf(s,"(signal) Out Coeff Index "); break;
		}
	}
}

void *coeffConv_new(t_symbol *s, int argc, t_atom *argv)
{
    t_coeffConv *x = (t_coeffConv *)newobject(coeffConv_class);
    dsp_setup((t_pxobject *)x,2);
    outlet_new((t_pxobject *)x, "signal");
	outlet_new((t_pxobject *)x, "signal");
	x->c_obj.z_misc |= Z_NO_INPLACE;
	
	//get arguments out of gimme list
	int order = atom_getintarg(0,argc,argv);
	t_symbol* sCoeffType = atom_getsymarg(1,argc,argv);
	t_symbol* tCoeffType = atom_getsymarg(2,argc,argv);
	
	//order bounds
	if (order < 1) {
		error("mbc.coeffConv~: mbc.coeffConv~ needs a positive integer value for order");
		order = 1;
	} else if (order > 199) {
		error("mbc.coeffConv~: max order is 199");
		order = 199;
	}
	
	//select in coeff type
	if (sCoeffType == gensym("parcor")) {
		x->c_sCoeffType = CT_PARCOR;
	} else if (sCoeffType == gensym("filter")) {
		x->c_sCoeffType = CT_FILTER;
	} else if (sCoeffType == gensym("area")) {
		x->c_sCoeffType = CT_AREA;
	} else {
		post("mbc.coeffConv~: no source coefficient type selected.  using default 'parcor' type");
		x->c_sCoeffType = CT_PARCOR;
	}
	
	//select out coeff type
	if (tCoeffType == gensym("parcor")) {
		x->c_tCoeffType = CT_PARCOR;
	} else if (tCoeffType == gensym("filter")) {
		x->c_tCoeffType = CT_FILTER;
	} else if (tCoeffType == gensym("area")) {
		x->c_tCoeffType = CT_AREA;
	} else {
		post("mbc.coeffConv~: no target coefficient type selected.  using default 'area' type");
		x->c_tCoeffType = CT_PARCOR;
	}
	
	//assign locals to globals and init
	x->c_order = order;
	coeffConv_init(x);
	
    return (x);
}



