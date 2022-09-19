/**
	@file
	mbc.allpole~ - an MSP object shell
	mark cartwright - mcartwright@gmail.com	

	@ingroup	lpcToolkit	
*/

#include "ext.h"							// standard Max include, always required (except in Jitter)
#include "ext_obex.h"						// required for new style objects
#include "z_dsp.h"							// required for MSP objects
#include <Accelerate/Accelerate.h>

#define MAX_ORDER 200
#define DEFAULT_INTERP 1.0
#define DEFAULT_DEEMPH 0
#define t_floatarg double
#define t_mbcfloat double
#define MBC_VDSP

////////////////////////// object struct
typedef struct _allpole 
{
	t_pxobject					ob;				// the object itself (t_pxobject in MSP)
	t_mbcfloat*					a_a;			//filer coefficients
	t_mbcfloat* 				a_aBuff;		//filter coeff input buffer
	t_mbcfloat* 				a_y;			//filter memory
	t_mbcfloat* 				a_tempVec;		//temporary buff for vector math
	double 						a_a1;			//deemph coeff
	double 						a_y1;			//deemph filter memory
	t_mbcfloat* 				a_Ar;			//tube areas
	t_mbcfloat**				a_A;			//for parcor to coeff conversion
	t_mbcfloat*					a_K;			//parcor coefficients
	t_mbcfloat*					a_interpCoeff;	//interpolated coefficients
	t_mbcfloat*					a_interpInc;	//interpolation increment
	t_mbcfloat 					a_G;
    float 						a_interp;		//interpolation time in ms
	long 						a_order;
	long 						a_deemph;
	int 						a_coeffType;	//type of coeffecients: CT_PARCOR,CT_FILTER,CT_AREA
	double 						a_fs;			//sampling rate
	long 						a_vsize;		//vector size
} t_allpole;

enum coeffTypes 
{
	CT_PARCOR,
	CT_FILTER,
	CT_AREA
};

///////////////////////// function prototypes
//// standard set
void *allpole_new(t_symbol *s, long argc, t_atom *argv);
void allpole_free(t_allpole *x);
void allpole_assist(t_allpole *x, void *b, long m, long a, char *s);

void allpole_float(t_allpole *x, double f);

void allpole_dsp64(t_allpole *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
//t_int *allpole_perf_filter(t_int *w);
//t_int *allpole_perf_area(t_int *w);
void allpole_perf_parcor(t_allpole *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void allpole_perf_parcorI(t_allpole *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
//t_int *allpole_perf_filterI(t_int *w);
//t_int *allpole_perf_areaI(t_int *w);
void allpole_interp(t_allpole *x, t_floatarg interp);
void allpole_order(t_allpole *x, int order);
void allpole_deemph(t_allpole *x, int deemph); 
void allpole_init(t_allpole *x);
void allpole_free(t_allpole *x);
void allpole_free_arrays(t_allpole *x);
void allpole_clear(t_allpole *x);
static inline void allpole_highOrdFilter(t_allpole* x, int N, int order, double* in, double* out);
static inline void allpole_solveForFiltCoefs(t_allpole* x, int order);
static inline void allpole_deemphFilter(t_allpole *x, int N, double* vec);

//////////////////////// global class pointer variable
void *allpole_class;


void ext_main(void *r)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	// OLD METHOD
	// setup((t_messlist **)&allpole_class, (method)allpole_new, (method)dsp_free, (short)sizeof(t_allpole), 0L, A_GIMME, 0);
	// addfloat((method)allpole_float);
	// you need this
	// addmess((method)allpole_dsp,				"dsp",			A_CANT, 0);
    // addmess((method)allpole_assist,			"assist",		A_CANT, 0);  
	// you need this
    // dsp_initclass();
	
	// NEW METHOD
	t_class *c;
	
	c = class_new("mbc.allpole~", (method)allpole_new, (method)allpole_free, (long)sizeof(t_allpole), 0L, A_GIMME, 0);

	class_addmethod(c, (method)allpole_dsp64, "dsp64", A_CANT, 0);
	class_addmethod(c, (method)allpole_interp,"interp",A_DEFFLOAT,0);
	class_addmethod(c, (method)allpole_order,"order",A_DEFLONG,0);
	class_addmethod(c, (method)allpole_clear,"clear",0);
	class_addmethod(c, (method)allpole_deemph,"deemph",A_DEFLONG,0);
	class_addmethod(c, (method)allpole_assist,"assist",A_CANT,0);
	
	class_dspinit(c);				// new style object version of dsp_initclass();
	class_register(CLASS_BOX, c);	// register class as a box class
	allpole_class = c;
}

//perform for FILTER COEFFICIENTS, INTERPOLATION OFF
/*t_int *allpole_perf_filter(t_int *w)
{
	t_float *in = (t_float *)(w[1]);
	t_float *coeffIn = (t_float *)(w[2]);
	t_float *coeffIdxIn = (t_float *)(w[3]);
	t_float *G_n = (t_float *)(w[4]);
	t_allpole *x = (t_allpole *)(w[5]);
	t_float *out = (t_float *)(w[6]);
	int n = (int)(w[7]);
	int order = x->a_order;
	int i, in_idx = 0, out_idx = 0;
	float val = 0.0;
	float sum = 0.0;
	
	while(n--) {
		if(IS_NAN_FLOAT(G_n[in_idx])) G_n[in_idx] = 0.0;
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[7]);
	
	//deemphasis?
	if (x->a_deemph) {
		while (n--) {
			in[in_idx] = in[in_idx] + x->a_a1*x->a_y1;
			x->a_y1 = in[in_idx];
		
			in_idx++;
		}
		in_idx = 0;
		n = (int)(w[7]);
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->a_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->a_a[i] = x->a_aBuff[i];
				}
				x->a_G = (float)(G_n[in_idx]); 				
			}
		}
		in_idx++;
	}
	
	n = (int)(w[7]);
	
	while (n--) { 
		vDSP_vmul(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sve(x->a_tempVec,1,&sum,order);
		val =  x->a_G * in[out_idx] + sum;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		out[out_idx] = val;
		out_idx++;
	}
	
	return (w+8);
}*/

//perform for PARCOR COEFFICIENTS, INTERPOLATION OFF
void allpole_perf_parcor(t_allpole *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	double *in = ins[0];
	double *coeffIn = ins[1];
	double *coeffIdxIn = ins[2];
	double *G_n = ins[3];
    
	double *out = outs[0];
    long N = sampleframes;
	int order = x->a_order;
	int i, n;
	
	for (n=0; n < N; n++)
	{
		if(IS_NAN_DOUBLE(G_n[n])) G_n[n] = 0.0;
		if(IS_NAN_DOUBLE(coeffIn[n])) coeffIn[n] = 0.0;
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	for (n=0; n < N; n++) 
	{
		if ((int)(coeffIdxIn[n]) > 0) 
		{
			x->a_aBuff[(int)(coeffIdxIn[n])-1] = (t_mbcfloat)coeffIn[n];
			if ((int)(coeffIdxIn[n]) == order) 
			{
				for(i = 0; i < order; i++) 
				{
					x->a_K[i+1] = x->a_aBuff[i];
				}
				x->a_G = (t_mbcfloat)(G_n[n]);
				
				allpole_solveForFiltCoefs(x, order);
			}
		}
	}
	
	allpole_highOrdFilter(x, N, order, in, out);
	
	//deemphasis?
	if (x->a_deemph) 
	{
		allpole_deemphFilter(x, N, out);
	}
}

void allpole_perf_parcorI(t_allpole *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	double *in = ins[0];
	double *coeffIn = ins[1];
	double *coeffIdxIn = ins[2];
	double *G_n = ins[3];
	
	double *out = outs[0];
	int N = sampleframes;
	int order = x->a_order;
	int i, n;
	double interpDiv; //interpolation denominator
	
	for (n=0; n < N; n++)
	{
		if(IS_NAN_DOUBLE(G_n[n])) G_n[n] = 0.0;
		if(IS_NAN_DOUBLE(coeffIn[n])) coeffIn[n] = 0.0;
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	for (n=0; n < N; n++) 
	{
		if ((int)(coeffIdxIn[n]) > 0) 
		{
			x->a_aBuff[(int)(coeffIdxIn[n])-1] = (t_mbcfloat)coeffIn[n];
			if ((int)(coeffIdxIn[n]) == order) 
			{
				interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
				for(i = 0; i < order; i++) 
				{
					x->a_K[i+1] = x->a_aBuff[i];
					//determine interpolation parameters
					x->a_interpInc[i+1] = (x->a_K[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
				}
				
				x->a_G = (t_mbcfloat)(G_n[n]);
			}
		}
	}
	
	//interpolate coefficients
	for (i=1; i <= order; i++) 
	{
		if (fabs(x->a_interpCoeff[i] - x->a_K[i]) > fabs(x->a_interpInc[i])) 
		{
			x->a_interpCoeff[i] += x->a_interpInc[i];
		} 
		else 
		{
			x->a_interpCoeff[i] = x->a_K[i];
		}
	}
	
	allpole_solveForFiltCoefs(x, order);

	allpole_highOrdFilter(x, N, order, in, out);
	
	//deemphasis?
	if (x->a_deemph) 
	{
		allpole_deemphFilter(x, N, out);
	}
}

/*
t_int *allpole_perf_parcor(t_int *w)
{
	t_float *in = (t_float *)(w[1]);
	t_float *coeffIn = (t_float *)(w[2]);
	t_float *coeffIdxIn = (t_float *)(w[3]);
	t_float *G_n = (t_float *)(w[4]);
	t_allpole *x = (t_allpole *)(w[5]);
	t_float *out = (t_float *)(w[6]);
	int N = (int)(w[7]);
	int order = x->a_order;
	int in_idx = 0;
	int out_idx = 0;
	float val = 0.0;
	float sum = 0.0;
	int n, i, j, i1, ji;
	
	for (n=N; n!= 0; n--) 
	{
		if(IS_NAN_FLOAT(G_n[in_idx])) G_n[in_idx] = 0.0;
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	// reset in_idx
	in_idx = 0;
	
	//deemphasis?
	if (x->a_deemph) 
	{
		for (n=N; n!=0; n--) 
		{
			in[in_idx] = in[in_idx] + x->a_a1*x->a_y1;
			x->a_y1 = in[in_idx];
		
			in_idx++;
		}
		
		// reset in_idx
		in_idx = 0;
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	for (n=N; n!=0; n--) 
	{
		if ((int)(coeffIdxIn[in_idx]) > 0) 
		{
			x->a_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) 
			{
				for(i = 0; i < order; i++) 
				{
					x->a_K[i+1] = x->a_aBuff[i];
				}
				x->a_G = (float)(G_n[in_idx]);
				
				//solve for filter coefficients
				for (i = 1; i <= order; i++) 
				{
					x->a_A[i][i] = x->a_K[i];
					i1 = i - 1;
					if (i1 >= 1) 
					{
						for (j = 1; j <=i1; j++) 
						{
							ji = i - j;
							x->a_A[j][i] = x->a_A[j][i1] - x->a_K[i] * x->a_A[ji][i1];
						}				
					}
				}
				
				for (j = 1; j <=order; j++) 
				{
					x->a_a[j-1] = x->a_A[j][order];
				}
			}
		}
		
		in_idx++;
	}
	
	for (n=N; n!=0; n--) 
	{ 
		vDSP_vmul(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sve(x->a_tempVec,1,&sum,order);
		val =  x->a_G * in[out_idx] + sum;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		out[out_idx] = (float)val;
		out_idx++;
	}
	
	return (w+8);
}*/

/*t_int *allpole_perf_area(t_int *w)
{
	t_float *in = (t_float *)(w[1]);
	t_float *coeffIn = (t_float *)(w[2]);
	t_float *coeffIdxIn = (t_float *)(w[3]);
	t_float *G_n = (t_float *)(w[4]);
	t_allpole *x = (t_allpole *)(w[5]);
	t_float *out = (t_float *)(w[6]);
	int n = (int)(w[7]);
	int order = x->a_order;
	//int orderp1 = order + 1;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	float val = 0.0;
	float sum = 0.0;
	float r;
	
	while(n--) {
		if(IS_NAN_FLOAT(G_n[in_idx])) G_n[in_idx] = 0.0;
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[7]);
	
	//deemphasis?
	if (x->a_deemph) {
		while (n--) {
			in[in_idx] = in[in_idx] + x->a_a1*x->a_y1;
			x->a_y1 = in[in_idx];
		
			in_idx++;
		}
		in_idx = 0;
		n = (int)(w[7]);
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->a_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->a_Ar[i+1] = x->a_aBuff[i];
					//cpost("\nAr(%d) = %f;",i+1,x->a_Ar[i+1]);
				}
				//x->a_G = (float)(G_n[in_idx]);
				x->a_G = 1.0; //using G from atov right now instead of input signal (recursively calculated below)
				
				//solve for filter coefficients from area coefficients
				for (i = 1; i <= order; i++) {
					//convert from area coefficients to reflection coefficients (remember reflection coefficients are -PARCOR);
					if (i != order)
						x->a_K[i] = -(x->a_Ar[i+1] - x->a_Ar[i])/(x->a_Ar[i+1] + x->a_Ar[i]);
					else
						x->a_K[i] = -0.7; //no losses at lips NOTE: need to make this changeable
					
					x->a_A[i][i] = x->a_K[i];
					x->a_G = x->a_G * (1.0 + -x->a_K[i]);
					
					i1 = i - 1;
					if (i1 >= 1) {
						for (j = 1; j <=i1; j++) {
							ji = i - j;
							x->a_A[j][i] = x->a_A[j][i1] - x->a_K[i] * x->a_A[ji][i1];
						}				
					}
				}
				
				for (j = 1; j <=order; j++) {
					x->a_a[j-1] = x->a_A[j][order];
				}
			}
		}
		in_idx++;
	}
	
	n = (int)(w[7]);
	
	while (n--) { 
		vDSP_vmul(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sve(x->a_tempVec,1,&sum,order);
		val =  x->a_G * in[out_idx] + sum;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		out[out_idx] = (float)val;
		out_idx++;
	}
	
	return (w+8);
}*/

/*t_int *allpole_perf_filterI(t_int *w) {
}*/

/*t_int *allpole_perf_parcorI(t_int *w)
{
	t_float *in = (t_float *)(w[1]);
	t_float *coeffIn = (t_float *)(w[2]);
	t_float *coeffIdxIn = (t_float *)(w[3]);
	t_float *G_n = (t_float *)(w[4]);
	t_allpole *x = (t_allpole *)(w[5]);
	t_float *out = (t_float *)(w[6]);
	int n = (int)(w[7]);
	int order = x->a_order;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	float val = 0.0;
	float sum = 0.0;
	float interpDiv; //interpolation denominator
	
	while(n--) {
		if(IS_NAN_FLOAT(G_n[in_idx])) G_n[in_idx] = 0.0;
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[7]);
	
	//deemphasis?
	if (x->a_deemph) {
		while (n--) {
			in[in_idx] = in[in_idx] + x->a_a1*x->a_y1;
			x->a_y1 = in[in_idx];
		
			in_idx++;
		}
		in_idx = 0;
		n = (int)(w[7]);
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->a_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
				for(i = 0; i < order; i++) {
					x->a_K[i+1] = x->a_aBuff[i];
					//determine interpolation parameters
					x->a_interpInc[i+1] = (x->a_K[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
				}
				x->a_G = (float)(G_n[in_idx]);
			}
		}
		in_idx++;
	}
	
	//interpolate coefficients
	for (i=1; i <= order; i++) {
		if (fabs(x->a_interpCoeff[i] - x->a_K[i]) > fabs(x->a_interpInc[i])) {
			x->a_interpCoeff[i] += x->a_interpInc[i];
		} else {
			x->a_interpCoeff[i] = x->a_K[i];
		}
	}
	
	//solve for filter coefficients
	for (i = 1; i <= order; i++) {
		x->a_A[i][i] = x->a_interpCoeff[i];
		i1 = i - 1;
		if (i1 >= 1) {
			for (j = 1; j <=i1; j++) {
				ji = i - j;
				x->a_A[j][i] = x->a_A[j][i1] - x->a_interpCoeff[i] * x->a_A[ji][i1];
			}				
		}
	}
	for (j = 1; j <=order; j++) {
		x->a_a[j-1] = x->a_A[j][order];
	}
	
	//perform filtering
	n = (int)(w[7]);
	
	while (n--) { 
		vDSP_vmul(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sve(x->a_tempVec,1,&sum,order);
		val =  x->a_G * in[out_idx] + sum;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		out[out_idx] = (float)val;
		out_idx++;
	}
	
	return (w+8);
}*/

/*t_int *allpole_perf_areaI(t_int *w) {
	
	t_float *in = (t_float *)(w[1]);
	t_float *coeffIn = (t_float *)(w[2]);
	t_float *coeffIdxIn = (t_float *)(w[3]);
	t_float *G_n = (t_float *)(w[4]);
	t_allpole *x = (t_allpole *)(w[5]);
	t_float *out = (t_float *)(w[6]);
	int n = (int)(w[7]);
	int order = x->a_order;
	int i, j, i1, ji, in_idx = 0, out_idx = 0;
	float val = 0.0;
	float sum = 0.0;
	float interpDiv; //interpolation denominator
	
	while(n--) {
		if(IS_NAN_FLOAT(G_n[in_idx])) G_n[in_idx] = 0.0;
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[7]);
	
	//deemphasis?
	if (x->a_deemph) {
		while (n--) {
			in[in_idx] = in[in_idx] + x->a_a1*x->a_y1;
			x->a_y1 = in[in_idx];
		
			in_idx++;
		}
		in_idx = 0;
		n = (int)(w[7]);
	}
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->a_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
				for(i = 0; i < order; i++) {
					x->a_Ar[i+1] = x->a_aBuff[i];
					//determine interpolation parameters
					x->a_interpInc[i+1] = (x->a_Ar[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
				}
			}
		}
		in_idx++;
	}
	
	//interpolate coefficients
	for (i=1; i <= order; i++) {
		if (fabs(x->a_interpCoeff[i] - x->a_Ar[i]) > fabs(x->a_interpInc[i])) {
			x->a_interpCoeff[i] += x->a_interpInc[i];
		} else {
			x->a_interpCoeff[i] = x->a_Ar[i];
		}
		//cpost("\ninterp[%d] = %f\tAr[%d]=%f\tinc[%d]=%f",i,x->a_interpCoeff[i],i,x->a_Ar[i],i,x->a_interpInc[i]);
	}
	
	x->a_G = 1.0; //recursively calculated G
	
	//solve for filter coefficients
	for (i = 1; i <= order; i++) {
		if (i != order)
			x->a_K[i] = -(x->a_interpCoeff[i+1] - x->a_interpCoeff[i])/(x->a_interpCoeff[i+1] + x->a_interpCoeff[i]);
		else
			x->a_K[i] = -0.7; //no losses at lips NOTE: need to make this changeable
			
		x->a_A[i][i] = x->a_K[i];
		x->a_G = x->a_G * (1.0 + -x->a_K[i]);
		i1 = i - 1;
		if (i1 >= 1) {
			for (j = 1; j <=i1; j++) {
				ji = i - j;
				x->a_A[j][i] = x->a_A[j][i1] - x->a_K[i] * x->a_A[ji][i1];
			}				
		}
	}
	for (j = 1; j <=order; j++) {
		x->a_a[j-1] = x->a_A[j][order];
	}
	
	//perform filtering
	n = (int)(w[7]);
	
	while (n--) { 
		vDSP_vmul(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sve(x->a_tempVec,1,&sum,order);
		val =  x->a_G * in[out_idx] + sum;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		out[out_idx] = (float)val;
		out_idx++;
	}
	
	return (w+8);
}*/

void allpole_init(t_allpole *x) {
	int i;
	
	x->a_a = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
	x->a_aBuff = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
	x->a_y = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
	x->a_tempVec = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
	x->a_Ar = (t_mbcfloat *) getbytes( MAX_ORDER * sizeof(t_mbcfloat));
	x->a_K = (t_mbcfloat *) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
	x->a_interpCoeff = (t_mbcfloat *) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
	x->a_interpInc = (t_mbcfloat *) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
	x->a_A = (t_mbcfloat **) getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat*));
	for(i=0; i<MAX_ORDER; i++) {
		x->a_A[i] = (t_mbcfloat *)getbytes( (MAX_ORDER + 1) * sizeof(t_mbcfloat));
	}
	//allpole_clear(x);
	
	x->a_a1 = 0.98;
}

void allpole_free(t_allpole *x) {
	dsp_free((t_pxobject *) x);
	allpole_free_arrays(x);
}

void allpole_free_arrays(t_allpole *x)
{
	int i;
	if (x->a_a) {
		freebytes((char *)x->a_a, MAX_ORDER * sizeof(t_mbcfloat));
		x->a_a = NULL;
	}
	if (x->a_aBuff) {
		freebytes(x->a_aBuff, MAX_ORDER * sizeof(t_mbcfloat));
		x->a_aBuff = NULL;
	}
	if (x->a_y) {
		freebytes((char *)x->a_y, MAX_ORDER * sizeof(t_mbcfloat));
		x->a_y = NULL;
	}
	if (x->a_tempVec) {
		freebytes((char *)x->a_tempVec, MAX_ORDER * sizeof(t_mbcfloat));
		x->a_tempVec = NULL;
	}
	if (x->a_Ar) {
		freebytes(x->a_Ar, MAX_ORDER * sizeof(t_mbcfloat));
		x->a_Ar = NULL;
	}
	if (x->a_K) {
		freebytes(x->a_K, (MAX_ORDER + 1) * sizeof(t_mbcfloat));
		x->a_K = NULL;
	}
	if (x->a_interpCoeff) {
		freebytes(x->a_interpCoeff, (MAX_ORDER + 1) * sizeof(t_mbcfloat));
		x->a_interpCoeff = NULL;
	}
	if (x->a_interpInc) {
		freebytes(x->a_interpInc, (MAX_ORDER + 1) * sizeof(t_mbcfloat));
		x->a_interpInc = NULL;
	}
	if (x->a_A) {
		for(i=0; i<MAX_ORDER; i++) {
			if (x->a_A[i]) {
				freebytes(x->a_A[i], (MAX_ORDER + 1) * sizeof(t_mbcfloat));
				x->a_A[i] = NULL;
			}
		}
		freebytes(x->a_A, MAX_ORDER * sizeof(t_mbcfloat*));
		x->a_A = NULL;
	}
}

void allpole_clear(t_allpole *x) {

	int i, j;
	for(i = 0; i < MAX_ORDER; i++) {
		x->a_a[i] = 0.0;
		x->a_aBuff[i] = 0.0;
		x->a_y[i] = 0.0;
		x->a_Ar[i] = 0.0;
		x->a_tempVec[i] = 0.0;
		for (j = 0; j < (MAX_ORDER+1); j++) {
			x->a_A[i][j] = 0.0;
		}
	}
	
	for(i = 0; i < (MAX_ORDER+1); i++) {
		x->a_K[i] = 0.0;
		x->a_interpCoeff[i] = 0.0;
		x->a_interpInc[i] = 0.0;
	}
	
	x->a_G = 0.0;
	x->a_y1 = 0.0;
}

// this function is called when the DAC is enabled, and "registers" a function
// for the signal chain. in this case, "allpole_perform"
void allpole_dsp64(t_allpole *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	//NOTE: need to specify parcor, filter, or area perform string !!!!!!
	x->a_fs = samplerate;
	x->a_vsize = maxvectorsize;
	
	allpole_clear(x);
	
	if (x->a_interp > 0.0) {
		switch (x->a_coeffType) {
			case CT_FILTER:
				//dsp_add(allpole_perf_filterI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
				//break;
			
			case CT_PARCOR:
                object_method(dsp64, gensym("dsp_add64"), x, allpole_perf_parcorI);
				break;
				
			case CT_AREA:
				//dsp_add(allpole_perf_areaI, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
				//break;
			
			default:
				error("mbc.allpole~: no coefficient type selected");
				break;
		}
	} else {
		switch (x->a_coeffType) {
			case CT_FILTER:
				//dsp_add(allpole_perf_filter, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
				//break;
			
			case CT_PARCOR:
                object_method(dsp64, gensym("dsp_add64"), x, allpole_perf_parcor);
				break;
				
			case CT_AREA:
				//dsp_add(allpole_perf_area, 7, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, x, sp[4]->s_vec, sp[0]->s_n);
				//break;
			
			default:
				error("mbc.allpole~: no coefficient type selected");
				break;
		}
	}	
}

void allpole_interp(t_allpole *x, t_floatarg interp) {

	x->a_interp = interp;
	int i;
	int order = x->a_order;
	double interpDiv;
	
	switch (x->a_coeffType) {
			case CT_FILTER:
				error("mbc.allpole~: interpolation is inactive with coefficient type: filter");
				break;
			
			case CT_PARCOR:
				interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
				for(i = 0; i < order; i++) {
					x->a_interpInc[i+1] = (x->a_K[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
				}
				break;
				
			case CT_AREA:
				interpDiv = x->a_interp * 0.001 * x->a_fs / x->a_vsize;
				for(i = 0; i < order; i++) {
					x->a_interpInc[i+1] = (x->a_Ar[i+1] - x->a_interpCoeff[i+1]) / interpDiv;
				}
				break;
			
			default:
				break;
		}
}

void allpole_order(t_allpole *x, int order) {
	if (order < 1) {
		error("mbc.allpole~: mbc.allpole~ needs a positive integer value for order");
		order = 1;
	} else if (order > 199) {
		error("mbc.allpole~: max order is 199");
		order = 199;
	}
	
	int i = x->a_order;
	x->a_order = order;
	for ( i = i+1; i < order; i++) {
		x->a_a[i] = 0.0;
		x->a_aBuff[i] = 0.0;
		x->a_y[i] = 0.0;
		x->a_Ar[i] = 0.0;
	}
}

void allpole_deemph(t_allpole *x, int deemph) {
	x->a_deemph = deemph;
}

void allpole_assist(t_allpole *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		switch (a) {
			case 0: sprintf(s,"(signal) Filter Input"); break;
			case 1: sprintf(s,"(signal) PARCOR Coefficients"); break;
			case 2: sprintf(s,"(signal) Coeff Index"); break;
			case 3: sprintf(s,"(signal) Filter Gain"); break;
		}
	}
	else {
		sprintf(s,"(signal) Filter Output");
	}
}

void *allpole_new(t_symbol *s, long argc, t_atom *argv)
{
	t_allpole *x = NULL;
	
	if (argc < 1)
	{
		error("mbc.allpole~: must specify filter order (this should match mbc.lpc~ and mbc.errfilt~)");
		return NULL;
	}
	
	if ((x = (t_allpole *)object_alloc(allpole_class))) {
		dsp_setup((t_pxobject *)x, 1);	// MSP inlets: arg is # of inlets and is REQUIRED! 
										// use 0 if you don't need inlets
		
		dsp_setup((t_pxobject *)x,4);
		outlet_new(x, "signal");
		
		//get arguments out of gimme list
		int order;
		float interp = DEFAULT_INTERP;
		int deemph = DEFAULT_DEEMPH;
		//t_symbol* coeffType = atom_getsymarg(3,argc,argv);
		
		switch (argc)
		{
			case 0:
				break;
			
			case 1:
				order = atom_getintarg(0,argc,argv);
				break;
				
			case 2:
				order = atom_getintarg(0,argc,argv);
				interp = atom_getfloatarg(1,argc,argv);
				break;
				
			case 3:
				order = atom_getintarg(0,argc,argv);
				interp = atom_getfloatarg(1,argc,argv);
				deemph = atom_getintarg(2,argc,argv);
				break;
				
			default:
				order = atom_getintarg(0,argc,argv);
				interp = atom_getfloatarg(1,argc,argv);
				deemph = atom_getintarg(2,argc,argv);
				error("mbc.allpole~: too many arguments");
		}
		
		//order bounds
		if (order < 1) {
			error("mbc.allpole~: mbc.allpole~ needs a positive integer value for order");
			order = 1;
		} else if (order > 199) {
			error("mbc.allpole~: max order is 199");
			order = 199;
		}
		
		//parse filter coeff type
#if 0
		if (coeffType == gensym("parcor")) {
			x->a_coeffType = CT_PARCOR;
		} else if (coeffType == gensym("filter")) {
			x->a_coeffType = CT_FILTER;
		} else if (coeffType == gensym("area")) {
			x->a_coeffType = CT_AREA;
		} else {
			//error("mbc.allpole~: coefficient must be specified as 'parcor', 'filter', or 'area'");
			x->a_coeffType = CT_FILTER;
		}
#else
		// FORCE TO PARCOR FOR NOW (FOR SIMPLICITY)
		x->a_coeffType = CT_PARCOR;
#endif
		
		//assign locals to globals and init
		x->a_order = order;
		x->a_interp = interp;
		x->a_deemph = deemph;
		allpole_init(x);
	}
	return (x);
}
inline void allpole_deemphFilter(t_allpole *x, int N, double* vec)
{
	int n;
	
	double a1 = x->a_a1;
	double y1 = x->a_y1;
	
	for (n=N; n!=0; n--) 
	{
		*vec = *vec + (a1 * y1);
		y1 = *vec++;
	}
	
	x->a_y1 = y1;
}

inline void allpole_solveForFiltCoefs(t_allpole* x, int order)
{
	int i, i1, ji, j;
	
	for (i = 1; i <= order; i++) 
	{
		x->a_A[i][i] = x->a_K[i];
		i1 = i - 1;
		
		if (i1 >= 1) 
		{
			for (j = 1; j <=i1; j++) 
			{
				ji = i - j;
				x->a_A[j][i] = x->a_A[j][i1] - x->a_K[i] * x->a_A[ji][i1];
			}				
		}
	}
	
	for (j = 1; j <=order; j++) 
	{
		x->a_a[j-1] = x->a_A[j][order];
	}
}

#if defined(MBC_VDSP)
inline void allpole_highOrdFilter(t_allpole* x, int N, int order, double* in, double* out)
{
	int n, i;
	t_mbcfloat val;
		
	for (n=N; n!=0; n--) 
	{ 
		vDSP_vmulD(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sveD(x->a_tempVec,1,&val,order);
		val =  x->a_G * (t_mbcfloat)(*(in++)) + val;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		*(out++) = (float)val;
	}
}
#else
inline void allpole_highOrdFilter(t_allpole* x, int N, int order, double* in, double* out)
{
	int n, i;
	t_mbcfloat val;
		
	for (n=N; n!=0; n--) 
	{ 
		val = 0;
		t_mbcfoat* pA = x->a_a;
		t_mbcfloat* pY = x->a_y;
		
		for (i=0; i < order; i++)
		{
			val += *pY * 
		}
		vDSP_vmulD(x->a_a,1,x->a_y,1,x->a_tempVec,1,order);
		vDSP_sveD(x->a_tempVec,1,&val,order);
		val =  x->a_G * (t_mbcfloat)(*(in++)) + val;
		for (i=order-1; i>0; i--) x->a_y[i] = x->a_y[i-1];
		x->a_y[0] = val;		
			
		*(out++) = (float)val;
	}
}
#endif // MBC_VDSP
