/**
	@file
	errfilt - an FIR filter for mbc.lpc~
	mark cartwright - mcartwright@gmail.com	

	@ingroup	lpcToolkit	
*/

#include "ext.h"							// standard Max include, always required (except in Jitter)
#include "ext_obex.h"						// required for new style objects
#include "z_dsp.h"							// required for MSP objects
#include <Accelerate/Accelerate.h>

#define MAX_ORDER 200

////////////////////////// object struct
typedef struct _errfilt 
{
	t_pxobject					ob;				// the object itself (t_pxobject in MSP)
	t_float* 					a_b;			//filer coefficients
	t_float*					a_bBuff;		//filter coeff input buffer
	t_float* 					a_x;			//filter memory
	t_float* 					a_A;			//tube areas
	t_float* 					a_tempVec;		//holds temporary values for vector math
    float 						a_interp;
	int 						a_order;
} t_errfilt;

///////////////////////// function prototypes
//// standard set
void *errfilt_new(int order, float interp);
void errfilt_free(t_errfilt *x);
void errfilt_assist(t_errfilt *x, void *b, long m, long a, char *s);


void errfilt_dsp64(t_errfilt *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void errfilt_interp(t_errfilt *x, float interp);
void errfilt_order(t_errfilt *x, int order); 
void errfilt_init(t_errfilt *x);
void errfilt_clear(t_errfilt *x);
void errfilt_perf_coeff(t_errfilt *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);

//////////////////////// global class pointer variable
void *errfilt_class;


void ext_main(void *r)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	// OLD METHOD
	// setup((t_messlist **)&errfilt_class, (method)errfilt_new, (method)dsp_free, (short)sizeof(t_errfilt), 0L, A_GIMME, 0);
	// addfloat((method)errfilt_float);
	// you need this
	// addmess((method)errfilt_dsp,				"dsp",			A_CANT, 0);
    // addmess((method)errfilt_assist,			"assist",		A_CANT, 0);  
	// you need this
    // dsp_initclass();
	
	// NEW METHOD
	t_class *c;
	
	c = class_new("mbc.errfilt~", (method)errfilt_new, (method)errfilt_free, (long)sizeof(t_errfilt), 0L, A_DEFLONG, A_DEFFLOAT, 0);
	
	class_addmethod(c, (method)errfilt_dsp64, "dsp64", A_CANT, 0);
	class_addmethod(c, (method)errfilt_interp,"interp",A_DEFFLOAT,0);
	class_addmethod(c, (method)errfilt_order,"order",A_LONG,0);
	class_addmethod(c, (method)errfilt_assist,"assist",A_CANT,0);
	
	class_dspinit(c);				// new style object version of dsp_initclass();
	class_register(CLASS_BOX, c);	// register class as a box class
	errfilt_class = c;
}

void errfilt_perf_coeff(t_errfilt *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	double *in = ins[0];
	double *coeffIn = ins[1];
	double *coeffIdxIn = ins[2];
    
	double *out = outs[0];
	int n = sampleframes;
	int order = x->a_order;
	int i, in_idx = 0, out_idx = 0;
	t_float val = 0.0;
	float sum = 0.0;
	
	while(n--) {
		//get rid of NANs
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = sampleframes;
	
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->a_bBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for(i = 0; i < order; i++) {
					x->a_b[i] = x->a_bBuff[i];
				}
			}
		}
		in_idx++;
	}
	
	n = sampleframes;
	
	while (n--) { 
		val = in[out_idx];
		vDSP_vmul(x->a_b,1,x->a_x,1,x->a_tempVec,1,order);
		vDSP_sve(x->a_tempVec,1,&sum,order);
		val -= sum;
		//for (i=0; i < order; i++) val -= x->a_b[i] * x->a_x[i];
		for (i=order-1; i>0; i--) x->a_x[i] = x->a_x[i-1];
		x->a_x[0] = in[out_idx];		
			
		out[out_idx] = (float)val;
		out_idx++;
	}
}


void errfilt_init(t_errfilt *x) {
	x->a_b = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->a_bBuff = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->a_x = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->a_A = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	x->a_tempVec = (t_float *) getbytes( MAX_ORDER * sizeof(t_float)); 
	errfilt_clear(x);
}

void errfilt_clear(t_errfilt *x) {
	int i;
	for(i = 0; i < MAX_ORDER; i++) {
		x->a_b[i] = 0.0;
		x->a_bBuff[i] = 0.0;
		x->a_x[i] = 0.0;
		x->a_A[i] = 0.0;
		x->a_tempVec[i] = 0.0;
	}
}

void errfilt_interp(t_errfilt *x, float interp) {
	x->a_interp = interp;
}

void errfilt_order(t_errfilt *x, int order) {
	if (order < 1) {
		error("mbc.errfilt~: mbc.errfilt~ needs a positive integer value for order");
		order = 1;
	} else if (order > 199) {
		error("mbc.errfilt~: max order is 199");
		order = 199;
	}
	
	x->a_order = order;
	errfilt_clear(x);
}

// this function is called when the DAC is enabled, and "registers" a function
// for the signal chain. in this case, "errfilt_perform"
void errfilt_dsp64(t_errfilt *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{

    // instead of calling dsp_add(), we send the "dsp_add64" message to the object representing the dsp chain
    // the arguments passed are:
    // 1: the dsp64 object passed-in by the calling
    // 2: the symbol of the "dsp_add64" message we are sending
    // 3: a pointer to your object
    // 4: a pointer to your 64-bit perform method
    // 5: flags to alter how the signal chain handles your object -- just pass 0
    // 6: a generic pointer that you can use to pass any additional data to your perform method
    
    object_method(dsp64, gensym("dsp_add64"), x, errfilt_perf_coeff, 0, NULL);
}

void errfilt_assist(t_errfilt *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		switch (a) {
			case 0: sprintf(s,"(signal) Filter Input"); break;
			case 1: sprintf(s,"(signal) Coefficients"); break;
			case 2: sprintf(s,"(signal) Coeff Index"); break;
		}
	}
	else {
		sprintf(s,"(signal) Filter Output");
	}
}

void errfilt_free(t_errfilt *x) {
	dsp_free((t_pxobject *) x);
	freebytes((char *)x->a_b, MAX_ORDER * sizeof(t_float));
	freebytes(x->a_bBuff, MAX_ORDER * sizeof(t_float));
	freebytes((char *)x->a_x, MAX_ORDER * sizeof(t_float));
	freebytes(x->a_A, MAX_ORDER * sizeof(t_float));
	freebytes((char *)x->a_tempVec, MAX_ORDER * sizeof(t_float));
}

void *errfilt_new(int order, float interp)
{
	t_errfilt *x = NULL;
	
	if (order == 0)
	{
		error("mbc.errfilt~: must specify filter order (this should match the order of mbc.lpc~ and mbc.allpole~)");
		return NULL;
	}
	
    if ((x = (t_errfilt *)object_alloc(errfilt_class))) 
	{
		dsp_setup((t_pxobject *)x,3);
		outlet_new(x, "signal");
		
		if (order < 1) {
			error("mbc.errfilt~: mbc.errfilt~ needs a positive integer value for order");
			order = 1;
		} else if (order > 199) {
			error("mbc.errfilt~: max order is 199");
			order = 199;
		}
		
		x->a_order = order;
		x->a_interp = interp;
		errfilt_init(x);
	}
	
    return (x);
}
