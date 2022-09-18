/**
	@file
	mbc.coeffDisp~ - an MSP object shell
	mark cartwright - mcartwright@gmail.com	

	@ingroup	lpcToolkit
*/

#include "ext.h"							// standard Max include, always required (except in Jitter)
#include "ext_obex.h"						// required for new style objects
#include "z_dsp.h"							// required for MSP objects

#define MAX_ORDER 200

////////////////////////// object struct
typedef struct _coeffDisp 
{
	t_pxobject					ob;			// the object itself (t_pxobject in MSP)
	t_atom* 					c_a;		//filter coefficients
	long						c_len;			// length of filter coefficients
	t_float* 					c_aBuff;	//filter coeff input buffer
	long 						c_order;
	int 						c_coeffType;//type of coeffecients: CT_FILTER,CT_AREA
	int 						c_outMax;
	void						*c_clock;
	void* 						c_out;
} t_coeffDisp;

enum coeffTypes {
	CT_FILTER,
	CT_AREA
};

///////////////////////// function prototypes
//// standard set
void *coeffDisp_new(t_symbol *s, long argc, t_atom *argv);
void coeffDisp_free(t_coeffDisp *x);
void coeffDisp_assist(t_coeffDisp *x, void *b, long m, long a, char *s);

void coeffDisp_dsp(t_coeffDisp *x, t_signal **sp, short *count);
t_int *coeffDisp_perf_filter(t_int *w);
t_int *coeffDisp_perf_area(t_int *w);
void coeffDisp_order(t_coeffDisp *x, int order);
void coeffDisp_tick(t_coeffDisp *x);
void coeffDisp_init(t_coeffDisp *x);
void coeffDisp_clear(t_coeffDisp *x);
//////////////////////// global class pointer variable
void *coeffDisp_class;


int main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	// OLD METHOD
	// setup((t_messlist **)&coeffDisp_class, (method)coeffDisp_new, (method)dsp_free, (short)sizeof(t_coeffDisp), 0L, A_GIMME, 0);
	// addfloat((method)coeffDisp_float);
	// you need this
	// addmess((method)coeffDisp_dsp,				"dsp",			A_CANT, 0);
    // addmess((method)coeffDisp_assist,			"assist",		A_CANT, 0);  
	// you need this
    // dsp_initclass();
	
	// NEW METHOD
	t_class *c;
	
	c = class_new("mbc.coeffDisp~", (method)coeffDisp_new, (method)coeffDisp_free, (long)sizeof(t_coeffDisp), 0L, A_GIMME, 0);
	
	class_addmethod(c, (method)coeffDisp_dsp,		"dsp",		A_CANT, 0);
	class_addmethod(c, (method)coeffDisp_assist,	"assist",	A_CANT, 0);
	class_addmethod(c, (method)coeffDisp_order,		"order",	A_DEFLONG,0);
	
	class_dspinit(c);				// new style object version of dsp_initclass();
	class_register(CLASS_BOX, c);	// register class as a box class
	coeffDisp_class = c;
	
	return 0;
}

t_int *coeffDisp_perf_filter(t_int *w)
{
	t_float *coeffIn = (t_float *)(w[1]);
	t_float *coeffIdxIn = (t_float *)(w[2]);
	t_float *G = (t_float *)(w[3]);
	t_coeffDisp *x = (t_coeffDisp *)(w[4]);
	int n = (int)(w[5]);
	int order = x->c_order;
	int outMax = x->c_outMax;
	int i, in_idx = 0;
	
	while(n--) {
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[5]);
		
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->c_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				(x->c_a[0]).a_type = A_FLOAT;
				(x->c_a[0]).a_w.w_float = G[in_idx];
				(x->c_a[outMax+1]).a_type = A_FLOAT;
				(x->c_a[outMax+1]).a_w.w_float = 1.0;
				for (i = 1; i < outMax; i++) {
					(x->c_a[i]).a_type = A_FLOAT;
					(x->c_a[i]).a_w.w_float = 0.0;
					(x->c_a[i+outMax+1]).a_type = A_FLOAT;
					(x->c_a[i+outMax+1]).a_w.w_float = -x->c_aBuff[i-1];
				}
				x->c_len = (outMax + 1) * 2;
				clock_delay(x->c_clock, 0);
	
			}
		}
		in_idx++;
	}
	
	return (w+6);
}

t_int *coeffDisp_perf_area(t_int *w) //TODO: complete correctly (right now it's the same as coeff)
{
	t_float *coeffIn = (t_float *)(w[1]);
	t_float *coeffIdxIn = (t_float *)(w[2]);
	t_coeffDisp *x = (t_coeffDisp *)(w[3]);
	int n = (int)(w[4]);
	int order = x->c_order;
	int i, in_idx = 0;
	
	while(n--) {
		if(IS_NAN_FLOAT(coeffIn[in_idx])) coeffIn[in_idx] = 0.0;
		in_idx++;
	}
	
	in_idx = 0;
	n = (int)(w[4]);
		
	//look at coefficient index, if not zeros, buffer in coefficients
	while (n--) {
		if ((int)(coeffIdxIn[in_idx]) > 0) {
			x->c_aBuff[(int)(coeffIdxIn[in_idx])-1] = coeffIn[in_idx];
			if ((int)(coeffIdxIn[in_idx]) == order) {
				for (i = 0; i < order; i++) {
					(x->c_a[i]).a_type = A_FLOAT;
					(x->c_a[i]).a_w.w_float = x->c_aBuff[i];
				}
				x->c_len = order;
				clock_delay(x->c_clock, 0);
			}
		}
		in_idx++;
	}
	
	return (w+5);
}

void coeffDisp_tick(t_coeffDisp *x)
{
	if (x->c_len < (MAX_ORDER * 2)) {
		outlet_list(x->c_out, 0L, x->c_len, x->c_a);
	}
}

void coeffDisp_init(t_coeffDisp *x) {
	x->c_a = (t_atom *) getbytes( MAX_ORDER * 2 * sizeof(t_atom));
	x->c_aBuff = (t_float *) getbytes( MAX_ORDER * sizeof(t_float));
	coeffDisp_clear(x);
}

void coeffDisp_free(t_coeffDisp *x) {
	dsp_free((t_pxobject *) x);
	object_free(x->c_clock);
	freebytes(x->c_a, MAX_ORDER * 2 * sizeof(t_atom));
	freebytes(x->c_aBuff, MAX_ORDER * sizeof(t_float));
}

void coeffDisp_clear(t_coeffDisp *x) {
	int i;
	for(i = 0; i < MAX_ORDER; i++) {
		x->c_aBuff[i] = 0.0;
	}
	for (i = 0; i < ((x->c_outMax + 1)*2); i++) {
		(x->c_a[i]).a_type = A_FLOAT;
		(x->c_a[i]).a_w.w_float = 0.0;
	}
}

void coeffDisp_dsp(t_coeffDisp *x, t_signal **sp, short *count)
{
	switch (x->c_coeffType) {
		case CT_FILTER:
			dsp_add(coeffDisp_perf_filter, 5, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, x, sp[0]->s_n);
			break;
			
		case CT_AREA:
			dsp_add(coeffDisp_perf_area, 4, sp[0]->s_vec, sp[1]->s_vec, x, sp[0]->s_n);
			break;
		
		default:
			error("mbc.coeffDisp~: no coefficient type selected");
	}
	
	coeffDisp_clear(x);
}

void coeffDisp_order(t_coeffDisp *x, int order) {
	if (order < 1) {
		error("mbc.coeffDisp~: mbc.coeffDisp~ needs a positive integer value for order");
		order = 1;
	} else if (order > 199) {
		error("mbc.coeffDisp~: max order is 199");
		order = 199;
	}
	
	x->c_outMax = order;
	if (order > 48) {
		//error("mbc.coeffDisp~: the maximum filter order than Max/MSP's filtergraph object can handle is 48.  Therefore, only the first 48 coefficients will be shown");
		x->c_outMax = 48;
	}
	
	int i = x->c_order;
	x->c_order = order;
	for ( i = i+1; i < order; i++) {
		x->c_aBuff[i] = 0.0;
	}
	
	for (i = 0; i < ((x->c_outMax + 1)*2); i++) {
		(x->c_a[i]).a_type = A_FLOAT;
		(x->c_a[i]).a_w.w_float = 0.0;
	}
}

void coeffDisp_assist(t_coeffDisp *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		switch (a) {
			case 0: sprintf(s,"(signal) Coefficients"); break;
			case 1: sprintf(s,"(signal) Coeff Index"); break;
			case 2: sprintf(s,"(signal) Filter Gain"); break;
		}
	}
	else {
		sprintf(s,"(list) Formatted Coefficient Output");
	}
}

void *coeffDisp_new(t_symbol *s, long argc, t_atom *argv)
{
	t_coeffDisp *x = NULL;
	
	if (argc < 1)
	{
		error("mbc.coeffDisp~: must specify filter order (this should match mbc.lpc~)");
		return NULL;
	}

	if (x = (t_coeffDisp *)object_alloc(coeffDisp_class)) {
		dsp_setup((t_pxobject *)x,3);
		x->c_out = listout(x);
		
		//get arguments out of gimme list
		int order = atom_getintarg(0,argc,argv);
		t_symbol* coeffType = atom_getsymarg(1,argc,argv);
		
		x->c_clock = clock_new(x, (method)coeffDisp_tick);
		
		//order bounds
		if (order < 1) {
			error("mbc.coeffDisp~: mbc.coeffDisp~ needs a positive integer value for order");
			order = 1;
		} else if (order > 199) {
			error("mbc.coeffDisp~: max order is 199");
			order = 199;
		}
		
		x->c_outMax = order;
		if (order > 48) {
			error("mbc.coeffDisp~: the maximum filter order than Max/MSP's filtergraph object can handle is 48.  Therefore, only the first 48 coefficients will be shown");
			x->c_outMax = 48;
		}
		
		//parse filter coeff type
		if (coeffType == gensym("filter")) {
			x->c_coeffType = CT_FILTER;
		} else if (coeffType == gensym("area")) {
			x->c_coeffType = CT_AREA;
		} else {
			x->c_coeffType = CT_FILTER;
		}
		
		//assign locals to globals and init
		x->c_order = order;
		coeffDisp_init(x);
	}
	return (x);
}
