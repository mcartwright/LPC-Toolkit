/*----------------------------------------------------------------------------------
Filename:		mbc.blit~.c
Project:		LPC Toolkit
Author:			Mark Cartwright
Created:		5/15/07
Updated:		10/21/10
Description:	band-limited impulse train generator external object for Max/MSP.  
				Uses the Sum of Windowed Sinc (SWS) method for blit generation.  See
				http://ccrma.stanford.edu/~stilti/papers/TimStilsonPhDThesis2006.pdf 
-------------------------------------------------------------------------------------*/


#include "ext.h"							// standard Max include, always required (except in Jitter)
#include "ext_obex.h"						// required for new style objects
#include "z_dsp.h"							// required for MSP objects

#define DEFAULT_PINC 440
#define DEFAULT_P 64;
#define DEFAULT_ZCS 1024;
#define DEFAULT_BANDLIMIT 1;

typedef struct _pulse 
{
	int	active;
	float phase;
} t_pulse;

////////////////////////// object struct
typedef struct _blit 
{
	t_pxobject					ob;				// the object itself (t_pxobject in MSP)
    float						b_pinc;			// == frequency
	double						b_phase;		// global phase
	t_float*					b_sincTable;	// wavetable
	t_pulse*					b_pulse;		// array of pulses
	double						b_fs;			// sampling rate
	int							b_zcs;			// samples per zero crossin in sinc table
	int							b_P;			// P in sinc calc... number of zero crossings
	long						b_length;		// table length
	int							b_bandlimit;	// bandlimiting on/off
} t_blit;

///////////////////////// function prototypes
//// standard set
void *blit_new(t_symbol *s, long argc, t_atom *argv);
void blit_free(t_blit *x);
void blit_assist(t_blit *x, void *b, long m, long a, char *s);

void blit_dsp(t_blit *x, t_signal **sp, short *count);
t_int *blit_sigperf(t_int *w);		//signal frequency input - bandlimited
t_int *blit_fltperf(t_int *w);		//float frequency input - bandlimited
t_int *blit_sigperf_a(t_int *w);	//signal frequency input - aliased
t_int *blit_fltperf_a(t_int *w);	//float frequency input - aliased
void blit_float(t_blit *x, double f);
void blit_int(t_blit *x, long n);
void blit_sincGen(t_blit *x);
void blit_winGen(double *win, long N);
//////////////////////// global class pointer variable
void *blit_class;


int main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	// OLD METHOD
	// setup((t_messlist **)&blit_class, (method)blit_new, (method)dsp_free, (short)sizeof(t_blit), 0L, A_GIMME, 0);
	// addfloat((method)blit_float);
	// you need this
	// addmess((method)blit_dsp,				"dsp",			A_CANT, 0);
    // addmess((method)blit_assist,			"assist",		A_CANT, 0);  
	// you need this
    // dsp_initclass();
	
	// NEW METHOD
	t_class *c;
	
	c = class_new("mbc.blit~", (method)blit_new, (method)blit_free, (long)sizeof(t_blit), 0L, A_GIMME, 0);
	
	class_addmethod(c, (method)blit_float,		"float",	A_FLOAT, 0);
	class_addmethod(c, (method)blit_int, "int", A_LONG, 0);
	class_addmethod(c, (method)blit_dsp,		"dsp",		A_CANT, 0);
	class_addmethod(c, (method)blit_assist,	"assist",	A_CANT, 0);
	
	class_dspinit(c);				// new style object version of dsp_initclass();
	class_register(CLASS_BOX, c);	// register class as a box class
	blit_class = c;
	
	return 0;
}

t_int *blit_sigperf(t_int *w) {
	int i;

	t_float *pinc = (t_float *)(w[1]);
	t_blit *x = (t_blit *)(w[2]);
	t_float *out = (t_float *)(w[3]);
	int n = (int)(w[4]);
	float thresh = x->b_fs;
	double phase = x->b_phase;
	float phasen1 = phase;
	long length = x->b_length;
	int zcs = x->b_zcs;
	float offset, step, idx, eta;
	int P = x->b_P;
	int PO2 = P / 2;
		
	while (n--) {
		phase += *pinc;
		if (phase >= thresh) {
			offset = 1.0 - ((thresh - phasen1)/(phase - phasen1));
			for (i = 0; i < PO2; i++) {
				if (!(x->b_pulse[i].active)) {
					x->b_pulse[i].active = 1;
					x->b_pulse[i].phase = offset;
					break;
				}
			}
			phase -= thresh;
		}
		phasen1 = phase;
		
		//render active pulses
		*out = 0.0;
		for (i = 0; i < PO2; i++) {
			if (x->b_pulse[i].active) {
				step = x->b_pulse[i].phase;
				idx = step * zcs;
				eta = idx - floor(idx);
				idx = floor(idx);
				*out += (float)(((1.0 - eta) * x->b_sincTable[(long)(idx)] + eta * x->b_sincTable[((long)(idx + 1.0)) % length]) * 0.89); //linear interpolation, the 0.89 is a scaling factor to keep peak below 1 (and therefore no aliasing... good up to 20k
				step += 1.0;
				if (step > (float)(P-1)) {
					x->b_pulse[i].active = 0;
					x->b_pulse[i].phase = 0.0;
				} else {
					x->b_pulse[i].phase = step;
				}
			}
		}
		out++;
		pinc++;
	}
	
	x->b_phase = phase;
	
	return (w+5);
}


t_int *blit_fltperf(t_int *w) {
	int i;

	t_blit *x = (t_blit *)(w[1]);
	t_float *out = (t_float *)(w[2]);
	int n = (int)(w[3]);
	float thresh = x->b_fs;
	double phase = x->b_phase;
	float phasen1 = phase;
	float pinc = x->b_pinc;
	long length = x->b_length;
	int zcs = x->b_zcs;
	float offset, step, idx, eta;
	int P = x->b_P;
	int PO2 = P / 2;
		
	while (n--) {
		phase += pinc;
		if (phase >= thresh) {
			//find sub-sample offset
			offset = 1.0 - ((thresh - phasen1)/(phase - phasen1));
			for (i = 0; i < PO2; i++) {
				if (!(x->b_pulse[i].active)) {
					x->b_pulse[i].active = 1;
					x->b_pulse[i].phase = offset;
					break;
				}
			}
			phase -= thresh;
		}
		phasen1 = phase;
		
		//render active pulses
		*out = 0.0;
		for (i = 0; i < PO2; i++) {
			if (x->b_pulse[i].active) {
				step = x->b_pulse[i].phase;
				idx = step * zcs;
				eta = idx - floor(idx);
				idx = floor(idx);
				*out += (float)(((1.0 - eta) * x->b_sincTable[(long)(idx)] + eta * x->b_sincTable[((long)(idx + 1.0)) % length]) * 0.89); //linear interpolation, the 0.89 is a scaling factor to keep peak below 1 (and therefore no aliasing... good up to 20k
				step += 1.0;
				if (step > (float)(P-1)) {
					x->b_pulse[i].active = 0;
					x->b_pulse[i].phase = 0.0;
				} else {
					x->b_pulse[i].phase = step;
				}
			}
		}
		out++;
	}
	
	x->b_phase = phase;
	
	return (w+4);
}

t_int *blit_sigperf_a(t_int *w) {
	t_float *pinc = (t_float *)(w[1]);
	t_blit *x = (t_blit *)(w[2]);
	t_float *out = (t_float *)(w[3]);
	int n = (int)(w[4]);
	float thresh = x->b_fs;
	double phase = x->b_phase;
		
	while (n--) {
		phase += *pinc;
		if (phase >= thresh) {
			*out = 1.0;
			phase -= thresh;
		} else {
			*out = 0.0;
		}

		out++;
		pinc++;
	}
	
	x->b_phase = phase;
	
	return (w+5);
}


t_int *blit_fltperf_a(t_int *w) {
	t_blit *x = (t_blit *)(w[1]);
	t_float *out = (t_float *)(w[2]);
	int n = (int)(w[3]);
	float thresh = x->b_fs;
	double phase = x->b_phase;
	float pinc = x->b_pinc;
		
	while (n--) {
		phase += pinc;
		if (phase >= thresh) {
			*out = 1.0;
			phase -= thresh;
		} else {
			*out = 0.0;
		}

		out++;
	}
	
	x->b_phase = phase;
	
	return (w+4);
}

void blit_sincGen(t_blit *x) {
	int i;
	double P = (double) x->b_P;
	double zcs = (double) x->b_zcs;
	double M = 2*floor(P/2) + 1;
	long length = x->b_length;
	double win[length];
	t_float* pSincTable = x->b_sincTable;
	double*	pWin = win;
	
	//generate blackman-harris window
	blit_winGen(win,length);
	
	//since this equation makes a zerophase window, we need to split it into 2 to make it causal again
	//first half
	for (i = (floor(length/2) + 1); i <= length; i++) {
		*pSincTable = sin(PI*(i/zcs)*(M/P))/(P*sin(PI*(i/zcs)/P)) * (*pWin);
		pSincTable++;
		pWin++;
	}
	//second half
	for (i = 1; i < (floor(length/2) + 1); i++) {
		*pSincTable = sin(PI*(i/zcs)*(M/P))/(P*sin(PI*(i/zcs)/P)) * (*pWin);
		pSincTable++;
		pWin++;
	}	
}

void blit_winGen(double *win, long N) {
	//this function generates a blackman harris window (low resolution/high dynamic range)
	int i;
	double a0,a1,a2,a3;
	
	a0 = 0.3635819;
	a1 = 0.4891775;
	a2 = 0.1365995;
	a3 = 0.0106411;
	
	for (i = 0; i < N; i++) {
		win[i] = a0 - a1*cos((2*PI*i)/(N-1)) + a2*cos((4*PI*i)/(N-1)) - a3*cos((6*PI*i)/(N-1));
	}
}

void blit_dsp(t_blit *x, t_signal **sp, short *count)
{
	x->b_fs = sys_getsr();
	x->b_phase = 0;
	
	if (count[0]) { // perform signal based frequency update
		if (x->b_bandlimit)
			dsp_add(blit_sigperf, 4, sp[0]->s_vec, x, sp[1]->s_vec, sp[0]->s_n);
		else
			dsp_add(blit_sigperf_a, 4, sp[0]->s_vec, x, sp[1]->s_vec, sp[0]->s_n);
	} else { //perform interrupt based frequency update
		if (x->b_bandlimit)
			dsp_add(blit_fltperf, 3, x, sp[1]->s_vec, sp[1]->s_n);
		else
			dsp_add(blit_fltperf_a, 3, x, sp[1]->s_vec, sp[1]->s_n);
	}
	
}

void blit_float(t_blit *x, double f)
{
    x->b_pinc = (float)f;
}

void blit_int(t_blit *x, long n)
{
	blit_float(x,(double)n);
}

void blit_assist(t_blit *x, void *b, long m, long a, char *s)
{
	if (m==ASSIST_INLET) {
		sprintf(s,"(signal/float) Frequency");
	} else {
		sprintf(s,"(signal) bandlimited impulse train");
	}
}

void blit_init(t_blit* x) {
	int i;
	
	for(i = 0; i < (x->b_P / 2); i++) {
		x->b_pulse[i].active = 0;
		x->b_pulse[i].phase = 0.0;
	}
	
	x->b_phase = 0.0;
}

void blit_free(t_blit* x) {
	dsp_free((t_pxobject *) x);
	sysmem_freeptr(x->b_sincTable);
	freebytes(x->b_pulse, (x->b_P / 2) * sizeof(t_pulse));
}

void *blit_new(t_symbol *s, long argc, t_atom *argv)
{
	t_blit *x = NULL;

	if (x = (t_blit *)object_alloc(blit_class)) {
		dsp_setup((t_pxobject *)x,1);
		outlet_new((t_pxobject *)x, "signal");
		
		x->b_pinc = DEFAULT_PINC;
		x->b_P = DEFAULT_P;
		x->b_zcs = DEFAULT_ZCS;
		x->b_bandlimit = DEFAULT_BANDLIMIT;
		//get arguments out of gimme list
		switch(argc) {
			case 0:
				break;
				
			case 1:
				x->b_pinc = atom_getfloatarg(0,argc,argv);
				break;
		
			case 2:
				x->b_pinc = atom_getfloatarg(0,argc,argv);
				x->b_P = atom_getintarg(1,argc,argv);
				break;
			
			case 3:
				x->b_pinc = atom_getfloatarg(0,argc,argv);
				x->b_P = atom_getintarg(1,argc,argv);
				x->b_zcs = atom_getintarg(2,argc,argv);
				break;
				
			case 4:
				x->b_pinc = atom_getfloatarg(0,argc,argv);
				x->b_P = atom_getintarg(1,argc,argv);
				x->b_zcs = atom_getintarg(2,argc,argv);
				x->b_bandlimit = atom_getintarg(3,argc,argv);
				break;
			
			default:
				x->b_pinc = atom_getfloatarg(0,argc,argv);
				x->b_P = atom_getintarg(1,argc,argv);
				x->b_zcs = atom_getintarg(2,argc,argv);
				x->b_bandlimit = atom_getintarg(3,argc,argv);
				error("mbc.blit~: too many arguments");
				break;
		}
		
		if ((int)(pow(2.0,round(log2((double)x->b_P)))) != x->b_P) {
			error("mbc.blit~: P argument must be a power of 2");
			x->b_P = pow(2.0,floor(log2(x->b_P)));
		}
	
		x->b_length = (long)(x->b_P * x->b_zcs - 1);
		x->b_sincTable = (t_float *) sysmem_newptrclear(x->b_length * sizeof(t_float));
		x->b_pulse = (t_pulse *) getbytes( (x->b_P / 2) * sizeof(t_pulse));
		blit_init(x);
		blit_sincGen(x);
	}
	return (x);
}
