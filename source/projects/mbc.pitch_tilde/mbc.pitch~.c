/*----------------------------------------------------------------------------------
Filename:		mbc.pitch~.c
Project:		LPC Toolkit
Author:			Mark Cartwright
Created:		5/22/07
Updated:		10/21/10
Description:	fundamental frequency estimator external object for Max/MSP.  
				Algorithm based of of IRCAM's Yin algorithm and McLeod's Tartini 
				algorithm.  See
				http://www.ircam.fr/pcm/cheveign/pss/2002_JASA_YIN.pdf
				and http://csweb.otago.ac.nz/tartini/papers/A_Smarter_Way_to_Find_Pitch.pdf
-------------------------------------------------------------------------------------*/

#include "ext.h"							// standard Max include, always required (except in Jitter)
#include "ext_obex.h"						// required for new style objects
#include "z_dsp.h"							// required for MSP objects
#include <math.h>
#include <Accelerate/Accelerate.h>

#define NLOG2(x) (ceil(log2(x)))
#define POW2(x) (1 << x);

#define DEFAULT_FRAME_RATE 30
#define DEFAULT_MINFREQ 60.0
#define DEFAULT_THRESH 0.1
#define DEFAULT_FS 44100
#define DEFAULT_V_SIZE 64

////////////////////////// object struct
typedef struct _pitch 
{
	t_pxobject					ob;					// the object itself (t_pxobject in MSP)
	t_float*					p_padframe_buff;	//padded input frame buffer
	t_float*					p_R;				//autocorrelation part
	t_float*					p_D;				//difference function
	t_float*					p_Dp;				//cumulative mean normalized difference function
	t_float*					p_M;				//the M in D = M - 2*R
	t_float*					p_tempVec;			//temporary storage for vector math
	float 						p_frame_rate;		//analysis frame rate
	int 						p_frame_size;		//analysis frame size, where fs = frame_rate * frame_size * 2
	int 						p_hop_size;			//hop_size = frame_size * 2 (b/c of overlap)
	int 						p_inframe_idx;		//current inframe buffer index
	int 						p_hop_idx;			//hop index
	long 						p_v_size;			//vector size
	float 						p_fs;				//sampling rate
	int 						p_maxnfft;			//fft length
	int 						p_log2nfft;			//log2(fft length)
	float 						p_thresh;			//detection threshold
	float 						p_minfreq;			//minimum detectable frequency
	int							p_overlap;			//samples of overlap (frame_size - hop_size)
	FFTSetup 					p_fftsetup;			//FFTSetup for vDSP FFT functions
	DSPSplitComplex 			p_fftSplitTemp;		//temp split complex vector structure
	float 						p_freq;				//estimated fundamental frequency
	float 						p_min;				//minimum at frequency location
	float 						p_clarity;			//clarity measure (~0-1... 1.0 for tonal) from Tartini
	void* 						p_freqOut;			//freq outlet
	void* 						p_clarityOut;		//clarity outlet
} t_pitch;

///////////////////////// function prototypes
//// standard set
void *pitch_new(t_symbol *s, long argc, t_atom *argv);
void pitch_free(t_pitch *x);
void pitch_free_arrays(t_pitch *x);
void pitch_assist(t_pitch *x, void *b, long m, long a, char *s);

void pitch_dsp(t_pitch *x, t_signal **sp, short *count);
t_int *pitch_perform(t_int *w);
void pitch_init(t_pitch *x);
void pitch_thresh(t_pitch *x, double thresh);
void pitch_frame_rate(t_pitch *x, double thresh);

//////////////////////// global class pointer variable
void *pitch_class;


int main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	// OLD METHOD
	// setup((t_messlist **)&pitch_class, (method)pitch_new, (method)dsp_free, (short)sizeof(t_pitch), 0L, A_GIMME, 0);
	// addfloat((method)pitch_float);
	// you need this
	// addmess((method)pitch_dsp,				"dsp",			A_CANT, 0);
    // addmess((method)pitch_assist,			"assist",		A_CANT, 0);  
	// you need this
    // dsp_initclass();
	
	// NEW METHOD
	t_class *c;
	
	c = class_new("mbc.pitch~", (method)pitch_new, (method)pitch_free, (long)sizeof(t_pitch), 0L, A_GIMME, 0);
	
	class_addmethod(c, (method)pitch_dsp,		"dsp",		A_CANT, 0);
	class_addmethod(c, (method)pitch_assist,	"assist",	A_CANT, 0);
	class_addmethod(c, (method)pitch_thresh,		"thresh",	A_FLOAT,0);
	class_addmethod(c, (method)pitch_frame_rate,	"frame_rate",A_FLOAT,0);
	
	class_dspinit(c);				// new style object version of dsp_initclass();
	class_register(CLASS_BOX, c);	// register class as a box class
	pitch_class = c;
	
	return 0;
}

t_int *pitch_perform(t_int *w) {	
	t_float *in = (t_float *)(w[1]);
	t_pitch *x = (t_pitch *) (w[2]);
	int n = (int)(w[3]);
	int length = x->p_frame_size;
	int hop_size = x->p_hop_size;
	int log2nfft = NLOG2(2*length);
	int nfft = POW2(log2nfft);
	int nfftO2 = nfft/2;
	float recipn = 1.0/nfft;
	int inframeidx = x->p_inframe_idx;
	int hopidx = x->p_hop_idx;
	float thresh = x->p_thresh;
	int overlap = x->p_overlap;
	float fs = x->p_fs;
	int i,j;
	float a,b,c,p,na,nb,nc,np;
	int in_idx = 0;
	float scale,sum;
	
	while (n--) {
		if (hopidx < hop_size) {
			hopidx++;
			in_idx++;
		} else {
			if (inframeidx < length) {
				//copy input into frame buff
				x->p_padframe_buff[inframeidx] = in[in_idx];

				inframeidx++;
				in_idx++;
			} else {
				//Step 1. Calculate Difference Function--------------------------------------------
				//1A) Autocorrelation via FFT------------------------------------------------------
							
				// zero pad
				vDSP_vclr(x->p_padframe_buff+length,1,nfft-length);
				//conver to split complex vector
				vDSP_ctoz( ( DSPComplex * ) x->p_padframe_buff, 2, &x->p_fftSplitTemp, 1, nfftO2);
				//perform forward in place fft
				vDSP_fft_zrip(x->p_fftsetup, &x->p_fftSplitTemp, 1, log2nfft, FFT_FORWARD);
				//scaling
				scale = 0.5;
				vDSP_vsmul(x->p_fftSplitTemp.realp,1,&scale,x->p_fftSplitTemp.realp,1,nfftO2);
				vDSP_vsmul(x->p_fftSplitTemp.imagp,1,&scale,x->p_fftSplitTemp.imagp,1,nfftO2);
				//compute PSD
				vDSP_zvmags(&x->p_fftSplitTemp,1,x->p_fftSplitTemp.realp,1,nfftO2);
				//clear imaginary part
				vDSP_vclr(x->p_fftSplitTemp.imagp,1,nfftO2);
				//perform inverse in place fft
				vDSP_fft_zrip(x->p_fftsetup, &x->p_fftSplitTemp, 1, log2nfft, FFT_INVERSE);
				//scaling
				vDSP_vsmul(x->p_fftSplitTemp.realp,1,&recipn,x->p_fftSplitTemp.realp,1,nfftO2);
				vDSP_vsmul(x->p_fftSplitTemp.imagp,1,&recipn,x->p_fftSplitTemp.imagp,1,nfftO2);
				//convert back to real number vector
				vDSP_ztoc(&x->p_fftSplitTemp, 1, (DSPComplex *)x->p_R, 2, nfftO2);
				
				//1B) Calculate M via recursion------------------------------------------------------ 
				x->p_M[0] = 2*x->p_R[0];
				for (i=1;i<length;i++) {
					x->p_M[i] = x->p_M[i-1] - (x->p_padframe_buff[length-i] * x->p_padframe_buff[length-i] + x->p_padframe_buff[i-1] * x->p_padframe_buff[i-1]);
				}
			
				//1C) Calculate D via M and R (D=M-2*R)----------------------------------------------
				//vDSP_vsub(x->p_M,1,x->p_R,1,x->p_D,1,length);
				//vDSP_vsub(x->p_D,1,x->p_R,1,x->p_D,1,length);
				for (i=0;i<length;i++) {
					x->p_D[i] = x->p_M[i] - x->p_R[i] - x->p_R[i];
				}
				
				//Step 2. -  Calculate cumulative mean norm difference function
				//dp = d/runningaverage(d)
				/*a = 1.0;
				b = 1.0;
				vDSP_vramp(&a,&b,x->p_tempVec,1,length);
				vDSP_vavlin(x->p_D,1,x->p_tempVec,x->p_Dp,1,length);
				vDSP_vdiv(x->p_D,1,x->p_Dp,1,x->p_Dp,1,length);
				x->p_Dp[0] = 1.0;*/
				for (i=1;i<length;i++) {
					sum = 0.0;
					for (j=0;j<i;j++) {
						sum += x->p_D[j];
					}
					sum = sum/i;
					x->p_Dp[i] = x->p_D[i]/sum;
				}
				x->p_Dp[0] = 1.0;
				
				//Step 3. find first (smallest) tau that is a local minimum below threshold-----------
				for (i=1;i<length;i++){
					if (x->p_Dp[i] < thresh) {
						if(x->p_Dp[i-1] > x->p_Dp[i]) {
							if(x->p_Dp[i] < x->p_Dp[i+1]) {
								//Step 4. - perform parabolic interpolation---------------------------
								a=x->p_Dp[i-1];
								b=x->p_Dp[i];
								c=x->p_Dp[i+1];
								p = 0.5*(c - a)/(2*b - a - c);
								x->p_freq = fs/(i+p);
								//x->p_min = b - 0.25*(a - c)*p;
								
								//NOTE: I probably don't need to do interpolation on this part...
								na=2*x->p_R[i-1]/x->p_M[i-1];
								nb=2*x->p_R[i]/x->p_M[i];
								nc=2*x->p_R[i+1]/x->p_M[i+1];
								np = 0.5*(nc - na)/(2*nb - na - nc);
								x->p_clarity = nb - 0.25*(na - nc)*np;
								//x->p_clarity = nb;
								
								outlet_float(x->p_freqOut,x->p_freq);
								outlet_float(x->p_clarityOut,x->p_clarity);
								break; 
							}
						}
					}
				}
				
				//copy right side to left side, move delayed input to output, this only works for 50% OL, change if we want more
				if (overlap > 0) {
					for (i=0; i < overlap; i++) {
						x->p_padframe_buff[i] = x->p_padframe_buff[i + x->p_hop_size];
					}
					inframeidx = overlap;
				} else {
					inframeidx = 0;
				}
				//hopidx = 0;
				hopidx = length;
				
				n++; //to avoid skipping a sample (since we already decremented
				while (n--) {
					if (hopidx < hop_size) {
						hopidx++;
						in_idx++;
					} else {
						//copy input into frame buff
						x->p_padframe_buff[inframeidx] = in[in_idx];
						
						inframeidx++;
						in_idx++;
					}
				}
				break;
			}
		}
	}
	
	x->p_inframe_idx = inframeidx;
	x->p_hop_idx = hopidx;

	return (w+4);
}

void pitch_dsp(t_pitch *x, t_signal **sp, short *count)
{
	if (sp[0]->s_n != x->p_v_size || sp[0]->s_sr != x->p_fs) {
		x->p_v_size = sp[0]->s_n;
		x->p_fs = sp[0]->s_sr;
		pitch_init(x);
	}
	
	dsp_add(pitch_perform, 3, sp[0]->s_vec, x, sp[0]->s_n);
}

void pitch_assist(t_pitch *x, void *b, long msg, long arg, char *dst)
{
	if (msg==ASSIST_INLET) {
		switch (arg) {
			case 0: sprintf(dst,"(signal) Input"); break;
		}
	}
	else {
		switch (arg) {
			case 0: sprintf(dst,"(Float) Detected Fundamental Frequency"); break;
			case 1: sprintf(dst,"(Float) Clarity"); break;
		}
	}
}

void pitch_thresh(t_pitch *x, double thresh) {
	x->p_thresh = (float)thresh;
	post("\nnew thresh = %f",x->p_thresh);
}

void pitch_frame_rate(t_pitch *x, double frame_rate) {
	x->p_frame_rate = (float)frame_rate;
	
	x->p_hop_size = (int) ceil(x->p_fs/frame_rate);
	
	if (x->p_v_size > x->p_hop_size) {
		error("mbc.pitch~: warning: frame_size is less than vector size.  this may cause errors, please reduce the frame rate or increase the vector size");
	}
	
	x->p_overlap = x->p_frame_size - x->p_hop_size;
	if (x->p_hop_idx > x->p_hop_size) x->p_hop_idx = 0;
	post("\nnew fr = %f",x->p_frame_rate);
}

void pitch_free(t_pitch *x) {
	dsp_free((t_pxobject *) x);
	pitch_free_arrays(x);
}

void pitch_free_arrays(t_pitch *x)
{
	if (x->p_padframe_buff) {
		sysmem_freeptr(x->p_padframe_buff);
		x->p_padframe_buff = NULL;
	}
	if (x->p_R) {
		sysmem_freeptr(x->p_R);
		x->p_R = NULL;
	}
	if (x->p_D) {
		sysmem_freeptr(x->p_D);
		x->p_D = NULL;
	}
	if (x->p_Dp) {
		sysmem_freeptr(x->p_Dp);
		x->p_Dp = NULL;
	}
	if (x->p_M) {
		sysmem_freeptr(x->p_M);
		x->p_M = NULL;
	}
	if (x->p_tempVec) {
		sysmem_freeptr(x->p_tempVec);
		x->p_tempVec = NULL;
	}
	if (x->p_fftSplitTemp.realp) {
		sysmem_freeptr(x->p_fftSplitTemp.realp);
		x->p_fftSplitTemp.realp = NULL;
	}
	if (x->p_fftSplitTemp.imagp) {
		sysmem_freeptr(x->p_fftSplitTemp.imagp);
		x->p_fftSplitTemp.imagp = NULL;
	}
	if (x->p_fftsetup) {
		vDSP_destroy_fftsetup(x->p_fftsetup);
		x->p_fftsetup = NULL;
	}
}

void pitch_init(t_pitch *x) {
	float fs = x->p_fs;
	
	x->p_frame_size = (int) ceil(fs/x->p_minfreq);
	x->p_hop_size = (int) ceil(fs/x->p_frame_rate);
	
	if (x->p_v_size > x->p_hop_size) {
		error("mbc.pitch~: warning: frame_size is less than vector size.  this may cause errors, please reduce the frame rate or increase the vector size");
	}
	
	x->p_overlap = x->p_frame_size - x->p_hop_size;
	x->p_inframe_idx = (x->p_overlap > 0) ? x->p_overlap : 0;
	x->p_hop_idx = 0;
	x->p_log2nfft = (int) NLOG2(2 * x->p_frame_size);
	x->p_maxnfft = POW2(x->p_log2nfft);
	
	// free memory if needed
	pitch_free_arrays(x);
	
	//allocate memory
	x->p_padframe_buff = (t_float *) sysmem_newptrclear( x->p_maxnfft * sizeof(t_float));
	x->p_R = (t_float *) sysmem_newptrclear( x->p_maxnfft * sizeof(t_float));
	x->p_D = (t_float *) sysmem_newptrclear( x->p_maxnfft * sizeof(t_float));
	x->p_Dp = (t_float *) sysmem_newptrclear( x->p_maxnfft * sizeof(t_float)); 
	x->p_M = (t_float *) sysmem_newptrclear( x->p_maxnfft * sizeof(t_float));
	x->p_tempVec = (t_float *) sysmem_newptrclear( x->p_maxnfft * sizeof(t_float));
	x->p_fftSplitTemp.realp = (float *)sysmem_newptrclear( (x->p_maxnfft / 2) * sizeof(float));
	x->p_fftSplitTemp.imagp = (float *)sysmem_newptrclear( (x->p_maxnfft / 2) * sizeof(float));
	
	//create FFTsetups
	x->p_fftsetup = vDSP_create_fftsetup(x->p_log2nfft,FFT_RADIX2);
	if (!x->p_fftsetup) error("mbc.pitch~: not enough available memory");
}

void *pitch_new(t_symbol *s, long argc, t_atom *argv)
{
	t_pitch *x = NULL;

	if (x = (t_pitch *)object_alloc(pitch_class)) {
		dsp_setup((t_pxobject *)x,1);
		x->p_clarityOut = floatout(x);
		x->p_freqOut = floatout(x);
		
		// defaults:
		x->p_frame_rate = DEFAULT_FRAME_RATE;
		x->p_minfreq = DEFAULT_MINFREQ;
		x->p_thresh = DEFAULT_THRESH;
		
		// from arguments:
		switch(argc) {
			case 0:
				break;
				
			case 1:
				x->p_frame_rate = atom_getfloatarg(0,argc,argv);
				break;
		
			case 2:
				x->p_frame_rate = atom_getfloatarg(0,argc,argv);
				x->p_minfreq = atom_getfloatarg(1,argc,argv);
				break;
			
			case 3:
				x->p_frame_rate = atom_getfloatarg(0,argc,argv);
				x->p_minfreq = atom_getfloatarg(1,argc,argv);
				x->p_thresh = atom_getfloatarg(2,argc,argv);
				break;
			
			default:
				break;
		}
	
		
		// system variables (note: these need to be checked again in pitch_dsp since they may differ if in a poly~ for instance)
		x->p_fs = sys_getsr();
		if (x->p_fs == 0)
		{
			x->p_fs = DEFAULT_FS;
		}
		
		x->p_v_size = sys_getblksize();
		if (x->p_v_size == 0)
		{
			x->p_v_size = DEFAULT_V_SIZE;
		}
		
		pitch_init(x);
	}
	return (x);
}
