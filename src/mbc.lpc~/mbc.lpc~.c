/**
	@file
	mbc.lpc~ - an MSP object shell
	mark cartwright - mcartwright@gmail.com	

	@ingroup	lpcToolkit	
*/

#include "ext.h"							// standard Max include, always required (except in Jitter)
#include "ext_obex.h"						// required for new style objects
#include "z_dsp.h"							// required for MSP objects
#include <math.h>
#include <Accelerate/Accelerate.h>

#define MAX_ORDER 200
#define NLOG2(x) (ceil(log2(x)))
#define POW2(x) (1 << x);
#define DEFAULT_FS 44100
#define DEFAULT_FRAMERATE 100
#define DEFAULT_V_SIZE 64

////////////////////////// object struct
typedef struct _lpc 
{
	t_pxobject					ob;					// the object itself (t_pxobject in MSP)
	t_float*					l_frame_buff;		//input frame buffer
	t_float*					l_winframe_buff;	//windowed input frame buffer
	t_float*					l_outCoeff_buff;	//coefficient signal
	t_float*					l_outParcor_buff; 	//PARCOR coeffs
	t_float*					l_outError_buff;	//error signal
	t_float*					l_win;				//analysis window
	t_float*					l_R;
	double*						l_W;
	double*						l_E;
	double*						l_K;
	double 						l_G;
	double**					l_A;
	double 						l_x1;				//last samples of pre-emphasis filter
	float 						l_b1;				//pre-emph coefficient
	int 						l_order;			//predictor order
	int 						l_order_max;		//max order according to fs = order * frame_rate
	int 						l_preemph;			//pre-epmhasis filter on/off
	int 						l_frame_rate;		//analysis frame rate
	int 						l_frame_size;		//analysis frame size, where fs = frame_rate * frame_size * 2
	int 						l_hop_size;			//hop_size = frame_size * 2 (b/c of overlap)
	int 						l_inframe_idx;		//current inframe buffer index
	int 						l_outframe_idx;		//current outframe buffer index
	long 						l_v_size;			//vector size
	float 						l_fs;				//sampling rate
	int 						l_maxnfft;			//fft length
	int 						l_log2nfft;			//log2(fft length)
	FFTSetup 					l_fftsetup;			//FFTSetup for vDSP FFT functions
	DSPSplitComplex 			l_fftSplitTemp;		//temp split complex vector structure
} t_lpc;

///////////////////////// function prototypes
//// standard set
void *lpc_new(long order, long framerate, long preemph);
void lpc_free(t_lpc *x);
void lpc_free_arrays(t_lpc *x);
void lpc_assist(t_lpc *x, void *b, long m, long a, char *s);

void lpc_dsp(t_lpc *x, t_signal **sp, short *count);
t_int *lpc_perform(t_int *w);

void lpc_order(t_lpc *x, int n);
void lpc_preemph(t_lpc *x, int n);
void lpc_init(t_lpc *x);
void lpc_hanning(t_lpc *x);
void lpc_hamming(t_lpc *x);
void lpc_bartlett(t_lpc *x);

//////////////////////// global class pointer variable
void *lpc_class;


int main(void)
{	
	// object initialization, note the use of dsp_free for the freemethod, which is required
	// unless you need to free allocated memory, in which case you should call dsp_free from
	// your custom free function.

	// OLD METHOD
	// setup((t_messlist **)&lpc_class, (method)lpc_new, (method)dsp_free, (short)sizeof(t_lpc), 0L, A_GIMME, 0);
	// addfloat((method)lpc_float);
	// you need this
	// addmess((method)lpc_dsp,				"dsp",			A_CANT, 0);
    // addmess((method)lpc_assist,			"assist",		A_CANT, 0);  
	// you need this
    // dsp_initclass();
	
	// NEW METHOD
	t_class *c;
	
	c = class_new("mbc.lpc~", (method)lpc_new, (method)lpc_free, (long)sizeof(t_lpc), 0L, A_DEFLONG, A_DEFLONG, A_DEFLONG, 0); //arglist: order, framerate, preemph
	
	class_addmethod(c, (method)lpc_dsp, "dsp", A_CANT, 0);
	class_addmethod(c, (method)lpc_order,"order",A_DEFLONG,0);
	class_addmethod(c, (method)lpc_preemph,"preemph",A_DEFLONG,0);
	class_addmethod(c, (method)lpc_assist,"assist",A_CANT,0);
	
	class_dspinit(c);				// new style object version of dsp_initclass();
	class_register(CLASS_BOX, c);	// register class as a box class
	lpc_class = c;
	
	return 0;
}

t_int *lpc_perform(t_int *w) 
{
	t_lpc *x = (t_lpc *) (w[1]);
	t_float *in = (t_float *)(w[2]);
	t_float *out_error = (t_float *)(w[3]);
	t_float *out_gain = (t_float *)(w[4]);
	t_float *out_coeff = (t_float *)(w[5]);
	t_float *out_parcor = (t_float *)(w[6]);
	t_float *out_index = (t_float *)(w[7]);
	int n = (int)(w[8]);
	int p = x->l_order;
	int length = x->l_frame_size;
	int log2nfft = NLOG2(length+p+1);
	int nfft = POW2(log2nfft);
	int nfftO2 = nfft/2;
	float recipn = 1.0/nfft;
	int inframeidx = x->l_inframe_idx;
	int outframeidx = x->l_outframe_idx;
	
	int i, j, i1, ji;
	int in_idx = 0, out_idx = 0;
	double val;
	float scale;
	
	if (x->l_preemph) 
	{
		while (n--) 
		{
			val = in[in_idx];
			in[in_idx] = val + x->l_b1 * x->l_x1;
			x->l_x1 = val;
			in_idx++;
		}
		n = (int)(w[8]);
		in_idx = 0;
	}
	
	while (n--) 
	{
		if (inframeidx < length) {
			//copy input into frame buff
			x->l_frame_buff[inframeidx] = in[in_idx];
			
			out_gain[out_idx] = x->l_G;
			out_error[out_idx] = x->l_outError_buff[outframeidx]; //for now
			if (outframeidx < x->l_order) {
				out_coeff[out_idx] = x->l_outCoeff_buff[outframeidx];
				out_parcor[out_idx] = x->l_outParcor_buff[outframeidx];
				out_index[out_idx] = outframeidx + 1;
			} else {
				out_coeff[out_idx] = 0.0;
				out_parcor[out_idx] = 0.0;
				out_index[out_idx] = 0;
			}
			
			inframeidx++;
			in_idx++;
			outframeidx++;
			out_idx++;
		} else {
			//perform durbin-levinson - for right now, just count to the order---------------
			//clear memory, is this necessary?
			for (i=0; i < p+1; i++){
				x->l_R[i] = 0.0;
				x->l_W[i] = 0.0;
				x->l_E[i] = 0.0;
				x->l_K[i] = 0.0;			
			}
			for(i=0; i<=p; i++) {
				for(j=0; j < p; j++) x->l_A[i][j] = 0.0;
			}
			//window frame buff
			vDSP_vmul(x->l_frame_buff, 1, x->l_win, 1, x->l_winframe_buff, 1, length);
			#ifdef DEBUG
				for(i=0;i<length;i++)cpost("\nwinframe(%d) = %g;",i+1,x->l_winframe_buff[i]);
			#endif
			
			//create r from auto correlation
			if ((2*nfft*log2(nfft)+nfft) > length*p) { //NOTE: change this to update only when order is changed!
			//time domain method
				for(i=0; i < p+1; i++) vDSP_dotpr(x->l_winframe_buff,1,x->l_winframe_buff+i,1,x->l_R+i,length-i);
			} else {
			//frequency domain method
				// zero pad
				vDSP_vclr(x->l_winframe_buff+length,1,nfft-length);
				//convert to split complex vector
				vDSP_ctoz( ( DSPComplex * ) x->l_winframe_buff, 2, &x->l_fftSplitTemp, 1, nfftO2);
				//perform forward in place fft
				vDSP_fft_zrip(x->l_fftsetup, &x->l_fftSplitTemp, 1, log2nfft, FFT_FORWARD);
				//scaling
				scale = 0.5;
				vDSP_vsmul(x->l_fftSplitTemp.realp,1,&scale,x->l_fftSplitTemp.realp,1,nfftO2);
				vDSP_vsmul(x->l_fftSplitTemp.imagp,1,&scale,x->l_fftSplitTemp.imagp,1,nfftO2);
				//compute PSD
				vDSP_zvmags(&x->l_fftSplitTemp,1,x->l_fftSplitTemp.realp,1,nfftO2);
				//clear imaginary part
				vDSP_vclr(x->l_fftSplitTemp.imagp,1,nfftO2);
				//perform inverse in place fft
				vDSP_fft_zrip(x->l_fftsetup, &x->l_fftSplitTemp, 1, log2nfft, FFT_INVERSE);
				//scaling
				vDSP_vsmul(x->l_fftSplitTemp.realp,1,&recipn,x->l_fftSplitTemp.realp,1,nfftO2);
				vDSP_vsmul(x->l_fftSplitTemp.imagp,1,&recipn,x->l_fftSplitTemp.imagp,1,nfftO2);
				//convert back to real number vector
				vDSP_ztoc(&x->l_fftSplitTemp, 1, (DSPComplex *)x->l_R, 2, nfftO2);	
			}
			
			
			/*for(i=0; i < p+1; i++) {
				x->l_R[i] = 0.0;
				for(j=0; j < length - i; j++) {
					x->l_R[i] += x->l_winframe_buff[j] * x->l_winframe_buff[i+j];
					//x->l_R[i] += x->l_win[j] * x->l_win[i+j];
				}
			}*/
			#ifdef DEBUG
					for(i=0;i< p+1; i++) cpost("\nR(%d) = %f;",i+1,x->l_R[i]);
			#endif
			
			x->l_W[0] = x->l_R[1];
			x->l_E[0] = x->l_R[0];
			
			for (i = 1; i <= p; i++) {
				x->l_K[i] = x->l_W[i-1] / x->l_E[i-1];
				
				x->l_A[i][i] = x->l_K[i];
				i1 = i - 1;
				if (i1 >= 1) {
					for (j = 1; j <=i1; j++) {
						ji = i - j;
						x->l_A[j][i] = x->l_A[j][i1] - x->l_K[i] * x->l_A[ji][i1];
					}				
				}
				
				x->l_E[i] = x->l_E[i-1] * (1.0 - x->l_K[i] * x->l_K[i]);
				
				if (i != p) {
					x->l_W[i] = x->l_R[i+1];
					for (j = 1; j <= i; j++) {
						x->l_W[i] -= x->l_A[j][i] * x->l_R[i-j+1];
					}
				}
			}
			
			x->l_G = sqrt(x->l_E[p]);
			for (i=0; i < p; i++) {
				x->l_outCoeff_buff[i] = (float)(x->l_A[i+1][p]);
				x->l_outParcor_buff[i] = (float)(x->l_K[i+1]);
				#ifdef DEBUG
					//cpost("\nParcor(%d) = %g;",i+1,x->l_K[i+1]);
					//cpost("\nCoeff(%d) = %g;",i+1,x->l_A[i+1][p]);
				#endif
			}
			
			//--------------------------------------------------------------------------------
			
			//copy right side to left side, move delayed input to output
			for (i=0; i < x->l_hop_size; i++) {
				x->l_outError_buff[i] = x->l_frame_buff[i];
				x->l_frame_buff[i] = x->l_frame_buff[i + x->l_hop_size];
			}
			
			inframeidx = x->l_hop_size;
			outframeidx = 0;
			n++; //to avoid skipping a sample (since we already decremented
			while (n--) {
				x->l_frame_buff[inframeidx] = in[in_idx];
				
				out_gain[out_idx] = (float)(x->l_G);
				out_error[out_idx] = x->l_outError_buff[outframeidx]; //for now
				if (outframeidx < x->l_order) {
					out_coeff[out_idx] = x->l_outCoeff_buff[outframeidx];
					out_parcor[out_idx] = x->l_outParcor_buff[outframeidx];
					out_index[out_idx] = (float)(outframeidx + 1);
				} else {
					out_coeff[out_idx] = 0.0;
					out_parcor[out_idx] = 0.0;
					out_index[out_idx] = 0;
				}
				
				inframeidx++;
				in_idx++;
				outframeidx++;
				out_idx++;
			}
			break;
		}
	}
	
	x->l_inframe_idx = inframeidx;
	x->l_outframe_idx = outframeidx;

	return (w+9);
}

void lpc_dsp(t_lpc *x, t_signal **sp, short *count)
{
	if (sp[0]->s_n != x->l_v_size || sp[0]->s_sr != x->l_fs) {
		x->l_v_size = sp[0]->s_n;
		x->l_fs = sp[0]->s_sr;
		lpc_init(x);
	}
	
	dsp_add(lpc_perform, 8, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[5]->s_vec, sp[0]->s_n);
}

void lpc_assist(t_lpc *x, void *b, long msg, long arg, char *dst)
{
	if (msg==ASSIST_INLET) {
		switch (arg) {
			case 0: sprintf(dst,"(signal) LPC Input"); break;
		}
	}
	else {
		switch (arg) {
			case 0: sprintf(dst,"(signal) Error Signal"); break;
			case 1: sprintf(dst,"(signal) Filter Gain"); break;
			case 2: sprintf(dst,"(signal) Filter Coefficients"); break;
			case 3: sprintf(dst,"(signal) PARCOR Coefficients"); break;
			case 4: sprintf(dst,"(signal) Coefficient Index"); break;
		}
	}
}

void lpc_free(t_lpc *x) {
	dsp_free((t_pxobject *) x);
	lpc_free_arrays(x);
}

void lpc_preemph(t_lpc *x, int n) {
	x->l_preemph = n;
}

void lpc_order(t_lpc *x, int n) {
	if (x->l_order_max < n) {
		error("mbc.lpc~: framerate * order must be less than or equal to sampling rate.  At the current setting maxorder is %d.  For a higher order, lower the framerate.", x->l_order_max);
		x->l_order = x->l_order_max;
	} else {
		x->l_order = n;
	}
	//lpc_init(x) ?
}

void lpc_init(t_lpc *x) {
	int i;
	
	if (x->l_fs < x->l_frame_rate * x->l_order) {
		x->l_frame_rate = (long)(x->l_fs / x->l_order);
		error("mbc.lpc~: framerate * order must be less than or equal to sampling rate.  framerate has been changed to %d", x->l_frame_rate);
	}
	
	x->l_hop_size = (int)(x->l_fs / x->l_frame_rate);
	
	if (x->l_v_size > x->l_hop_size) {
		error("mbc.lpc~: warning: frame_size is less than vector size.  this may cause errors, please reduce the frame rate or increase the vector size");
	}
	
	x->l_order_max = x->l_fs / x->l_frame_rate;
	x->l_frame_size = x->l_hop_size * 2;
	x->l_inframe_idx = x->l_hop_size;
	x->l_outframe_idx = 0;
	x->l_log2nfft = (int) NLOG2(x->l_frame_size+MAX_ORDER+1);
	x->l_maxnfft = POW2(x->l_log2nfft);
	
	// free memory if needed
	lpc_free_arrays(x);
	
	//allocate memory
	x->l_frame_buff = (t_float *) sysmem_newptrclear( x->l_frame_size * sizeof(t_float));
	x->l_winframe_buff = (t_float *) sysmem_newptrclear( x->l_maxnfft * sizeof(t_float));
	x->l_outCoeff_buff = (t_float *) sysmem_newptrclear(MAX_ORDER * sizeof(t_float));
	x->l_outParcor_buff = (t_float *) sysmem_newptrclear(MAX_ORDER * sizeof(t_float));
	x->l_outError_buff = (t_float *) sysmem_newptrclear( x->l_hop_size * sizeof(t_float));
	x->l_win = (t_float *) sysmem_newptrclear( x->l_frame_size * sizeof(t_float));
	x->l_R = (t_float *) sysmem_newptrclear( x->l_maxnfft * sizeof(t_float));
	x->l_W = (double *) sysmem_newptrclear( (MAX_ORDER + 1) * sizeof(double));
	x->l_E = (double *) sysmem_newptrclear( (MAX_ORDER + 1) * sizeof(double));
	x->l_K = (double *) sysmem_newptrclear( (MAX_ORDER + 1) * sizeof(double));
	x->l_A = (double **) sysmem_newptrclear( (MAX_ORDER + 1) * sizeof(double*));
	for(i=0; i<MAX_ORDER; i++) {
		x->l_A[i] = (double *)sysmem_newptrclear( (MAX_ORDER + 1) * sizeof(double));
	}
	x->l_fftSplitTemp.realp = (float *)sysmem_newptrclear( (x->l_maxnfft / 2) * sizeof(float));
	x->l_fftSplitTemp.imagp = (float *)sysmem_newptrclear( (x->l_maxnfft / 2) * sizeof(float));
		
	x->l_x1 = 0;
	x->l_G = 0.0;
	
	x->l_b1 = -0.98;
	
	//calculate window
	lpc_hamming(x);
	//create FFTsetups
	x->l_fftsetup = vDSP_create_fftsetup(x->l_log2nfft,FFT_RADIX2);
	if (!x->l_fftsetup) error("mbc.lpc~: not enough available memory");
	
}

void lpc_free_arrays(t_lpc *x)
{
	int i;
	
	if (x->l_frame_buff) {
		sysmem_freeptr(x->l_frame_buff);
		x->l_frame_buff = NULL;
	}
	if (x->l_winframe_buff) {
		sysmem_freeptr(x->l_winframe_buff);
		x->l_winframe_buff = NULL;
	}
	if (x->l_outCoeff_buff) {
		sysmem_freeptr(x->l_outCoeff_buff);
		x->l_outCoeff_buff = NULL;
	}
	if (x->l_outError_buff) {
		sysmem_freeptr(x->l_outError_buff);
		x->l_outError_buff = NULL;
	}
	if (x->l_outParcor_buff) {
		sysmem_freeptr(x->l_outParcor_buff);
		x->l_outParcor_buff = NULL;
	}
	if (x->l_win) {
		sysmem_freeptr(x->l_win);
		x->l_win = NULL;
	}
	if (x->l_R) {
		sysmem_freeptr(x->l_R);
		x->l_R = NULL;
	}
	if (x->l_W) {
		sysmem_freeptr(x->l_W);
		x->l_W = NULL;
	}
	if (x->l_E) {
		sysmem_freeptr(x->l_E);
		x->l_E = NULL;
	}
	if (x->l_K) {
		sysmem_freeptr(x->l_K);
		x->l_K = NULL;
	}
	if (x->l_A) {
		for(i=0; i<MAX_ORDER; i++) {
			if (x->l_A[i]) {
				sysmem_freeptr(x->l_A[i]);
				x->l_A[i] = NULL;
			}
		}
		sysmem_freeptr(x->l_A);
		x->l_A = NULL;
	}
	if (x->l_fftSplitTemp.realp) {
		sysmem_freeptr(x->l_fftSplitTemp.realp);
		x->l_fftSplitTemp.realp = NULL;
	}
	if (x->l_fftSplitTemp.imagp) {
		sysmem_freeptr(x->l_fftSplitTemp.imagp);
		x->l_fftSplitTemp.imagp = NULL;
	}

	if (x->l_fftsetup) {
		vDSP_destroy_fftsetup(x->l_fftsetup);
		x->l_fftsetup = NULL;
	}
}

void lpc_hanning(t_lpc *x) {
	int N = x->l_frame_size;
	int n;
	
	for (n=0; n<N; n++) {
		x->l_win[n] = 0.5 * (1.0 - cos((TWOPI*(n+1))/(N+1)));
	}
}

void lpc_hamming(t_lpc *x) {
	int N = x->l_frame_size;
	int n;
	
	for (n=0; n<N; n++) {
		x->l_win[n] = (0.54 - 0.46*cos((TWOPI*n)/(N-1)));
	}
}

void lpc_bartlett(t_lpc *x) {
	int N = x->l_frame_size;
	int n;
	
	for (n=0; n<N; n++) {
		x->l_win[n] = 2.0/(N - 1) * ((N-1)/2 - fabs(n - (N-1)/2.0));
	}
}

void *lpc_new(long order, long framerate, long preemph)
{
	t_lpc *x = NULL;
	
	if (order == 0)
	{
		error("mbc.lpc~: must specify the lpc order");
		return NULL;
	}
	
	if (framerate == 0)
	{
		framerate = DEFAULT_FRAMERATE;
	}
	
	if (x = (t_lpc *)object_alloc(lpc_class)) {
		dsp_setup((t_pxobject *)x,1);
		outlet_new((t_pxobject *)x, "signal");
		outlet_new((t_pxobject *)x, "signal");
		outlet_new((t_pxobject *)x, "signal");
		outlet_new((t_pxobject *)x, "signal");
		outlet_new((t_pxobject *)x, "signal");
		
		// from arguments:
		x->l_frame_rate = framerate;
		x->l_order = order;
		x->l_preemph = preemph;
		
		// system variables (note: these need to be checked again in lpc_dsp since they may differ if in a poly~ for instance)
		x->l_fs = sys_getsr();
		if (x->l_fs == 0)
		{
			x->l_fs = DEFAULT_FS;
		}
		
		x->l_v_size = sys_getblksize();
		if (x->l_v_size == 0)
		{
			x->l_v_size = DEFAULT_V_SIZE;
		}
		
		lpc_init(x);
	}
	return (x);
}
