/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Fousek, FEE CTU Prague

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/

#ifndef opts_h
#define opts_h

class opts {

	// class opts
	//
	// maintains all configuration information,
	// loads configuration from config file and from command line
	// checks for correctness

	int	parse(char*,char*); 		// parse one cmdline option (left, [right] side)
	int check_config();				// check settings after parsing all options
	void settings();            	// dump settings
	enum IOendians { big, little };
	void usage();					// print usage text
	void set_preset();				// set presets of options
	
  public:

	enum NRwhen { beforeFB, afterFB };
    
    // IO
  	char format_in[39];  				// input file format
  	char format_out[39];  				// output file format
	IOendians endian_in;  				// input endian type
	IOendians endian_out; 				// output endian type
	bool pipe_in;						// online input
	bool pipe_out;						// online output
	float preem;						// preemphasis coeff (0 = none)
	int fs;           					// sampling frequency in Hz
 	int nfeacoefs;                      // # of input coeffs
	double dither;    					// additive dithering constant (0 = none)
	bool remove_dc;                     // remove DC from frame (right before FFT)
        char list[999];       				// list of files to be processed (in file)
    char in[999];       				// input file (in single file mode or online mode)
    char out[999];       				// output file (in single file mode or online mode)
    char pfilename[999];     			// name of output pfile (when applicable)
    char arkfilename[999];     			// name of output ark kaldi format(when applicable)
	
	// Segmentation
	double window_ms;						// window length in milliseconds (by user)
	double wshift_ms;       				// window shift in milliseconds
		
	// Filter Bank	
	char fb_scale[39];					// filter bank freq scale type
 	char fb_shape[39];					// shape of one filter from filter bank
 	bool fb_norm;                       // normalize filter area (before eqld)
	bool fb_power;						// power spectrum to FB?
 	bool fb_eqld;                       // apply Equal Loudness?
 	bool fb_inld;                       // apply Intensity-Loudness Power Law?
	char fb_definition[9999]; 		    // filterbank description string
	bool fb_printself;                  // undocumented option which prints out the filter bank in ASCII

	// Noise Reduction    
	char vadmode[39];					// VAD mode (class NR is responsible)
 	char filevad[999];    				// optional VAD filename
 	char nr_mode[39];					// SS mode (exten,ss2fw,... class NR responsible)
 	double nr_p;         				// integrator constant
	double nr_q;         				// threshold constant (used for VAD)
	double nr_a;         				// spectrum powering constant
	double nr_b;         				// spectrum oversubtraction constant
 	int nr_initsegs;                    // # of initial speechless frames in every file
 	char nr_rasta[999];  				// RASTA filter file
 	bool rasta;                         // perform RASTA?
 	NRwhen nr_when;                     // when to do NR
 	
 	// Parametrization
 	char fea_kind[39];                  // feature extraction method (spec/lpa/lpc/dctc)
 	int fea_lporder;                    // LP order
 	int fea_ncepcoefs;                  // # of cepstral coeffs
 	bool fea_c0;                        // write 0th ceps coeff flag
 	bool fea_E;                         // write Energy flag
 	bool fea_rawenergy;                 // energy should be computed from input segment
	int fea_lifter;       				// cepstrum liftering coefficient
 	char preset[39];                    // output features presets (mfcc, plpc....)
        int fea_trapdct_traplen;            // length of the TRAP to which the DCT is applied
        int fea_trapdct_ndct;               // number of the first DCT coeffs to output (per band)
 
        // Voice Activity Detection
        char vad_apply_mode[39];            // turn on VAD and silence/drop non-speech segments
        char vad_out_mode[39];              // turn on VAD and save VAD output file

        char vad_out[999];                  // VAD output file for single-file mode

        char vad_cri_mode[39];              // VAD criterion function
        char vad_thr_mode[39];              // VAD thresholding method

        bool vad_energy_db;                 // normalize energy in decibels

        char vad_cepdist_mode[39];          // cepstral distance mode
        double vad_cepdist_p;               // cepstral background estimation smoothness
        int vad_cepdist_init;               // number of segments to force updating cepstral background estimation

        int vad_lpc_coefs;                  // number of LPC cepstral coefficients

        double vad_absolute_thr;            // absolute threshold placement
        int vad_perc_init;                  // number of initial segments when criterial value is not confident
        double vad_perc_thr;                // percentual threshold placement

        int vad_adapt_init;                 // number of voiceless init segments for adaptive threshold
        double vad_adapt_q;                 // threshold remembering coefficient q
        double vad_adapt_za;                // speech to background distance parameter za

        int vad_dyn_init;
        double vad_dyn_perc;                // vad dynamic threshold percentrual placement
        double vad_dyn_min;
        double vad_dyn_qmaxinc;
        double vad_dyn_qmaxdec;
        double vad_dyn_qmindec;
        double vad_dyn_qmininc;

        int vad_filter_order;               // vad output median filter order

 	// Misc.
	int window;							// window length in samples (will be computed)
	int wshift;       					// window shift in samples
	int wfft;         					// FFT size
	int wfftby2;                        // non-redundant FFT size (half fft + 1)
 	bool verbose;   					// verbose mode switch
	bool quiet;       					// quiet mode switch
	bool info;                          // switch for program config printout
 	char config[999];                   // config file name
 	bool natural_little;       			// natural byte order (autodetected)
	bool swap_in; 						// byte order swapping (autodetected from -endian_Xs)
	bool swap_out; 						// byte order swapping (autodetected from -endian_Xs)
	char *version;    					// program version
	bool phase_needed;                  // flag whether to compute phase

        float fea_Z_exp;                    // write p flag ->  CMN with exponential moving  average
        float fea_Z_block;                  // write t flag -> length of window in milliseconds 
        int length_b;
	bool remove_dc1;                     // remove DC from frame (right before FFT)
        bool stat_cmvn;
        char fcmvn_stat_out[999];
        bool apply_cmvn;
        char fcmvn_stat_in[999];
        int d_win;          // delta window size
        int a_win;          // delta-delta window size
        int t_win;          // delta-delta-delta window size
        int delta_w;
        bool fea_delta;
        int n_order;
        bool fea_trap;
        int trap_win;
        char ffilters[999];
        float weight_of_td_iir_mfcc_bank;

	opts (int&, char**&, char*);     	// opts(argc,argv,version)
	~opts();
  };


#endif //opts
