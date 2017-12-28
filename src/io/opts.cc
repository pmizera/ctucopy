/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Fousek, FEE CTU Prague
 Copyright 2015 Petr Mizera, FEE CTU Prague

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

#include "stdafx.h"
#include "opts.h"
#include <fstream>
#include <math.h>


using namespace std;

opts::opts(int& argc, char **&argv, char *_version) {

	using namespace std;

     // first defaults are defined:

    // IO
   	strcpy(format_in,"");   // user has to choose
  	strcpy(format_out,"");  // user has to choose
	preem = 0.0;            // no preem
	fs = 0;
	dither = 0.0;			// add dither
        strcpy(list,"");
   	strcpy(in,"");
  	strcpy(out,"");
        strcpy(pfilename,"");
        strcpy(arkfilename,"");
	endian_in = little;
	endian_out = little;
	pipe_in = false;
	pipe_out = false;
	remove_dc = true;		// remove offset
	remove_dc1 = false;		// remove offset

	// Segmentation
	window_ms = 25.;
	wshift_ms = 10.;

	// Filter Bank
	strcpy(fb_scale,"mel");
 	strcpy(fb_shape,"triang");
	fb_power = true;
	fb_norm  = true;
 	fb_eqld = true;
 	fb_inld = true;
	strcpy(fb_definition,"26filters");
	fb_printself = false;



	// Noise Reduction
	strcpy(vadmode,"none");
 	strcpy(filevad,"");
 	strcpy(nr_mode,"none");
 	nr_p = 0.95;
	nr_q = 0.99;
	nr_a = 1.;
	nr_b = 1.;
 	nr_initsegs = 10;
 	strcpy(nr_rasta, "");
 	rasta = false;
 	nr_when = beforeFB;

 	// Parametrization
 	strcpy(fea_kind,"lpc");
 	fea_lporder = 12;
 	fea_ncepcoefs = 12;
 	fea_c0 = true;
 	fea_E = false;
        fea_Z_exp = -1;
        fea_Z_block = -1;
        stat_cmvn = false;
   	strcpy(fcmvn_stat_out,"");
        apply_cmvn = false;
   	strcpy(fcmvn_stat_in,"");
        d_win=2;
        a_win=2;
        t_win=2;
        fea_delta=false;
        n_order=0;
        fea_trap=false;
        trap_win=5;
        nfeacoefs = 13;
   	strcpy(ffilters,"");
 	fea_rawenergy = false;
	fea_lifter = 22;
 	strcpy(preset,"user");
        weight_of_td_iir_mfcc_bank=2.026;

	// Voice Activity Detection
        strcpy(vad_apply_mode, "none");
        strcpy(vad_out_mode, "none");

        strcpy(vad_out, "");

	strcpy(vad_cri_mode, "energy");
	strcpy(vad_thr_mode, "perc");

        vad_energy_db = true;

        strcpy(vad_cepdist_mode, "lpc");
        vad_cepdist_p = 0.8;
        vad_cepdist_init = 4;

        vad_lpc_coefs = 14;

	vad_absolute_thr = 1.0; // TODO default value
        vad_perc_init = 10;
	vad_perc_thr = 50.0;

        vad_adapt_init = 20;
        vad_adapt_q = 0.9;
        vad_adapt_za = 2.0;

        vad_dyn_init = 5;
        vad_dyn_perc = 50.0;
        vad_dyn_min = 1.0; // TODO default values
        vad_dyn_qmaxinc = 0.8;
        vad_dyn_qmaxdec = 0.995;
        vad_dyn_qmindec = 0.8;
        vad_dyn_qmininc = 0.9999;

        vad_filter_order = 3;

 	// Misc.
 	verbose = false;
	quiet = false;
	info = false;
 	strcpy(config,"");
 	version = _version;

  	// natural endianness test
	short word = 0x4321;
	if((*(char *)& word) != 0x21) natural_little = false;
	else natural_little = true;

	// Parse command line
	bool conf = false; // do we have to parse config file?
	if (argc == 1) {usage(); throw "OPTS: No command line options!"; }

	// look for a "-C " option in cmdline and assign filename
	for ( int j = 1 ; j < argc ; j++ ) {
		if (0 == strcmp(argv[j],"-C")) {
			strcpy(config,argv[++j]);
			conf = true;
		}
	}

	// 1st, parse config file if needed
	if (conf) {
		ifstream cfg(config); // open config
		if (!cfg) { 
			throw "OPTS: Cannot open config file!";
		} 
		char line[9999];	char *l,*r;
		//		std::string line;	char *l,*r;
		while (cfg.getline(line,9999)) {
			char *pchr = strstr (line,"#");
			if (pchr) line[pchr-line] = '\0'; // discard anything after char "#"
			l = strtok(line," \t");
			if (l) { 
				r = strtok (NULL, " \t"); 
				parse(l,r);
			}
		}
	}

	// command line arguments processing
	for (int j = 1; j < argc; j++) { 
		if (0==strncmp(argv[j],"-",1)) {
			if ((j+1)==argc) parse(argv[j],0);
			else
				if (0==strncmp(argv[j+1],"-",1)) parse(argv[j],0);
				else parse(argv[j],argv[j+1]);
		}
	}
	check_config();
};

void opts::set_preset() {
	// MFCC
	if (0==strcmp(preset,"mfcc")) {
		strcpy(fb_scale, "mel");
		strcpy(fb_shape, "triang");
		fb_power = true;
		strcpy(fb_definition,"1-26/26filters");
		strcpy(nr_mode, "none");
		rasta = false;
		fb_eqld = false;
		fb_inld = false;
		strcpy(fea_kind,"dctc");
		fea_ncepcoefs = 12;
		fea_c0 = true;
		fea_E = false;
		fea_lifter = 22;
		fea_rawenergy = false;
	} else
	// PLP
	if (0==strcmp(preset,"plpc")) {
		strcpy(fb_scale, "bark");
		strcpy(fb_shape, "trapez");
		fb_power = true;
		strcpy(fb_definition,"1-15/15filters");
		strcpy(nr_mode, "none");
		rasta = false;
		fb_eqld = true;
		fb_inld = true;
		strcpy(fea_kind,"lpc");
		fea_lporder = 12;
		fea_ncepcoefs = 12;
		fea_c0 = true;
		fea_E = false;
		fea_lifter = 22;
		fea_rawenergy = false;
	} else 

	// EXTEN enhancer
	if (0==strcmp(preset,"exten")) {
		window_ms = 32.; // segmentation 32/16 ms
		wshift_ms = 16.;
		strcpy(fb_definition,"none");
		strcpy(fb_scale, "none");
		strcpy(fb_shape, "none");
		nr_a = 2.;
		fb_eqld = false;
		fb_inld = false;
		fb_power = false;
		fb_norm = false;
		strcpy(nr_mode, "exten");
		strcpy(fea_kind,"none");
		fea_c0 = false;
		fea_E = false;
		fea_lifter = 0;
		fea_rawenergy = false;
	} else 
	throw ("OPTS: Unknown preset!");
}

int opts::check_config() {

	if (fs == 0) throw "OPTS: Please specify sampling rate!";

	window = (int) floor(.5 + window_ms/1000. * (double)fs);
	wshift = (int) floor(.5 + wshift_ms/1000. * (double)fs);

        if((stat_cmvn == true || apply_cmvn == true) && (fea_Z_block!=-1 || fea_Z_exp!=-1)){
          cerr << "OPTS: CMN or CMVN can not be set together with CMS option (exp or block CMS normalisations)!" << endl;
        }

        if(fea_Z_block!=-1 && fea_Z_exp!=-1){
          cerr << "OPTS: Please set exp CMS or block CMS normalisation method!" << endl;
        }

        if(fea_Z_block!=-1)
  	   length_b = (int) floor((fea_Z_block-window_ms)/wshift_ms)+1;

        if(fea_Z_exp!=-1)
            fea_Z_exp = (float)1-(2*wshift_ms)/fea_Z_exp;

	// compute fft size
	for (int i = 1048576; i > 4; i /= 2) {
		if ((window/i) == 1)
			wfft = i*(1 + ((window % i) != 0));
	}

	// half FFT size
	wfftby2 = int(wfft/2) + 1;

	// set byte swapping
	swap_in  = (endian_in == little)  ^ natural_little; // XOR
	swap_out = (endian_out == little) ^ natural_little;

	// decide whether we need to compute phase
	if (0==strcmp(format_out,"raw")||0==strcmp(format_out,"wave")) 
		phase_needed = true; else phase_needed = false;

	// Brutal hack for Burg VAD (see nr.cc for info)
	if (0==strcmp(vadmode,"burg")) phase_needed = true;

	if (preem >= 1.0 || preem < 0.0) throw ("OPTS: Preemphasis not in range <0,1)!");

//	if ((strcmp(in,"")==0)^(0==strcmp(out,""))) 
//		throw "OPTS: Both -i and -o options have to be set in single file mode!";

	if ((pipe_in || 0!=strcmp(in,"")) ^	(pipe_out || 0!=strcmp(out,""))) 
		throw "OPTS: Single file mode has to be set at both sides (input and output)!";

	if (pipe_out && (0==strcmp(format_out,"wave") || 0==strcmp(format_out,"pfile")))
		throw "OPTS: Online output available only for raw and htk formats!";

	//	if (mod == exten && vadmod != not_specified) // VAD specified even when exten
//		if (!nowarning)
//			cerr << "\tWARNING: Using `exten' processing mode, VAD overriden!" << endl;

	bool warn_about_power_off = false; // printout needs to be AFTER settings()
	if ((0==strcmp(format_out,"raw") || 0==strcmp(format_out,"wave"))
		&& fb_power) {
		fb_power = false; // otherwise in would produce power spectrum
		warn_about_power_off = true;
	}

	// print configuration information
	if (info) settings();

	if (warn_about_power_off && verbose) 
		cerr << "OPTS: Output = signal -> -fb_power forced to off!" << endl;

	return 0;
}

void opts::settings() {

    string _IOendians [] = { "big", "little" };
    string _NRwhen[] = { "before filter bank", "after filter bank" };
	string _bool[] = { "no", "yes" };

	cerr << "\nCtuCopy configuration:\n"
		 << "\tprogram version                      = " << version << endl
		 << "\tconfiguration file                   = ";
	if (0!=strcmp(config,"")) cerr << config << endl; 
	else cerr << "-not used-" << endl;
	cerr << "- IO:" << endl
	     << "\tinput format                         = " << format_in << endl
	     << "\toutput format                        = " << format_out << endl;
	if (0==strcmp(format_out,"pfile"))
		cerr
		 << "\tpfile name                           = " << pfilename << endl;
	if (0==strcmp(format_out,"ark"))
		cerr
		 << "\tarkfile name                           = " << arkfilename << endl;
            cerr << "\tinput byte order                     = " << _IOendians[endian_in] 
		 << " endian" << endl
		 << "\toutput byte order                    = " << _IOendians[endian_out] 
		 << " endian" << endl
		 << "\t(natural byte order                  = " << _IOendians[natural_little] 
		 << " endian)" << endl
		 << "\tonline input                         = " << _bool[pipe_in] << endl
		 << "\tonline output                        = " << _bool[pipe_out] << endl
		 << "\tlist of files                        = ";
	if (0!=strcmp(list,"")) cerr << list << endl; 
	else cerr << "-not used-" << endl;
	cerr << "\tsingle input file                    = ";
	if (0!=strcmp(in,"")) cerr << in << endl; 
	else cerr << "-not used-" << endl;
	cerr << "\tsingle output file                   = ";
	if (0!=strcmp(out,"")) cerr << out << endl; 
	else cerr << "-not used-" << endl;
	cerr << "\tpreemphasis                          = " << preem << endl
		 << "\tsampling rate                        = " << fs << endl
		 << "\tdithering constant                   = " << dither << endl
		 << "\tremove DC                            = " << _bool[remove_dc] << endl
		 << "\tremove DC_1 v. 3.2.2                  = " << _bool[remove_dc1] << endl
		 
		 << "- Segmentation:" << endl
		 << "\twindow size                          = " << window_ms 
		 << " ms (" << window << " samples)" << endl
		 << "\twindow shift (samples)               = " << wshift_ms
		 << " ms (" << wshift << " samples)" << endl
		 << "\t(FFT size                            = " << wfft << ")"<< endl
		 
		 << "- Filter Bank:" << endl
		 << "\tfrequency scale                      = " << fb_scale << endl
		 << "\tfilter shape                         = " << fb_shape << endl
		 << "\tnormalize filter area                = " << _bool[fb_norm] << endl
		 << "\tuse power spectrum                   = " << _bool[fb_power] << endl
		 << "\tapply Equal Loudness                 = " << _bool[fb_eqld] << endl
		 << "\tapply Intensity Loudness Power Law   = " << _bool[fb_inld] << endl
		 << "\tfilter bank spec.                    = " << fb_definition << endl
		 
		 << "- Noise Reduction:" << endl
		 << "\tVAD mode                             = " << vadmode << endl
		 << "\tfile with VAD                        = ";
	if (0!=strcmp(filevad,"")) cerr << filevad << endl; 
	else cerr << "-not used-" << endl;
	cerr << "\tnoise reduction                      = " << nr_mode << endl
		 << "\tintegrator const. p                  = " << nr_p << endl	 
		 << "\tVAD threshold q                      = " << nr_q << endl	 
		 << "\tspectrum power const. a              = " << nr_a << endl	 
		 << "\tspectrum oversubtraction const. b    = " << nr_b << endl	 
		 << "\tinitial speechless segments          = " << nr_initsegs << endl;
//		 << "\tperform RASTA                        = " << _bool[rasta] << endl;
	if (rasta) cerr 
		 << "\tRASTA specification file             = " << nr_rasta << endl;
	cerr << "\twhen to do Noise Reduction           = " << _NRwhen[nr_when] << endl
	
		 << "- Feature Extraction:" << endl
		 << "\tfeature kind                         = " << fea_kind << endl
		 << "\tLP order                             = " << fea_lporder << endl
		 << "\t# of cepstral coeffs                 = " << fea_ncepcoefs << endl
		 << "\tinclude 0th cepstral coeff           = " << _bool[fea_c0] << endl
		 << "\tinclude Energy                       = " << _bool[fea_E] << endl
		 << "\tcompute raw energy                   = " << _bool[fea_rawenergy] << endl
		 << "\tliftering const.                     = " << fea_lifter << endl
		 << "\tpreset                               = " << preset << endl
                 << "\tCMS fea_Z_exp   (off = -1)           = " << fea_Z_exp << endl
                 << "\tCMS fea_Z_block (off = -1)           = " << fea_Z_block << endl
                 << "\tCMVN stat_cmvn                       = " << _bool[stat_cmvn] << endl
                 << "\tCMVN apply_cmvn                      = " << _bool[apply_cmvn] << endl
		 << "- Voice Activity Detection:" << endl
		 << "\tVAD apply mode                       = " << vad_apply_mode << endl
		 << "\tVAD out mode                         = " << vad_out_mode << endl
		 << "\tVAD output file (single-file mode)   = " << vad_out << endl
		 << "\t-- VAD criterial function            = " << vad_cri_mode << endl
		 << "\tnormalize energy in decibels         = " << _bool[vad_energy_db] << endl
                 << "\tcepstral distance mode               = " << vad_cepdist_mode << endl
		 << "\tsmoothness of cepstral background estimation p = " << vad_cepdist_p << endl
		 << "\tinit voiceless segments for cepstral background estimation = " << vad_cepdist_init << endl
		 << "\tnumber of LPC cepstral coefficients  = " << vad_lpc_coefs << endl
		 << "\t-- VAD thresholding method           = " << vad_thr_mode << endl
		 << "\tabsolute threshold placement         = " << vad_absolute_thr << endl
		 << "\tinit segments for perc. threshold    = " << vad_perc_init << endl
		 << "\tpercentual threshold placement       = " << vad_perc_thr << endl
		 << "\tadaptive threshold init (voiceless) segments = " << vad_adapt_init << endl
		 << "\tadaptive threshold smoothness q      = " << vad_adapt_q << endl
		 << "\tspeech-to-background parameter z_alpha = " << vad_adapt_za << endl
		 << "\tdynamic threshold init (voiceless) segments = " << vad_dyn_init << endl
		 << "\tthreshold placement in dynamic range   = " << vad_dyn_perc << endl
		 << "\tdynamic range minimum                    = " << vad_dyn_min << endl
		 << "\tsmoothness of increasing dynamic maximum = " << vad_dyn_qmaxinc << endl
		 << "\tsmoothness of decreasing dynamic maximum = " << vad_dyn_qmaxdec << endl
		 << "\tsmoothness of decreasing dynamic minimum = " << vad_dyn_qmindec << endl
		 << "\tsmoothness of increasing dynamic minimum = " << vad_dyn_qmininc << endl
                 << "\tvad output median filter order           = " << vad_filter_order << endl
		 << "- Miscellaneous:" << endl
		 << "\tverbose                              = " << _bool[verbose] << endl
		 << "\tquiet mode                           = " << _bool[quiet] << endl
		 << "\tprint info                           = " << _bool[info] << endl << endl;
};

void opts::usage() {	
	 
 	cerr << endl
 		 << "CtuCopy - universal feature extractor and speech enhancer." << endl
		 << "Version " << version << ", (C) FEE CTU Prague 2003-2014." << endl;

	cerr << "\n  Usage: ctucopy <options> \n" << endl

	     << "  Optional Parameters                                                         <defaults>" << endl
	     << "  ......................................................................................" << endl
	     << "  INPUT/OUTPUT:                                                                         " << endl
	     << "    -format_in <fmt>        - input file format                                      <->" << endl
		 << "                fmt    = raw     - raw file (sequence of 16b integers)                  " << endl
		 << "                         alaw    - raw file, A-law encoded (sequence of 8b integers)    " << endl
		 << "                         mulaw   - raw file, mu-law encoded (sequence of 8b integers)   " << endl
		 << "                         wave    - MS wave file, PCM 16b, mono                          " << endl
                 << "                         htk     - HTK file                                             " << endl
		 << "    -nfeacoefs  <int>            - number of input coeffs                           <13>" << endl
		 << "    -filters    <filename>       - single input file with coefs of filters.        <none>"<< endl
     << "    -format_out <fmt>       - output file format                                   <htk>" << endl
		 << "                fmt    = raw     - raw file (sequence of 16b integers)                  " << endl
		 << "                         wave    - MS wave file, PCM 16b, mono                          " << endl
		 << "                         htk     - HTK file                                             " << endl
		 << "                         pfile=<filename>   - ICSI pfile format                         " << endl
		 << "                         ark=<filename>     - KALDI ark,scp format                      " << endl
		 << "    -endian_in <little|big> - input byte order (not for MS wave)                <little>" << endl
		 << "    -endian_out <little|big>- output byte order (consider big for HTK)          <little>" << endl
		 << "    -online_in              - input from STDIN                                 <offline>" << endl
		 << "    -online_out             - output to STDOUT (raw or htk only)               <offline>" << endl
		 << "    -S <filename>           - list of files to be processed: source and target    <none>" << endl
         << "                              pair per line, separated by tab(s) or space(s).           " << endl
         << "                              When format_out = pfile, then target is discarded.        " << endl
         << "    -i <filename>           - single input file (overrides -S)                    <none>" << endl
         << "    -o <filename>           - single output file (overrides -S)                   <none>" << endl
         << "    -preem <float>          - preemphasis in range <0.0-1.0), 0 = no preemphasis   <0.0>" << endl
         << "    -fs <int>               - sampling frequency in Hz                                  " << endl
		 << "    -dither <float>         - input dither constant, 0 = no dither                 <1.0>" << endl
		 << "    -remove_dc <on|off>     - remove DC from frame                                  <on>" << endl
		 << "    -remove_dc1 <on|off>    - remove DC from frame                                  <off>" << endl
		 << "                                                                                        " << endl
         << "  SEGMENTATION:                                                                         " << endl
		 << "    -w <float>              - window size (ms)                                      <25>" << endl
		 << "    -s <float>              - window shift size (ms)                                <10>" << endl
		 << "                                                                                        " << endl
         << "  FILTER BANK:                                                                          " << endl
		 << "    -fb_scale <scale>       - type of frequency scale warping                      <mel>" << endl
		 << "               scale   = mel     - melodic scale                                        " << endl
		 << "                         bark    - Bark scale                                           " << endl
		 << "                         lin     - linear scale                                         " << endl
		 << "                         expolog - Hansen's scale                                       " << endl
		 << "    -fb_shape <type>        - shape of filters                                  <triang>" << endl
		 << "               type    = triang  - triangular with 50% overlap                          " << endl
		 << "                         rect    - rectangular without overlap                          " << endl
		 << "                         trapez  - critical band-like shape (Bark scale only)           " << endl
	     << "               NOTE: for trapez case the option \"-fb_definition\" is discarded,        " << endl
	     << "                     and a fixed special filter bank defined for PLP analysis is used.  " << endl
         << "    -fb_norm <on|off>       - normalize filter to unit area                         <on>" << endl
         << "    -fb_power <on|off>      - use power spectrum as input to filter bank            <on>" << endl
	     << "    -fb_eqld <on|off>       - apply Equal Loudness curve                            <on>" << endl
	     << "    -fb_inld <on|off>       - apply Intensity Loudness Power Law                    <on>" << endl
         << "    -fb_definition <string> - definition of filter bank.                     <26filters>" << endl
         << "                              It is a combination of tokens, each of which specifying   " << endl
         << "                              one or more filter(s) to be added to the filter bank.     " << endl
         << "               string = token[,token]...                                                " << endl
         << "               token   = [[X-YHz:]K-L/]Nfilters                                         " << endl
	     << "                           X, Y ... optional frequency limits of the freq. scale        " << endl
	     << "                                    to be used for filter(s) design [floats, Hz]        " << endl
	     << "                                    X - Y defaults to 0 - fs/2 Hz                       " << endl
	     << "                           N    ... # of equidistant filters on fb_scale                " << endl
	     << "                           K, L ... optionally specifies a subset of the above filters  " << endl
	     << "               * see man page for more information and examples                         " << endl
		 << "                                                                                        " << endl
         << "  NOISE REDUCTION:                                                                      " << endl
	     << "    -nr_mode <mode>         - algorithm of additive noise suppression             <none>" << endl
		 << "               mode    = exten   - Extended Spectral Subtraction                        " << endl
		 << "                         hwss    - half-wave rectified spectral subtraction with VAD    " << endl
		 << "                         fwss    - full-wave rectified spectral subtraction with VAD    " << endl
		 << "                         2fwss   - 2-pass full-wave rectified s.s. with VAD             " << endl
		 << "                         none    - no noise reduction                                   " << endl
	     << "    -nr_p <double>          - integrator constant for Spectral Subtr.             <0.95>" << endl
	     << "    -nr_q <double>          - VAD threshold                                       <0.99>" << endl
	     << "    -nr_a <double>          - magnitude spectrum power constant for Spectr. Subtr. <1.0>" << endl
	     << "    -nr_b <double>          - noise oversubtraction factor (for modes *ss)         <1.0>" << endl
	     << "    -nr_initsegs            - # of initial frames without speech (for modes *ss)    <10>" << endl
	     << "    -vad <string>           - voice activity detector specification               <none>" << endl
		 << "               string  = burg    - Burg's built-in detector                             " << endl
		 << "                         file=<name>                                                    " << endl
		 << "                                 - VAD specified in a file (simple text file containing " << endl
		 << "                                   sequence of numbers 0 or 1 corresponding to frames)  " << endl
//		 << "    -rasta <filename>       - perform RASTA filtering with impulse responses         <->" << endl
//	     << "                              loaded from a file                                        " << endl
	     << "    -nr_when <beforeFB,afterFB>                                                         " << endl
	     << "                            - apply noise reduction before/after filter bank  <beforeFB>" << endl
		 << "                                                                                        " << endl
         << "  FEATURE EXTRACTION:                                                                   " << endl
	     << "    -fea_kind <kind>        - output feature kind                                  <lpc>" << endl
		 << "               kind    = spec    - magnitude spectrum                                   " << endl
		 << "                         logspec - log magnitude spectrum                               " << endl
		 << "                         lpa     - linear predictive coeffs                             " << endl
		 << "                         lpc     - LP cepstral coeffs (as in PLPC)                      " << endl
		 << "                         dctc    - DCT cepstral coeffs (as in MFCC)                     " << endl
		 << "                         trapdct,<int>,<int> - TRAP-DCT cepstral coeffs calculated      " << endl
		 << "                                   from log mel spectrum; first arg = <TRAP length>,    " << endl
                 << "                                   second arg = <# of first DCTs to output per band>    " << endl
		 << "                         td-iir-mfcc  - td-iir-mfcc features                            " << endl
	     << "    -fea_delta <kind>    - delta, acceleration and third differential coefficients   <d_a>" << endl
	     << "                     kind = d      - delta diff. coefficients                               " << endl
	     << "                     kind = d_a    - delta, acceleration diff. coefficients                 " << endl
	     << "                     kind = d_a_t  - delta , acc., and third diff. coefficients             " << endl
	     << "    -d_win    <int>         - Delta window size                                            <2>" << endl
	     << "    -a_win    <int>         - Acceleration diff. coeff. window size                        <2>" << endl
	     << "    -t_win    <int>         - third diff. coeff. window size                               <2>" << endl
	     << "    -fea_trap <int>         - Create context information from features defined -fea_kind.  <5>" << endl
             << "                            - (Length is number of frame and must be odd)               " <<endl
	     << "    -fea_lporder   <int>    - LPC order                                             <12>" << endl
	     << "    -fea_ncepcoefs <int>    - number of cepstral coeffs (without c0)                <12>" << endl
	     << "    -fea_c0 <on|off>        - add zeroth cepstral coeff to feature vector           <on>" << endl
	     << "    -fea_E  <on|off>        - add frame energy to feature vector                   <off>" << endl
	     << "    -fea_rawenergy <on|off> - compute energy directly from input signal frame      <off>" << endl
	     << "    -fea_lifter <int>       - cepstrum liftering constant (1 = off)                 <22>" << endl
	     << "    -fea_Z_exp <int>        - CMN, average cepstrum is computed as exponential moving average (off = -1)   <2000>" << endl
             << "                              Value of fea_Z_exp correspons to time constant (window size) across which average cepstrum is computed." << endl
	     << "    -fea_Z_block <int>      - CMN, average cepstrum is computed as moving average (off = -1)               <2000>" << endl
             << "                              Value of fea_Z_block correspons to time constant (window size) across which average cepstrum is computed." << endl
	     << "    -stat_cmvn  <filename>  - single output file with stat. of CMVN              <none>" << endl
	     << "    -apply_cmvn <filename>  - apply cmvn, single input file with stat. of CMVN   <none>" << endl
	     << "    -weight_of_td_iir_mfcc_bank <double>  - weight of td-iir-mfcc bank          <2.026>" << endl
             << endl
             << "  VOICE ACTIVITY DETECTION:                                                         " << endl
             << "    -vad_apply_mode <mode>  - how to apply VAD on output                            <none>" << endl
             << "               mode    = none     - do not apply VAD on output                      " << endl
             << "                         silence  - silence non-speech segments                     " << endl
             << "                         drop     - drop non-speech segments                        " << endl
             << endl
             << "    -vad_out_mode <mode>    - how to safe VAD output                                <none>" << endl
             << "               mode    = none     - do not save VAD output file                     " << endl
             << "                         vad      - only save VAD output file                       " << endl
             << "                         debug    - save VAD output file and debug files            " << endl
             << endl
             << "    -vad_out  <filename>    - VAD output file for single-file ctucopy mode          <none>"<< endl
             << endl
             << "    -vad_cri_mode <mode>          - VAD criterion mode                              <energy>" << endl
             << "               mode    = energy   - energetic detector                              " << endl
             << "                         cepdist  - cepstral distance                               " << endl
             << endl
             << "    -vad_energy_db <on|off> - normalize energy in decibels                          <on>" << endl
             << endl
             << "    -vad_cepdist_mode <mode> - cepstral distance from...                            <lpc>" << endl
             << "               mode    = lpc - VAD internal LPC cepstrum estimator                  " << endl
             << "                         fea - features computed inside fea module                  " << endl
             << "                         in  - input features (if input is HTK feature file)        " << endl
             << "    -vad_cepdist_p      <double> - smoothness of cepstral background estimation     <0.8>" << endl
             << "    -vad_cepdist_init   <int>    - init segments for cepstral background            <4>" << endl
             << "    -vad_lpc_coefs <int>     - number of LPC cepstral coefficients                  <14>" << endl
             << endl
             << "    -vad_thr_mode <mode>    - VAD threshold-setting method                          <perc>" << endl
             << "                mode   = absolute - set threshold by absolute value                 " << endl
             << "                         perc     - percentual threshold placed in criterial value range" << endl
             << "                         adapt    - adaptive threshold (Harrison detector)          " << endl
             << "                         dyn      - adaptive threshold based on criterial function dynamic " << endl
             << endl
             << "    -vad_absolute_thr <double>  - absolute threshold placement                      <TODO>" << endl
             << "    -vad_perc_init <int>        - init segments for perc. threshold                 <10.0>" << endl
             << "    -vad_perc_thr  <double>     - percentual threshold placement                    <50.0>" << endl
             << endl
             << "    -vad_adapt_init  <int>  - voiceless init segments for adaptive threshold        <20>" << endl
             << "    -vad_adapt_q  <double>  - adaptive threshold smoothness q                       <0.9>" << endl
             << "    -vad_adapt_za <double>  - speech-to-background parameter z_alpha                <2.0>" << endl
             << endl
             << "    -vad_dyn_init    <int>    - voiceless init segments for dynamic threshold       <5>" << endl
             << "    -vad_dyn_perc    <double> - threshold placement in dynamic range (percentual)   <50.0>" << endl
             << "    -vad_dyn_min     <double> - dynamic range minimum                               <TODO>" << endl
             << "    -vad_dyn_qmaxinc <double> - smoothness of increasing dynamic maximum            <0.8>" << endl
             << "    -vad_dyn_qmaxdec <double> - smoothness of decreasing dynamic maximum            <0.995>" << endl
             << "    -vad_dyn_qmindec <double> - smoothness of decreasing dynamic minimum            <0.8>" << endl
             << "    -vad_dyn_qmininc <double> - smoothness of increasing dynamic minimum            <0.9999>" << endl
             << endl
             << "    -vad_filter_order <int>   - vad output filter order                             <3>" << endl
             << endl
         << "  PRESETS:                                                                              " << endl
	     << "    -preset <type>          - apply a preset to the above options. Use -v option  <user>" << endl
	     << "                              to list exact settings when used. See man page for        " << endl
	     << "                              full macro definitions.                                   " << endl
		 << "               preset  = mfcc    - MFC cepstral coeffs                                  " << endl
		 << "                         plpc    - PLP cepstral coeffs                                  " << endl
		 << "                         exten   - Speech enhancement using Extended Spectral Subtr.    " << endl
		 << "                                                                                        " << endl
         << "  MISCELLANEOUS:                                                                        " << endl
	     << "    -verbose or -v          - verbose mode, prints all warnings and info           <off>" << endl 
	     << "    -quiet                  - suppress console output                              <off>" << endl
	     << "    -info                   - print CtuCopy settings                               <off>" << endl
         << "    -C <filename>           - load configuration file                             <none>" << endl
         << "                      SYNTAX: same as command line, one option per line                 " << endl  
         << "                              comments beginning with #, whitespace allowed.            " << endl  
	     << "    -h or --help            - print this help and exit" << endl;
};


int	opts::parse(char*l,char*r) {

	if (0==strcmp(l,"-S"))    { if(r) strcpy(list,r); }
	else if (!strcmp(l,"-i")) { if(r) strcpy(in,r); }
	else if (!strcmp(l,"-o")) { if(r) strcpy(out,r); }
	else if (!strcmp(l,"-format_in")) 	{ if(r) strcpy(format_in,r); }
	else if (!strcmp(l,"-format_out") && r)	{
		if (strstr(r,"pfile=") != NULL) {
			char*next = strchr(r,'=');
			r[next-r]=0;
			strcpy(pfilename,next+1);
			strcpy(format_out,"pfile");
		}else if (strstr(r,"ark=") != NULL) {
			char*next = strchr(r,'=');
			r[next-r]=0;
			strcpy(arkfilename,next+1);
			strcpy(format_out,"ark");
		}
		else strcpy(format_out,r);
	}
	else if (!strcmp(l,"-endian_in") && r) { 
		if (!strcmp(r,"big")) endian_in = big; 
		else if (!strcmp(r,"little")) endian_in = little;
	}
	else if (!strcmp(l,"-endian_out") && r) { 
		if (!strcmp(r,"big")) endian_out = big; 
		else if (!strcmp(r,"little")) endian_out = little;
	}	
	else if (!strcmp(l,"-online_in"))     { pipe_in = true; }
	else if (!strcmp(l,"-online_out"))     { pipe_out = true; }
	else if (!strcmp(l,"-fb_printself"))     { fb_printself = true; }
	else if (!strcmp(l,"-preem"))		{ if(r) preem = (float) atof(r); }
	else if (!strcmp(l,"-fea_Z_exp"))		{ if(r) fea_Z_exp = (float) atof(r); }
	else if (!strcmp(l,"-fea_Z_block"))		{ if(r) fea_Z_block = (float) atof(r); }
	else if (!strcmp(l,"-stat_cmvn") && r){
               strcpy(fcmvn_stat_out,r);
               stat_cmvn=true;
        }
	else if (!strcmp(l,"-apply_cmvn") && r){
               strcpy(fcmvn_stat_in,r);
               apply_cmvn=true;
        }
	else if (!strcmp(l,"-fea_delta") && r){
                fea_delta=true;
                fea_trap=false;
		if      (!strcmp(r,"d"))     n_order=1;
		else if (!strcmp(r,"d_a"))   n_order=2;
		else if (!strcmp(r,"d_a_t")) n_order=3;
                else fea_delta=false;
        }
	else if (!strcmp(l,"-fea_trap") && r){
             if(!fea_delta){
                fea_trap=true;
                if(r){
                     trap_win = atoi(r);
                     fea_delta=true;
                     n_order=1;
                     d_win=(trap_win-1)/2;
                }
             }
        }
	else if (!strcmp(l,"-filters") && r){
               strcpy(ffilters,r);
        }
	else if (!strcmp(l,"-fs")) 			{ if(r) fs = atoi(r); }
	else if (!strcmp(l,"-dither"))  	{ if(r) dither = atof(r); }
	else if (!strcmp(l,"-remove_dc") && r) {
		if (!strcmp(r,"on")) remove_dc = true;
		else if (!strcmp(r,"off")) remove_dc = false;
	}
	else if (!strcmp(l,"-remove_dc1") && r) {
		if (!strcmp(r,"on")) remove_dc1 = true;
		else if (!strcmp(r,"off")) remove_dc1 = false;
	}
	else if (!strcmp(l,"-w")) 			{ if(r) window_ms = atof(r); }
	else if (!strcmp(l,"-s")) 			{ if(r) wshift_ms = atof(r); }
	
	else if (!strcmp(l,"-fb_scale")) 	{ if(r) strcpy(fb_scale,r); }
	else if (!strcmp(l,"-fb_shape")) 	{ if(r) strcpy(fb_shape,r); }
	else if (!strcmp(l,"-fb_norm") && r) {
		if (!strcmp(r,"on")) fb_norm = true;
		else if (!strcmp(r,"off")) fb_norm = false;
	}
	else if (!strcmp(l,"-fb_power") && r) {
		if (!strcmp(r,"on")) fb_power = true;
		else if (!strcmp(r,"off")) fb_power = false;
	}
	else if (!strcmp(l,"-fb_eqld") && r) { 
		if (!strcmp(r,"on")) fb_eqld = true; 
		else if (!strcmp(r,"off")) fb_eqld = false;
	}
	else if (!strcmp(l,"-fb_inld") && r) { 
		if (!strcmp(r,"on")) fb_inld = true; 
		else if (!strcmp(r,"off")) fb_inld = false;
	}
	else if (!strcmp(l,"-fb_definition")) { if(r) strcpy(fb_definition,r); }
	else if (!strcmp(l,"-vad") && r) {
		if (!strcmp(r, "burg")) strcpy(vadmode,"burg");
		else if (strstr(r,"file=") != NULL) {
			char*next = strchr(r,'='); 
			r[next-r]=0;
			strcpy(filevad,next+1); 
			strcpy(vadmode,"file");
		}
		else throw "OPTS: Syntax error in option -vad !";
	}
	else if (!strcmp(l,"-nr_mode")) 	{ if(r) strcpy(nr_mode,r); }
	else if (!strcmp(l,"-nr_p")) 		{ if(r) nr_p = atof(r); }
	else if (!strcmp(l,"-nr_q")) 		{ if(r) nr_q = atof(r); }
	else if (!strcmp(l,"-nr_a")) 		{ if(r) nr_a = atof(r); }
	else if (!strcmp(l,"-nr_b")) 		{ if(r) nr_b = atof(r); }
	else if (!strcmp(l,"-nr_initsegs")) { if(r) nr_initsegs = atoi(r); }
	else if (!strcmp(l,"-nr_rasta")) 	{
		rasta = true;
		strcpy(nr_rasta,r);
	}
	else if (!strcmp(l,"-nr_when") && r) { 
		if (!strcmp(r,"beforeFB")) nr_when = beforeFB; 
		else if (!strcmp(r,"afterFB")) nr_when = afterFB;
	}
	else if (!strcmp(l,"-fea_kind")){
           if(!r)
              throw "OPTS: Missing argument to '-fea_kind' option!";
           if(strstr(r,"trapdct") != NULL){
              char*next = strchr(r,',');
              if (next==NULL)
                throw "OPTS: Syntax error in option -fea_kind! (should be -fea_kind trapdct,<X>,<Y>)";
              r[next-r]=0;
              strcpy(fea_kind,r);
              fea_trapdct_traplen = atoi(next+1);
              next = strchr(next+1,',');
              if(next==NULL)
                 throw "OPTS: Syntax error in option -fea_kind! (should be -fea_kind trapdct,<X>,<Y>)";
              fea_trapdct_ndct = atoi(next+1);
           }else
              strcpy(fea_kind,r);
        }
	else if (!strcmp(l,"-d_win"))     { if(r) d_win       = atoi(r); }
	else if (!strcmp(l,"-a_win"))     { if(r) a_win       = atoi(r); }
	else if (!strcmp(l,"-t_win"))     { if(r) t_win       = atoi(r); }
	else if (!strcmp(l,"-fea_lporder"))   { if(r) fea_lporder   = atoi(r); }
	else if (!strcmp(l,"-fea_ncepcoefs")) { if(r) fea_ncepcoefs = atoi(r); }
	else if (!strcmp(l,"-nfeacoefs")) { if(r) nfeacoefs = atoi(r); }
	else if (!strcmp(l,"-fea_c0") && r) {
               if (!strcmp(r,"on"))       fea_c0 = true;
               else if (!strcmp(r,"off")) fea_c0 = false;
        }
	else if (!strcmp(l,"-fea_E") && r) { 
		if (!strcmp(r,"on")) fea_E = true; 
		else if (!strcmp(r,"off")) fea_E = false;
	}
	else if (!strcmp(l,"-fea_rawenergy") && r) { 
		if (!strcmp(r,"on")) fea_rawenergy = true; 
		else if (!strcmp(r,"off")) fea_rawenergy = false;
	}
	else if (!strcmp(l,"-weight_of_td_iir_mfcc_bank"))  { if(r) weight_of_td_iir_mfcc_bank = (float) atof(r);}
	else if (!strcmp(l,"-fea_lifter"))  { if(r) fea_lifter = atoi(r); }

// VAD -------------------------------------------------------------------------
	else if (!strcmp(l,"-vad_apply_mode"))  { if(r) strcpy(vad_apply_mode,r); }
	else if (!strcmp(l,"-vad_out_mode"))    { if(r) strcpy(vad_out_mode,r); }
	else if (!strcmp(l,"-vad_out"))         { if(r) strcpy(vad_out,r); }
	else if (!strcmp(l,"-vad_cri_mode"))    { if(r) strcpy(vad_cri_mode,r); }
	else if (!strcmp(l,"-vad_thr_mode"))    { if(r) strcpy(vad_thr_mode,r); }
	else if (!strcmp(l,"-vad_energy_db") && r) {
		if (!strcmp(r,"on")) vad_energy_db = true;
		else if (!strcmp(r,"off")) vad_energy_db = false;
	}
        else if (!strcmp(l,"-vad_cepdist_mode")) { if(r) strcpy(vad_cepdist_mode,r); }
	else if (!strcmp(l,"-vad_cepdist_p"))   { if(r) vad_cepdist_p = atof(r); }
	else if (!strcmp(l,"-vad_cepdist_init")) { if(r) vad_cepdist_init = atoi(r); }
	else if (!strcmp(l,"-vad_lpc_coefs"))   { if(r) vad_lpc_coefs = atoi(r); }
	else if (!strcmp(l,"-vad_absolute_thr"))  { if(r) vad_absolute_thr = atof(r); }
	else if (!strcmp(l,"-vad_perc_init"))  { if(r) vad_perc_init = atoi(r); }
	else if (!strcmp(l,"-vad_perc_thr"))  { if(r) vad_perc_thr = atof(r); }
	else if (!strcmp(l,"-vad_adapt_init"))  { if(r) vad_adapt_init = atoi(r); }
	else if (!strcmp(l,"-vad_adapt_q"))     { if(r) vad_adapt_q = atof(r); }
	else if (!strcmp(l,"-vad_adapt_za"))    { if(r) vad_adapt_za = atof(r); }
	else if (!strcmp(l,"-vad_dyn_init"))    { if(r) vad_dyn_init = atoi(r); }
	else if (!strcmp(l,"-vad_dyn_perc"))    { if(r) vad_dyn_perc = atof(r); }
	else if (!strcmp(l,"-vad_dyn_min"))     { if(r) vad_dyn_min = atof(r); }
	else if (!strcmp(l,"-vad_dyn_qmaxinc")) { if(r) vad_dyn_qmaxinc = atof(r); }
	else if (!strcmp(l,"-vad_dyn_qmaxdec")) { if(r) vad_dyn_qmaxdec = atof(r); }
	else if (!strcmp(l,"-vad_dyn_qmindec")) { if(r) vad_dyn_qmindec = atof(r); }
	else if (!strcmp(l,"-vad_dyn_qmininc")) { if(r) vad_dyn_qmininc = atof(r); }
	else if (!strcmp(l,"-vad_filter_order")){ if(r) vad_filter_order = atoi(r); }
// VAD END ---------------------------------------------------------------------

	else if (!strcmp(l,"-preset")) 		{ if(r) {strcpy(preset,r);set_preset();}}
	else if (!strcmp(l,"-verbose"))     { verbose = true; quiet = false; info = true; }
	else if (!strcmp(l,"-v"))   		{ verbose = true; quiet = false; info = true; }
	else if (!strcmp(l,"-quiet")) 		{ quiet = true; verbose = false; info = false; }
	else if (!strcmp(l,"-info")) 		{ info = true; quiet = false; }
	else if (!strcmp(l,"-C")) 			{ if(r) strcpy(config,r); }
	else if (!strcmp(l,"-h")) 			{ usage(); exit(0); }
	else if (!strcmp(l,"--help")) 		{ usage(); exit(0); }
	else if (r) {
		cerr << "OPTS: Syntax error in option \"" <<l<<" "<<r<<"\"."<<endl; throw "";
	} else {
		cerr << "OPTS: Syntax error in option \"" <<l<<"\"."<<endl; throw "";
	}
	return 0;
}


opts::~opts() {
	//	cerr << "opts DELETED.." << endl;
};
