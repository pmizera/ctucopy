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

#include "stdafx.h"
#include "in.h"
#include "amulaw.h"
#include <fstream>
//#include "base/types.h"

//////////////////////////////////////////////////////////////
//////////// DECLARATION OF SPECIALIZED CLASSES //////////////
//////////////////////////////////////////////////////////////

// ---------------- RAW (linear 16b PCM) -----------------
class rawIN : public _IN {
	// specialization for raw signal as input format

	protected:
		double * cbuffer; 	       // circ. buffer for input samples [window]
		Sint * cache;   	       // cache for input reading
		Mat<double> * _fbuffer;	       // circ. buffer for input filtered samples [window]

		int start; 		       // pointer to cbuffer start
		int ostart; 		       // pointer to obuffer start
		double * fft_in;
		double * fft_out;
                ifstream f_filters;            // load coefs. of filters
                Mat<double> *_coefs;           // all coefs. of filts [number of filters num. of coenfs]
                Mat<double> *_state;           // state of filter for sec. canon. structure  IIR
		double * Ebuffer; 	       // buffer with energy inband
                double * wdct;
                double normcoef;

        fftw_plan plan;
        FILE *fin;
        double preemtmp;    	// last sample from previous window to compute the right preemph

        virtual bool loadframe(int,int); 	// position in cbuffer, how many samples to load
        virtual bool compute_td_iir_mfcc(int,int);

	public:
		rawIN (opts *);
		virtual ~rawIN ();
		virtual bool new_file(char*);
		virtual bool loadf_filters(char*);
		virtual bool get_frame();
};

// ---------------- A-law -----------------
class alawIN : public rawIN {
	// specialization for A-law signal as input format
    protected:
        virtual bool loadframe(int,int); 	// position in cbuffer, how many samples to load
		char * lcache;   	        		// cache for input reading (char! <- Alaw)
	public:
		alawIN (opts *);
		virtual ~alawIN ();
		virtual bool new_file(char*);
		virtual bool get_frame();
};

// ---------------- mu-law -----------------
class mulawIN : public alawIN {
	// specialization for A-law signal as input format
    protected:
        virtual bool loadframe(int,int); 	// position in cbuffer, how many samples to load
	public:
		mulawIN (opts *o) : alawIN(o) {};
		virtual ~mulawIN () {};
		virtual bool new_file(char*);
		virtual bool get_frame();
};

// ---------------- MS WAVE file -----------------
class waveIN : public rawIN {
	// specialization for MS Wave signal as input format
    protected:
        uLint samples_left;
        virtual bool loadframe(int,int); 	// position in cbuffer, how many samples to load
		char * lcache;   	           		// cache for input reading (char!)
	public:
		waveIN (opts* o) : rawIN(o) {};
		virtual ~waveIN () {};
		virtual bool new_file(char*);
		virtual bool get_frame();
};

// ---------------- HTK file -----------------
class htkIN : public  rawIN {

       protected:
                 uSint fea_size_Bytes;             // of bytes per feature vector
                 int fea_size;                     // feature vector size
                 float * cache;                  // input buffer
          public:
                 htkIN (opts *);
                 virtual ~htkIN ();
                 virtual bool new_file(char*);
		 virtual bool get_frame();
};

// ---------------- CMVN file ------------------------------
class cmvnIN : public rawIN {

        protected:
                int nfea;  // input vector size
                int num_spk;
                int len_spk;
                char **list_ID_spk;
                ifstream fin_cmvn;
        public:
                cmvnIN(opts *);
                virtual ~cmvnIN ();
                virtual bool new_file(char *);
                virtual int count_speakers(char *);
//                virtual void save_frame(int num_spk);
};



//////////////////////////////////////////////////////////////
//////////// DEFINITIONS OF ALL CLASSES //////////////////////
//////////////////////////////////////////////////////////////

// -----   implementation of methods common to all classes   --------
void _IN::hamming (double alpha) {
	// init Hamming window of size o->window in W[i]
	double pi = 2.*asin(1.) ;
	for (int j = 0 ; j < o->window ; j++ )
	 W[j] = alpha - (1-alpha) * cos(2*pi*j/(o->window-1.)) ;
};

void _IN::ByteSwap16 (Sint *nValue) {
    uSint *x = (uSint*) nValue;
    *nValue = (uSint)(((*x >> 8)) | (*x << 8));
}

void _IN::ByteSwap32 (uLint *x) {
        *x = (((*x&0x000000FF)<<24)+((*x&0x0000FF00)<<8)+((*x&0x00FF0000)>>8)+((*x&0xFF000000)>>24));
};

void _IN::ByteSwap16u (uSint* x) {
        *x = (((*x) >> 8) | (*x << 8));
};

void _IN::SwapFloat (float* x) {
        uLint *f = (uLint*) x;
        *f = (((*f&0x000000FF)<<24)+((*f&0x0000FF00)<<8)+((*f&0x00FF0000)>>8)+((*f&0xFF000000)>>24));
};

// prototype for signal-like input classes (raw, a-law,...)
_IN::_IN(opts* o_in) {
    o = o_in;
    try {
      // grab space
      if(0==strcmp(o->fea_kind,"td-iir-mfcc")){  // do td-iir-mfcc
         _Xsabs = new Vec<double>(o->fea_ncepcoefs+1);
      }else{
        _Xsabs = new Vec<double>(o->wfftby2);			 // Xabs
      }

      W = new double[o->window]; // Hamming win (will be allocated further)
      if (o->phase_needed) _Xsph  = new Vec<double>(o->wfftby2); // Xph
    }
    catch (...) { throw "IO: Not enough memory!";};
}

_IN::~_IN() {
   	if (o->phase_needed) delete _Xsph;
	delete [] W;
	delete _Xsabs;
}

double _IN::c_ph (double re, double im) {
	// compute phase from Re,Im representation
	static const double hpi = 1.57079632679490;
	static const double pi  = 3.14159265358979;
	if (re == 0.0) {
		if (im > 0.0) return hpi;
		return -hpi;
	} else {
		double y = atan (im/re);
		if (re < 0.0 && im >= 0.0 )  y += pi;
		if (re < 0.0 && im  < 0.0 )  y -= pi;
		return y;
	}
}


// ******** INPUT = RAW ***************************************
rawIN::rawIN(opts* o) : _IN(o) {
	srand(1); // for dither
	hamming(0.54);
    try {
    	cbuffer = new double [o->window];
    	cache   = new Sint [o->window];

    	_fbuffer = new Mat<double>(24,o->window);
        _coefs   = new Mat<double>(24,10);
        _state   = new Mat<double>(24,5);
    	Ebuffer  = new double [24];
        wdct = new double [4*24];
    }
    catch (...) {
    	throw ("IN: Not enough memory!");
    }
    start = 0;	// pointer to circullar buffer start
    ostart = 0;	// pointer to circullar buffer start

    fft_in  = (double*) fftw_malloc(sizeof(double) * o->wfft);
    fft_out = (double*) fftw_malloc(sizeof(double) * o->wfft);

    if (fft_in == 0 || fft_out == 0) throw "IN: Not enough memory!";

    //     init FFT (fftw library)
    plan = fftw_plan_r2r_1d(o->wfft, fft_in, fft_out, FFTW_R2HC, FFTW_MEASURE);

    fin = NULL;

    // DCT coefficients init:
    // w[i] = cos (2*pi*i/4N)
    for (int i=0; i < 4*24; i++)
       wdct[i] = cos(3.14159265358979*(double)i/(2*24));

    normcoef = sqrt(2.0/24);
};


bool rawIN::loadf_filters(char * filename) {
   Mat<double> &coefs = *_coefs;

   f_filters.open(filename);
   if(!f_filters) return false;
   char line[1000];
   int numf = 0;
   int numc = 10;

   // loop over all filters
   while(f_filters.getline(line,1000)){
          coefs(numf,0) = atof(strtok(line,"\t"));
          for(int i=1; i<numc; i++){
             coefs(numf,i) = atof(strtok(NULL,"\t"));
          }
          numf++;
   }
   f_filters.close();

return true;
};

bool rawIN::new_file(char * filename) {
   // close previous file, open new file and load initial (window-wshift) samples
   if (fin) fclose(fin);

   if (o->pipe_in) fin = fdopen (0, "rb");		// STDIN
   else fin = fopen (filename,"rb");	                // filename

   if (!fin) throw "IN: Cannot open data file!";

   start = 0;
   preemtmp = 0.;

    // fill-in (window-wshift) samples in cbuffer
    if (!loadframe(start, o->window - o->wshift)) throw "IO: Signal shorter than one frame!";
    return true;
};

bool rawIN::compute_td_iir_mfcc(int pos, int num) {

    Mat<double> &coefs   = *_coefs;
    Mat<double> &state   = *_state;
    Mat<double> &fbuffer = *_fbuffer;

    for (int ff = 0; ff < 24; ff++){
       for (int i = 0; i < num; i++){
          state(ff,4) = coefs(ff,5) * cache[i];
          for (int a = 6; a < 10; a++){
             state(ff,4) -= coefs(ff,a) * state(ff,4-(a-5));
          }
          fbuffer(ff,(pos + i)%o->window) = coefs(ff,0) * state(ff,4);
          for (int b = 0; b < 4; b++){
             fbuffer(ff,(pos + i)%o->window) += coefs(ff,4-b) * state(ff,b);
          }
          for (int c = 0; c < 4; c++)
             state(ff,c)=state(ff,c+1);
          fbuffer(ff,(pos + i)%o->window) *= W[(pos + i)%o->window];
       }
    }
    return true;
};

bool rawIN::get_frame() {
        // Re+Im <-> Abs+Ang - yes, costly, but...
        // 0'th and wlen/2'th coeffs need to be processed separately (pure reals)
        static Vec<double> &Xabs = *_Xsabs;  // just a pointer
        static Vec<double> &Xph = *_Xsph;    // just a pointer
        Mat<double> &fbuffer = *_fbuffer;

	// load frame, fill in fft_in, possibly apply preemphasis/Hamming, run fft, compute E
	// load wshift samples at position = (start+(win-shift))%win
	if (!loadframe((start + o->window - o->wshift) % o->window, o->wshift)) return false;
//-------------------------------------------------------------------------------------------------------------------
// do TD-IIR-MFCC
        if(0==strcmp(o->fea_kind,"td-iir-mfcc")){
          for (int ff = 0; ff < 24; ff++){
           Ebuffer[23-ff] = 0;
           for (int i = 0; i < o->window; i++){
             Ebuffer[23-ff] += fbuffer(ff,i) * fbuffer(ff,i);
           }
           Ebuffer[23-ff] = log((o->window*o->window*Ebuffer[23-ff]/o->window)/o->weight_of_td_iir_mfcc_bank); // take logarithm of input fvec
          }

          // DCT projection
          for (int i = 0; i < 13; i++) {
             Xabs[i] = 0;
             for (int k = 1; k <= 24; k++) { // indexation from 1
                int id = (2*k - 1)*i % (4*24);
               Xabs[i] += Ebuffer[k - 1]*wdct[id];
                 // c[i] = sum_k=1^N [ Xk * w_(i*(2k-1)) ]
                 //      = sum..     [ Xk * cos(2*pi*i*(2k-1)/4N) ]
                 //      = sum..     [ Kx * cos(pi*i*(k-0.5)/N) ]
             }
             Xabs[i] *= normcoef;
          }
          start = (start + o->wshift) % o->window;//  - move ptr forward by wshift
          return true;
        }
//-------------------------------------------------------------------------------------------------------------------
	// DC offset removal from frame
	if (o->remove_dc1) {
		double offset = 0;
		for (int i = 0; i < o->window; i++)
			offset += cbuffer[i];
		offset /= o->window;
		for (int i = 0; i < o->window; i++)
			cbuffer[i] -= offset;
	}

	// compute raw energy
	E = -1.;
	if (o->fea_E && o->fea_rawenergy) {
		// compute it only when user requires raw energy
		E = 0;
		// E = log(sum(sig[x] ^ 2))
		for (int i = 1; i < o->window; i++)
			E += cbuffer[(start+i)%o->window] * cbuffer[(start+i)%o->window];
		E = log(E);
	}

	// fill in fft_in, apply preemphasis and Hamming
	if (o->preem>0.) {
		fft_in[0] = W[0] * (cbuffer[start] - o->preem * preemtmp); // transition sample from prev. frame
    	for (int i = 1; i < o->window; i++)
    		// carefully read circullar buffer, apply preemphasis
    		fft_in[i] = W[i] * (cbuffer[(start+i)%o->window] - o->preem*cbuffer[(start+i-1)%o->window]);
	} else
		for (int i = 0; i < o->window; i++)
			// without preemph
    		fft_in[i] = W[i] * cbuffer[(start+i)%o->window];

	// DC offset removal from frame
	if (o->remove_dc) {
		double offset = 0;
		for (int i = 0; i < o->window; i++)
			offset += fft_in[i];
		offset /= o->window;
		for (int i = 0; i < o->window; i++)
			fft_in[i] -= offset;
	}

	preemtmp = cbuffer[(start+o->wshift-1)%o->window]; // remember last sample for future
	start = (start + o->wshift) % o->window;//  - move ptr forward by wshift

    for (int i=o->window; i<o->wfft; i++) fft_in[i] = 0.; // pad with zeros before FFT
    fftw_execute(plan); // FFT

    if (o->remove_dc) Xabs[0] = 1e-10; // fixed floor
	else Xabs[0] = fft_out[0]*fft_out[0];
    Xabs[o->wfftby2-1] = fft_out[o->wfftby2-1]*fft_out[o->wfftby2-1];
 	for (int i=1; i<o->wfftby2-1; i++)
		// Re^2 + Im^2
		Xabs[i] = fft_out[i]*fft_out[i] + fft_out[o->wfft-i]*fft_out[o->wfft-i];
	if (o->phase_needed) {
		// this is costly so we only compute phase when it is needed
		Xph[0] = 0;
		if (fft_out[o->wfftby2-1] >= 0) Xph[o->wfftby2-1] = 0; else Xph[o->wfftby2-1] = 3.14159265358979;
		for (int i=1; i<o->wfftby2-1; i++)
			Xph[i] = c_ph(fft_out[i], fft_out[o->wfft-i]);
	}
	// now we are in power domain
	// compute energy after preemphasis, dither, Hamming
	if (o->fea_E && !o->fea_rawenergy) {
		// E = sum(X_i^2) over symmetr. spec == 2xsum over half spec + X_0^2 + X_n/2^2.
		// Instead of multiplying X(1..n/2-1) by two we divide X0 and Xn/2+1 by two
		// and after all multiply by two (little faster).
		E = Xabs[0]/2. + Xabs[o->wfftby2 - 1]/2.;
		for (int i=1; i<o->wfftby2 - 1; i++)
			E += Xabs[i];
		E = log(E*2.);
	}

	if (!o->fb_power)  // no one wants power, let's sqrt :-/
		for (int i=0; i<o->wfftby2; i++)
			Xabs[i] = sqrt(Xabs[i]);
	return true;
};

rawIN::~rawIN() {
   fftw_destroy_plan(plan);
   fftw_free(fft_in);
   fftw_free(fft_out);
   delete [] cache;
   delete [] cbuffer;
   delete [] Ebuffer;
   delete [] wdct;
   delete _fbuffer;
   delete _coefs;
   delete _state;
};

bool rawIN::loadframe(int pos,int num) {
   // loads <num> data into cbuffer at <pos>ition (via cache)

   // read raw data to cache
   if (num != (int) fread (cache, sizeof(Sint), num, fin )) return false;

   // swap byte order if needed
   if (o->swap_in) for (int i=0; i<num; i++) ByteSwap16(&cache[i]);

//---------------------------------------
//do td-iir-mfcc
       if(0==strcmp(o->fea_kind,"td-iir-mfcc")){
         compute_td_iir_mfcc(pos,num);
         return true;
       }
//---------------------------------------

   // add dither, move the data to cbuffer
   if (o->dither != 0.)
     for (int i=0; i<num; i++)  // dithered version
       cbuffer[(pos+i)%o->window] = cache[i] + (2.*(double)rand()/(double)RAND_MAX - 1.)*o->dither;
   else
     for (int i=0; i<num; i++)
       cbuffer[(pos+i)%o->window] = cache[i];

    return true;
};


// ******** INPUT = A-law ***************************************
alawIN::alawIN(opts* o) : rawIN(o) {
    try { lcache   = new char [o->window]; }
    catch (...) { throw ("IN: Not enough memory!"); }
};

bool alawIN::new_file(char * filename) {
	return rawIN::new_file(filename);
};

bool alawIN::get_frame() {
	return rawIN::get_frame();
}

alawIN::~alawIN() {
	delete [] lcache;
};

bool alawIN::loadframe(int pos,int num) {
	// see the same for rawIN

	if (num != (int) fread (lcache, sizeof(char), num, fin )) return false;

	// alaw to linear
	for (int i=0; i<num; i++)
		alaw2lin(lcache[i],cache[i],1); // trailing 1 means "alaw mode"

	// swap byte order if needed
// if (o->swap_in) for (int i=0; i<num; i++) ByteSwap16(&cache[i]);
// ...no point swapping bytes when samples are 1-byte long

//---------------------------------------
//do td-iir-mfcc
       if(0==strcmp(o->fea_kind,"td-iir-mfcc")){
         rawIN::compute_td_iir_mfcc(pos,num);
         return true;
       }
//---------------------------------------

    // add dither, move the data to cbuffer
    if (o->dither)
    	for (int i=0; i<num; i++)  // dithered version
        	cbuffer[(pos+i)%o->window] = cache[i] + (2.*(double)rand()/(double)RAND_MAX - 1.)*o->dither;
    else
        for (int i=0; i<num; i++) cbuffer[(pos+i)%o->window] = cache[i];

    return true;
};


// ******** INPUT = mu-law ***************************************

bool mulawIN::new_file(char * filename) {
	return alawIN::new_file(filename);
};

bool mulawIN::get_frame() {
	return alawIN::get_frame();
}

bool mulawIN::loadframe(int pos,int num) {
	// see the same for rawIN

	// read raw data
	if (num != (int) fread (lcache, sizeof(char), num, fin )) return false;

	// mulaw to linear
	for (int i=0; i<num; i++)
		alaw2lin(lcache[i],cache[i],0); // trailing 0 means "mulaw mode"

	// swap byte order if needed
//    if (o->swap_in) for (int i=0; i<num; i++) ByteSwap16(&cache[i]);
// ...no point swapping bytes when samples are 1-byte long

    // add dither, move the data to cbuffer
    if (o->dither)
    	for (int i=0; i<num; i++)  // dithered version
        	cbuffer[(pos+i)%o->window] = cache[i] + (2.*(double)rand()/(double)RAND_MAX - 1.)*o->dither;
    else
        for (int i=0; i<num; i++) cbuffer[(pos+i)%o->window] = cache[i];

    return true;
};


// ******** INPUT = MS Wave ***************************************

bool waveIN::new_file(char * filename) {
	// open file and load initial (window-wshift) samples

	if (fin) fclose (fin);
	if (o->pipe_in) fin = fdopen (0, "rb");	// STDIN
	else fin = fopen (filename,"rb");		// filename
	if (!fin) throw ("IN: Cannot open file!");

    char id[5] = ""; // four bytes+1 to hold text IDs
    uLint size; // 32 bit value to hold file size
    Sint compression, channels, block_align, bits_per_sample; // 16b values
    uLint this_chunk_size, sample_rate, avg_bytes_sec, data_size; // 32b bit values

	// RIFF chunk ("RIFFxxxxWAVE")
	fread(id, 1, 4, fin); // first four bytes have to be "RIFF"
    if (0!=strcmp(id,"RIFF")) throw ("IN: No RIFF header in file!");
    fread(&size, 4, 1, fin); // read in 32bit size value
    fread(id, 1, 4, fin); // these 4 bytes need to be "WAVE"
    if (0!=strcmp(id,"WAVE")) throw ("IN: Not a WAVE file!");

	// FORMAT chunk
    fread(id, 1, 4, fin); // these bytes should be "fmt ", not checking
    fread(&this_chunk_size, 4, 1, fin); // not checking as for PCM always 16
    fread(&compression, 2, 1, fin);
    fread(&channels, 2, 1, fin); // 1=mono, 2=stereo
    fread(&sample_rate, 4, 1, fin); // fs in Hz
    fread(&avg_bytes_sec, 4, 1, fin);
    fread(&block_align, 2, 1, fin);
    fread(&bits_per_sample, 2, 1, fin); // 8/16 bps
    fread(id, 1, 4, fin); // these bytes should be "data", not checking
    fread(&data_size, 4, 1, fin); // # of bytes of sound data in file

	// checks
	if (compression != 1) throw ("IN: Not a PCM WAVE file!");
	if (long(sample_rate) != o->fs)
		throw ("IN: WAVE file reports different sampling rate than specified!");

	if (channels != 1) throw ("IN: Input WAVE file is not mono!");
    if (bits_per_sample != 16) throw ("IN: Not 16 bits per sample!");

	samples_left = data_size/2;
    start = 0;
    preemtmp = 0;

	// fill-in (window-wshift) samples in cbuffer
	if (!loadframe(start, o->window - o->wshift)) throw "IO: Signal shorter than one frame!";
	return true;
};

bool waveIN::get_frame() {
	return rawIN::get_frame();
}

bool waveIN::loadframe(int pos,int num) {
	// see rawIN::loadframe for comments
	if (samples_left < (uLint)num) return false;
	samples_left -= num;

	// read raw data
	if (num != (int) fread (cache, sizeof(Sint), num, fin )) return false;

    // add dither, move the data to cbuffer
    if (o->dither)
    	for (int i=0; i<num; i++)  // dithered version
        	cbuffer[(pos+i)%o->window] = cache[i] + (2.*(double)rand()/(double)RAND_MAX - 1.)*o->dither;
    else
        for (int i=0; i<num; i++) cbuffer[(pos+i)%o->window] = cache[i];

    return true;
};


// ---------------- HTK file input ------------------------------
htkIN::htkIN(opts* o) :  rawIN(o) {
   fea_size = o->nfeacoefs;
   _fvec   = new Vec<double>(fea_size);

};

bool htkIN::new_file(char * filename) {
   // close previous file, open new file and load initial (window-wshift) samples
   if (fin)fclose(fin);

   if (o->pipe_in) fin = fdopen (0, "rb");		// STDIN
   else fin = fopen (filename,"rb");	                // filename

   if (!fin) throw "IN: Cannot open data file!";

   uLint frames;  // frames couter - for HTK header
   uLint period;
   uSint kind;
   uSint fea_size_BytesSWP; // number of bytes per feature ( 4-byte floats )

   // *** read HTK file header ***
   if (1 != fread(&frames, 4, 1, fin)) throw ("OUT: Error in stream writing!\n");
   if (1 != fread(&period, 4, 1, fin)) throw ("OUT: Error in stream writing!\n");
   if (1 != fread(&fea_size_BytesSWP, 2, 1, fin)) throw ("OUT: Error in stream writing!\n");
   if (1 != fread(&kind, 2, 1, fin)) throw ("OUT: Error in stream writing!\n");

   if (o->swap_in) {
      ByteSwap16u(&fea_size_BytesSWP);
      ByteSwap16u(&kind);
      ByteSwap32(&period);
      ByteSwap32(&frames);
   }
   fea_size=fea_size_BytesSWP/4;
//   cout<<fea_size<<endl;

/*  if(o->fea_c0) kind = kind | 020000;            // set c0 bit
  if(o->fea_E)  kind = kind | 000100;            // set E bit

if(020000==(11)^b) kind = 11; // LP cepstral coeffs (as in PLPC)
if(020000
        if (0==strcmp(o->fea_kind,"lpc"))          kind = 11; // LP cepstral coeffs (as in PLPC)
        else if (0==strcmp(o->fea_kind,"dctc"))    kind = 6;  // DCT cepstral coeffs (as in MFCC)
        else if (0==strcmp(o->fea_kind,"trapdct")) kind = 9;  // TRAP-DCT cepstral coeffs calculated
        else if (0==strcmp(o->fea_kind,"spec"))    kind = 8;  // spec    - linear mel-filter bank
        else if (0==strcmp(o->fea_kind,"logspec")) kind = 7;  // logspec - log mel-filter bank
        else kind = 6;

int a,b;
b=8198;
a=020000^b;
cout<<a<<endl;
*/

//--------------------------------------- View HTK header
//    cout << "kind: "<<kind << endl;
//    cout << "period: "<< period << endl;
//    cout << "frames: "<< frames << endl;
//    cout << "fea_size_BytesSWP: "<< fea_size_BytesSWP << endl;
//---------------------------------------
   try {
    	 cache = new float [fea_size];
   }
   catch (...) {
   	throw ("IN: Not enough memory!");
   }
   return true;
};

bool htkIN::get_frame() {
    static Vec<double> &fvec = *_fvec; // just a pointer

        if(fea_size != (int) fread(cache, sizeof(float), fea_size, fin)) return false;

        if(o->swap_in)
          for (int i=0; i<fea_size; i++)
            SwapFloat(&cache[i]);

        for(int i=0; i<fea_size; i++)
           fvec[i]=(double)cache[i];

//---------------------------------------- HList - View features
//        for(int i=0;i<fea_size;i++)
//           cout<< fvec[i] <<" ";
//           cout << endl;
//----------------------------------------
	return true;
};

htkIN::~htkIN() {
	delete [] cache;
};

// ---------------- cmvn file input ------------------------------
cmvnIN::cmvnIN(opts* o): rawIN(o) {
    num_spk=0;
    nfea = o->fea_ncepcoefs+1;
    num_spk=count_speakers(o->fcmvn_stat_in);
    if(num_spk){
      try {
           _M_cmean  = new Vec<double>*[num_spk];
           _M_cvar   = new Vec<double>*[num_spk];
           _list_ID_spk = new char*[num_spk];
      }
      catch (...) {
           throw ("IN: Not enough memory!");
      }
      new_file(o->fcmvn_stat_in);
    }else{
      _list_ID_spk = NULL;
    }
}

bool cmvnIN::new_file(char * filename){
   fin_cmvn.open(filename);
   if(!fin_cmvn) return false;
   char line[1000];
   int j=0;
   int new_spk=1;
   char *spk;
   // loop over all speakers
   while(fin_cmvn.getline(line,1000)){
        spk  = strtok(line," \t");      // input speaker name
        if((0==strcmp(spk,"mean") || 0==strcmp(spk,"var")) && new_spk){
          try{
             _M_cmean[j]     = new Vec<double>(nfea);
             _M_cvar[j]      = new Vec<double>(nfea);
             _list_ID_spk[j] = new char[(int)strlen(spk)];
          }
          catch(...) { throw ("IN: Not enough memory!"); }
          strcpy(_list_ID_spk[j],spk);
          new_spk=0;
        }
        if(!strcmp(spk,"mean")){
          Vec<double> & X = *_M_cmean[j];
//          cout<<"mean:";
          for(int i=0; i<nfea; i++){
             X[i] = atof(strtok(NULL," \t"));
//             cout<<" "<<X[i];
          }
//          cout<<endl;
        }else if(!strcmp(spk,"var")){
          Vec<double> & Y = *_M_cvar[j];
//          cout<<"var:";
          for(int i=0; i<nfea; i++){
             Y[i] = atof(strtok(NULL," \t"));
//             cout<<" "<<Y[i];
          }
          j++;
          new_spk=1;
//          cout<<endl;
        }
   }
   return true;
}

int cmvnIN::count_speakers(char *filename){
  int num_spk = 0;
  char line[1000];
  char *spk;
  ifstream f(filename);
  if(!f) return 0;
  while(f.getline(line,1000)){
       spk  = strtok(line," \t");      // input file name
       if(0==strcmp(spk,"mean") || 0==strcmp(spk,"var"))
         num_spk++;
  }
  f.close();
 return num_spk/2;
}

cmvnIN::~cmvnIN() {
        fin_cmvn.close();
};

// --------- public class IN -----------
IN::IN (opts* o_in) {
	o = o_in;

	if (0==strcmp(o->format_in,"raw"))    in = new rawIN(o);
	else if (0==strcmp(o->format_in,"alaw"))   in = new alawIN(o);
	else if (0==strcmp(o->format_in,"mulaw"))  in = new mulawIN(o);
	else if (0==strcmp(o->format_in,"wave"))   in = new waveIN(o);
	else if (0==strcmp(o->format_in,"htk"))    in = new htkIN(o);
	else if (0==strcmp(o->format_in,"cmvn"))   in = new cmvnIN(o);
	else throw ("IN: Unknown input file format!");

	_Xsabs = in->_Xsabs;
	_Xsph  = in->_Xsph;
        _fvec  = in->_fvec;
	E = &in->E;
        _M_cmean = in->_M_cmean;
        _M_cvar  = in->_M_cvar;
        _list_ID_spk = in->_list_ID_spk;
};

bool IN::new_file(char *filename) {
	return in->new_file(filename);  	// new signal file initializer
}

bool IN::loadf_filters(char *filename) {
	return in->loadf_filters(filename);  	// load coefs of filters from file
}

bool IN::get_frame() {
	return in->get_frame();
}

IN::~IN() {
	delete in;
}
