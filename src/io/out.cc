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
#include <iostream>
#include "out.h"
#include "pfile.h"
#include <math.h>

_OUT::_OUT(opts* o_in, Vec<double> *Xabs, Vec<double> *Xph, double * _E) {
	o = o_in;
	E = _E;
	size = Xabs->get_size();
	_Xsabs = Xabs;
	_Xsph  = Xph;

        if(sizeof(float) != 4)		// if not, all file outputs would be messed up
          throw "OUT: TOTAL PANIC! Size of float is not four bytes! :-O";
};
_OUT::_OUT(opts* o_in, Vec<double> **Xabs, Vec<double> **Xph,char **list_ID_spk) {
        _M_cmean =Xabs;
        _M_cvar  =Xph;
        _list_ID_spk =list_ID_spk;
};

void _OUT::ByteSwap16 (uSint* x) {
        *x = (((*x) >> 8) | (*x << 8));
};

void _OUT::ByteSwap16s (Sint *x) {
        *x = (((*x) >> 8) | (*x << 8));
};

void _OUT::ByteSwap32 (uLint *x) {
        *x = (((*x&0x000000FF)<<24)+((*x&0x0000FF00)<<8)+((*x&0x00FF0000)>>8)+((*x&0xFF000000)>>24));
};

void _OUT::SwapFloat (float* x) {
        uLint *f = (uLint*) x;
        *f = (((*f&0x000000FF)<<24)+((*f&0x0000FF00)<<8)+((*f&0x00FF0000)>>8)+((*f&0xFF000000)>>24));
};



// ---------------- HTK file output ------------------------------
class htkOUT : public _OUT {

	protected:
        FILE *fout;
        uLint frames;			// frames couter - for HTK header
		uSint fea_size_Bytes;	// # of bytes per feature vector
		int fea_size;						// feature vector size
		int Xsize;							// Xabs size (either spectrum or fea vector size without c0, E)
		float * obuffer;					// output buffer
		unsigned int get_fea_size();	// returns overall HTK feature size (including E, c0)
		void close();						// close one file
 	public:
		htkOUT (opts *, Vec<double>*, Vec<double>*, double *);
		virtual ~htkOUT ();
		virtual void new_file(char *);
		virtual void save_frame();
};

htkOUT::htkOUT(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph, double * E): _OUT(o,_Xsabs,_Xsph, E) {

    fout = NULL;
    // get feature size
	if (get_fea_size() > 32767) // HTK stores this in two bytes signed
		throw "OUT: HTK format does not support more than 32767 features!";
	fea_size = get_fea_size();     // feature vector size
	fea_size_Bytes = 4 * fea_size; // in bytes
        Xsize = _Xsabs->get_size();     // size of input Xabs vector (disregarding E, c0)
	try {
		obuffer = new float[fea_size];
    }
    catch (...) { throw ("OUT: Not enough memory!"); }

};

unsigned int htkOUT::get_fea_size() {
	// compute true output feature size
	if ((0==strcmp(o->fea_kind, "lpa")      ||
	     0==strcmp(o->fea_kind, "spec")     ||
	     0==strcmp(o->fea_kind, "logspec")) && o->fea_c0) {

          o->fea_c0 = false;
	  if (o->verbose)
            cerr << "OUT: WARNING: feature c0 not available for this fea_kind!" << endl;
	}

	unsigned int size = _Xsabs->get_size(); // in "spec" case we are done

	if (0==strcmp(o->fea_kind, "lpa")) size -= 1; // a0 is not written
	if (0==strcmp(o->fea_kind, "lpc")  && !o->fea_c0) size --;
	if (0==strcmp(o->fea_kind, "dctc") && !o->fea_c0) size --;
	if (o->fea_E) size++; // energy can be added even to spectrum
	return size;
};

void htkOUT::close() {	// updates HTK header and close
	if (!o->pipe_out) {
		// in output to pipe case no header is written
		if (0 != fseek(fout, 0L, SEEK_SET))
			throw "OUT: Cannot rewind output file for HTK header update!";
		if (o->swap_out) ByteSwap32(&frames);
		if (1 != fwrite (&frames, 4, 1, fout))
			throw "OUT: Error in stream writing!";
	}
	frames = 0; // discard the possibly byteswapped rubbish
	fclose(fout);
};

void htkOUT::new_file(char * filename) {
	// close last file, open new one
 	if (fout) close();

 	if (o->pipe_out) {			// output goes to STDOUT, no HTK header!
 		fout = fdopen(1,"wb");
	    if (!fout) throw "OUT: Cannot open standard output!";
	    frames = 0;
	    return;
 	}

 	// output to file
 	fout = fopen (filename,"wb");		// file
        if (!fout) throw "OUT: Cannot create output file!";

        frames = 0;

	// *** write HTK file header ***
        uLint period = (uLint)floor(.5 + 10000000.*o->wshift/(double)o->fs); // 100ns units
        uSint kind;
        if (0==strcmp(o->fea_kind,"lpc"))          kind = 11; // LP cepstral coeffs (as in PLPC)
        else if (0==strcmp(o->fea_kind,"dctc"))    kind = 6;  // DCT cepstral coeffs (as in MFCC)
        else if (0==strcmp(o->fea_kind,"trapdct")) kind = 9;  // TRAP-DCT cepstral coeffs calculated
        else if (0==strcmp(o->fea_kind,"spec"))    kind = 8;  // spec    - linear mel-filter bank
        else if (0==strcmp(o->fea_kind,"logspec")) kind = 7;  // logspec - log mel-filter bank
        else kind = 9;

	if (o->fea_c0) kind = kind | 020000;       	// set c0 bit
	if (o->fea_E)  kind = kind | 000100;   		// set E bit
	if (o->fea_delta && o->n_order>=1)  kind = kind | 000400;   		// set D bit
	if (o->fea_delta && o->n_order>=2)  kind = kind | 001000;   		// set A bit
	if (o->fea_delta && o->n_order==3)  kind = kind | 100000;   		// set T bit
        uSint fea_size_BytesSWP = fea_size_Bytes;

	if (o->swap_out) {
    	  ByteSwap16(&fea_size_BytesSWP);
          ByteSwap16(&kind);
          ByteSwap32(&period);
        }
	// # of frames will be updated at the end of file
	if (1 != fwrite (&frames, 4, 1, fout)) throw ("OUT: Error in stream writing!\n");
	if (1 != fwrite (&period, 4, 1, fout)) throw ("OUT: Error in stream writing!\n");
        if (1 != fwrite (&fea_size_BytesSWP, 2, 1, fout)) throw ("OUT: Error in stream writing!\n");
        if (1 != fwrite (&kind, 2, 1, fout)) throw ("OUT: Error in stream writing!\n");
};

void htkOUT::save_frame() {
   Vec<double> & X = *_Xsabs;

   if(0==strcmp(o->format_in, "htk")){
     for(int i=0; i<fea_size; i++)
       obuffer[i] = (float)X[i];
   }else{
        // fill in output buffer
        if(o->fea_trap) strcpy(o->fea_kind,"spec");
        if(0==strcmp(o->fea_kind, "spec") || 0==strcmp(o->fea_kind, "logspec") || 0==strcmp(o->fea_kind, "trapdct")) {       // spectrum output
            for(int i=0; i<Xsize; i++)
               obuffer[i] = (float) X[i];
            if(o->fea_E)
               obuffer[Xsize] = (float) *E;
	}else{
           // for a[k] or c[k], feature vector starts from index 1 not 0 !
           for(int j=0;j<=o->n_order;j++){
             for(int i=(((o->fea_ncepcoefs+1)*j)+1);i<((o->fea_ncepcoefs+1)*(j+1));i++){
                obuffer[i - 1]   = (float) X[i];
             }
             if(o->fea_c0){
                obuffer  [(o->fea_ncepcoefs)*(j+1)+j] = (float) X[(o->fea_ncepcoefs+1)*j]; // obuffer  [Xsize - 1] = (float) X[0];
                if(o->fea_E){
                  obuffer[Xsize] = (float) *E; //                 obuffer[(o->fea_ncepcoefs+1)*(j+1)] = (float) *E;
                }
             }else if (o->fea_E)
               obuffer[(o->fea_ncepcoefs)*(j+1)+j] = (float) *E;
           }
        }
   }

   if(o->swap_out)
     for(int i=0; i<fea_size; i++)
       SwapFloat(&obuffer[i]);

   // write it
   if(fea_size != (int) fwrite(obuffer, sizeof(float), fea_size, fout))
     throw "OUT: Error in stream writing!";
   frames++;
};

htkOUT::~htkOUT() {
	close();
	delete [] obuffer;
};

// ---------------- pfile output ------------------------------

class pfileOUT : public _OUT {

	protected:
 		PFile pf;
		int fea_size;		// full feature vector size
		int Xsize;			// feature vector size without E, c0
		float * obuffer;	// buffer

		bool firstframe;	// a little workaround (needed for pfile.h class)
		unsigned int *lab;	// in PFILE format also labels can be stored, this is a dummy ptr to labels
		uSint get_fea_size();	// includes also E, c0

 	public:
		pfileOUT (opts *, Vec<double>*, Vec<double>*, double *);
		virtual ~pfileOUT ();
		virtual void new_file(char *);
		virtual void save_frame();
};

pfileOUT::pfileOUT(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph, double * E) : _OUT(o,_Xsabs,_Xsph, E) {

	// get feature size
	fea_size = get_fea_size();    	// feature vector size
        Xsize = _Xsabs->get_size();    	// size of input Xabs vector (disregarding E, c0)
	firstframe = true;		// only in case of the very 1st frame we don't call pfile.startNewSent
	try {
		obuffer = new float[fea_size];	// data
		lab = new unsigned int;		// this needs to be available to pfile class, but will not be used
        }
        catch (...) { throw ("OUT: Not enough memory!"); }
        if (!pf.open(o->pfilename,"w", size, 0))
         throw "OUT: Cannot create output file!";
};

uSint pfileOUT::get_fea_size() {
	// compute output feature size
	if ((0==strcmp(o->fea_kind, "lpa") ||
		0==strcmp(o->fea_kind, "spec") ||
		0==strcmp(o->fea_kind, "logspec")) &&
		o->fea_c0) {
		o->fea_c0 = false;
		if (o->verbose) cerr << "OUT: WARNING: feature c0 not available for this fea_kind!" << endl;
	}
	uSint size = _Xsabs->get_size(); // in "spec" case we are done
	if (0==strcmp(o->fea_kind, "lpa")) size -= 1; // a0 is not written
	if (0==strcmp(o->fea_kind, "lpc") && !o->fea_c0) size --;
	if (0==strcmp(o->fea_kind, "dctc") && !o->fea_c0) size --;
	if (o->fea_E) size++; // energy can be added even to spectrum
	return size;
};

void pfileOUT::new_file(char * filename) {
	// one pfile contains all "speech files" inside called sentences,
	// so we only put "line break" here
	if (firstframe) { firstframe = false; }
	else { pf.startNewSent(); }
};

void pfileOUT::save_frame() {
    Vec<double> & X = *_Xsabs;

	// fill in output buffer
        if(o->fea_trap) strcpy(o->fea_kind,"spec");
        if(0==strcmp(o->fea_kind, "spec") || 0==strcmp(o->fea_kind, "logspec") || 0==strcmp(o->fea_kind, "trapdct")){
            for(int i=0; i<Xsize; i++)
               obuffer[i] = (float) X[i];
            if (o->fea_E)
               obuffer[Xsize] = (float) *E;
        }else{
           // for a[k] or c[k] features out vector starts with i=1..N
         for(int j=0;j<=o->n_order;j++){
          for(int i=(((o->fea_ncepcoefs+1)*j)+1);i<((o->fea_ncepcoefs+1)*(j+1));i++){
              obuffer[i - 1] = (float) X[i];
           if(o->fea_c0){
              obuffer  [(o->fea_ncepcoefs)*(j+1)+j] = (float) X[(o->fea_ncepcoefs+1)*j]; // obuffer  [Xsize - 1] = (float) X[0];
             if(o->fea_E)
                obuffer[Xsize] = (float) *E;
           }else if(o->fea_E)
              obuffer[(o->fea_ncepcoefs)*(j+1)+j] = (float) *E; //  obuffer[Xsize - 1] = (float) *E;
	  }
         }
        }
	// write it
	pf.write(obuffer, lab);
};

pfileOUT::~pfileOUT() {
	pf.startNewSent();	// last sent has to be closed properly
	pf.close();			// close pfile
	delete lab;
	delete [] obuffer;
};


// ---------------- GENERAL signal output ------------------------------

class sigOUT : public _OUT {

	protected:
 		double * cbuffer; 	    // circ. buffer for output samples [window]
		Sint * cache;       	// cache for file output [window] - mostly only [wstep] will be used
		int start; 				// pointer to cbuffer
		int w,s;				// window and shift sizes
		double * fft_in;
		double * fft_out;
        fftw_plan plan;
        FILE *fout;
        double correction;      // OLA correction
        Vec<double> *Win_corr;             ///< OLA ripple correction window

		virtual void open(char *) {};	// these guys will be implemented by specialized classes
		virtual void write(int) {};
		virtual void close() {};

    	void fill_cache(int N);		// move N samples from cbuffer to cache (also float -> int)

 	public:
		sigOUT (opts *, Vec<double>*, Vec<double>*, double *);
		virtual void new_file(char *);
		virtual void save_frame();
		virtual ~sigOUT ();
};


sigOUT::sigOUT(opts* o, Vec<double>*Xabs, Vec<double>*Xph, double *E) : _OUT(o,Xabs,Xph,E) {

	/* ---- get OLA window correction factor (depends on w, s) -------
	 * IDEA: When using general integer values for w and s and not the nice values when
	 * w/s=2 or w/s=4, we have to compute the ripple by going through all possible window
	 * placements and for every point compute the sum of all windows contribution. Among
	 * these contributions we remember the extremes and finally the ratio (max/min - 1)*100%
	 * gives the ripple.
	 */
	s = o->wshift;
	w = o->window;
	double alpha = 0.54;
	double pi = 2.*asin(1.);
	double min = 999.; // for OLA ripple computation
	correction = 0.;
	for (int i=0; i<s; i++) {
		// go through all different start offsets (repeats with period = s)
		int x = i; // initial offset of 1st window
		double y = 0.;
		while (x<w) {
			// while we are within the region of one window length,
			// add together contributions of all windows that affect value at this position
        	y += alpha - (1-alpha) * cos(2*pi*(double)x/(w - 1.)) ;
            x+=s;
        }
        if (y > correction) correction = y;	// max
        if (y < min) min = y;				// min
	}
	double ripple = floor(.5 + (correction/min - 1.)*10000.)/100.; // in %
	if (w/(double)s != 4 && w/(double)s != 2 && ripple>1.0 && o->verbose) {
		cerr << "OUT: WARNING! Non-optimal overlap of windows -> envelope ripple = "
			 << ripple << " %." << endl;
		cerr << "              Consider shift=window/2^N or other local minima." << endl;
	}

	// allocations
	try {
    	cbuffer = new double [o->window];
    	cache   = new Sint [o->window];
    }
    catch (...) { throw ("OUT: Not enough memory!"); }
	fft_in  = (double*) fftw_malloc(sizeof(double) * o->wfft);
    fft_out = (double*) fftw_malloc(sizeof(double) * o->wfft);
    if (fft_in == 0 || fft_out == 0) throw "OUT: Not enough memory!";

    plan = fftw_plan_r2r_1d(o->wfft, fft_in, fft_out, FFTW_HC2R, FFTW_MEASURE); // plan FFT
    fout = NULL;
	start = 0;
};

void sigOUT::new_file(char * filename) {
	// open file
	if (fout) close();	// close last one first
	open(filename);
	start = 0;
	for (int i=0; i<w; i++)
		cbuffer[i] = 0.;
};

void sigOUT::save_frame() {
	// when advance to the next frame, zero-fill the "forgotten" part of input buffer
	for (int i=0; i<s; i++)
		cbuffer[(start+w-s+i)%w] = 0.;

	Vec<double> & Xa = *_Xsabs;

	if (o->fb_power)  // in case of power input we need to sqrt Xabs :-/
		for (int i=0; i<size; i++)
			Xa[i] = sqrt(Xa[i]);

    // fill fft_in (abs,phi -> Re,Im), 1/N
    Vec<double> & Xp = *_Xsph;
    fft_in[0] = Xa[0]/(double)o->wfft; // DC
    fft_in[size-1] = Xa[size-1]/(double)o->wfft; if(Xp[size-1]!=0) Xa[size-1]=-Xa[size-1];
	for (int i=1; i<size-1; i++) {
		double ampl = Xa[i]/(double)o->wfft;
		fft_in[i] = ampl * cos(Xp[i]);
		fft_in[o->wfft-i] = ampl * sin(Xp[i]);
	}
	fftw_execute(plan); // iFFT

	// OLA
	for (int i=0; i<w; i++)
		cbuffer[(start+i)%w] += fft_out[i]; // discarding last wlen-window samples

	fill_cache(s);
	write(s);
	start += s;	// move circ. buffer start
};

void sigOUT::fill_cache(int n) {
	// move output signal to integer values and possibly clip, fill in cache
	// n is the # of samples to write (differs for the 1st frame and others)
	for (int i=0; i<n; i++) {
		int value = (int) floor(cbuffer[(start+i)%w] / correction);
		if (fabs((float)value) > 32767) {
			cache[i] = (value<0) ? -32767 : 32767;
			if (o->verbose) cerr << "OUT: WARNING: Clipped output sample (" << value
				 << ") to fit 16b dynamic range!" << endl;
		} else
			cache[i] = value;
		if (fabs((float)value) > 2*32768 && !o->quiet)
			cerr << "OUT: WARNING: Available dynamic range = 16bits, "
				 << "output sample even bigger than 17bits!" << endl;
	}
};

sigOUT::~sigOUT() {
	fftw_destroy_plan(plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
	delete [] cache;
	delete [] cbuffer;
};

// ---------------- raw signal output ------------------------------

class rawOUT : public sigOUT {

	protected:
		virtual void open(char *);
		virtual void write(int);
		virtual void close();

 	public:
		rawOUT (opts *, Vec<double>*, Vec<double>*, double *);
		virtual ~rawOUT() { close (); };
};

rawOUT::rawOUT(opts* o, Vec<double>*Xabs, Vec<double>*Xph, double *E) : sigOUT(o,Xabs,Xph,E) {};

void rawOUT::open(char * filename) {
	if (o->pipe_out) fout = fdopen(1,"wb");		        // STDOUT
	else fout = fopen (filename,"wb");			// file
	if (!fout) throw "OUT: Cannot open output stream!";
};

void rawOUT::close() {
    fill_cache(w-s);
    write(w-s);
	fclose(fout);
};

void rawOUT::write(int n) {

	// swap byte order if needed
	if (o->swap_out) for (int i=0; i<n; i++) ByteSwap16s(&cache[i]);

	// save raw data
	if (n != (int) fwrite(cache, sizeof(Sint), n, fout ))
		throw "OUT: Error in stream writing!";
};


// ---------------- wave signal output ------------------------------

class waveOUT : public sigOUT {

	protected:
		uLint samples;

		virtual void open(char *);
		virtual void write(int);
		virtual void close();

 	public:
		waveOUT (opts *, Vec<double>*, Vec<double>*, double *);
		virtual ~waveOUT() { close(); }
};


waveOUT::waveOUT(opts* o, Vec<double>*Xabs, Vec<double>*Xph, double *E) : sigOUT(o,Xabs,Xph,E) {};

void waveOUT::open(char * filename) {
	fout = fopen (filename,"wb");
	if (!fout) throw "OUT: Cannot open data file!";

	// write WAVE RIFF header
	fwrite("RIFF", 1, 4, fout);
	uLint size = 0; fwrite(&size, 4, 1, fout);
	fwrite("WAVE", sizeof(char), 4, fout);
	fwrite("fmt ", sizeof(char), 4, fout);
	Lint pcm=16; 	fwrite(&pcm, 4,1,fout);
	Sint cmpr=1;	fwrite(&cmpr, 2, 1, fout);
	Sint mono=1;	fwrite(&mono, 2, 1, fout);
	Lint fs=o->fs;  fwrite(&fs, 4, 1, fout);
	Lint Bps=fs*2;  fwrite(&Bps, 4, 1, fout);
	Sint Bpsamp=2;	fwrite(&Bpsamp, 2, 1, fout);
	Sint bpsamp=16;fwrite(&bpsamp, 2, 1, fout);
	fwrite("data", 1, 4, fout);
	fwrite(&size, 4, 1, fout);

	samples = 0;
};

void waveOUT::close() {
    fill_cache(w-s);
    write(w-s);

	// update WAVE header - need to fill in overall # of samples
	uLint data_size = 2*samples;
	uLint whole_size = data_size + 36; // 36+8 bytes for header
	if (0 != fseek(fout, 4L, SEEK_SET)) 
		throw "OUT: Cannot rewind output file for RIFF WAVE header update!";
	if (1 != fwrite (&whole_size, 4, 1, fout)) 
		throw "OUT: Error in stream writing!";
	if (0 != fseek(fout, 40L, SEEK_SET)) 
		throw "OUT: Cannot rewind output file for RIFF WAVE header update!";
	if (1 != fwrite (&data_size, 4, 1, fout)) 
		throw "OUT: Error in stream writing!";
	fclose(fout);
};

void waveOUT::write(int n) {
	// save raw data
	if (n != (int) fwrite(cache, 2, n, fout ))
		throw "OUT: Error in stream writing!";
	samples += n;	// count samples (for WAVE header)
};

// ---------------- CMVN file output ------------------------------
class cmvnOUT : public _OUT {

	protected:
        FILE *fout;
 	public:
		cmvnOUT (opts *, Vec<double>**, Vec<double>**,char **);
		virtual ~cmvnOUT ();
		virtual void new_file(char *);
		virtual void save_frame(int num_spk);
};


cmvnOUT::cmvnOUT(opts* o, Vec<double> **_Xsabs, Vec<double> **_Xsph,char **list_ID_spk): _OUT(o,_Xsabs,_Xsph,list_ID_spk) {
    fout = NULL;
}

void cmvnOUT::new_file(char * filename) {
	// close last file, open new one
 	if (fout) fclose(fout);
 	// output to file
 	fout = fopen (filename,"wt"); // file
        if (!fout) throw "OUT: Cannot create output file with stat. of cmvn!";
};

void cmvnOUT::save_frame(int num_spk){
 for(int j=0; j<num_spk; j++){
    Vec<double> & X = *_M_cmean[j];
    Vec<double> & Y = *_M_cvar[j];
    char & Z = *_list_ID_spk[j];
    char x[999];
    int Xsize = _M_cmean[0]->get_size()-1;
    fprintf(fout, "%s\nmean\t",&Z);
    for (int i = 0; i < Xsize*2; i++){
      if(i<Xsize-1){
        sprintf(x, "%f ",X[i] );
      }else if(i<Xsize){
        sprintf(x, "%f",X[i] );
        fprintf(fout, "%s",x);
        fprintf(fout, "%s","\nvar\t");
        continue;
      }else if(i==Xsize*2-1){
        sprintf(x, "%f\n",Y[i-Xsize]);
      }else if(i>=Xsize){
        sprintf(x, "%f ",Y[i-Xsize]);
      }
      fprintf(fout, "%s",x);
    }
 }
};

cmvnOUT::~cmvnOUT() {
        fclose(fout);
};

// ---------------- ark file output (KALDI format)------------------------
class arkOUT : public _OUT {

	protected:
                  FILE *fark;
                  FILE *fscp;
                  uLint frames;           // frames couter - for ark header
//                  long int frames;           // frames couter - for ark header
                  int fea_size;           // full feature vector size
                  int Xsize;              // feature vector size without E, c0
                  float * obuffer;        // buffer

                  bool firstframe;        // a little workaround (needed for pfile.h class)
                  uSint get_fea_size();   // includes also E, c0
                  __off64_t  ark_index;
                  long int token_index;
                  char scpfilename[999];  // name of output ark kaldi format(when applicable)
                  long int frame_len_index;
                  void update_header();   // update ark header
           public:
                  arkOUT (opts *, Vec<double>*, Vec<double>*,double *);
                  virtual ~arkOUT ();
                  virtual bool rename_path_ark_to_scp(char *string);
                  virtual void new_file(char *);
                  virtual void save_frame();
};

arkOUT::arkOUT(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph, double * E) : _OUT(o,_Xsabs,_Xsph, E) {
	// get feature size
	fea_size = get_fea_size();    	// feature vector size
        Xsize = _Xsabs->get_size();    	// size of input Xabs vector (disregarding E, c0)
	firstframe = true;		// only in case of the very 1st frame we don't call pfile.startNewSent
        frames=0;
        frame_len_index=0;
	try {
		obuffer = new float[fea_size];	// data
        }
        catch (...) { throw ("OUT: Not enough memory!"); }
        // output to file
        fark = fopen64(o->arkfilename,"wb");           // file
        if (!fark) throw "OUT: Cannot create output ark file!";

        strcpy(scpfilename, o->arkfilename);
        rename_path_ark_to_scp(scpfilename);
        if ((fscp = fopen64(scpfilename, "wt")) == NULL) throw "Cannot open output scp file for writing!";
};

void arkOUT::update_header() {  // updates ark header
     if(0 != fseeko64(fark, (__off64_t) frame_len_index, SEEK_SET))
       throw "OUT: Cannot rewind output file for ark header update!";
     if(o->swap_out)
       ByteSwap32(&frames);
     if(1 != fwrite (&frames, 4, 1, fark))
       throw "OUT: Error in stream writing!";
     fseeko64(fark, (__off64_t)0, SEEK_END);
     frames =0;
     firstframe = true;
};

void arkOUT::new_file(char * filename) {
        if(frames){
          update_header();
        }
	if(firstframe){
          firstframe = false;
          if(fprintf(fark, "%s %cBFM %c", filename, 0, 4) < 0){
            printf("Couldn't write header to %s, quitting!\n", filename);
          }else{
            frame_len_index = ftello64(fark);
          }
          fwrite(&frames, 4, 1, fark);
          fprintf(fark, "%c", 4);
          fwrite(&fea_size, 4, 1, fark);
        }
        if(frames == 0){

          ark_index = ftello64(fark) -15;   // fast access index
          fprintf(fscp, "%s %s:%llu\n", filename, o->arkfilename, ark_index);
        }
};

uSint arkOUT::get_fea_size() {
	// compute output feature size
	if ((0==strcmp(o->fea_kind, "lpa") ||
		0==strcmp(o->fea_kind, "spec") ||
		0==strcmp(o->fea_kind, "logspec")) &&
		o->fea_c0) {
		o->fea_c0 = false;
		if (o->verbose) cerr << "OUT: WARNING: feature c0 not available for this fea_kind!" << endl;
	}
	uSint size = _Xsabs->get_size(); // in "spec" case we are done
	if (0==strcmp(o->fea_kind, "lpa")) size -= 1; // a0 is not written
	if (0==strcmp(o->fea_kind, "lpc") && !o->fea_c0) size --;
	if (0==strcmp(o->fea_kind, "dctc") && !o->fea_c0) size --;
	if (o->fea_E) size++; // energy can be added even to spectrum

	return size;
};

void arkOUT::save_frame() {
    Vec<double> & X = *_Xsabs;

    if(0==strcmp(o->format_in, "htk")){
          for(int i=0; i<fea_size; i++)
        obuffer[i] = (float)X[i];
    }else{
        // fill in output buffer
        if(o->fea_trap) strcpy(o->fea_kind,"spec");
        if(0==strcmp(o->fea_kind, "spec") || 0==strcmp(o->fea_kind, "logspec") || 0==strcmp(o->fea_kind, "trapdct")) {       // spectrum output
            for(int i=0; i<Xsize; i++)
               obuffer[i] = (float) X[i];
            if (o->fea_E)
               obuffer[Xsize] = (float) *E;
        }else{
           // for a[k] or c[k] features out vector starts with i=1..N
         for(int j=0;j<=o->n_order;j++){
          for(int i=(((o->fea_ncepcoefs+1)*j)+1);i<((o->fea_ncepcoefs+1)*(j+1));i++){
              obuffer[i - 1] = (float) X[i];
           if(o->fea_c0){
              obuffer  [(o->fea_ncepcoefs)*(j+1)+j] = (float) X[(o->fea_ncepcoefs+1)*j]; // obuffer  [Xsize - 1] = (float) X[0];
             if(o->fea_E)
                obuffer[Xsize] = (float) *E;
           }else if(o->fea_E)
              obuffer[(o->fea_ncepcoefs)*(j+1)+j] = (float) *E; //  obuffer[Xsize - 1] = (float) *E;
	  }
         }
        }
    }

    // write it
    if(fea_size != (int) fwrite(obuffer, sizeof(float), fea_size, fark))
      throw "OUT: Error in stream writing!";
    frames++;
};

bool arkOUT::rename_path_ark_to_scp(char *string){
         char line[1000]={'\0'};;
         char *to;
         to = strtok(string,".");
         while (to != NULL){
              if(strcmp(to,"ark")){
                strcat(line,to);
                strcat(line,".");
              }else{
                strcpy(string,line);
                strcat(string,"scp");
                return true;
              }
              to = strtok(NULL,".");
        }
        strcpy(string,line);
        strcat(string,"scp");
        return false;
}

arkOUT::~arkOUT() {
      update_header();
      fclose(fark);
      fclose(fscp);
      delete [] obuffer;
};

// --------- public class OUT -----------
OUT::OUT (opts* o_in, Vec<double> * Xabs, Vec<double> * Xph, double *_E) {
	o = o_in;
	_Xsabs = Xabs;
	_Xsph  = Xph;
	E = _E;
	if      (0==strcmp(o->format_out,"htk"))    out = new htkOUT(o,Xabs,Xph,E);
	else if (0==strcmp(o->format_out,"pfile"))  out = new pfileOUT(o,Xabs,Xph,E);
	else if (0==strcmp(o->format_out,"ark"))    out = new arkOUT(o,Xabs,Xph,E);
	else if (0==strcmp(o->format_out,"raw"))    out = new rawOUT(o,Xabs,Xph,E);
	else if (0==strcmp(o->format_out,"wave"))   out = new waveOUT(o,Xabs,Xph,E);
	else throw ("OUT: Unknown output file format!");
};
OUT::OUT (opts* o_in, Vec<double> ** M_cmean, Vec<double> ** M_cvar,char **list_ID_spk) {
	o = o_in;
	if      (0==strcmp(o->format_out,"cmvn"))    out = new cmvnOUT(o,M_cmean,M_cvar,list_ID_spk);
	else throw ("OUT: Unknown output file format!");
};

void OUT::new_file(char * filename) {
	out->new_file(filename);
}

void OUT::save_frame() {
	out->save_frame();
}
void OUT::save_frame(int num_spk) {
	out->save_frame(num_spk);
}

OUT::~OUT() {
	delete out;
}
