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

#ifndef IN_H_
#define IN_H_

#include "opts.h"
#include "../base/types.h"
#include <fftw3.h>
#include <math.h>

class _IN {
    // PRIVATE CLASS - prototype inherited by all speech-enhancing classes
	// one file processing utility
	// opening, loading, buffering, segmentation..

  protected:
	opts *o;
	double *W;              // Hamming window coefficiens

	void ByteSwap16 (Sint*);            // byte swapping utility
	void ByteSwap16u (uSint*);            // byte swapping utility
        void ByteSwap32 (uLint *);
        void SwapFloat (float *);
	void hamming (double alpha); 		// initializes Hamming window
	double c_abs (double *re, double *im) { return sqrt (*re * *re + *im * *im); } // abs(Re,Im)
	double c_ph  (double, double);		// angle(Re,Im) (rad)
        int do_td_iir_mfcc;
  public:
	Vec<double> *_Xsabs;     // buffer for magnitude spectral samples
	Vec<double> *_Xsph;      // buffer for phase spectral samples
        Vec<double> *_fvec;
	double E;				// frame energy. Should not be public but...
        Vec<double> ** _M_cmean;
        Vec<double> ** _M_cvar;
        char ** _list_ID_spk;

	_IN(opts *);
	virtual bool new_file(char*) {return 0;};  // init new file
        virtual bool loadf_filters(char*) {return 0;};
	virtual bool get_frame() {return 0;};      // loads one frame, sets Xabs,Xph,E
	virtual ~_IN();
//        virtual Vec<double>* get_fvec() {return _fvec;}
};

class IN {
    // ------------ USE THIS ONE -----------------
    // this is the only class intended to be public
    // it implements all needed methods using subclasses <x>IN

  protected:
	opts *o;				// options
	_IN* in;

  public:
	Vec<double> *_Xsabs;	 // ptr to buffer for magnitude spectral samples
	Vec<double> *_Xsph;      // ptr to buffer for phase spectral samples
        Vec<double> *_fvec;
	double * E;		 // ptr to frame energy
        Vec<double> ** _M_cmean;
        Vec<double> ** _M_cvar;
        char ** _list_ID_spk;

	IN(opts *);
//	double get_E();
	bool new_file(char*);  	// init new file
        bool loadf_filters(char*);
	bool get_frame();    	// loads one frame, sets Xabs,Xph,E
	~IN();
};

#endif /*IN_H_*/

