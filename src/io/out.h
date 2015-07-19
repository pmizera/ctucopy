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

#ifndef OUT_H_
#define OUT_H_

#include "opts.h"
#include "../base/types.h"
#include <math.h>
#include <fftw3.h>

class _OUT {
    // prototype class inherited by all classes

  protected:
	opts *o;
	Vec<double> * _Xsabs;
	Vec<double> * _Xsph;
	Vec<double> ** _M_cmean;
	Vec<double> ** _M_cvar;
        char ** _list_ID_spk;
	double *E;
	int size;

    void ByteSwap16 (uSint *);
    void ByteSwap16s (Sint *);
    void ByteSwap32 (uLint *);
    void SwapFloat (float*);

  public:
	_OUT(opts *, Vec<double> *, Vec<double> *, double *);
	_OUT(opts *, Vec<double> **, Vec<double> **,char **list_ID_spk);
	virtual void new_file(char * filename) {};
	virtual void save_frame() {};
	virtual void save_frame(int num_spk) {};
	virtual ~_OUT() {};
};

//---------------------------------------------------

class OUT {
    // this is the only class intended to be public
    // it implements all needed methods using subclasses <x>OUT

  protected:
  	opts* o;
	_OUT* out;
	Vec<double> * _Xsabs;
	Vec<double> * _Xsph;
	Vec<double> ** _M_cmean;
	Vec<double> ** _M_cvar;
        char ** list_ID_spk;
	double *E;

  public:
	OUT(opts *, Vec<double> *, Vec<double> *, double *);
	OUT(opts *, Vec<double> **, Vec<double> **,char **);
	// OUT(options, ptr to Xabs, ptr to Xphase, ptr to Energy)
	// note: Xabs is generally ptr to output vector - can be features,
	//		 Xph in that case is disregarded
	void new_file(char * filename);		// announces new file
	void save_frame(int num_spk);		// saves one frame of signal/features
	void save_frame();                      // saves one frame of signal/features
	~OUT();
};

#endif /*OUT_H_*/

