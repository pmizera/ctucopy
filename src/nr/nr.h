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

#ifndef NR_H_
#define NR_H_

#include "../io/opts.h"
#include "../base/types.h"
#include <math.h>

class _NR {
    // --- INTERNAL CLASS, please use the one below! -------
    // prototype class inherited by all speech-enhancing classes

  protected:
	opts *o;
	int size;
	Vec<double> * _Xsabs;
	Vec<double> * _Xsph;
	
	void compute_E();

  public:
	double E;

	_NR(opts *, Vec<double> *, Vec<double> *);
	virtual bool new_file() {return 0;};
	virtual bool process_frame() {return 0;};
	virtual ~_NR() {};
};

//---------------------------------------------------

class NR {
    // this is the only class intended to be public
    // it implements all needed methods using subclasses <x>NR

  protected:
  	opts* o;
	_NR* nr;
	Vec<double> * _Xsabs;
	Vec<double> * _Xsph;

  public:
  	double * E;		// ptr to frame energy
  
	NR(opts *, Vec<double> *, Vec<double> *); 	// NR(options, ptr to Xabs, ptr to Xphase)
	bool new_file();							// announces new file
	bool process_frame();						// suppress noise from one frame
	~NR();
};

#endif /*NR_H_*/

