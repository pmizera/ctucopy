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

#ifndef FB_H_
#define FB_H_

#include "../io/stdafx.h"
#include "../io/opts.h"
#include "../base/types.h"
#include <math.h>

class FB {
    // Filter Bank class. Implements nonlinear projection of input amplitude
    // spectrum to output spectrum.

  protected:
  	opts* o;
	double *hz_axis, *warp_axis; 	// linear and nonlinear scales
	double warp_low, warp_high;		// freq boundaries on warped scale
	bool PLP; 						// PLP filterbank flag (PLP has special V.I.P. settings)
	int subbanks;                   // number of subbanks (FB::parse will set this)
	struct subbank{					// struct for one subbank
		double f_start;
		double f_stop; 
		int bands;
		int band_first;
		int band_last;
	};
	subbank* bank;					// this stores complete FB definition
			
	void parse(); 		// parses user fb_definition
	void check_and_design();      // checks subbank boundaries and designs the bank
	void plp_design();  // designs the one and only PLP filter bank with trapezoidal filters
	void init_scale();  // initializes warped freq. scale
	void get_filter(double*, double, double, int, int);	// design only one filter
	// get_filter(ptr to output coeffs vector,freq_low, freq_high, #filters there, which one);
	void optimize();	// to be run after all filters have been designed
	
  public:

	double **mat;		// actual filters in rows of a matrix
	int size;			// number of output filters
	Vec<double> * _Xs;	// ptr to input ampl. spectrum 
	Vec<double> * _Y;	// ptr to output spectrum

	FB(opts *, Vec<double> *);
	bool new_file();	// annnounce new file - most likely unused here
	bool project_frame();	// projection X -> Y
	~FB();
};

#endif /*FB_H_*/

