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
#include "../io/stdafx.h"

#include <math.h>
#include "../vdet/CepstralDet.h"
#include <fftw3.h>
#include "nr.h"

using namespace Voice; // for detectors
using namespace DSP;   // for detectors


_NR::_NR(opts* o_in, Vec<double> *Xabs, Vec<double> *Xph) {
	o = o_in;
	size = Xabs->get_size();
	_Xsabs = Xabs;
	_Xsph  = Xph;
	E = -1.;
};

void _NR::compute_E() {
	// compute frame energy
    Vec<double> & X = *_Xsabs;
 	if (o->fea_E && !o->fea_rawenergy) {
		E = X[0]*X[0]/2. + X[size - 1]*X[size - 1]/2.;
		for (int i=1; i<size - 1; i++)
			E += X[i]*X[i];
		E = log(E*2.);
	}
}


// ---------------- no noise reduction     -----------------

class noNR : public _NR {
	public:
		noNR (opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph) : _NR(o,_Xsabs,_Xsph) {};
		virtual ~noNR () {};
		virtual bool new_file() { return true; };
		virtual bool process_frame() { compute_E(); return true; };
};

// ---------------- Extended Spectral Subtraction -----------------

class extenNR : public _NR {

	protected:
		double a, p;
		double *N, *Navg, *H, *Yavg;
		// noise est, smoothed noise est, Wiener transfer func, smoothed speech estim
	public:
		extenNR (opts *, Vec<double>*, Vec<double>*);
		virtual ~extenNR ();
		virtual bool new_file();                 // new file initializer, return =
		virtual bool process_frame();            // loads frame, returns ?
};

extenNR::extenNR(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph) : _NR(o,_Xsabs,_Xsph) {
	a = o->nr_a;	// powering factor of the spectrum
	p = o->nr_p;	// integrator (smoothing) factor

    try {
        N = new double [size];
        Navg = new double [size];
        H    = new double [size];
        Yavg = new double [size];
    }
    catch (...) { throw ("NR: Not enough memory!"); }
};

bool extenNR::new_file() {
    // initialization
    for (int i = 0; i < size; i++) {
	    Navg[i] = 0.95;	// assume at the beginning 95% of input is noise
		Yavg[i] = 0.05;	// so that we get H = 0.95
   	}
	return true;
};

bool extenNR::process_frame() {
    Vec<double> & X = *_Xsabs;

        // initialization
    if(Navg[0]==-1){
    for (int i = 0; i < size; i++) {
            Navg[i] = X[i];     // assume at the beginning 95% of input is noise
            Yavg[i] = 0; // so that we get H = 0.95
        }
    }
	// adaptation of Wiener filter transfer function
    if (a == 1.0)		// amplitude mode
    	for (int i = 0; i < size; i++)
        	H[i] = Navg[i] / (Navg[i] + Yavg[i]);
    else if (a == 2.0)	// power mode
    	for (int i = 0; i < size; i++)
            H[i] = Navg[i] / sqrt(Navg[i]*Navg[i] + Yavg[i]*Yavg[i]);
    else				// general rational power - slow
        for (int i = 0; i < size; i++)
        	H[i] = Navg[i] / pow(pow(Navg[i],a)+pow(Yavg[i],a), 1./a);

	// current noise estimation
	for (int i = 0; i < size; i++)
		N[i] = H[i] * X[i];

	// noise averaging
	for (int i = 0; i < size; i++)
    	Navg[i] = p * Navg[i] + (1 - p) * N[i];

	// spectral subtraction with FWR
   	for (int i = 0; i < size; i++) {
    	if (X[i] > Navg[i])
    		Yavg[i] = X[i] - Navg[i];
	    else
	    	Yavg[i] = Navg[i] - X[i];
   	}

    // final noise subtraction
    for (int i = 0; i < size; i++)
//  		X[i] -= Yavg[i] ;
  		X[i] -= N[i] ;

	compute_E();

    return true;
};

extenNR::~extenNR() {
	delete [] H;
	delete [] Navg;
	delete [] N;
	delete [] Yavg;
};


// ---------------- HWSS - Standard Spectral Subtraction with VAD -----------------

class hwssNR : public _NR {

	protected:
		bool vad;
		double a,b,p;
		int ninitsegs;
		double *Navg;

		CepstralDetector<BurgCepstrumEstimator> *cdet;
		double * fft_in, * fft_out; fftw_plan plan; // BRUTAL HACK for Burg (see below)
		FILE *fvad;

	public:
		hwssNR (opts *, Vec<double>*, Vec<double>*);
		virtual ~hwssNR ();
		virtual bool new_file();                 // new file initializer, return =
		virtual bool process_frame();            // loads frame, returns ?
		bool vad_init();
		bool vad_get_frame();
};

/*
 * Hack for Burg: Normally, Burg detector expects signal frame on it's input which
 * it itself converts to cepstrum. Since I did not have time to go through the code
 * to understand how to use directly spec/ceps at the input, I do the ugly hack that
 * I go back to time domain using iFFT and give the signal to Burg.
 * Definitely could be fixed...
*/

hwssNR::hwssNR(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph) : _NR(o,_Xsabs,_Xsph) {
	ninitsegs = o->nr_initsegs;
	a = o->nr_a;	// spectrum powering const
	b = o->nr_b;	// spectrum oversubtraction factor (<1=less >1=more)
	p = o->nr_p;	// smoothing const
    try {
        Navg = new double [size];
    }
    catch (...) { throw ("NR: Not enough memory!"); }

	// BRUTAL HACK HERE---------
    // when `Burg' VAD, then iFFT needed to get time signal back
	if (0==strcmp(o->vadmode,"burg")) {
	    if (_Xsph==NULL || _Xsabs->get_size() != _Xsph->get_size())
    		throw "NR: Cannot use Burg detector after filter bank!";
		fft_in  = (double*) fftw_malloc(sizeof(double) * o->wfft);
    	fft_out = (double*) fftw_malloc(sizeof(double) * o->wfft);
    	if (fft_in == 0 || fft_out == 0) throw "IN: Not enough memory!";
    	plan = fftw_plan_r2r_1d(o->wfft, fft_in, fft_out, FFTW_HC2R, FFTW_MEASURE);
	}
	// BRUTAL HACK end
	cdet = NULL;

	// if VAD from file, open it
	fvad = NULL;
	if (strcmp(o->vadmode,"file")==0) {
		fvad = fopen ( o->filevad, "r" );
		if (fvad == NULL) throw "NR: Unable to open VAD file!\n";
	}
};

bool hwssNR::new_file() {
    // initialization
	ninitsegs = o->nr_initsegs;
	vad_init();
    Vec<double> & X = *_Xsabs;
    for (int i = 0; i < size; i++) {
        Navg[i] = pow(X[i],a);
		X[i] *= 0.1;	// init phase -> out = 0.1 * in
    }
	return true;
};

bool hwssNR::process_frame() {
    Vec<double> & X = *_Xsabs;
	ninitsegs--;

	// dynamic expansion
    if (a == 2.0)
    	for (int i = 0; i < size; i++)
            X[i] *= X[i];
    else
    	if (a != 1.0)
        	for (int i = 0; i < size; i++)
        		X[i] = pow(X[i], a);

	// get VAD
	vad = vad_get_frame();

    // noise estimate smoothing:
	if (vad == 0 || ninitsegs > 0)
    	for (int i = 0; i < size; i++)
	    	Navg[i] = p*Navg[i] + (1 - p)*X[i];

    // half-wave rectification spectral subtraction:
  	for (int i = 0; i < size; i++) {
		X[i] -= b*Navg[i];
        if (X[i] < 0.) X[i] = 0.;
  	}

	// dynamic compression
	if (a == 2.)
	  	for (int i = 0; i < size; i++)
			X[i] = sqrt(X[i]);
	else if (a != 1.)
	  	for (int i = 0; i < size; i++)
			X[i] = pow(X[i], 1.0/a);

	compute_E();
    return true;
};

bool hwssNR::vad_init() {
	// if VAD is Burg's method, then initialize external Burg's VAD
	if (0==strcmp(o->vadmode,"burg")) {
		typedef CepstralDetector<BurgCepstrumEstimator> CDBurg;
		int ncep = o->fea_ncepcoefs; double q = o->nr_q;
        CDBurg::Options opt(ninitsegs, ncep, p, q);
        if (cdet != NULL) delete cdet;
        cdet = new CDBurg(o->window, opt);
	}
	else if (strcmp(o->vadmode,"file")==0) {
		// no separation of files in VAD file (yes, it is dangerous, but...)
	} else if (strcmp(o->vadmode,"none")==0) {
		throw "NR: Please specify Voice Activity Detector!";
	} else
		throw "NR: Unknown VAD mode!";
	return true;
};

bool hwssNR::vad_get_frame() {
	if (0==strcmp(o->vadmode,"burg")) {
	    // BRUTAL HACK start
		Vec<double> & Xa = *_Xsabs;
		Vec<double> & Xp = *_Xsph;

		for (int i=0; i<o->wfftby2; i++)
			fft_in[i] = Xa[i] * cos(Xp[i]);
		for (int i=o->wfft-1;i>=o->wfftby2; i--)
			fft_in[i] = Xa[o->wfft - i] * sin(Xp[o->wfft - i]);
		fftw_execute(plan); // iFFT
		vad = (int) cdet->Process(&fft_out[0], &fft_out[0] + o->window);
//                cout << vad << endl;
		return bool(vad);
		// BRUTAL HACK end
	}
	else if (strcmp(o->vadmode,"file")==0) {
		// process VAD file
		char vad = fgetc(fvad);
		if (vad != EOF) return bool(vad);
		else throw "NR: Unexpected end of VAD file!";
	} else
		throw "NR: Unknown VAD mode!";
	return true;
};

hwssNR::~hwssNR() {
	if (0==strcmp(o->vadmode,"burg")) {
		delete cdet;
		fftw_destroy_plan(plan); // BRUTAL HACK
		fftw_free(fft_in); // BRUTAL HACK
		fftw_free(fft_out); // BRUTAL HACK
	}
	if (fvad != NULL) fclose(fvad);
	delete [] Navg;
};

// ---------------- FWSS - Standard Spectral Subtraction with VAD -----------------

class fwssNR : public hwssNR {

	protected:

	public:
		fwssNR (opts *, Vec<double>*, Vec<double>*);
		virtual bool process_frame();            // loads frame, returns ?
};

fwssNR::fwssNR(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph) : hwssNR(o,_Xsabs,_Xsph) {};

bool fwssNR::process_frame() {
    Vec<double> & X = *_Xsabs;

	// dynamic expansion
    if (a == 2.0)
    	for (int i = 0; i < size; i++)
            X[i] *= X[i];
    else
    	if (a != 1.0)
        	for (int i = 0; i < size; i++)
        		X[i] = pow(X[i], a);

	// get VAD
	vad = vad_get_frame();

    // noise estimate smoothing:
	if (vad == 0 || ninitsegs > 0)
    	for (int i = 0; i < size; i++)
	    	Navg[i] = p*Navg[i] + (1 - p)*X[i];

    // full-wave rectification spectral subtraction:
  	for (int i = 0; i < size; i++) {
		X[i] -= b*Navg[i];
        if (X[i] < 0.) X[i] = -X[i];
  	}

	// dynamic compression
	if (a == 2.)
	  	for (int i = 0; i < size; i++)
			X[i] = sqrt(X[i]);
	else if (a != 1.)
	  	for (int i = 0; i < size; i++)
			X[i] = pow(X[i], 1.0/a);

	compute_E();

	ninitsegs--;
    return true;
};


// ---------------- 2-step FWSS - Standard Spectral Subtraction with VAD -----------------
//             (dfwss = double full wave SS)

class dfwssNR : public hwssNR {

	protected:
		double *Nravg;		// residual noise after 1st SS step

	public:
		dfwssNR (opts *, Vec<double>*, Vec<double>*);
		virtual ~dfwssNR ();
		virtual bool new_file();
		virtual bool process_frame();
};


dfwssNR::dfwssNR(opts* o, Vec<double> *_Xsabs, Vec<double> *_Xsph) : hwssNR(o,_Xsabs,_Xsph) {
	if (o->nr_a != 1. && o->verbose)
		cerr << "NR: WARNING! Power constant nr_a not used in 2-pass FW SS mode!" << endl;
	try {
        Nravg = new double [size];
    }
    catch (...) { throw ("NR: Not enough memory!"); }
};

bool dfwssNR::new_file() {
    // initialization
	ninitsegs = o->nr_initsegs;
	vad_init();
	Vec<double> & X = *_Xsabs;
    for (int i = 0; i < size; i++) {
        Navg[i] = X[i];
        Nravg[i] = 0.;
		X[i] *= 0.1;
    }
	return true;
};

bool dfwssNR::process_frame() {
    Vec<double> & X = *_Xsabs;

	// get VAD
	vad = vad_get_frame();

    // rough noise estimate smoothing:
	if (vad == 0 || ninitsegs > 0)
    	for (int i = 0; i < size; i++)
	    	Navg[i] = p*Navg[i] + (1 - p)*X[i];

    // full-wave rectification spectral subtraction:
  	for (int i = 0; i < size; i++) {
		X[i] -= Navg[i];
        if (X[i] < 0.) X[i] = -X[i];
  	}

    // residual noise estimate smoothing:
	if (vad == 0 || ninitsegs > 0)
    	for (int i = 0; i < size; i++)
	    	Nravg[i] = p*Nravg[i] + (1 - p)*X[i];

    // full-wave rectification spectral subtraction:
  	for (int i = 0; i < size; i++) {
		X[i] -= Nravg[i];
        if (X[i] < 0.) X[i] = -X[i];
  	}

	compute_E();

	ninitsegs--;
    return true;
};

dfwssNR::~dfwssNR() {
	delete [] Nravg;
};

// --------- public class NR -----------

NR::NR (opts* o_in, Vec<double> * Xabs, Vec<double> * Xph) {
	o = o_in;
	_Xsabs = Xabs;
	_Xsph  = Xph;

	if (0==strcmp(o->nr_mode,"exten")) nr = new extenNR(o,Xabs,Xph);
	else if (0==strcmp(o->nr_mode,"hwss")) nr = new hwssNR(o,Xabs,Xph);
	else if (0==strcmp(o->nr_mode,"fwss")) nr = new fwssNR(o,Xabs,Xph);
	else if (0==strcmp(o->nr_mode,"2fwss")) nr = new dfwssNR(o,Xabs,Xph);
	else if (0==strcmp(o->nr_mode,"none")) nr = new noNR(o,Xabs,Xph);
	else throw ("NR: Unknown noise reduction mode!");
	E = &nr->E;
};

bool NR::new_file() {
	return nr->new_file();
}

bool NR::process_frame() {
	return nr->process_frame();
}

NR::~NR() {
	delete nr;
}
