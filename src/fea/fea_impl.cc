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

#include "fea_impl.h"
#include <fstream>

_FEA::_FEA(opts* o_in, Vec<double> *X) {
	o = o_in;
	Nin = X->get_size();
	Nout = 0;
	_Xs = X;
};

// ---------------- dummy class for spectrum (almost no function) -------------

specFEA::specFEA(opts* o, Vec<double> *_Xs) : _FEA(o,_Xs) {
   try {
        _fvec = new Vec<double>(Nin);
    }
    catch (...) { throw ("FEA: Not enough memory!"); }
};

bool specFEA::process_frame() {
    Vec<double> & X = *_Xs;
    Vec<double> & F = *_fvec;

    for (int i = 0; i < Nin; i++)
    	F[i] = X[i];

	// compute frame energy
 	if (o->fea_E && !o->fea_rawenergy) {
		E = X[0]*X[0]/2. + X[Nin - 1]*X[Nin - 1]/2.;
		for (int i=1; i<Nin - 1; i++)
			E += X[i]*X[i];
		E = log(E*2.);
	}

    return true;
};

specFEA::~specFEA() {
	delete _fvec;
};

// ---------------- class for log spectrum  -------------

bool logspecFEA::process_frame() {
    Vec<double> & X = *_Xs;
    Vec<double> & F = *_fvec;

    for (int i = 0; i < Nin; i++)
    	F[i] = log(X[i]);

	// compute frame energy
 	if (o->fea_E && !o->fea_rawenergy) {
		E = X[0]*X[0]/2. + X[Nin - 1]*X[Nin - 1]/2.;
		for (int i=1; i<Nin - 1; i++)
			E += X[i]*X[i];
		E = log(E*2.);
	}

    return true;
};

// ---------------- cepstral coeffs via DCT -----------------

dctcFEA::dctcFEA(opts* o, Vec<double> *_Xs) : _FEA(o,_Xs) {
	Nout = o->fea_ncepcoefs + 1; // plus c0
    try {
       _fvec = new Vec<double>(Nout);
       wdct = new double [4*Nin];
       lift= new double [Nout - 1];
    }
    catch (...) { throw ("FEA: Not enough memory!"); }

	// DCT coefficients init:
	// w[i] = cos (2*pi*i/4N)
	for (int i=0; i < 4*Nin; i++)
		wdct[i] = cos(3.1415926535898*(double)i/(2*Nin));

	// lifter coefficients init
	for (int n = 0; n < Nout - 1; n++)
		lift[n] = 1 + ((double)o->fea_lifter)/2*sin(3.141592653589793*(n + 1.)/((double)o->fea_lifter));

	// HTK compatibility coefficient
	normcoef = sqrt(2.0/Nin);

};

bool dctcFEA::process_frame() {
    Vec<double> & X = *_Xs;
    Vec<double> & c = *_fvec;

    // take logarithm of input spectrum
    for (int i = 0; i < Nin; i++) X[i] = log(X[i]);

    // DCT projection
    for (int i = 0; i < Nout; i++) {
		c[i] = 0;
		for (int k = 1; k <= Nin; k++) { // indexation from 1
			int id = (2*k - 1)*i % (4*Nin);
			c[i] += X[k - 1]*wdct[id];
			// c[i] = sum_k=1^N [ Xk * w_(i*(2k-1)) ]
			//      = sum..     [ Xk * cos(2*pi*i*(2k-1)/4N) ]
			//      = sum..     [ Kx * cos(pi*i*(k-0.5)/N) ]
		}
		c[i] *= normcoef; // for HTK compatibility (ck = sqrt(2/N) * sum ...)
    }

    // liftering
    if (o->fea_lifter > 1) {
      for (int n = 1; n < Nout; n++)
         c[n] *= lift[n - 1];
    }

    return true;
};

dctcFEA::~dctcFEA() {
	delete _fvec;
	delete [] wdct;
	delete [] lift;
};

// ---------------- linear prediction coefficients a[k] -----------------

lpaFEA::lpaFEA(opts* o, Vec<double> *_Xs) : _FEA(o,_Xs) {
	lporder = o->fea_lporder;
	Nout = lporder + 1; // plus energy
	Ninfull = (Nin - 1)*2;
	try {
       	_fvec = new Vec<double>(Nout);
		RRe = new double [lporder + 1];
		rc  = new double [lporder + 1];
		a   = new double [lporder + 1];
		aa  = new double [lporder + 1];
		WRe = new double [Ninfull];
		P   = new double [lporder + 1];
    }
    catch (...) { throw ("FEA: Not enough memory!"); }

   	// iDFT coefficients init:
	for (int i = 0; i < Ninfull; i++) WRe[i] = cos(2*3.141592653589793*i/Ninfull);
	// WIm = sin(...), sin = shifted cos
	//	for (int i = Xsize/4; i < Xsize; i++)	WIm[i] = WRe[i - Xsize/4];
	//	for (int i = 0; i < Xsize/4; i++) WIm[i] = WRe[i + 3*Xsize/4];
};

bool lpaFEA::process_frame() {
	// taking power only when IN-LD not used (power: MFCC=yes, PLP=no)
	if (!o->fb_inld) {
		Vec<double> & X = *_Xs;
		for (int i=0; i<Nin; i++)
			X[i] *= X[i];
	}

	idft();
	LevDurb();

	Vec<double> & fea = *_fvec;
	for (int i=0; i<Nout; i++) fea[i] = a[i];

	E = log(RRe[0]);
	return true;
};

void lpaFEA::idft() {
	Vec<double> & X = *_Xs;

	// --------- iDFT  ----------  ( Power Spect -> Autocorr ) ---------------------
	// R[k] = 1/N * sum_n=0:N [ Xn * cos (2*pi*n*k/N) ] N=full size spectrum (2*(Nin-1))
	// summation over full spect === sum over half with two times taken values for n=1:Nin-1
	// speedup: adding X0 and XNin separately plus 2 times the left half
	int id = 0;
	for (int k = 0; k <= lporder; k++) {
		RRe[k] = X[0]/2.;                       // taking only half - will be compensated at the end
		for (int n = 1; n < Nin - 1; n++) {       // sum over half a spectrum
			id = (n*k) % Ninfull;
			RRe[k] += X[n] * WRe[id];              // should be *2, will be done at the end
		}
		RRe[k] += (1 - 2*(k % 2)) * X[Nin - 1]/2.; // also only half....
		RRe[k] /= double(Ninfull)/2;                // iDFT = 1/N * sum(...), we compensate -> 2/N
	}
};

void lpaFEA::LevDurb() {
	// ---------- Levinson-Durbin ---- (Autocorr -> LP coeffs) ---------------------

	P[0]  = RRe[0];
	rc[1] = -RRe[1]/RRe[0];
	P[1]  = P[0]*(1 - rc[1]*rc[1]);
	aa[1] = rc[1];
	a[0]  = aa[0] = 1;              // first LPCoefficient == 1
	int n = 0;
	for (int ik = 2; ik <= lporder; ik++) {
		double dm = RRe[ik];
		for (n = 1; n <= ik - 1; n++) dm += aa[n]*RRe[ik - n];
		rc[ik] = - dm/P[ik-1];        // last reflection coefficient

		a[ik] = rc[ik];               // last LPC == last r.c.
		// Levinson's recursion
		for (n = 1; n <= ik - 1; n++) a[n] = aa[n] + rc[ik]*aa[ik-n];
		// a -> aa
		for (n = 1; n <= ik; n++) aa[n] = a[n];

		P[ik] = P[ik-1]*(1 - rc[ik]*rc[ik]);
	}
};

lpaFEA::~lpaFEA() {
	delete [] P;
	delete [] WRe;
	delete [] aa;
	delete [] a;
	delete [] rc;
	delete [] RRe;
	delete _fvec;
};


// ---------------- linear prediction cepstral coefficients c[k] -----------------

lpcFEA::lpcFEA(opts* o, Vec<double> *_Xs) : lpaFEA(o,_Xs) {
	Nout 	= o->fea_ncepcoefs + 1; // plus c0
	try { 
		delete _fvec;
       	_fvec = new Vec<double>(Nout);
		lift= new double [Nout - 1];
    }
    catch (...) { throw ("FEA: Not enough memory!"); }

	// lifter coefficients init
	for (int n = 0; n < Nout - 1; n++)
		lift[n] = 1 + ((double)o->fea_lifter)/2*sin(3.141592653589793*(n + 1.)/((double)o->fea_lifter));
};

bool lpcFEA::process_frame() {

	lpaFEA::process_frame();
	a2c();

	// liftering
	if (o->fea_lifter > 1) {
		Vec<double> & c = *_fvec;
		for (int n = 1; n < Nout; n++)
			c[n] *= lift[n - 1];
	}

	return true;
};

void lpcFEA::a2c() {
	Vec<double> & c = *_fvec;

	int k = 0;
	double sum;
	c[0] = log(P[lporder]);
	// " c0 = 2*log( sqrt(alpha) )" ; alpha = Pred.error = P = gain^2 = G^2
	for (int n = 1; n < Nout; n++) {
		sum = 0;
		if (n <= lporder) {
			for (k = 1; k <= n - 1; k++) sum += (n - k)*c[n - k]*a[k];
			c[n] = -a[n] - sum/n;
		}
		else {
			for (k = 1; k <= lporder; k++) sum += (n - k)*c[n - k]*a[k];
			c[n] = -sum/n;
		}
	}
};

lpcFEA::~lpcFEA() {
	delete [] lift;
	delete _fvec;
	_fvec = NULL;
};
