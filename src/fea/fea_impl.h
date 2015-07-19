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

#ifndef FEA_IMPL_H_
#define FEA_IMPL_H_

#include "../io/stdafx.h"
#include "../io/opts.h"
#include "../base/types.h"
#include <math.h>

class _FEA {
    // prototype class inherited by all subclasses

  protected:
	opts *o;
	Vec<double> * _Xs;
	Vec<double> * _fvec;
	int Nin, Nout;

  public:
	double E;

	_FEA(opts *, Vec<double> *);
	virtual bool process_frame() {return 0;};
	virtual bool flush_frame() {return 0;};
	virtual bool new_file() {return 0;};
        virtual Vec<double>* get_fvec() {return _fvec;}
	virtual ~_FEA() {};
};

// ---------------- dummy class for spectrum -------------
class specFEA : public _FEA {
	public:
		specFEA (opts *, Vec<double>*);
		virtual ~specFEA ();
		virtual bool process_frame();
};

// ---------------- class for log spectrum  -------------
class logspecFEA : public specFEA {
	public:
		logspecFEA(opts *o, Vec<double>*X):specFEA(o,X) {};
		virtual bool process_frame();
};

// ---------------- cepstral coeffs via DCT -----------------
class dctcFEA : public _FEA {
	protected:
		double * wdct;
		double * lift;
		double normcoef;
	public:
		dctcFEA (opts *, Vec<double>*);
		virtual ~dctcFEA ();
		virtual bool process_frame();
};

// ---------------- linear prediction coefficients a[k] -----------------
class lpaFEA : public _FEA {
	protected:
		int lporder, Ninfull; // size of full symm spec
		double *RRe, *rc, *aa, *WRe, *P, *a;
	public:
		lpaFEA (opts *, Vec<double>*);
		virtual ~lpaFEA ();
		virtual bool process_frame();
		void idft();
		void LevDurb();
};

// ---------------- linear prediction cepstral coefficients c[k] -----------------
class lpcFEA : public lpaFEA {
		double * lift;
	public:
		lpcFEA (opts *, Vec<double>*);
		virtual ~lpcFEA ();
		virtual bool process_frame();
		void a2c();
};


#endif /*FEA_IMPL_H_*/

