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

#include "fea_trap.h"

trapdctFEA::trapdctFEA(opts* o, Vec<double> *_Xs) : logspecFEA(o,_Xs) {
    traplen = o->fea_trapdct_traplen;
    ndct    = o->fea_trapdct_ndct;
    nfea    = _Xs->get_size();
    firstframe = true;
    htraplen= (traplen+1)/2;
    avail   = 0;
    start   = 0;

    if (traplen%2==0)
        throw("FEA: TRAP length must be odd!");
    if (ndct>=traplen)
        throw("FEA: Number of DCT coeffs must be less than TRAP length (c0 is not output)!");

    try {
        _fvec = new Vec<double>(ndct*nfea);
        buf   = new Mat<double>(traplen,nfea);
        hamm  = new double[traplen];
    }
    catch (...) { throw ("FEA: Not enough memory!"); }

	// init Hamming window
	for (int i=0; i<traplen;i++)
		hamm[i] = 0.54 - (1-0.54) * cos(2*3.14159265359*i/(traplen - 1.));

    // FFTW init
	fftw_in  = (double*) fftw_malloc(sizeof(double) * traplen);
	fftw_out = (double*) fftw_malloc(sizeof(double) * traplen);
	if (fftw_in == 0 || fftw_out == 0) throw "FEA: Out of memory!";
	plan = fftw_plan_r2r_1d(traplen, fftw_in, fftw_out, FFTW_REDFT10, FFTW_MEASURE);
};


bool trapdctFEA::process_frame() {

    logspecFEA::process_frame();

    static Vec<double> & Xin = *logspecFEA::_fvec;
    static Mat<double> &b = *buf;

    // copy fea vector to buffer
    for (int i=0;i<nfea;i++)
        b(start,i) = Xin[i];

    if (firstframe) {
        // fill the left context with the current input
        for (int i=1;i<htraplen;i++)
            for (int j=0;j<nfea;j++)
                b((start-i+traplen)%traplen,j) = Xin[j];
        firstframe = false;
    }

    start++; start%=traplen; // move on

    if (avail<htraplen-1) {
        avail++;
        return false;
    }

    core();    // calculate TRAP-DCTs

    return true;
};

void trapdctFEA::core() {
    static Vec<double> & fvec = *_fvec;
    static Mat<double> &b = *buf;
    for (int i=0;i<nfea;i++) {
        double sum=0;

        // form TRAP
        for (int j=0;j<traplen;j++) {
            fftw_in[j] = b((start+j)%traplen,i);
            sum += fftw_in[j];
        }

        sum /= traplen;

        // remove mean, apply Hamming window
        for (int j=0;j<traplen;j++)
            fftw_in[j] = (fftw_in[j]-sum) * hamm[j];

        // DCT, store without the zeroth coeff
        fftw_execute(plan);
        for (int j=0;j<ndct;j++) {
            fvec[i*ndct+j] = fftw_out[j+1];
        }
    }
};


bool trapdctFEA::flush_frame() {

    static Mat<double> &b = *buf;
    if (avail) {
        // replicate last frame
        for (int i=0;i<nfea;i++)
            b(start,i) = b((start-1+traplen)%traplen,i);

        start++; start%=traplen; // move on
        avail--;

        core();    // calculate TRAP-DCTs

        return true;
    }
    return false;
};


trapdctFEA::~trapdctFEA() {
	fftw_destroy_plan(plan);
	fftw_free(fftw_in);
	fftw_free(fftw_out);
    delete [] hamm;
    delete _fvec;
    delete buf;
};
