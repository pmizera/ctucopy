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

#ifndef TRAP_H_
#define TRAP_H_

#include "../io/stdafx.h"
#include "fea_impl.h"
#include <fftw3.h>


class trapdctFEA : public logspecFEA {
    int traplen, ndct, htraplen; // htraplen = (traplen+1)/2
    int start;                   // circular buffer start ptr
    int avail;                   // number of buffered frames
    int nfea;                    // input vector size
    bool firstframe;             // upon process_frame call the frame is repeated to create artificial context
    Mat<double> * buf;           // buffer for input frames
    double *hamm;                // Hamming window
  	double *fftw_in, *fftw_out;
	fftw_plan plan;
    Vec<double> * _fvec;         // output vector

    virtual void core();

public:
    trapdctFEA(opts *o, Vec<double>*X);
    virtual bool process_frame();
    virtual bool flush_frame();
    virtual bool new_file() {firstframe=true; start=0; return 0;};
    virtual Vec<double>* get_fvec() {return _fvec;}
    virtual ~trapdctFEA();
};


#endif /*TRAP_H_*/

