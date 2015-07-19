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

#ifndef FEA_H_
#define FEA_H_

#include "../io/stdafx.h"
#include "../io/opts.h"
#include "../base/types.h"
#include "math.h"
#include <stdio.h>
#include "fea_impl.h"
#include "fea_trap.h"
#include "fea_delta.h"

class FEA {
    // this is the only class intended to be public
    // it implements all needed methods using subclasses <x>FEA

    // !!! Methods may DESTROY input vectors Xabs, Xph !!!!!!

  protected:
  	opts* o;
	_FEA* fea;
	Vec<double> * _Xs;

  public:
	Vec<double> * _fvec;
	double * E;

	FEA(opts *, Vec<double> *);
	bool new_file();
	bool process_frame();
	bool flush_frame();
	~FEA();
};

#endif /*FEA_H_*/

