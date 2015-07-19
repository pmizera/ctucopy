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

#include "fea.h"

FEA::FEA (opts* o_in, Vec<double> * X) {
	o = o_in;
	_Xs = X;

	if (0==strcmp(o->fea_kind,"spec"))         fea = new specFEA(o,X);
	else if (0==strcmp(o->fea_kind,"logspec")) fea = new logspecFEA(o,X);
        else if (0==strcmp(o->fea_kind,"dctc"))    fea = new dctcFEA(o,X);
	else if (0==strcmp(o->fea_kind,"lpa"))     fea = new lpaFEA(o,X);
	else if (0==strcmp(o->fea_kind,"lpc"))     fea = new lpcFEA(o,X);
	else if (0==strcmp(o->fea_kind,"trapdct")) fea = new trapdctFEA(o,X);
	else throw ("FEA: Unknown feature kind!");

        _fvec  = fea->get_fvec();
	E = &fea->E;
};

bool FEA::new_file() {
	return fea->new_file();
}

bool FEA::process_frame() {
	return fea->process_frame();
}

bool FEA::flush_frame() {
	return fea->flush_frame();
}

FEA::~FEA() {
	delete fea;
}

