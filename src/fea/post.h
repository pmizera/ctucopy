/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Mizera, FEE CTU Prague

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

#ifndef POST_H_
#define POST_H_

#include "../io/stdafx.h"
#include "../io/opts.h"
#include "../base/types.h"
#include "math.h"
#include <stdio.h>
#include "post_impl.h"

class POST {

  protected:
        opts* o;
        _POST* post;
        Vec<double> * _fvec;

  public:
        Vec<double> ** _mvec_cmean;
        Vec<double> ** _mvec_cvar;
        char ** _list_ID_spk;

        POST(opts *, Vec<double> *);
        POST(opts *, Vec<double> *, Vec<double> **, Vec<double> **, char **);
        bool sum_fea();
        bool sum_cv();
        bool stat_cm();
        bool stat_cv();
        bool cmvn();
        int add_spk(char* id_spk);
        ~POST();
};

#endif /*POST_H_*/

