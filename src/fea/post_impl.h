
/*
 CtuCopy 3 - Universal feature extractor and speech enhancer.
 Copyright 2014 Petr Mizera, FEE CTU Prague

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

#ifndef POST_IMPL_H_
#define POST_IMPL_H_

#include "../io/stdafx.h"
#include "../io/opts.h"
#include "../base/types.h"
#include <math.h>

class _POST {
  // prototype class inherited by all subclasses
  protected:
        opts *o;
        Vec<double> * _fvec;
        Vec<double> ** _mvec_cmean;
        Vec<double> ** _mvec_cvar;
        char ** _list_ID_spk;
        virtual bool new_spk_1(char *ID_spk){return 0;};
  public:
        _POST(opts *, Vec<double> *);
        _POST(opts *, Vec<double> *, Vec<double> **, Vec<double> **, char **);
        virtual bool sum_fea() {return 0;};
        virtual bool sum_cv() {return 0;};
        virtual bool stat_cm() {return 0;};
        virtual bool stat_cv() {return 0;};
        virtual void clean(){};
        virtual bool process_frame() {return 0;};
        virtual int add_spk(char *id_spk) {return 0;};

        virtual Vec<double>**  get_fvec_cmean() {return _mvec_cmean;}
        virtual Vec<double>**  get_fvec_cvar()  {return _mvec_cvar;}
        virtual Vec<double>*   get_fvec()       {return _fvec;}
        virtual char **        get_list_ID_spk(){return _list_ID_spk;}

        virtual ~_POST() {};
};

class cmvn_POST : public _POST{

    int n_frame;
    int nfea;                    // input vector size
    int num_spk;
    int id_spk;
    virtual bool new_spk_1(char *ID_spk);

public:
    cmvn_POST(opts* o, Vec<double> *);
    cmvn_POST(opts* o, Vec<double> *, Vec<double> **, Vec<double> **,char **);
    virtual bool sum_fea();
    virtual bool sum_cv();
    virtual bool stat_cm();
    virtual bool stat_cv();
    virtual bool process_frame();
    virtual int add_spk(char *id_spk);

    virtual ~cmvn_POST();
};

class cms_POST : public _POST{

    int nfea;                    // input vector size
    string cms_type;
    int x_cbuffer;
    float * sumM;
    float ** Z_obufferr;
    int num_frame;

public:
    cms_POST(opts* o, Vec<double> *);
    virtual bool process_frame();
    virtual void clean();
    virtual ~cms_POST();
};


#endif /*POST_IMPL_H_*/
