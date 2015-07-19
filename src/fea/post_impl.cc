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

#include "post_impl.h"

_POST::_POST(opts* o_in, Vec<double> *fvec){
        o = o_in;
        _fvec = fvec;
};

_POST::_POST(opts* o_in, Vec<double> * fvec, Vec<double> ** m_cmean, Vec<double> ** m_cvar, char ** list_ID_spk){
        o = o_in;
        _fvec = fvec;
};

cmvn_POST::cmvn_POST(opts* o, Vec<double> * _fvec): _POST(o,_fvec){
    nfea    = _fvec->get_size()+1;
    num_spk =0;
    id_spk =0;
    try {
         _mvec_cmean = new Vec<double>*[5000];
         _mvec_cvar  = new Vec<double>*[5000];
         _list_ID_spk = new char*[5000];
    }
    catch (...) { throw ("FEA: Not enough memory!"); }
};

cmvn_POST::cmvn_POST(opts* o, Vec<double> * _fvec, Vec<double> ** _m_cmean, Vec<double> ** _m_cvar, char ** _l_ID_spk): _POST(o, _fvec, _mvec_cmean, _mvec_cvar, _list_ID_spk){
        nfea    = _fvec->get_size()+1;
        _mvec_cmean  = _m_cmean;
        _mvec_cvar   = _m_cvar;
        _list_ID_spk = _l_ID_spk;
        num_spk = sizeof (_list_ID_spk)/sizeof(_list_ID_spk[0]);
        id_spk=0;
};

bool cmvn_POST::sum_fea(){
    Vec<double> & F = *_fvec;
    Vec<double> & M_cmean = *_mvec_cmean[id_spk];

    if(0==strcmp(o->format_in, "htk")){
      for (int i = 0; i < nfea-1; i++)
        M_cmean[i] += F[i];
    }else{
      for (int i = 1; i < nfea-1; i++)
        M_cmean[i-1] += F[i];
      M_cmean[nfea-2] += F[0];
    }

    M_cmean[nfea-1]++;
    return true;
};

bool cmvn_POST::stat_cm(){
    for(int j=0; j<num_spk; j++){
      Vec<double> & M_cmean = *_mvec_cmean[j];
      for (int i = 0; i < nfea-1; i++)
        M_cmean[i] /= M_cmean[nfea-1];
    }
    return true;
};

bool cmvn_POST::sum_cv(){
    Vec<double> & F = *_fvec;
    Vec<double> & M_cmean = *_mvec_cmean[id_spk];
    Vec<double> & M_cvar  = *_mvec_cvar[id_spk];

    if(0==strcmp(o->format_in, "htk")){
      for (int i = 0; i < nfea-1; i++)
        M_cvar[i] += ( F[i] - M_cmean[i] ) * ( F[i] - M_cmean[i] );
    }else{
      for (int i = 1; i < nfea-1; i++)
        M_cvar[i-1] += ( F[i] - M_cmean[i-1] ) * ( F[i] - M_cmean[i-1] );
      M_cvar[nfea-2] += ( F[0] - M_cmean[nfea-2] ) * ( F[0] - M_cmean[nfea-2] );
    }
    M_cvar[nfea-1]++;
    return true;
};

bool cmvn_POST::stat_cv(){
    for(int j=0; j<num_spk; j++){
      Vec<double> & M_cvar = *_mvec_cvar[j];
      for (int i = 0; i < nfea-1; i++)
        M_cvar[i] /= (M_cvar[nfea-1]-1);
   }

   return true;
};

bool cmvn_POST::cmvn(){
    Vec<double> & F = *_fvec;
    Vec<double> & M_cmean = *_mvec_cmean[id_spk];
    Vec<double> & M_cvar  = *_mvec_cvar[id_spk];

    if(0==strcmp(o->format_in,"htk")){
      for (int i = 0; i < nfea-1; i++)
        F[i] =  (F[i] - M_cmean[i]) / M_cvar[i];
    }else{
      for (int i = 1; i < nfea-1; i++)
        F[i] =  (F[i] - M_cmean[i-1]) / M_cvar[i-1];
      F[0] =  (F[0] - M_cmean[nfea-2]) / M_cvar[nfea-2];
    }
    return true;
}

int cmvn_POST::add_spk(char *ID_spk){
    int new_spk=0;
    if(!num_spk){                         //add the first speaker
      new_spk_1(ID_spk);
      num_spk++;
    }
    for(int ii=0;ii<num_spk;ii++){        //add new speaker
       if(strcmp(_list_ID_spk[ii],ID_spk)){
         new_spk=1;
         id_spk=ii;
       }else{
         id_spk=ii;
         new_spk=0;
         break;
       }
    }
    if(new_spk){
      id_spk=num_spk;
      new_spk_1(ID_spk);
      num_spk++;
    }
 return num_spk;
}

bool cmvn_POST::new_spk_1(char *ID_spk){
    try {
        _mvec_cmean[id_spk]  = new Vec<double>(nfea);
        _mvec_cvar[id_spk]   = new Vec<double>(nfea);
        _list_ID_spk[id_spk] = new char[(int)strlen(ID_spk)];
    }
    catch (...) { throw ("FEA: Not enough memory!"); }
    strcpy(_list_ID_spk[id_spk],ID_spk);
    return true;
};

cmvn_POST::~cmvn_POST() {
    delete _mvec_cmean;
    delete _mvec_cvar;
};
