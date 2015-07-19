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

#include "post.h"

POST::POST (opts* o_in, Vec<double> * fvec) {
        o = o_in;
        _fvec = fvec;

        if(o->stat_cmvn||o->apply_cmvn)  post = new cmvn_POST(o,fvec);
        else throw ("POST: Unknown post-processing!");

        _mvec_cmean  = post->get_fvec_cmean();
        _mvec_cvar   = post->get_fvec_cvar();
        _list_ID_spk = post->get_list_ID_spk();
        _fvec  = post->get_fvec();
};

POST::POST (opts* o_in, Vec<double> * fvec, Vec<double> ** m_cmean, Vec<double> ** m_cvar,char ** l_ID_spk) {
        o = o_in;
        _fvec = fvec;
        _mvec_cmean  = m_cmean;
        _mvec_cvar   = m_cvar;
        _list_ID_spk = l_ID_spk;

        if(o->apply_cmvn)  post = new cmvn_POST(o,fvec,m_cmean,m_cvar,l_ID_spk);
        else throw ("POST: Unknown post-processing!");

        _fvec  = post->get_fvec();
};

bool POST::sum_fea(){
        return post->sum_fea();
}

bool POST::sum_cv(){
        return post->sum_cv();
}

bool POST::stat_cm(){
        return post->stat_cm();
}

bool POST::stat_cv(){
        return post->stat_cv();
}

bool POST::cmvn(){
        return post->cmvn();
}

int POST::add_spk(char* id_spk){
        return post->add_spk(id_spk);
}

POST::~POST() {
        delete post;
}
