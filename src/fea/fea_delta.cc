/*
 CtuCopy 4 - Universal feature extractor and speech enhancer.
 Copyright 2015 Petr Mizera, FEE CTU Prague

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

#include "fea_delta.h"

deltaFEA::deltaFEA(opts* o_in, Vec<double> *_in_fvec,int _n_order,int _delta_w){
    o = o_in;    // pointer to options
    in_fvec=_in_fvec;
    fea_c = o->fea_ncepcoefs+1;

    delta_w = _delta_w;
    dwlen   = 2*delta_w+1;
    n_order = _n_order;
    nfea    = fea_c*(n_order-1);

    start   = 0;
    end     = 0;
    endframe = true;
    frame   =1;

    //inicialization
    avail  = delta_w+1;
    index  = avail-1;
    num_c  = 1;
    den    = 0;
    if(o->fea_trap)n_order =dwlen;

    if(delta_w < 1){
       if(o->fea_trap) throw("FEA: Trap window size must be >= 3!");
       else throw("FEA: Delta window size must be > 1!");
    }
    try {
        _fvec   = new Vec<double>(fea_c*n_order);
        buf     = new Mat<double>(dwlen, nfea);
        buf_d   = new Mat<double>(dwlen, nfea);
    }
    catch (...) { throw ("FEA: Not enough memory!"); }

   for(int i=1;i<=delta_w;i++){
      den += i*i;
   }
   den=2*den;
};

bool deltaFEA::new_file(){
  start   = 0;
  end     = 0;
  endframe = true;
  frame   =1;

  avail  = delta_w+1;
  index  = avail-1;
  num_c  = 1;
  return false;
}
bool deltaFEA::process_frame() {
    Vec<double> & ofvec = *_fvec;
    Vec<double> & fvec = *in_fvec;
    Mat<double> &b = *buf;

    if(avail!=0){
      if(frame==1)      num_c = delta_w;    // repeat the first frame
      else if(frame==2) num_c = 2;          // repeat the second frame
      else              num_c = 1;          // buffering of frames

      init_cbuffer();
      avail--;
      if(!avail){
        if(o->fea_trap)trap();
        else delta();
        for(int j=0; j<nfea; j++){
           b(dwlen-1,j)=b(dwlen-2,j);
        }
        for (int i=0;i<nfea;i++){
           ofvec[i] = b(delta_w,i);
        }
        start++;
        start%=dwlen; // move on
        frame++;
        return true;
      }
      start++;
      start%=dwlen; // move on
      frame++;
      return false;
    }

    // copy fea vector to buffer
    int aa=start;
    if(aa==0)aa=dwlen-1;
    else aa--;

    aa%=dwlen;
    for(int i=0;i<nfea;i++){
        b(aa,i) = fvec[i];
    }

    if(o->fea_trap) trap();
    else delta();     // calculate delta Coefficients

    int qqq=aa-delta_w;
    if(qqq<0)
      qqq=dwlen+qqq;

    qqq%=dwlen;
    if(!o->fea_trap){
      for (int i=0;i<nfea;i++){
        ofvec[i] = b(qqq,i);
       }
    }
    end=start;
    start++;start%=dwlen; // move on

    frame++;
    return true;
};

bool deltaFEA::init_cbuffer() {
    Vec<double> & f = *in_fvec;
    Mat<double> &b = *buf;

    for(int ii=0; ii<num_c; ii++){
      for(int j=0; j<nfea; j++){
         b(index,j) = f[j];
      }
    index++;
    index %=dwlen;
    }
    return true;
};

bool deltaFEA::delta() {
    Vec<double> & fvec = *_fvec;
    Mat<double> &b     = *buf;
    Mat<double> &d     = *buf_d;

    for(int i=0; i<dwlen; i++){
      for(int j=0; j<nfea; j++){
        d(i,j) = b((start+i)%dwlen,j);
      }
    }
    for(int j=fea_c*(n_order-2);j<fea_c*(n_order-1);j++){
     double x=0;
     for(int i=1;i<=delta_w;i++){
        x +=  i*(d(delta_w+i,j) - d(delta_w-i,j));
     }
     fvec[j+fea_c] = x/den;
    }
    return true;
};

bool deltaFEA::trap() {
    Vec<double> & fvec = *_fvec;
    Mat<double> &b     = *buf;

    for(int i=0; i<fea_c; i++){
      for(int j=0; j<dwlen; j++){
        fvec[i*dwlen+j] = b((start+j)%dwlen,i);
  //      cout<<" "<< b((start+j)%dwlen,i);
  //      cout<<" "<< (start+j)%dwlen;
      }
    }
//cout<<endl;

    return true;
};

bool deltaFEA::flush_frame() {
    Vec<double> & ofvec = *_fvec;
    Mat<double> &b = *buf;

    if(endframe){
      index=end;
      endframe = false;
    }

    if(avail < delta_w){
      init_cbuffer();

      if(o->fea_trap) trap();
      else delta();               // calculate delta Coefficients

      int qqq=start-delta_w-1;
      if(qqq<0)
        qqq=dwlen+qqq;

      qqq%=dwlen;
      for (int i=0;i<nfea;i++){
          ofvec[i] = b(qqq,i);
      }
      start++; start%=dwlen; // move on
      avail++;
      return true;
    }
    return false;
};


deltaFEA::~deltaFEA() {
    delete _fvec;
    delete buf;
    delete buf_d;
};
