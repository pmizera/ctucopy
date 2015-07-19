#ifndef DELTA_H_
#define DELTA_H_

//#include <math.h>

#include "../io/stdafx.h"
//#include "fea_impl.h"
#include "../base/types.h"
#include <fftw3.h>
#include "../io/opts.h"

class deltaFEA{
    opts  *o;               // options class
    int delta_w;                 // Delta window size
    int dwlen;                   // dwlen=2*delta_w+1;
    int start;                   // circular buffer start ptr
    int end;                     // index of end frame in circular buffer
    int nfea;                    // input vector size
    int frame;
    bool endframe;               // create artificial context
    int fea_c;

    //inicialization
    int avail;                   // number of buffered frames (avail=delta_w+1;)
    int index;                   // for initialization circular buffet ( index=delta_w+1 )
    int num_c;                   // number of cycles for repeat one frames
    int den;
    int n_order;

    virtual bool delta();
    virtual bool trap();
    virtual bool init_cbuffer();

protected:
    Vec<double> *in_fvec;         // output vector
    Mat<double> *buf;           // buffer for input frames
    Mat<double> *buf_d;           // buffer for delta

public:
    Vec<double> *_fvec;         // output vector

    deltaFEA(opts *o, Vec<double>*X,int _n_order,int _delta_w);
    virtual bool new_file();
    virtual bool process_frame();
    virtual bool flush_frame();
    virtual ~deltaFEA();
};

#endif /*DELTA_H_*/
