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

#include <math.h>
#include <fftw3.h>
#include <limits>
#include <linux/limits.h>

#include "../io/stdafx.h"

#include "vad.h"

#include "../vdet/CepstralDet.h"
using namespace Voice;

#define  MAX(a,b) (a>b ? a : b)
#define  MIN(a,b) (a<b ? a : b)


#define  CRI_MODE_ENERGY(opt) !strcmp(opt, "energy")
#define  CRI_MODE_CEPDIST(opt) !strcmp(opt, "cepdist")

#define  THR_MODE_ABSOLUTE(opt) !strcmp(opt, "absolute")
#define  THR_MODE_PERC(opt) !strcmp(opt, "perc")
#define  THR_MODE_ADAPT(opt) !strcmp(opt, "adapt")
#define  THR_MODE_DYN(opt) !strcmp(opt, "dyn")

#define  APPLY_MODE_NONE(opt) !strcmp(opt, "none")
#define  APPLY_MODE_SILENCE(opt) !strcmp(opt, "silence")
#define  APPLY_MODE_DROP(opt) !strcmp(opt, "drop")

#define  OUT_MODE_NONE(opt) !strcmp(opt, "none")
#define  OUT_MODE_VAD(opt) !strcmp(opt, "vad")
#define  OUT_MODE_DEBUG(opt) !strcmp(opt, "debug")

#define  CEPDIST_MODE_LPC(opt) !strcmp(opt, "lpc")
#define  CEPDIST_MODE_FEA(opt) !strcmp(opt, "fea")
#define  CEPDIST_MODE_IN(opt) !strcmp(opt, "in")

inline double DB(double x)
{
    return 10.0d*log10(std::numeric_limits<double>::min()+x);
}

// -----------------------------------------------------------------------------
// ---------------- energetic VAD criterial ------------------------------------
class VADcri_energy : public VADcri {

private:        
                double energy;
                FileWriter *file_energy;
                bool opt_db;
                int size;

public:
                VADcri_energy(const opts*, Vec<double>*, Vec<double>*, Vec<double>*, Vec<double>*);
                virtual ~VADcri_energy();

                virtual void new_file(const char* filename);
                virtual double process_frame(int frame_index);
                virtual void save_frame(int frame_index);
                virtual void consume_vad(int frame_index, bool vad);
};

VADcri_energy::VADcri_energy(const opts* options, Vec<double> *Xsabs, Vec<double> *Xsph, Vec<double>* InFvec, Vec<double>* FeaFvec)
: VADcri(options, Xsabs, Xsph, InFvec, FeaFvec),
  opt_db(options->vad_energy_db)
{
    file_energy = NULL; 
    size = xsabsvec->get_size();
};

VADcri_energy::~VADcri_energy() {
    DELETE(file_energy);
};

void VADcri_energy::new_file(const char* filename) {
    DELETE(file_energy);
    file_energy = new FileWriter(filename, "energy");
}

double VADcri_energy::process_frame(int frame_index) {
    Vec<double>& X = *xsabsvec;

    energy = 0.0d;
    for(int i=0; i<size; i++)
        energy += X[i]*X[i];
    
    if(opt_db)
        energy = DB(energy);

    return energy;
};

void VADcri_energy::save_frame(int frame_index) {
    file_energy->write(energy);
}

void VADcri_energy::consume_vad(int frame_index, bool vad)
{ return; };

// -----------------------------------------------------------------------------
// ---------------- cepstral distance criterial class --------------------------
class VADcri_cepdist : public VADcri {
    private:
        const char *opt_mode;
        const double opt_p;
        const int opt_init;
        const int opt_lpc_coefs;
        const int opt_window; // segment length
        const int opt_wfft, opt_wfftby2;

        double *c0;
        double *ci;
        int size;
                
        double dist;
        FileWriter *file_cepdist;
        FileWriter *file_c0init;

        BurgCepstrumEstimator *lpc;
        double *fft_in, *fft_out;
        fftw_plan plan;

    public:
        VADcri_cepdist(const opts*, Vec<double>*, Vec<double>*, Vec<double>*, Vec<double> *);
        virtual ~VADcri_cepdist();

        virtual void new_file(const char* filename);
        virtual double process_frame(int frame_index);
        virtual void save_frame(int frame_index);
        virtual void consume_vad(int frame_index, bool vad);
};

VADcri_cepdist::VADcri_cepdist(const opts* options, Vec<double> *Xsabs, Vec<double> *Xsph, Vec<double> *InFvec, Vec<double> *FeaFvec)
: VADcri(options, Xsabs, Xsph, InFvec, FeaFvec),
  opt_mode(options->vad_cepdist_mode),
  opt_p(options->vad_cepdist_p),
  opt_init(options->vad_cepdist_init),
  opt_lpc_coefs(options->vad_lpc_coefs),
  opt_window(options->window),
  opt_wfft(options->wfft),
  opt_wfftby2(options->wfftby2)
{
    file_cepdist = NULL;
    file_c0init = NULL;

    if(CEPDIST_MODE_LPC(opt_mode))
    {
        size = opt_lpc_coefs;
        lpc = new BurgCepstrumEstimator(opt_window, size);

        if((xsphvec == NULL) || (xsabsvec->get_size() != xsphvec->get_size()))
            throw("VADcri_cepdist: cannot perform iFFT!");

        fft_in  = (double*) fftw_malloc(sizeof(double) * opt_wfft);
        fft_out = (double*) fftw_malloc(sizeof(double) * opt_wfft);

        if(fft_in == 0 || fft_out == 0)
            throw "VADcri_cepdist: not enough memory!";

        plan = fftw_plan_r2r_1d(opt_wfft, fft_in, fft_out, FFTW_HC2R, FFTW_MEASURE);
    }
    else if(CEPDIST_MODE_FEA(opt_mode))
    {
        if(!feafvec)
            throw("VADcri_cepdist: vad_cepdist_mode=fea and fea->_fvec is not available!");
        size = feafvec->get_size();
    }
    else if(CEPDIST_MODE_IN(opt_mode))
    {
        if(!infvec)
            throw("VADcri_cepdist: vad_cepdist_mode=in and in->_fvec is not available!");
        size = infvec->get_size();
    }
    else
        throw("VADcri_cepdist: unknown vad_cepdist_mode!");

    c0 = new double[size];
    ci = new double[size];
};

VADcri_cepdist::~VADcri_cepdist()
{
    if(CEPDIST_MODE_LPC(opt_mode))
    {
        fftw_destroy_plan(plan);
        fftw_free(fft_in);
        fftw_free(fft_out);
    }

    DELETE(file_cepdist);
    DELETE(file_c0init);
    delete[] c0;
    delete[] ci;
};

void VADcri_cepdist::new_file(const char* filename)
{
    DELETE(file_cepdist);
    file_cepdist = new FileWriter(filename, "cepdist");
    DELETE(file_c0init);
    file_c0init = new FileWriter(filename, "c0init");
};

double VADcri_cepdist::process_frame(int frame_index)
{
    if(CEPDIST_MODE_LPC(opt_mode))
    {
        Vec<double> &Xa = *xsabsvec;
        Vec<double> &Xp = *xsphvec;

        for (int i=0; i < opt_wfftby2; i++)
            fft_in[i] = Xa[i] * cos(Xp[i]);
        for (int i=opt_wfft-1; i >= opt_wfftby2; i--)
            fft_in[i] = Xa[opt_wfft - i] * sin(Xp[opt_wfft - i]);

        fftw_execute(plan);
        lpc->Process(&fft_out[0], &fft_out[0] + opt_window);

        for(int i=0; i<size; i++)
            ci[i] = (*lpc)[i];
    }
    else if(CEPDIST_MODE_FEA(opt_mode))
    {
        Vec<double>& F = *feafvec;
        for(int i=0; i<size; i++)
            ci[i] = F[i];
    }
    else if(CEPDIST_MODE_IN(opt_mode))
    {
        Vec<double>& F = *infvec;
        for(int i=0; i<size; i++)
            ci[i] = F[i];
    }

    if(frame_index == 0)
    {
        for (int i = 0; i<size; i++)
            c0[i] = ci[i];
        dist = 0.0;
    }
    else if(frame_index == 1)
    {
        for (int i = 0; i<size; i++)
            c0[i] = (c0[i] + ci[i]) / 2.0;

        double sum = 0.0;
        for(int i=1; i<size; i++)
            sum += (ci[i]-c0[i])*(ci[i]-c0[i]);
        dist = 4.3429 * sqrt(2*sum);
    }
    else
    {
        double sum = 0.0;
        for(int i=1; i<size; i++)
            sum += (ci[i]-c0[i])*(ci[i]-c0[i]);
        dist = 4.3429 * sqrt(2*sum);
    }

    return dist;
};

void VADcri_cepdist::save_frame(int frame_index)
{
    file_cepdist->write(dist);

    if(frame_index <= opt_init)
        file_c0init->write(true);
    else
        file_c0init->write(false);
};

void VADcri_cepdist::consume_vad(int frame_index, bool vad) {
    if(vad && (frame_index > opt_init))
        return;

    for (int i = 0; i<size; i++)
      c0[i] = opt_p*c0[i] + (1.0-opt_p)*ci[i];
};

// -----------------------------------------------------------------------------
// ---------------- absolute threshold VAD class -------------------------------
class VADthr_absolute : public VADthr {
       private:
            const double opt_thr;

            FileWriter *file_thr;

       public:
            VADthr_absolute(const opts*);
            virtual ~VADthr_absolute();
            virtual void new_file(const char* filename);
            virtual bool process_frame(int frame_index, double cri);
            virtual void save_frame(int frame_index);
};

VADthr_absolute::VADthr_absolute(const opts* options)
: VADthr(options),
  opt_thr(options->vad_absolute_thr)
{
    file_thr = NULL;
};

VADthr_absolute::~VADthr_absolute()
{
    DELETE(file_thr);
};

void VADthr_absolute::new_file(const char* filename) {
    DELETE(file_thr);
    file_thr = new FileWriter(filename, "thr");
};

bool VADthr_absolute::process_frame(int frame_index, double cri) {
    return (cri >= opt_thr);
};

void VADthr_absolute::save_frame(int frame_index) {
    file_thr->write(opt_thr);
};

// -----------------------------------------------------------------------------
// ---------------- percentual threshold VAD class -----------------------------
class VADthr_perc : public VADthr {
       private:
            const double opt_thr;
            const double opt_init;
            double crimin, crimax;
            double thr;

            FileWriter *file_crimin;
            FileWriter *file_crimax;
            FileWriter *file_thr;

       public:
            VADthr_perc (const opts*);
            virtual ~VADthr_perc();
            virtual void new_file(const char* filename);
            virtual bool process_frame(int frame_index, double cri);
            virtual void save_frame(int frame_index);
};

VADthr_perc::VADthr_perc(const opts* options)
: VADthr(options),
  opt_thr(options->vad_perc_thr),
  opt_init(options->vad_perc_init)
{
    file_crimin = NULL;
    file_crimax = NULL;
    file_thr = NULL;
};

VADthr_perc::~VADthr_perc()
{
    DELETE(file_crimin);
    DELETE(file_crimax);
    DELETE(file_thr);
};

void VADthr_perc::new_file(const char* filename) {
    DELETE(file_crimin);
    file_crimin = new FileWriter(filename, "crimin");
    DELETE(file_crimax);
    file_crimax = new FileWriter(filename, "crimax");
    DELETE(file_thr);
    file_thr = new FileWriter(filename, "thr");
};

bool VADthr_perc::process_frame(int frame_index, double cri) {
    if(frame_index == 0 || frame_index < opt_init)
    {
      crimin = cri;
      crimax = cri;
    }
    else
    {
      crimin = MIN(cri, crimin);
      crimax = MAX(cri, crimax);
    }

    thr = crimin + (opt_thr/100.0)*(crimax - crimin);
    return (cri >= thr);
};

void VADthr_perc::save_frame(int frame_index) {
    file_crimin->write(crimin);
    file_crimax->write(crimax);
    file_thr->write(thr);
};

// -----------------------------------------------------------------------------
// ---------------- adaptive threshold VAD class -------------------------------
class VADthr_adapt : public VADthr {
  private:
    const int opt_init;
    const double opt_q;
    const double opt_za;

    double crimean;
    double crimean2;
    double crivar;
    double thr;
    bool vad;

    FileWriter *file_init;
    FileWriter *file_crimean;
    FileWriter *file_crimean2;
    FileWriter *file_crivar;
    FileWriter *file_thr;

  public:
    VADthr_adapt (const opts* options);
    virtual ~VADthr_adapt ();
        
    virtual void new_file(const char* filename);
    virtual bool process_frame(int frame_index, double cri);
    virtual void save_frame(int frame_index);
};

VADthr_adapt::VADthr_adapt(const opts* options)
: VADthr(options),
  opt_init(options->vad_adapt_init),
  opt_q(options->vad_adapt_q),
  opt_za(options->vad_adapt_za)
{
    file_init = NULL;
    file_crimean = NULL;
    file_crimean2 = NULL;
    file_crivar = NULL;
    file_thr = NULL;
};

VADthr_adapt::~VADthr_adapt() {
    DELETE(file_init);
    DELETE(file_crimean);
    DELETE(file_crimean2);
    DELETE(file_crivar);
    DELETE(file_thr);
};

void VADthr_adapt::new_file(const char* filename) {
    DELETE(file_init);
    file_init = new FileWriter(filename, "init");
    DELETE(file_crimean);
    file_crimean = new FileWriter(filename, "crimean");
    DELETE(file_crimean2);
    file_crimean2 = new FileWriter(filename, "crimean2");
    DELETE(file_crivar);
    file_crivar = new FileWriter(filename, "crivar");
    DELETE(file_thr);
    file_thr = new FileWriter(filename, "thr");
};

bool VADthr_adapt::process_frame(int frame_index, double cri)
{
    if(frame_index == 0)
    {
      thr = cri;
      crimean = cri;
      crimean2 = cri*cri;
      crivar = 0.0;
      vad = false;
    }
    else
    {
        thr = crimean + opt_za * sqrt(crivar);

        if((cri < thr) || (frame_index <= opt_init))
        {
            crimean = opt_q*crimean + (1.0 - opt_q)*cri;
            crimean2 = opt_q*crimean2 + (1.0 - opt_q)*cri*cri;
            crivar = crimean2 - crimean*crimean;
            vad = false;
        }
        else
            vad = true;
    }

    return vad;
};

void VADthr_adapt::save_frame(int frame_index) {
    if(frame_index <= opt_init)
        file_init->write(true);
    else
        file_init->write(false);
    
    file_crimean->write(crimean);
    file_crimean2->write(crimean2);
    file_crivar->write(crivar);
    file_thr->write(thr);
}

// -----------------------------------------------------------------------------
// ---------------- adaptive threshold based on dynamics VAD class -------------
class VADthr_dyn : public VADthr {
       private:
               const int opt_init;
               const double opt_perc;
               const double opt_dynmin;
               const double opt_qmaxinc;
               const double opt_qmaxdec;
               const double opt_qmindec;
               const double opt_qmininc;

               double dmin, dmax;
               double dyn;
               double thr;

               FileWriter *file_dmin;
               FileWriter *file_dmax;
               FileWriter *file_dyn;
               FileWriter *file_dynmin;
               FileWriter *file_thr;

       public:
               VADthr_dyn (const opts* options);
               virtual ~VADthr_dyn ();
               
               virtual void new_file(const char *filename);
               virtual bool process_frame(int frame_index, double cri);
               virtual void save_frame(int frame_index);
};

VADthr_dyn::VADthr_dyn(const opts* options)
: VADthr(options),
    opt_init(options->vad_dyn_init),
    opt_perc(options->vad_dyn_perc),
    opt_dynmin(options->vad_dyn_min),
    opt_qmaxinc(options->vad_dyn_qmaxinc),
    opt_qmaxdec(options->vad_dyn_qmaxdec),
    opt_qmindec(options->vad_dyn_qmindec),
    opt_qmininc(options->vad_dyn_qmininc)
{
    file_dmin = NULL;
    file_dmax = NULL;
    file_dyn = NULL;
    file_dynmin = NULL;
    file_thr = NULL;
};

VADthr_dyn::~VADthr_dyn() {
    DELETE(file_dmin);
    DELETE(file_dmax);
    DELETE(file_dyn);
    DELETE(file_dynmin);
    DELETE(file_thr);
};

void VADthr_dyn::new_file(const char *filename) {
    DELETE(file_dmin);
    file_dmin = new FileWriter(filename, "dmin");
    DELETE(file_dmax);
    file_dmax = new FileWriter(filename, "dmax");
    DELETE(file_dyn);
    file_dyn = new FileWriter(filename, "dyn");
    DELETE(file_dynmin);
    file_dynmin = new FileWriter(filename, "dynmin");
    DELETE(file_thr);
    file_thr = new FileWriter(filename, "thr");
} 

bool VADthr_dyn::process_frame(int frame_index, double cri) {
    bool vad;

    if(frame_index < MAX(1, opt_init))
    {
        dmax = cri;
        dmin = cri;
        dyn = 0.0;
        thr = cri;
        
        vad = false;
    }
    else if(frame_index == MAX(1, opt_init))
    {
        dmax = MAX(dmax, cri);
        dmax += opt_dynmin/10.0;

        dmin = MIN(dmin, cri);
        dmin -= opt_dynmin/10.0;

        dyn = dmax - dmin;
        thr = cri;

        vad = false;
    }
    else
    {
        if(dmax < cri)
            dmax = opt_qmaxinc*dmax + (1.0-opt_qmaxinc)*cri; // increase maximum
        else
            dmax = opt_qmaxdec*dmax + (1.0-opt_qmaxdec)*cri; // decrease maximum
    
        if(dmin > cri)
            dmin = opt_qmindec*dmin + (1.0-opt_qmindec)*cri; // decrease minimum
        else
            dmin = opt_qmininc*dmin + (1.0-opt_qmininc)*cri; // increase minumum

        dyn = dmax - dmin;
        thr = dmin + (opt_perc/100.0)*dyn;

        if((cri > thr) && (dyn > opt_dynmin))
            vad = true;
        else
            vad = false;
    }

    return vad;
};

void VADthr_dyn::save_frame(int frame_index)
{
    file_dmin->write(dmin);
    file_dmax->write(dmax);
    file_dyn->write(dyn);
    file_dynmin->write(opt_dynmin);
    file_thr->write(thr);
};

// -----------------------------------------------------------------------------
// ---------------- public class VAD -------------------------------------------

VAD::VAD(const opts* options, Vec<double> *Xsabs, Vec<double> *Xsph, Vec<double> *InFvec, Vec<double> *FeaFvec, Vec<double> **FeaFvecOUT)
:   opt_apply_mode(options->vad_apply_mode),
    opt_out_mode(options->vad_out_mode),
    opt_filter_order(options->vad_filter_order),
    xsabsvec(Xsabs)
{
    file_vad0 = NULL;
    file_vad = NULL;

    if(CRI_MODE_ENERGY(options->vad_cri_mode)) vadcri = new VADcri_energy(options, Xsabs, Xsph, InFvec, FeaFvec);
    else if (CRI_MODE_CEPDIST(options->vad_cri_mode)) vadcri = new VADcri_cepdist(options, Xsabs, Xsph, InFvec, FeaFvec);
    else throw("VAD: unknown vad_cri_mode!");

    if(THR_MODE_ABSOLUTE(options->vad_thr_mode)) vadthr = new VADthr_absolute(options);
    else if(THR_MODE_PERC(options->vad_thr_mode)) vadthr = new VADthr_perc(options);
    else if(THR_MODE_ADAPT(options->vad_thr_mode)) vadthr = new VADthr_adapt(options);
    else if(THR_MODE_DYN(options->vad_thr_mode)) vadthr = new VADthr_dyn(options);
    else throw("VAD: unknown vad_thr_mode!");

    vadfilter = new medianFilter(opt_filter_order, FeaFvec, FeaFvecOUT);
    vad_ready = false;
};

VAD::~VAD() {
    DELETE(file_vad0)
    DELETE(file_vad)
    delete vadcri;
    delete vadthr;
    DELETE(vadfilter);
};

void VAD::new_file(const char *filename) {
    frame_index = 0;

    if(OUT_MODE_NONE(opt_out_mode))
        return;

    if(!strlen(filename))
        throw("VAD::new_file(): invalid filename!");

    DELETE(file_vad);
    file_vad = new FileWriter(filename);   

    if(OUT_MODE_DEBUG(opt_out_mode))
    {
        DELETE(file_vad0);
        file_vad0 = new FileWriter(filename, "vad0");

        vadcri->new_file(filename);
        vadthr->new_file(filename);
    }
};

bool VAD::process_frame() {
    double cri;

    cri = vadcri->process_frame(frame_index);
    vad0 = vadthr->process_frame(frame_index, cri);
    vadcri->consume_vad(frame_index, vad0);
    vad = vadfilter->push(vad0);
    vad_ready = vadfilter->ready;

    frame_index++;
    return vad;
};

void VAD::clean() {
  vadfilter->cleanFilter();
  vad=false;
}

void VAD::save_frame() {
    
    if(OUT_MODE_NONE(opt_out_mode))
        return;

    if(vad_ready){
      file_vad->write(vad);

      if(OUT_MODE_DEBUG(opt_out_mode))
      {
        file_vad0->write(vad0);
        vadcri->save_frame(frame_index);
        vadthr->save_frame(frame_index);
      }
    }
};

void VAD::silence_frame() {
    if(vad)
        return;
    if(!APPLY_MODE_SILENCE(opt_apply_mode))
        return;

    int size = xsabsvec->get_size();
    for(int i=0; i<size; i++)
        (*xsabsvec)[i] = 0.0;
};

bool VAD::drop_frame() {
    return (!vad && APPLY_MODE_DROP(opt_apply_mode));
};

bool VAD::flush_frame(){
    vad = vadfilter->flush_frame();
    return vadfilter->ready;
}
