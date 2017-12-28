/*
 CtuCopy 4 - Universal feature extractor and speech enhancer.
 Copyright 2012 Petr Fousek, FEE CTU Prague
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

#include <fstream>
#include "batch.h"
#include "../base/types.h"
#include "string.h"

BATCH::BATCH (opts *o_in) {
    o = o_in;    // pointer to options

    // decide on using classes FB and FEA
    if( (0==strcmp(o->format_out,"htk") || 0==strcmp(o->format_out,"pfile") || 0==strcmp(o->format_out,"ark")) && 0!=strcmp(o->format_in,"htk") && 0!=strcmp(o->fea_kind,"td-iir-mfcc")){
         do_fea = true;
    }else{
         do_fea = false;
    }

    if(strcmp(o->vad_apply_mode, "none") || strcmp(o->vad_out_mode, "none")) {
        do_vad = true;
    } else {
        do_vad = false;
    }

    in = new IN(o);
    out_fea=1;
    process_stage=0;
    out = NULL;out_stat_cmvn = NULL;fea = NULL;post = NULL;fb = NULL;nr = NULL;
    in_stat_cmvn = NULL;fea_d = NULL;fea_d_d = NULL;fea_d_d_d = NULL;

    if(do_fea){
          if(o->nr_when == o->afterFB) { // NR performed on FB output
               fb = new FB(o, in->_Xsabs);
               nr = new NR(o, fb->_Y, in->_Xsph);
          }else{                                               // first NR, then FB
               nr = new NR(o, in->_Xsabs, in->_Xsph);
               fb = new FB(o, in->_Xsabs);
          }
          fea = new FEA(o, fb->_Y);
          init_delta(fea->_fvec,in->_Xsph);
          init_cmvn(fea->_fvec);
    }else{  // signal output -> no FB, no FEA or input FEA
         if(0==strcmp(o->format_in,"htk")) {// input is HTK features -> DELTA, STACK, CMVN  or  RAW -> VAD, td-iir-mfcc, NR
               init_delta(in->_fvec,in->_Xsph);
               init_cmvn(in->_fvec);
         }else if(0==strcmp(o->fea_kind,"td-iir-mfcc")) {// input is raw and features are based on td-iir-mfcc
               out = new OUT(o, in->_Xsabs, in->_Xsph, in->E);
         }else{  // input is raw; only speech enhancement
               nr  = new NR(o, in->_Xsabs, in->_Xsph);
               out = new OUT(o, in->_Xsabs, in->_Xsph, in->E);
         }
    }

};
void BATCH::init_out(Vec<double> *Xabs, Vec<double> * Xph) {
	// properly initializes output class
	// output type = features (htk/pfile/ibm)

        if(do_vad){
          Vec<double> *FeaFvecOut;
          vad = new VAD(o, in->_Xsabs, in->_Xsph, Xabs, Xabs, &FeaFvecOut);
	  if(0==strcmp(o->fea_kind,"dctc")) { // DCT-derived features, c0 represents energy after NR,FB
	    if(o->nr_when == o->afterFB){
               // we'd like to compute E after NR but in case FB precedes NR it is
	       // impossible -> get normal energy before FB instead
               out = new OUT(o, FeaFvecOut, Xph, in->E);	// IN decides how to compute E
            }else{
               // this is better, get E after NR right before FB
                  out = new OUT(o, FeaFvecOut, Xph, nr->E);
            }
          }else{
            // LPC-derived fea -> E derived from 1st autocorr coeff
            // or Spectral features -> E computed from input spectrum (after FB,NR)
            if(0==strcmp(o->format_in,"htk")) {
              out = new OUT(o, FeaFvecOut, Xph, in->E);
            }else{
              out = new OUT(o, FeaFvecOut, Xph, fea->E);
            }
          }
          return;
        }

	if(o->fea_rawenergy){
          out = new OUT(o, Xabs, Xph, in->E); 		// E computed at the very beginning, IN knows how
	}else{
	  if(0==strcmp(o->fea_kind,"dctc")) { // DCT-derived features, c0 represents energy after NR,FB
	    if(o->nr_when == o->afterFB){
               // we'd like to compute E after NR but in case FB precedes NR it is
	       // impossible -> get normal energy before FB instead
               out = new OUT(o, Xabs, Xph, in->E);	// IN decides how to compute E
            }else{
               // this is better, get E after NR right before FB
                  out = new OUT(o, Xabs, Xph, nr->E);
            }
          }else{
            // LPC-derived fea -> E derived from 1st autocorr coeff
            // or Spectral features -> E computed from input spectrum (after FB,NR)
            if(0==strcmp(o->format_in,"htk")) {
              out = new OUT(o, Xabs, Xph, in->E);
            }else{
              out = new OUT(o, Xabs, Xph, fea->E);
            }
          }
        }
};

void BATCH::init_delta(Vec<double> * in_fvec, Vec<double> * in_Xsph){
          if(o->n_order>=1) fea_d     = new deltaFEA(o,in_fvec,2,o->d_win);
          if(o->n_order>=2) fea_d_d   = new deltaFEA(o,fea_d->_fvec,3,o->a_win);
          if(o->n_order==3) fea_d_d_d = new deltaFEA(o,fea_d_d->_fvec,4,o->t_win);
          if(o->n_order==0)       init_out(in_fvec, in_Xsph);
           else if(o->n_order==1) init_out(fea_d->_fvec, in_Xsph);
           else if(o->n_order==2) init_out(fea_d_d->_fvec, in_Xsph);
           else if(o->n_order==3) init_out(fea_d_d_d->_fvec, in_Xsph);
};
void BATCH::init_cmvn(Vec<double> * in_fvec){
          if(o->apply_cmvn){
              char tmp[(int)strlen(o->format_in)] ;
              strcpy(tmp,o->format_in);
              strcpy(o->format_in,"cmvn");
              in_stat_cmvn = new IN(o);
              strcpy(o->format_in,tmp);
            if(in_stat_cmvn->_list_ID_spk!=NULL){
              list_ID_spk=in_stat_cmvn->_list_ID_spk;
              if(o->n_order==0) post = new POST(o, in_fvec,in_stat_cmvn->_M_cmean,in_stat_cmvn->_M_cvar, in_stat_cmvn->_list_ID_spk);
              if(o->n_order==1) post = new POST(o, fea_d->_fvec,in_stat_cmvn->_M_cmean,in_stat_cmvn->_M_cvar, in_stat_cmvn->_list_ID_spk);
              if(o->n_order==2) post = new POST(o, fea_d_d->_fvec,in_stat_cmvn->_M_cmean,in_stat_cmvn->_M_cvar, in_stat_cmvn->_list_ID_spk);
              if(o->n_order==3) post = new POST(o, fea_d_d_d->_fvec,in_stat_cmvn->_M_cmean,in_stat_cmvn->_M_cvar, in_stat_cmvn->_list_ID_spk);
              o->stat_cmvn=false;
            }else{
              cout<<"IN: Cannot open stat. cmvn file!"<<endl;
              cout<<"IN: Stat. cmvn file is being created: "<< o->fcmvn_stat_in <<endl;
              strcpy(o->fcmvn_stat_out,o->fcmvn_stat_in);
              o->stat_cmvn=true;
              init_post(in_fvec);
            }
          }else if(o->fea_Z_exp>0 || o->fea_Z_block>0){
            init_post(in_fvec);
          }else if(o->stat_cmvn){
            out=NULL; out_fea=0; // only Statistics for CMVN
            init_post(in_fvec);
          }
};
void BATCH::init_post(Vec<double> * in_fvec) {
             if(o->n_order==0) post = new POST(o, in_fvec);
             if(o->n_order==1) post = new POST(o, fea_d->_fvec);
             if(o->n_order==2) post = new POST(o, fea_d_d->_fvec);
             if(o->n_order==3) post = new POST(o, fea_d_d_d->_fvec);
             if(!(o->fea_Z_exp>0 || o->fea_Z_block>0)){
               char tmp[(int)strlen(o->format_out)] ;
               strcpy(tmp,o->format_out);
               strcpy(o->format_out,"cmvn");
               out_stat_cmvn  = new OUT(o, post->_mvec_cmean, post->_mvec_cvar,post->_list_ID_spk);
               strcpy(o->format_out,tmp);
             }
};
void BATCH::fea_delta(){//-------------------------- Compute differencial features
     if(o->n_order>=1){
       if(fea_d->process_frame()){
         if(o->n_order>=2){
           if(fea_d_d->process_frame()){
             if(o->n_order==3){
               if(fea_d_d_d->process_frame()){
                 cmvn_stat();
               }
             }else{
                 cmvn_stat();  // default
             }
           }
         }else{
              cmvn_stat();
         }
       }
     }else{
       cmvn_stat();
     }
}
void BATCH::cmvn_stat(){//-------------------------- Accumulate statistic for CMVN
       if(!process_stage && o->stat_cmvn){
         post->sum_fea();
       }else if(process_stage==1 && o->stat_cmvn){
         post->sum_cv();
       }else if((!process_stage || process_stage==2) && o->apply_cmvn || (o->fea_Z_exp>0 || o->fea_Z_block>0)){
         post->process_frame();
         save_frame();
       }else{
         save_frame();// fea can now be delayed
       }
}
void BATCH::process_frame() {
       if(do_fea){//-------------------------------- Do feature extraction ?
         if(o->nr_when == o->afterFB){
           fb->project_frame();
           nr->process_frame();
         }else {
           nr->process_frame();
           fb->project_frame();
         }
         if(fea->process_frame()){
           fea_delta();
         }
       }else if(o->fea_delta || o->fea_trap){//------- Input is features vector or does trap ?
         fea_delta();
       }else if(o->stat_cmvn || o->apply_cmvn){//----- Do CMVN (both compute statistic and aplly CMVN on features) ?
         cmvn_stat();
       }else if(0==strcmp(o->fea_kind,"td-iir-mfcc")) {//---- Do td-iir-mfcc fea
         save_frame();
       }else{ //-------------------------------------- Do only noise reduction?
         if(nr != NULL)//------only for create ark file for KALDI from HTK input features
             nr->process_frame();
         save_frame();
       }
}

void BATCH::save_frame() {
    if(do_vad){
        vad->process_frame();
        vad->silence_frame();
        vad->save_frame();
        if(!vad->vad_ready)
          return;
        if(vad->drop_frame())
            return;
    }
    out->save_frame();
}

void BATCH::flush_vad() {
    while(vad->flush_frame()){
        vad->save_frame();
        if(!vad->drop_frame())
            out->save_frame();
    }
}

void BATCH::flush_fea() {
       if(do_fea || (o->fea_delta||o->fea_trap)){
         if(!((o->fea_delta||o->fea_trap))){
            while(fea->flush_frame()){ // fea->flush_frame()
                  cmvn_stat();
            }
         }
         if(o->n_order>=1){
           while(fea_d->flush_frame()){
                 if(o->n_order>=2){
                   if(fea_d_d->process_frame()){
                     if(o->n_order>=3){
                       if(fea_d_d_d->process_frame()){
                         cmvn_stat();
                       }
                     }else{
                       cmvn_stat();
                     }
                   }
                 }else{
                   cmvn_stat();
                 }
           }
         }
         if(o->n_order>=2){
           while(fea_d_d->flush_frame()){
                if(o->n_order>=3){
                  if(fea_d_d_d->process_frame()){
                     cmvn_stat();
                  }
                }else{
                  cmvn_stat();
                }
           }
         }
         if(o->n_order>=3){
           while(fea_d_d_d->flush_frame()){
                cmvn_stat();
           }
         }
       }
       if(do_vad){
         flush_vad();
         vad->clean();
       }
}
int BATCH::process() {
    // ------- SINGLE FILE MODE -------
    if (o->pipe_in || 0!=strcmp(o->in,"")) {	// output is also single mode (ensured by OPTS)
      if (o->verbose) cerr << "processing: " << o->in << " ";
      in->new_file(o->in);		// in knows when to use o->in
      nr->new_file();

      if(do_fea) {
        fb->new_file();
        fea->new_file();
      }

      out->new_file(o->out);		// out knows when to use o->out
      vad->new_file(o->vad_out);

      // loop over segments
      int segs = 0; // frames counter
      while (in->get_frame()) {
        process_frame();
      }

      // possibly flush 'fea'
      flush_fea();

      if (o->verbose) cerr << "- " << segs << " frames." << endl;
      return 0;
    }

    // ------- BATCH MODE -------
    char line[1000]; 		// to store one pair <input file - output file>
    char *fin, *fout, *ID_spk, *foutvad;

    int num_spk=0;
    int post_CN_stage;
    if(o->stat_cmvn){
      post_CN_stage=1;              // Statistics for cmvn will computed.
      if(o->apply_cmvn)
        post_CN_stage=2;            // CMVN will applied to features.
    }else{
      post_CN_stage=0;              // Standard processing
    };

    for(int stage=0; stage<=post_CN_stage; stage++){
      process_stage = stage;
      ifstream list;
      if (o->list == NULL) throw "BATCH: Nothing to do!";
      list.open(o->list); // open list
      if(!list) throw "BATCH: Cannot open list file!";

      if(0==strcmp(o->fea_kind,"td-iir-mfcc"))
        in->loadf_filters(o->ffilters);

      while(list.getline(line,1000)){      // loop over all files
        fin    = strtok(line," \t");        // input file name
        fout   = strtok (NULL," \t");       // output file name (in case of pfile output discarded)
        ID_spk = strtok(NULL," \t");        // ID speaker for CMVN
        foutvad = strtok (NULL," \t");      // VAD output file name

        if(!fin || !fout) throw "BATCH: Bad list format!";
        if(do_vad && !foutvad) throw "BATCH: Bad list format!";
        if(o->verbose) cerr << "processing: " << fin << " ";

        if(o->stat_cmvn){
          if(!ID_spk){
             ID_spk=fout;
             if(!fin || !ID_spk) throw "BATCH: Bad list format for computed stat. cmvn!";
          }
          num_spk=post->add_spk(ID_spk);
        }else if(o->apply_cmvn){
          if(!fin || !ID_spk || !fout) throw " Bad list format for applying cmvn!";
          num_spk=post->add_spk(ID_spk);
        }

        // announce new file, let objects close the previous and open a new one
	in->new_file(fin);
	if(do_fea){
	  fb->new_file();
	  fea->new_file();
	  if(o->n_order>=1) fea_d->new_file();
	  if(o->n_order>=2) fea_d_d->new_file();
	  if(o->n_order==3) fea_d_d_d->new_file();
	}else if(o->fea_delta||o->fea_trap){
	  if(o->n_order>=1) fea_d->new_file();
	  if(o->n_order>=2) fea_d_d->new_file();
	  if(o->n_order==3) fea_d_d_d->new_file();
        }

        if(strcmp(o->format_in,"htk") && strcmp(o->fea_kind,"td-iir-mfcc")){
          nr->new_file();
        }

        if(out_fea){
          if(stage==post_CN_stage){
            out->new_file(fout);
            if(o->fea_Z_exp>0 || o->fea_Z_block>0)
              post->clean();
          }
          if( 0==strcmp(o->format_out,"ark") && ((stage==0||stage==1 ) && o->stat_cmvn) )
            out_fea=0;
        }

        if(do_vad)
          vad->new_file(foutvad);

        int segs = 0; // frames counter
        while(in->get_frame()){  //  loop over all segments in the file
          process_frame();
          segs++;
        }
        // possibly flush 'fea' ???
        flush_fea();

        if(o->verbose) cerr << "- " << segs << " frames." << endl;
      }// END - loop over all files

      if(!stage && o->stat_cmvn){
        post->stat_cm();
      }else if(stage==1 && o->stat_cmvn){ // Save statistics
        post->stat_cv();
        out_stat_cmvn->new_file(o->fcmvn_stat_out);
        out_stat_cmvn->save_frame(num_spk);
        out_fea=1;
      }
   }
   return 0;
}

BATCH::~BATCH () {

    if(in != NULL)
      delete in;
    if(out != NULL){
      if(out_fea)
       delete out;
    }
    if(do_fea)
      delete fea;
    if(post!= NULL)
      delete post;
    if(fb != NULL)
      delete fb;
    if(nr != NULL)
      delete nr;
    if(in_stat_cmvn != NULL)
      delete in_stat_cmvn;
    if(fea_d!= NULL)
      delete fea_d;
    if(fea_d_d!= NULL)
      delete fea_d_d;
    if(fea_d_d_d!= NULL)
      delete fea_d_d_d;
    if(out_stat_cmvn != NULL)
      delete out_stat_cmvn;
    if(do_vad)
      delete vad;
};

