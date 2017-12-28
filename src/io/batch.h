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

#ifndef batch_h
#define batch_h

#include "stdafx.h"
#include "opts.h"
#include "in.h"
#include "../nr/nr.h"
#include "../fea/fb.h"
#include "../fea/fea.h"
#include "out.h"
#include "../fea/post.h"
#include "../vad/vad.h"

class BATCH {
	//
	// class BATCH
	//
	// This is the main processing class. It constructs the processing chain,
	// initializes all needed objects, processes all input files in a loop
	// and then closes all objects.

	bool do_fea;	// flag whether to include FB and FEA objects
        bool do_vad;
        char **list_ID_spk;
        int out_fea;
        int process_stage;
	opts  *o;		// options class
	IN    *in;		// input class
	NR    *nr;		// noise reduction class
	FB    *fb;		// filter bank class
	FEA   *fea;		// feature extraction class
        VAD   *vad;             // voice activity detection class
	deltaFEA   *fea_d;	// delta
	deltaFEA   *fea_d_d;    // delta-delta
	deltaFEA   *fea_d_d_d;  // delta-delta-delta
	POST  *post;		// Post-processing with features as a CMVN, delta, delta-delta,...
	OUT   *out;		// output class
	OUT   *out_stat_cmvn;	// output object for save stat. cmvn
        IN    *in_stat_cmvn;    // input object for load stat. cmvn
	void init_out(Vec<double> *, Vec<double> *);
  void process_frame();
  void save_frame();
  void cmvn_stat();
  void fea_delta();
  void flush_fea();
  void flush_vad();
  void init_cmvn(Vec<double> *);
  void init_delta(Vec<double> *, Vec<double> *);
  void init_post(Vec<double> *);
public:
	BATCH (opts*);
	int process();
	~BATCH ();
};

#endif

