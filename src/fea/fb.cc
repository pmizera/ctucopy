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

#include "fb.h"

FB::FB (opts* o_in, Vec<double> * X) {
	o = o_in;
	_Xs = X;
	size = 0;

	PLP = (0==strcmp(o->fb_shape,"trapez")); // PLP design only if shape=trapez !

	try {
		// grab space
		hz_axis   = new double [o->wfftby2];	// lin axis in Hz
		warp_axis = new double [o->wfftby2];	// the same but in warped units
		mat  = new double*[999]; // hope 999 filters is enough
		// matrix of <Nfilters> x <input_spectr_size + 2>
		// last 2 columns hold last and first nonzero element of rows, respect.
		for (int i = 0; i < 999; i++) 
			mat[i] = new double [o->wfftby2 + 2];
		bank = new subbank[999];
		subbanks = 0;
	}
	catch (...) { throw "FB: Not enough memory!"; }
		for (int i = 0; i < 999; i++) // clean up before design
		for (int j = 0; j < o->wfftby2 + 2; j++) 
			mat[i][j] = 0.;

	if (PLP) {
		if (0!=strcmp(o->fb_scale, "bark") && o->verbose) 
			cerr << "FB: WARNING: -fb_shape set to trapez => Using Bark scale!" << endl;
		if (o->verbose)
			cerr << "FB: WARNING: PLP filter bank => discarding FB definition,"
				 << "             enabling Intensity-Loudness Power Law and Equal Loudness!" 
				 << endl;
		strcpy(o->fb_scale, "bark");
		o->fb_inld = true;
		o->fb_eqld = true;
	}

	init_scale();	// design linear Hz axis and it's warped version
	if (PLP) plp_design();	// treat PLP in special  way
	else {
		parse();			// else parse filter bank specification
		check_and_design();
	}
	optimize();				// finally optimize filter bank (see below)

//	for (int i = 0; i < size; i++) {
//        for (int j = 0; j < o->wfftby2; j++)
//        	cout << mat[i][j] << "\t";
//        cout << endl;
//    }
//    throw;

	try { _Y = new Vec<double>(size); } // alloc output spect.
	catch (...) { throw "FB: Not enough memory!"; }
};

bool FB::new_file() {
	return true;
}

bool FB::project_frame() {
	// the main thing here ;) - project spectrum X -> Y
	static Vec<double> &X = *_Xs;	// set pointers
	static Vec<double> &Y = *_Y;
	for (int i=0; i<size; i++) {
		Y[i] = 0;
        for (int k = (int) mat[i][o->wfftby2+1]; k <= (int) mat[i][o->wfftby2]; k++)
        	// multiply only where mat[i] has nonzero values, EQLD already in mat
			Y[i] += X[k]*mat[i][k];
		if (o->fb_inld)
			// nonlinear step, has to be performed separately
			Y[i] = pow(Y[i], 0.33);
 	}
	return true;
}


FB::~FB() {
	delete _Y;
	for (int i = 0; i < 999; i++)
		delete [] mat[i];
	delete [] mat;
	delete [] bank;
	delete [] hz_axis;
	delete [] warp_axis;
}


void FB::init_scale() {
	// precomputing the warped "xscale-frequency" axis, xscale - lin, mel, bark or expolog
	for (int i=0; i<o->wfftby2; i++)
		hz_axis[i] = (double)i*o->fs/(double)o->wfft; // linear
		
	// linear	
	if (0==strcmp(o->fb_scale,"lin")) {
		for (int i=0; i<o->wfftby2; i++)
			warp_axis[i] = hz_axis[i];
	}
	
	// Bark
	else if (0==strcmp(o->fb_scale,"bark")) {
		for (int i=0; i<o->wfftby2; i++) 
			warp_axis[i] = 6.*log(hz_axis[i]/600. + sqrt((hz_axis[i]/600.)*(hz_axis[i]/600.) + 1.));
	}
	
	// expolog
	else if (0==strcmp(o->fb_scale,"expolog")) {
		for (int i=0; i<o->wfftby2; i++) {
			if (hz_axis[i] <= 2000) 
				warp_axis[i] = 700.*(pow(10., hz_axis[i]/3988.) - 1.); // exp
			else            
				warp_axis[i] = 2595.*log10(1. + hz_axis[i]/700.); // mel ~ log
		}
	}
	
	// melodic
	else if (0==strcmp(o->fb_scale,"mel")) {
		for (int i=0; i<o->wfftby2; i++) 
			warp_axis[i] = 2595*log10(1. + hz_axis[i]/700.);
	}
}	

void FB::plp_design() {
	// in PLP, many things are fixed such as # of bands, EQLD=yes, INLD=yes,
	// that is why it has it's own routine. Only here are trapez filters.
	
	// maximum of Bark scale
	double maxBark = 6*log(o->fs/1200. + sqrt((o->fs/1200.)*(o->fs/1200.) + 1.));

	// number of Bark's intervals
	int nBark = (int) (floor(maxBark + .5));
	int plpsize = nBark - 1; // every filter "lies around" it's center frequency,
	                         // 1st is around 1Bark, last is around nBark-1 Bark

	// real size of one step in Bark (FB needs to end exactly at fs/2)
	double Barkstep = maxBark/(double)nBark;
	
	if (o->verbose) cerr << "FB: Filter bank size = " << plpsize << endl;

	for (int i = 0; i < plpsize; i++) {
		// loop over all filters
		double Om = (i + 1)*Barkstep;             		// current center position [Bark]
		double om = 3.1415926535898*1200*sinh(Om/6);   	// Bark to rad/s ( = 2*pi*Hz)

		// get equal loudness coefficient
		double eqnum = om*om*om*om*(om*om + 5.68e7);
		double eqden = 0.;
		if (o->fs <= 10000)
			eqden = (om*om + 6.3e6)*(om*om + 6.3e6)*(om*om + 3.8e8);
		else
			eqden = (om*om + 6.3e6)*(om*om + 6.3e6)*(om*om + 3.8e8)
				*(om*om*om*om*om*om + 9.58e26);
		double eqloud = eqnum/eqden;

		for (int k=0; k<o->wfftby2; k++) {
			// compute this filter shape over Hz axis using trapez (see Hermansky's PLP paper)
			double diff = warp_axis[k] - Om;
			if (diff >= -1.3 && diff <= -.5)
				mat[i][k] = pow(10., 2.5*(0.5 + diff));
			else if (fabs(diff) < 0.5)
				mat[i][k] = 1;
			else if (diff >= 0.5 && diff <= 2.5)
				mat[i][k] = pow(10., 0.5 - diff);
			else 
				mat[i][k] = 0;

			if (o->fb_eqld) // apply the equal loudness funtion
				mat[i][k] *= eqloud;
		}
		// btw, we do not normalize filter area here...
	}
	size = plpsize;
};

void FB::parse() {
	// ***** parse filterbank specification string ******
	const char delim[] = ",";
	char *copy = o->fb_definition;
	char *subfilt = strtok(copy,delim);
	int nfilt = 0;
	subbanks = 0;
	
	while (subfilt != NULL) { 
		// every while cycle parses parses one token of type
		// token = [[X-YHz:]K-L/]Nfilters
		// specifying K-L neighboring filters out of N filters placed equidistantly 
		// on warped axis from f_low Hz to f_high Hz 

		int offset = 0;				// byte offset within token
		double flow = 0.;			// low cutoff - default
		double fhigh = o->fs/2.;	// high cutoff - default
		int start = 1;				// index of 1st filter - indexed from 1
		int stop = 0;				// last one - this means no design
		
		// Now, parse the token.
		// get 1st number 
		double x = atof(subfilt + offset); offset += strspn(subfilt+offset,".1234567890");
		
		if (0==strncmp(subfilt+offset,"-",1)) {
			// next char == "-"
			offset++;
			// get 2nd number
			double y = atof(subfilt + offset); offset += strspn(subfilt+offset,".1234567890");

			// is this freq specification?
			if (0==strncmp(subfilt+offset,"Hz:",3)) {
				// yes, then we can assign all
				offset += 3;
				flow = x;
				fhigh = y;
				start = atoi(subfilt + offset); offset += strspn(subfilt+offset,"1234567890");
				if (0!=strncmp(subfilt+offset,"-",1)) throw "FB: Filter bank specification parse error!";
				stop  = atoi(subfilt + ++offset); offset += strspn(subfilt+offset,"1234567890");
				if (0!=strncmp(subfilt+offset,"/",1)) throw "FB: Filter bank specification parse error!";
				nfilt = atoi(subfilt + ++offset); offset += strspn(subfilt+offset,"1234567890");
				if (0!=strncmp(subfilt+offset,"filters",7)) throw "FB: Filter bank specification parse error!";
			} 
			else if (0==strncmp(subfilt+offset,"/",1)) {
				// OK, it was not freq. spec. but filter range specification
				start = (int) x;
				stop = (int) y;
				offset++;
				nfilt = atoi(subfilt + offset); offset += strspn(subfilt+offset,"1234567890");
				if (0!=strncmp(subfilt+offset,"filters",7)) throw "FB: Filter bank specification parse error!";
			} else throw "FB: Filter bank specification parse error!";
		} 
		else if (0==strncmp(subfilt+offset,"filters",7)) {
			// user specified only # of filters
			stop = (int) x;
			nfilt = stop;
		} else 
			throw "FB: Filter bank specification parse error!";
		subfilt = strtok(NULL,delim);
	
		bank[subbanks].f_start    = flow;
		bank[subbanks].f_stop     = fhigh;
		bank[subbanks].bands      = nfilt;
		bank[subbanks].band_first = start;
		bank[subbanks].band_last  = stop;
		subbanks++;
	}
};

void FB::check_and_design() {
		
	if (0==strcmp(o->fb_shape,"rect")) { // rectangular shape only !
		// *** First, check if any two subbanks are placed next to each
		// other. If yes, modify boundaries so that there is no overlap.
	
		double df = o->fs/(double)o->wfft;
		// frequency difference (distance between two bins on frq axis)
	
		for (int i=0; i<subbanks; i++) {
			// for each subbank we have to check endings
			bool no_join = true; // unless we find something that joins..
			for (int j=0; j<subbanks; j++) 
				if (bank[i].f_stop == bank[j].f_start) no_join = false;

			if (o->verbose) { 
				cerr << "FB: Adding filters: "<<bank[i].f_start<<" - "<<bank[i].f_stop<<" Hz, "
				 <<bank[i].band_first<<" - "<<bank[i].band_last<<" out of "
				 <<bank[i].bands<<" filters, last bin ";
				 if (no_join) cerr <<"included."<< endl; else cerr << "excluded."<<endl;
			}
				
			// if nothing joins, make sure f_stop belongs to the bank:
			if (no_join) bank[i].f_stop += df;
		}
	} else { // any other shape
		if (o->verbose) 
			for (int i=0; i<subbanks; i++)
				cerr << "FB: Adding filters: "<<bank[i].f_start<<" - "<<bank[i].f_stop<<" Hz, "
				 <<bank[i].band_first<<" - "<<bank[i].band_last<<" out of "
				 <<bank[i].bands<<" filters."<<endl;
	}		
		
	// *** Second, design the bank subbank by subbank, filter by filter.
	size = 0;
	for (int sb=0; sb<subbanks; sb++) {
		// go through all subbanks
		for (int b=bank[sb].band_first; b<=bank[sb].band_last; b++) { 
			// call design routine for every single filter
			
//            cerr << "get filter ("<<mat[size]<<","<<bank[sb].f_start<<","<<bank[sb].f_stop<<","
//                                       <<bank[sb].bands <<","<<b<<");"<<endl;
			
			get_filter(mat[size], bank[sb].f_start, bank[sb].f_stop, bank[sb].bands, b);
			size++;
			if (size == 999) throw "FB: Too many filters in FB!";
		}
	}
};


void FB::get_filter(double *vec, double f_start, double f_stop, int bands, int band_index) {
	
	double w_high, w_low;
	// here we will add (band_index)th filter of filter subset containing (bands) filters

	// ---- 1st, find frequency boundaries of subset in warped units --------

	// linear	
	if (0==strcmp(o->fb_scale,"lin")) {
		w_high = f_stop;
		w_low  = f_start;
	}
	
	// Bark
	else if (0==strcmp(o->fb_scale,"bark")) {
		w_high = 6.*log(f_stop/600. + sqrt((f_stop/600.)*(f_stop/600.) + 1.));
		w_low  = 6.*log(f_start/600. + sqrt((f_start/600.)*(f_start/600.) + 1.));
	}
	
	// expolog
	else if (0==strcmp(o->fb_scale,"expolog")) {
		if (f_stop <= 2000) w_high = 700.*(pow(10., f_stop/3988.) - 1.);
		else                w_high = 2595.*log10(1. + f_stop/700.);
		if (f_start <= 2000) w_low = 700.*(pow(10., f_start/3988.) - 1.);
		else                 w_low = 2595.*log10(1. + f_start/700.0);
	}
	
	// melodic
	else if (0==strcmp(o->fb_scale,"mel")) {
		w_high = 2595*log10(1 + f_stop/700.0);
		w_low  = 2595*log10(1 + f_start/700.0);
	}
	else throw "FB: Unknown frequency scale!";
	
	// ----------- design a filter based on the desired shape ----------
	
	// start and stop position (in warped units) of the required filter
	double w_start, w_end;
	if (0==strcmp(o->fb_shape,"rect")) { // overlap fixed at 0% !
		w_start = w_low + (band_index - 1.)*(w_high - w_low)/(double)(bands);
		w_end   = w_low + (band_index + 0.)*(w_high - w_low)/(double)(bands);
	} else if (0==strcmp(o->fb_shape,"triang")) { // overlap fixed at 50% !
		w_start = w_low + (band_index - 1.)*(w_high - w_low)/(double)(bands + 1);
		w_end   = w_low + (band_index + 1.)*(w_high - w_low)/(double)(bands + 1);
	} else throw "FB: Unknown filter shape!";
	
	// compute Equal Loudness coefficient && filter area (one for each band)
	double eqloud = 1.;
	if (o->fb_eqld) {

		// get approx. location of window center on warped scale -> hz,
		// EQLD will be approximated by that point
		double w_mid = w_start + (w_end - w_start)/2.;
		double f_mid = 0;
		if (0==strcmp(o->fb_scale,"lin")) // linear
			f_mid = w_mid;
		if (0==strcmp(o->fb_scale,"bark")) // Bark
			f_mid = 600*sinh(w_mid/6.);
		if (0==strcmp(o->fb_scale,"expolog")) {
			if (w_mid <= 1521.4) {
				f_mid = 3988.*log10(1. + (w_mid/700.));
			}else{
				f_mid = 700.*(pow(10., w_mid/2595.) - 1);
			}
		}
		if (0==strcmp(o->fb_scale,"mel")) 
			f_mid = 700.*(pow(10., w_mid/2595.) - 1.);
					
		// get equal loudness coefficient from om = 2*pi*f_mid
		double om = 2 * 3.141592653589793 * f_mid;
		double eqnum = om*om*om*om*(om*om + 5.68e7);
		double eqden = 0.;
		if (o->fs <= 10000)
			eqden = (om*om + 6.3e6)*(om*om + 6.3e6)*(om*om + 3.8e8);
		else
			eqden = (om*om + 6.3e6)*(om*om + 6.3e6)*(om*om + 3.8e8)
				*(om*om*om*om*om*om + 9.58e26);
		eqloud = eqnum/eqden;
	}

	//	-------- now the real design ----------
	// rect shape
	if (0==strcmp(o->fb_shape,"rect")) {
		double area = 0; // filter area
		for (int i=0; i<o->wfftby2; i++) {
			// loop over fr. axis
			if (warp_axis[i] >= w_start && warp_axis[i] < w_end) {
				// this point belongs to current rectangle
				vec[i] = 1;
				area++;
			} else {
				// point out of current rectangle
				vec[i] = 0;
			}
		}
		// normalize filter area, apply EQ-LD
		if (o->fb_norm)
			for (int i=0; i<o->wfftby2; i++) vec[i] *= eqloud/area;
		else
			for (int i=0; i<o->wfftby2; i++) vec[i] *= eqloud;
	}
	// triang shape
	else if (0==strcmp(o->fb_shape,"triang")) {
		double area = 0; // filter area
		for (int i=0; i<o->wfftby2; i++) {
			if (warp_axis[i] < w_start || warp_axis[i] > w_end)
				// this point is not in current triangle
				vec[i] = 0;
			else {
				// write one triangle point height into mat
				double w_mid = w_start + (w_end - w_start)/2.;
				vec[i] = 1. - 2.*fabs(w_mid - warp_axis[i])/(w_end - w_start);
				area += vec[i];
			}
		}
		// normalize filter area, apply EQ-LD
		if (o->fb_norm)
			for (int i=0; i<o->wfftby2; i++) vec[i] *= eqloud/area;
		else
			for (int i=0; i<o->wfftby2; i++) vec[i] *= eqloud;
	}

	else throw "FB: Unknown filter shape!";
};


void FB::optimize() {
	// This speeds up further processing stage. Instead of computing
	// projection of input spectrum to the filter over whole frequency axis
	// we multiply-add only at filter's nonzero positions.
	// So here we find left and right boundaries and write them 
	// at the end of the vector.
	
	for (int band=0; band<size; band++) {
		int	k = 0;
        // write index of the first non-zero element to last position!!!
        while (mat[band][k++] == 0);
        mat[band][o->wfftby2 + 1] = k - 1;
        // write index of the last non-zero element to (last - 1) position!!!
        while (mat[band][k++] != 0);
        mat[band][o->wfftby2] = k - 2;
	}
	
	// this is a small hack for printing out the FB in ASCII
	if (o->fb_printself) {
		for (int band=0; band<size; band++) {
			for (int i=0; i<o->wfftby2; i++)
				cerr << mat[band][i] << "\t";
			cerr << endl;
		}
	}
};
