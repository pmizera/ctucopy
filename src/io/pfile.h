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

// This file is a courtesy of Petr Schwarz, Speech@FIT Brno

// This class reads and writes pfiles

// Modifications:
// 24. 4. 2004 by Petr Schwarz
// - createSentIndexFromData is called if the index does not exist and "rw" mode
// - goTo changes the file pointer in the "rw" mode
// - goTo changes the actual sentence number for writing
// - write moves the file pointer to the next frame

#ifndef _PFILE_H
#define _PFILE_H

#define PFILE_DEFAULT_HDR_SIZE 32768
#define PFILE_SENT_TABLE_SIZE   1000

#include "stdafx.h"

#define FOPEN fopen64

class PFile
{
	protected:
		FILE *fp;
		bool alredy_exists;

		// header informations
		struct 
		{
			unsigned int version;
			__off64_t size;
			unsigned int nsents;
			unsigned int nframes;
			unsigned int first_fea_col;
			unsigned int nfea;
			unsigned int first_lbl;
			unsigned int nlbl;
			__off64_t data_size;
			__off64_t data_offset;
			unsigned int data_dim;
			unsigned int data_nrows;
			unsigned int data_ncols;
			__off64_t sentt_size;
			__off64_t sentt_offset;
			unsigned int sentt_dim;
		} hdr;

		// Opening file mode
		char file_mode[4];
		// 
		unsigned int sent_id;
		unsigned int p_sent_id;
		unsigned int frame_id;
		unsigned int abs_frame_id;
		unsigned int act_sent_id;

		// Absolute numbers of beginning frames for each sentence
		unsigned int *sent_idx;
		__off64_t sent_idx_max_size;

		// Reset eos flag till next call of the read function
		bool reset_eos_flag;

		void nullVars();
		void error(const char *msg);
		void readHeader(char *file);
		void writeHeader();
		unsigned int swapBytes(unsigned int v);
		void vctSwapBytes(unsigned int *vct, int len);
		void loadSentIndex(char *file);
		void createEmptySentIndex();
		void createSentIndexFromData();
		void writeSentIndex();
	public:
		PFile();
		~PFile();
		// mode: r - read, w - write, rw - read/write
		bool open(char *file, char *mode, unsigned int fea_num = 0, unsigned int lbl_num = 0);
		void close();
		// End of file
		bool eof();
		// End of sentence
		bool eos();
		// Reset eos flag till the next call of read function
		// - this enables reading of pfile in two 'while' cycles
		// else 'while' and 'do' have to be used
		void resetEos() {reset_eos_flag = true;};
		unsigned int getFeaNum() {return hdr.nfea;};
		unsigned int getLblNum() {return hdr.nlbl;};
		// These two functions return actual sentence ID and frame ID  
		// (for frame which will be read in next call of the 'read' function)
		// Is valid if eof() == false
		unsigned int getSentID() {return sent_id;};
		unsigned int getFrameID() {return frame_id;};
		// Returns the farme ID which is related to begin of file instead of sentence
		unsigned int getAbsFrameId() {return abs_frame_id;};		
		unsigned int getFramesNum(unsigned int sid);
		unsigned int getFramesNum() {return getFramesNum(getSentID());};
		// Returns number of frames for the longest sentence
		unsigned int getMaxFramesNum();
		unsigned int getAbsFramesNum() {return hdr.nframes;};
		unsigned int getSentencesNum() {return hdr.nsents;};

		void goTo(unsigned int sid, unsigned int fid);
		void rewind() {goTo(0, 0);};
		bool read(float *features, unsigned int *labels);
		void write(float *features, unsigned int *labels);
		void startNewSent();
};
#endif

