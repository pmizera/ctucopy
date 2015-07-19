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


#include "stdafx.h" 
#include "pfile.h"

PFile::PFile()
{
	nullVars();
}

PFile::~PFile()
{
	if(fp)
		close();
	if(sent_idx)
		delete [] sent_idx;
}

void PFile::nullVars()
{
	fp = 0;
	sent_idx = 0;
	sent_id = 0;
	p_sent_id = 0;
	frame_id = 0;
	abs_frame_id = 0;
	memset(&hdr, 0, sizeof(hdr));
	hdr.version = 0xFFFF;   // Version 0 is the first version
        reset_eos_flag = 0;
	sent_idx_max_size = (__off64_t)0;
	act_sent_id = 0;
}

void PFile::error(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
	exit(1);
}

bool PFile::open(char *file, char *mode, unsigned int fea_num, unsigned int lbl_num)
{
	if(fp)
		error("ERROR: PFile:open - pfile is alredy open");
	
	nullVars();
		
	// check opening modes
	if(strcmp(mode, "r") != 0 && strcmp(mode, "w") != 0 && strcmp(mode, "rw") != 0)
		error("ERROR: PFile::open - invalid opening mode. Use 'r', 'w' or 'rw'.");

	// check whater the file alredy exists or not 
	alredy_exists = false;
	fp = FOPEN(file, "rb");
	if(fp)
	{
		alredy_exists = true;
		fclose(fp);
	}

	// Select appropriate opening mode
	strcpy(file_mode, "rb");

	if(strcmp(mode, "rw") == 0)
	{
		if(alredy_exists)
			strcpy(file_mode, "r+b");
		else
			strcpy(file_mode, "wb");
	}
	
	if(strcmp(mode, "w") == 0)
		strcpy(file_mode, "wb");

	// Open the pfile
	fp = FOPEN(file, file_mode);
	if(!fp)
	{
		char msg[1024];
		sprintf(msg, "ERROR: Can not open the pfile: %s.", file);
		error(msg);
	}

	if(alredy_exists && strcmp(mode, "w") != 0)
	{
		readHeader(file);

		if(hdr.sentt_size == (__off64_t)0)
			createSentIndexFromData();
		else
			loadSentIndex(file);
		rewind();
	}
	else
	{
		hdr.version = 0;
		hdr.size = (__off64_t)PFILE_DEFAULT_HDR_SIZE;
		hdr.nsents = 0;
		hdr.nframes = 0;
		hdr.first_fea_col = 2;
		hdr.nfea = fea_num;
		hdr.first_lbl = hdr.first_fea_col + hdr.nfea;
		hdr.nlbl = lbl_num;
		hdr.data_size = (__off64_t)0;   
		hdr.data_offset = (__off64_t)0;
		hdr.data_dim = 2;
		hdr.data_nrows = 0;
		hdr.data_ncols = fea_num + lbl_num + 2;
		hdr.sentt_size = (__off64_t)0;
		hdr.sentt_offset = (__off64_t)0;
		hdr.sentt_dim = 1;
		createEmptySentIndex();
		writeHeader();
		rewind();
	}
	return true;
}

void PFile::close()
{
	if(strcmp(file_mode, "rb") != 0)
	{
		writeSentIndex();
		writeHeader();
	}
	fclose(fp);
	fp = 0;
}

bool PFile::eof()
{
	return abs_frame_id >= hdr.nframes;
}

bool PFile::eos()
{
	return eof() || (p_sent_id != sent_id && !reset_eos_flag);	
}

void PFile::readHeader(char *file)
{
	// Check if the file is pfile
	fseeko64(fp, (__off64_t)0, SEEK_SET);
	char key[256];
	fscanf(fp, "%255s", key);
	if(strcmp(key, "-pfile_header") != 0)
	{
		char msg[1024];
		sprintf(msg, "ERROR: The file is not pfile: %s.", file);
		error(msg);
	}
	
	// Parse the header
	fseeko64(fp, (__off64_t)0, SEEK_SET);
	while(fscanf(fp, "%255s", key) && strcmp(key, "-end") != 0)
	{
		if(strcmp(key, "-pfile_header") == 0)
		{
			fscanf(fp, " version %u size %Lu", &hdr.version, &hdr.size);
                        //printf(">> %u %d %Lu\n", hdr.version, sizeof(hdr.size), hdr.size);
			if(hdr.version != 0)
			{
				char msg[1024];
				sprintf(msg, "ERROR: Only pfiles version 0 are supported - the %s file.", file);
				error(msg);
			}
		}
		else if(strcmp(key, "-num_sentences") == 0)
		{
			fscanf(fp, "%u", &hdr.nsents);
			//printf(">> %u\n", hdr.nsents);
		}
		else if(strcmp(key, "-num_frames") == 0)
		{
			fscanf(fp, "%u", &hdr.nframes);
			//printf(">> %u\n", hdr.nframes);
		}
		else if(strcmp(key, "-first_feature_column") == 0)
		{
			fscanf(fp, "%u", &hdr.first_fea_col);
			//printf(">> %u\n", hdr.first_fea_col);
		}
		else if(strcmp(key, "-num_features") == 0)
		{
			fscanf(fp, "%u", &hdr.nfea);
			//printf(">> %u\n", hdr.nfea);
		}
		else if(strcmp(key, "-first_label_column") == 0)
		{
			fscanf(fp, "%u", &hdr.first_lbl);
			//printf(">> %u\n", hdr.first_lbl);
		}
		else if(strcmp(key, "-num_labels") == 0)
		{
			fscanf(fp, "%u", &hdr.nlbl);
			//printf(">> %u\n", hdr.nlbl);
		}
		else if(strcmp(key, "-format") == 0)
		{
			// skip the formating string
			int ch;
			do
			{
				ch = fgetc(fp);
			}while(ch != '-');
			ungetc(ch, fp);
		}
		else if(strcmp(key, "-data") == 0)
		{
			fscanf(fp, " size %Lu offset %Lu ndim %u nrow %u ncol %u", &hdr.data_size, &hdr.data_offset, &hdr.data_dim, &hdr.data_nrows, &hdr.data_ncols);
			//printf(">> %u %u %u %u %u\n", hdr.data_size, hdr.data_offset, hdr.data_dim, hdr.data_nrows, hdr.data_ncols);
		}
		else if(strcmp(key, "-sent_table_data") == 0)
		{
			fscanf(fp, " size %Lu offset %Lu ndim %u", &hdr.sentt_size, &hdr.sentt_offset, &hdr.sentt_dim);
			//printf(">> %u %u %u\n", hdr.sentt_size, hdr.sentt_offset, hdr.sentt_dim);
		}
		//puts(key);
	}
}

bool PFile::read(float *features, unsigned int *labels)
{
	reset_eos_flag = 0;
	
	if(eof())
		return false;

	unsigned int r = fread(features, sizeof(float), hdr.nfea, fp);
	if(r != hdr.nfea)
		error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        vctSwapBytes((unsigned int *)features, hdr.nfea);
	r = fread(labels, sizeof(unsigned int), hdr.nlbl, fp);
	if(r != hdr.nlbl)
		error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        vctSwapBytes(labels, hdr.nlbl);

	abs_frame_id++;
	if(abs_frame_id < hdr.nframes)
	{
		p_sent_id = sent_id;
		//printf("zzz %u %u\n", p_sent_id, sent_id);
		r = fread(&sent_id, sizeof(unsigned int), 1, fp);
		if(r != 1)
			error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        	sent_id = swapBytes(sent_id);
		r = fread(&frame_id, sizeof(unsigned int), 1, fp);
		if(r != 1)
			error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        	frame_id = swapBytes(frame_id);
		//printf(">>>> %u %u\n", sent_id, frame_id);
		return true;
	}
	else
	{
		p_sent_id = sent_id;
		sent_id++; 
	}
	
        return false;
}

unsigned int PFile::swapBytes(unsigned int v)
{
	unsigned int b0, b1, b2, b3;
	b0 = v >> 24;
	b1 = (v >> 8) & 0x0000FF00;
	b2 = (v << 8) & 0x00FF0000;
		b3 = v << 24;
	return b0 | b1 | b2 | b3;
}

void PFile::vctSwapBytes(unsigned int *vct, int len)
{
        int i;
        for(i = 0; i < len; i++)
                vct[i] = swapBytes(vct[i]);
}

// Loads a sentence table from the pfile
// The sentence table starts at byte header_size + sentence_table_offset * sizeof(unsigned int)
// and it is single column of unsigned ints. The first value is zero and each other adds
// number of sentence frames  
void PFile::loadSentIndex(char *file)
{
	if(sent_idx)
		delete [] sent_idx;

	sent_idx_max_size = hdr.sentt_size;
	sent_idx = (unsigned int *)new unsigned int[hdr.sentt_size];

	fseeko64(fp, hdr.sentt_offset * (__off64_t)sizeof(unsigned int) + hdr.size, SEEK_SET);

	unsigned int r = fread(sent_idx, sizeof(unsigned int), (size_t)hdr.sentt_size, fp);
	if(r != (size_t)hdr.sentt_size)
	{
		char msg[1024];
		sprintf(msg, "ERROR: Can not read from the pfile. Probably corrupted pfile - %s.", file); 
		error(msg);
	}

	vctSwapBytes(sent_idx, (size_t)hdr.sentt_size);

/*	unsigned int i;
	for(i = 0; i < (size_t)hdr.sentt_size; i++)
		printf("%u %u\n", i, sent_idx[i]);*/
}

// Create the sentence table directly from the data. This function goes through the pfile and 
// if the sentence id changes, it add one entry to the sentence table.
void PFile::createSentIndexFromData()
{
	if(sent_idx)
		delete [] sent_idx;


	sent_idx_max_size = (__off64_t)hdr.nsents + 1;
	sent_idx = (unsigned int *)new unsigned int [sent_idx_max_size];

	unsigned int sid, psid = 0xFFFF, fid;
	unsigned int i, j = 0;

	for(i = 0; i < hdr.nframes; i++)
	{
		fseeko64(fp, hdr.data_offset + (__off64_t)hdr.size + (__off64_t)i * (__off64_t)(sizeof(unsigned int) * hdr.data_ncols), SEEK_SET);
		unsigned int r = fread(&sid, sizeof(unsigned int), 1, fp);   // Preload sentID and frameID
		if(r != 1)
			error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        	sid = swapBytes(sid);
		r = fread(&fid, sizeof(unsigned int), 1, fp);
		if(r != 1)
			error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        	fid = swapBytes(fid);
		if(sid != psid)
		{
			//printf("%d %d\n", (int)j, (int)i);

			sent_idx[j] = i;
			j++;
		}
		psid = sid;
	}	
	sent_idx[j] = i;
	hdr.sentt_size = (__off64_t)(hdr.nsents + 1);
	hdr.sentt_dim = 1;
}

void PFile::createEmptySentIndex()
{
	if(sent_idx)
		delete [] sent_idx;

	sent_idx_max_size = PFILE_SENT_TABLE_SIZE;
	sent_idx = (unsigned int *)new unsigned int [(size_t)sent_idx_max_size];

	sent_idx[0] = 0;
	sent_idx[1] = 0;	

	hdr.sentt_size = (__off64_t)(hdr.nsents + 1);
}

unsigned int PFile::getFramesNum(unsigned int sid)
{
	if(sid > hdr.nsents)
		error("ERROR: PFile::getFramesNum - The function is executed for nonexist sentence.");
	else if(sid == hdr.nsents)
		return 0;
	return (sent_idx[sid + 1] - sent_idx[sid]);
}

void PFile::goTo(unsigned int sid, unsigned int fid)
{
//	printf("%d %d\n", sid, fid);
	if(sid > hdr.nsents)
		error("ERROR: PFile::goTo - The sentence doe not exist.");
	if(fid > getFramesNum(sid) || (sid == hdr.nsents - 1 && fid == getFramesNum(sid) 
                                               && strcmp(file_mode, "rb") == 0))
		error("ERROR: PFile::goTo - The frame does not exist.");

	unsigned long fr = sent_idx[sid] + fid;

	fseeko64(fp, hdr.data_offset + (__off64_t)hdr.size + (__off64_t)fr * (__off64_t)(sizeof(unsigned int) * hdr.data_ncols), SEEK_SET);

	if(fid < getFramesNum(sid))
	{
		unsigned int r = fread(&sent_id, sizeof(unsigned int), 1, fp);   // Preload sentID and frameID
		if(r != 1)
			error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        	sent_id = swapBytes(sent_id);
		r = fread(&frame_id, sizeof(unsigned int), 1, fp);
		if(r != 1)
			error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        	frame_id = swapBytes(frame_id);
	}
	else
	{
		sent_id = sid;
		frame_id = fid;
	}

	act_sent_id = sid;
	abs_frame_id = sent_idx[sid] + frame_id;
}

unsigned int PFile::getMaxFramesNum()
{
	unsigned int i;
	unsigned int n, max = 0;
	for(i = 0; i < hdr.nsents; i++)
	{
		n = getFramesNum(i);
		if(n > max)
			max = n;
	}
	return max;
}

void PFile::writeHeader()
{
	// reserve enough space for header
	fseeko64(fp, (__off64_t)0, SEEK_SET);
	size_t i;
	for(i = 0; i < (size_t)hdr.size; i++)
		fputc(0, fp);

	// set remained fields
	hdr.data_nrows = hdr.nframes;
	hdr.data_size = (__off64_t)(hdr.nfea + hdr.nlbl + 2) * (__off64_t)hdr.nframes;   // +2 = sent_id & frame_id
		
	// write all fields
	fseeko64(fp, (__off64_t)0, SEEK_SET);

	fprintf(fp, "-pfile_header version %u size %Lu\n", hdr.version, hdr.size);
	fprintf(fp, "-num_sentences %u\n", hdr.nsents);
	fprintf(fp, "-num_frames %u\n", hdr.nframes);
	fprintf(fp, "-first_feature_column %u\n", hdr.first_fea_col);
	fprintf(fp, "-num_features %u\n", hdr.nfea);
	fprintf(fp, "-first_label_column %u\n", hdr.first_lbl);
	fprintf(fp, "-num_labels %u\n", hdr.nlbl);
	fprintf(fp, "-format dd");
	for(i = 0; i < hdr.nfea; i++)
		fprintf(fp, "f");
	for(i = 0; i < hdr.nlbl; i++)
		fprintf(fp, "d");
	fprintf(fp, "\n");
	fprintf(fp, "-data size %Lu offset %Lu ndim %u nrow %u ncol %u\n", hdr.data_size, hdr.data_offset, hdr.data_dim, hdr.data_nrows, hdr.data_ncols);
	fprintf(fp, "-sent_table_data size %Lu offset %Lu ndim %u\n", hdr.sentt_size, hdr.sentt_offset, hdr.sentt_dim);
	fprintf(fp, "-end\n");

//	printf("!!! %u\n" , hdr.nframes);
}

void PFile::write(float *features, unsigned int *labels)
{
	if(abs_frame_id != hdr.nframes && abs_frame_id >= sent_idx[sent_id + 1])
	{
		fprintf(stderr, "ERROR: Writing behind end of sentence\n");
		exit(1);
	}
	
	fseeko64(fp, hdr.data_offset + (__off64_t)hdr.size + (__off64_t)abs_frame_id * (__off64_t)(sizeof(unsigned int) * hdr.data_ncols), SEEK_SET);
	sent_id = act_sent_id;
	frame_id = abs_frame_id - sent_idx[sent_id]; 

	unsigned int sent_id_sw = swapBytes(sent_id);
	unsigned int frame_id_sw = swapBytes(frame_id);

	//printf("- %u %u\n", sent_id, frame_id);
	
	unsigned int w;
	w = fwrite(&sent_id_sw, sizeof(unsigned int), 1, fp);
	if(w != 1)
		error("ERROR: Can not write to the pfile.");

	w = fwrite(&frame_id_sw, sizeof(unsigned int), 1, fp);
	if(w != 1)
		error("ERROR: Can not write to the pfile.");

        vctSwapBytes((unsigned int *)features, hdr.nfea);
        vctSwapBytes((unsigned int *)labels, hdr.nlbl);

	w = fwrite(features, sizeof(float), hdr.nfea, fp);
	if(w != hdr.nfea)
		error("ERROR: Can not write to the pfile.");
	w = fwrite(labels, sizeof(unsigned int), hdr.nlbl, fp);
	if(w != hdr.nlbl)
		error("ERROR: Can not write to the pfile.");

	// Swap it back just for the case anybody wants to use it 
        vctSwapBytes((unsigned int *)features, hdr.nfea);
        vctSwapBytes((unsigned int *)labels, hdr.nlbl);	

	if(abs_frame_id == hdr.nframes)                   // if there it is a new frame in the file
	{
		hdr.nframes++;
		sent_idx[sent_id + 1]++;
	}
		
	abs_frame_id++;
	
	if(strcmp(file_mode, "r+b") == 0)
	{
		if(abs_frame_id < hdr.nframes)
		{
			p_sent_id = sent_id;
			//printf("zzz %u %u\n", p_sent_id, sent_id);
			unsigned int r = fread(&sent_id, sizeof(unsigned int), 1, fp);
			if(r != 1)
				error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        		sent_id = swapBytes(sent_id);
			r = fread(&frame_id, sizeof(unsigned int), 1, fp);
			if(r != 1)
				error("ERROR: Can not read from the pfile. Probably corrupted pfile.");
        		frame_id = swapBytes(frame_id);
		}
		else
		{
			p_sent_id = sent_id;
			sent_id++; 
		}
	}
}

void PFile::startNewSent()
{
//	printf("nsent %d\n", hdr.nsents);
	if(act_sent_id == hdr.nsents)
	{
		if((size_t)sent_idx_max_size <= hdr.nsents + 1)   // If the sentence index table is too short, extend it
		{
			int new_size = (size_t)sent_idx_max_size + PFILE_SENT_TABLE_SIZE;
			unsigned int *new_sent_idx = new unsigned int [new_size];
			int i;
			for(unsigned int i = 0; i < (size_t)sent_idx_max_size; i++)
				new_sent_idx[i] = sent_idx[i];
			for(i = (size_t)sent_idx_max_size; i < new_size; i++)
				new_sent_idx[i] = 0;
			sent_idx_max_size = (__off64_t)new_size;
			unsigned int *tmp = sent_idx;
			sent_idx = new_sent_idx;
			delete [] tmp;
		}
		sent_idx[act_sent_id + 1] = hdr.nframes;
		hdr.nsents++;
		hdr.sentt_size = (__off64_t)(hdr.nsents + 1);
		act_sent_id++;
	}
	else
	{
		act_sent_id++;
		goTo(act_sent_id, 0);
	}
	
}

void PFile::writeSentIndex()
{
	if(sent_idx == 0)
		free(sent_idx);

	hdr.sentt_offset = (__off64_t)(hdr.nfea + hdr.nlbl + 2) * (__off64_t)hdr.nframes; // right behind data
	fseeko64(fp, hdr.sentt_offset * (__off64_t)sizeof(unsigned int) + (__off64_t)hdr.size, SEEK_SET);

	vctSwapBytes(sent_idx, (size_t)hdr.sentt_size);	
	
	unsigned int w = fwrite(sent_idx, sizeof(unsigned int), (size_t)hdr.sentt_size, fp);
	if(w != (size_t)hdr.sentt_size)
	{
		char msg[1024];
		sprintf(msg, "ERROR: Can not write sentence index to the pfile."); 
		error(msg);
	}

	vctSwapBytes(sent_idx, (size_t)hdr.sentt_size);
}
