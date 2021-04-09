/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//most obvious input / output operations

#ifndef _IO_h
#define _IO_h


#include "containers.h"
#include<cstdio>

typedef std::map<std::string, shared_ptr<DNA>> DNAmap;


//void openOutFiles(string files, string fmt,string );
//void prepareOutFiles(OptContainer& cmdArgs);
//void read_fastq(OptContainer& cmdArgs, OutputStreamer* MD,string fileS);
bool read_paired(OptContainer& cmdArgs, shared_ptr<OutputStreamer> MD, 
	shared_ptr<InputStreamer>,bool , int Nthreads);
bool read_paired_DNAready(vector< shared_ptr<DNA>> tdn, bool MIDuse,
	shared_ptr<OutputStreamer> MD, int curThread);
bool read_paired_STRready(vector<vector< string >> tdn, bool MIDuse,
	shared_ptr<OutputStreamer> MD, int curThread, bool keepPairHd, qual_score FastqVer);
//shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2,shared_ptr<DNA> MIDseq,

//bool read_tripple(OptContainer& cmdArgs, OutputStreamer* MD, InputStreamer*);

void separateByFile(Filters* mainFilter, OptContainer& cmdArgs);

void threadAnalyzeDNA(shared_ptr<DNA> tdn, shared_ptr<OutputStreamer> MD,int thrCnt);
//void trippleThreadAnalyzeDNA(shared_ptr<OutputStreamer> MD, shared_ptr<DNA> tdn,shared_ptr<DNA>dnaTemp2,shared_ptr<DNA> MIDseq,bool changePHead);//,int thrCnt=0);

void read_single(OptContainer& cmdArgs, shared_ptr<OutputStreamer> MD, shared_ptr<InputStreamer> IS);

bool readCmdArgs(int argc, char* argv[],OptContainer& cmdArgs);




//specialized functions .. end sdm after execution
void rewriteNumbers(OptContainer& cmdArgs);


void Announce_sdm();
void help_head();
void general_help();
void printCmdsHelp();
void printOptionHelp();
void printMapHelp();
void printVersion();


//bool readCmdArgs(int argc, char* argv[],base_map<char*, char*, lstr>& cmdArgs);
#endif