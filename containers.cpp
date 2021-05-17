/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand
email: Falk.Hildebrand@gmail.com

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


#include <mutex>
#include "containers.h"

using namespace std;


void trim(string& str){
	// trim trailing spaces
	size_t endpos = str.find_last_not_of(" \t");
	if( string::npos != endpos )
	{
	    str = str.substr( 0, endpos+1 );
	}

	// trim leading spaces
	size_t startpos = str.find_first_not_of(" \t");
	if( string::npos != startpos )
	{
	    str = str.substr( startpos );
	}
}
//from http://stackoverflow.com/questions/8888748/how-to-check-if-given-c-string-or-char-contains-only-digits
bool is_digits(const std::string &str)
{
	return std::all_of(str.begin(), str.end(), ::isdigit); // C++11
}


//



ReadSubset::ReadSubset(const string inf, const string default_outfile):
RemainderStrPos(-1), newHD(0), outFiles(0), outFilesIdx(0) {
	string line;
	ifstream in(inf.c_str());
	if (!in){
		cerr << "Could not find " << inf << " read subset file. Exiting.\n"; exit(90);
	}
	int ini_ColPerRow(0), cnt(0), skips(0);

	//check read subset format
	while (!safeGetline(in, line).eof()) {
		if (line.substr(0, 1) == "#"){ skips++; continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, '\t')) {
			ColsPerRow++;
		}

		if (cnt == 0){
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow){
				cerr << "Number of columns in read subset file on line " << cnt + skips << " is " << ColsPerRow << ". Expected " << ini_ColPerRow << " columns.\n";
				exit(91);
			}
			if (cnt > 1000){
				break;
			}
		}
		cnt++;
	}
	if (ini_ColPerRow == 0){
		cerr << "Read Subset File exists, but appears to be badly formated (0 columns detected). Exiting\n"; exit(92);
	}
	if (cnt == 0){
		cerr << "Read Subset File exists, but appears to be badly formated (0 lines detected). Exiting\n"; exit(92);
	}
	in.clear();
	in.seekg(0, ios::beg);
	//extract read subset content
	cnt = 0;
	map<string, int> tmpFiles;
	map<string, int>::iterator tmpFilIT;
	unordered_map <string, int>::iterator TarIT;
	//parameters were set for in matrix, now read line by line
	while (!safeGetline(in, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (cnt != 0 && line.substr(0, 1) == "#"){ continue; }
		if (line.length() < 5){ continue; }
		stringstream ss; string segments; int cnt2 = 0;
		ss << line;
		while (getline(ss, segments, '\t')) {
			if (cnt2 == 0){
				//free target RD from paired end info
				remove_paired_info(segments);
				TarIT = Targets.find(segments);
				if (TarIT == Targets.end()){
					Targets[segments] = cnt;
				} else {
					break;
				}
			}
			else if (cnt2 == 1){
				newHD.push_back(segments);
			}
			else if (cnt2 == 2){
				uint Idx(0);
				//1: test if outfile already exists
				tmpFilIT = tmpFiles.find(segments);
				if (tmpFilIT != tmpFiles.end()){
					Idx = (*tmpFilIT).second;
				}
				else {//create new key
					Idx = (uint) outFiles.size();
					outFiles.push_back(segments);
					tmpFiles[segments] = Idx;
				}
				//2: add the index to outfile to this position
				outFilesIdx.push_back(Idx);
			}
			cnt2++;
		}
		if (cnt2 == 1) {//no specific header
			newHD.push_back("");
			tmpFilIT = tmpFiles.find("Default");
			uint Idx(0);
			if (tmpFilIT != tmpFiles.end()) {
				Idx = (*tmpFilIT).second;
			} else {//create new key
				Idx = (uint)outFiles.size();
				outFiles.push_back("Default");
				tmpFiles["Default"] = Idx;
			}

			outFilesIdx.push_back(Idx);
		}
		if (cnt2>0) { cnt++; }
	}
	//is extra file info in read subset file?
	if (ini_ColPerRow < 3){
		outFilesIdx.resize(Targets.size(), 0);
		outFiles.resize(1, default_outfile);
	}

}

void ReadSubset::findMatches(shared_ptr<InputStreamer> IS, shared_ptr<OutputStreamer> MD, bool mocatFix) {
	vector<shared_ptr<DNA>> match; //shared_ptr<DNA> match2(NULL);
	bool cont(true), cont2(true);
	int idx(0);
	int pairs = IS->pairNum();
	unordered_map <string, int>::iterator SEEK;
	bool b_doHD = newHD.size() > 0;
	//bool sync(false);//meaningless placeholder
	while (cont) {
		match = IS->getDNAMC();
		if (match[0] == NULL) { cont = false; break; }
		//if (pairs > 1) {			match2 = IS->getDNA(cont2, 1);		}
		string curID = match[0]->getPositionFreeId();
		if (mocatFix && curID.length()>5 ) {
			curID = header_stem(curID);
			//cerr << curID << endl;
			//exit(0);
		} 
		//cerr << curID << endl;
		SEEK = Targets.find(curID);
		if (SEEK == Targets.end()) {//no hit
			if (RemainderStrPos != -1) {
				MD->writeSelectiveStream(match[0], 0, RemainderStrPos);
				if (pairs > 1) {
					MD->writeSelectiveStream(match[1], 1, RemainderStrPos);
				}
			} else {// nothing written, but still need to delete
//				delete match; if (pairs > 1) {delete match2;}
			}
			continue;
		} 
		//serious work 
		//cerr << "H";
		idx = SEEK->second;
		if (b_doHD && newHD[idx] != "" && pairs>1) {
            match[0]->setHeader(newHD[idx] + "#1:0");
            match[1]->setHeader(newHD[idx] + "#2:0");
		}
		MD->writeSelectiveStream(match[0], 0, outFilesIdx[idx]);
		if (pairs > 1) {
			MD->writeSelectiveStream(match[1], 1, outFilesIdx[idx]);
		}

		//finished, clean up to reduce search space (faster??)
		Targets.erase(SEEK);
		if (Targets.size() == 0 && RemainderStrPos==-1) {
			return;
		}
	}
	cerr << Targets.size() << " seqs remaining (not found in current file)\n";
}





OutputStreamer::OutputStreamer(Filters* fil, OptContainer& cmdArgs,
        std::ios_base::openmode writeStatus, 
		shared_ptr<ReadSubset> RDSset,
		int numThreads, string fileExt, int forceFmt) :
        MFil(fil), subFilter(0), //DNAsP1(0), DNAsP2(0), DNAsS1(0), DNAsS2(0),
        //DNAsNoHead(0), DNAsP1_alt(0), DNAsP2_alt(0), DNAsS1_alt(0), DNAsS2_alt(0),
        suppressOutWrite(0), write2File(true),// mem_used(false),
        DNAinMem(0), writeThreadStatus(0),
        fastQver(fil->getuserReqFastqVer()),
        fastQoutVer(fil->getuserReqFastqOutVer()), BWriteQual(false),
        BWriteFastQ(false), b_multiOutStream(false), pairedSeq(-1),
		b_changeFQheadVer(false),
		b_checkedHeaderChange(false),
        b_oneLinerFasta(false), b_doDereplicate(false), b_writeGreenQual(true), b_writeYellowQual(true),
        maxReadsPerOFile(fil->maxReadsOutput()), 
		demultiBPperSR(fil->getDemultiBPperSR()), 
		BPwrittenInSR(0), BPwrittenInSRmerg(0),
		ReadsWritten(fil->writtenReads()),
        maxRdsOut(-1), stopAll(false),
        leadingOutf(""), locCmdArgs(cmdArgs), dereplicator(NULL), cntDerep(0), wrMode(ios::out),
        sFile(0), qFile(0), fqFile(0),
		of_merged_fq(NULL),
        //sFileStr(0), qFileStr(0), fqFileStr(0), fqNoBCFile(0), 
		totalFileStrms(0),
        bDoDemultiplexIntoFiles(false), demultiSinglFiles(0), //demultiSinglFilesF(0),
		demultiMergeFiles(0),
        onlyCompletePairsDemulti(false),
		b_merge_pairs_demulti_(false), b_merge_pairs_filter_(false),
		b_merge_pairs_(false), merger(1,nullptr),
		Nthrds(numThreads)


/*	qFilePos(0), sFilePos(0), fqFilePos(0),
	qFile2Pos(0), sFile2Pos(0), fqFile2Pos(0),//second pair_
	qFileSPos(0), sFileSPos(0), fqFileSPos(0),//singleton
	qFileS2Pos(0), sFileS2Pos(0), fqFileS2Pos(0)//singleton*/
{
	
 	pairedSeq = MFil->isPaired();
	if (cmdArgs.find("-suppressOutput") != cmdArgs.end()) {
		suppressOutWrite = atoi( cmdArgs["-suppressOutput"].c_str() );
	}
	if (cmdArgs.find("-pairedDemulti") != cmdArgs.end() && cmdArgs["-pairedDemulti"] == "1") { //only demultiplex proper pairs
		onlyCompletePairsDemulti = true;
	}
	if (cmdArgs.find("-oneLineFastaFormat") != cmdArgs.end() && cmdArgs["-oneLineFastaFormat"] == "1") {
		b_oneLinerFasta = true;
	}

	
	if (suppressOutWrite == 3 || suppressOutWrite == 1){
		b_writeGreenQual = false;
	}
	if (suppressOutWrite == 3 || suppressOutWrite == 2){
		b_writeYellowQual = false;
	}
	if (fil->doSubselReads() && RDSset != nullptr ) {
		this->openSeveralOutstreams(locCmdArgs,RDSset, writeStatus);
	} else {//standard one file output stream
		this->openOutStreams(locCmdArgs, MFil->getFileIncrementor(), 
			writeStatus, fileExt, forceFmt);
	}
	if (cmdArgs["-o_fastq_noBC"] != "") {
		openNoBCoutstrean(cmdArgs["-o_fastq_noBC"]);
	}
	if (cmdArgs.find("-o_demultiplex") != cmdArgs.end()) {//demulti: always write mode out
		generateDemultiOutFiles(cmdArgs["-o_demultiplex"],fil, ios::out);
	}
	if (cmdArgs.find("-merge_pairs_filter") != cmdArgs.end()
		&& cmdArgs["-merge_pairs_filter"] == "1") {
		b_merge_pairs_filter_ = true;
	}
	if (cmdArgs.find("-merge_pairs_demulti") != cmdArgs.end() 
		&& cmdArgs["-merge_pairs_demulti"] == "1") {
		b_merge_pairs_demulti_ = true;
	}
	if (b_merge_pairs_demulti_ || b_merge_pairs_filter_) {
		b_merge_pairs_ = true;
	}
	setSubfilters(Nthrds);

	//threads = futures(num_threads);
	//if (pairedSeq){cerr<<"paired MFil\n";}
}


OutputStreamer::~OutputStreamer(){
		//delete MFil;
	cdbg("Destr OutputStreamer" );
	//delAllDNAvectors();
	//cdbg(".. done\n" );
	dmltMTX.lock();
	this->closeOutStreams(true);
	closeOutFilesDemulti();
	for (uint i=0; i<subFilter.size();i++){
		delete subFilter[i];
	}
	dmltMTX.unlock();
	cdbg("Subfilters deleted, streams closed\n" );
	//delete optim;
}
void OutputStreamer::delAllDNAvectors() {
#ifdef DEBUG
	cerr << "cleaning MD..";
#endif
	//hard delete all saved data, no write to HDD
	for (size_t i = 0; i < sFile.size(); i++) {
		for (size_t j = 0; j < sFile[i].size(); j++) {
			if (sFile[i][j] != nullptr) {
				sFile[i][j]->emptyStream();
			}
		}
	}
	for (size_t i = 0; i < qFile.size(); i++) {
		for (size_t j = 0; j < qFile[i].size(); j++) {
			if (qFile[i][j] != nullptr) {
				qFile[i][j]->emptyStream();
			}
		}
	}
	for (size_t i = 0; i < fqFile.size(); i++) {
		for (size_t j = 0; j < fqFile[i].size(); j++) {
			if (fqFile[i][j] != nullptr) {
				fqFile[i][j]->emptyStream();
			}
		}
	}
	/*
	DNAsP1.resize(0);
	DNAsP2.resize(0);
	DNAsS1.resize(0);
	DNAsS2.resize(0);
	DNAsNoHead.resize(0);
	DNAsP1_alt.resize(0);
	DNAsP2_alt.resize(0);
	DNAsS1_alt.resize(0);
	DNAsS2_alt.resize(0);
#ifdef DEBUG
	cerr << " finished\n";
#endif

*/
}

void OutputStreamer::analyzeDNA(shared_ptr<DNA> d, int FilterUse, int pair, int& idx, int thr) {
	if (!d){
		return ;
	}
	Filters * curFil = this->getFilters(thr);

	if ( !curFil->doFilterAtAll() ) {
		bool isP1 = max(0, pair) == 0;
		if (idx < 0 && !isP1 && !curFil->doubleBarcodes()) {
			;
		} else if (idx < 0) {
            idx = curFil->detectCutBC(d, isP1); //still need to check for BC
		}
		if (idx >= 0) { //prevent second read pair_ from being flagged as true
			d->setPassed(true);
		}
		return;
	}

	if (curFil->secondaryOutput() ){
		//pricipally safe to call from thread
        curFil->checkYellowAndGreen(d, pair, idx);
	} else if (FilterUse==-1){
		curFil->check(d, false, pair, idx);
	} else {
		cerr << "Invalid control path in analyzeDNA\n"; exit(55);
		//Filters* f = subFilter[FilterUse];
		//subFilter[FilterUse].check(d, false, pair_, idx);
        //f->check(d, false, pair, idx);
        
        
    }
	if (idx == -1 ) {
		d->setBarcodeDetected(false); 
	}

	//count this as failure if BC was present
	//d->prepareWrite(fastQoutVer);
}

vector<bool> OutputStreamer::analyzeDNA(shared_ptr<DNA> p1, shared_ptr<DNA> p2, shared_ptr<DNA> mid,
	bool changePHead, int FilterUse) {
	cerr << "deprecated analuze DNA"; exit(2323);
	vector<bool> ret(2, true);
	//1st: check if DNA pointer valid
	if (p1 == NULL){
		ret[0] = false;
	} 
	if (p2 == NULL){
		ret[1] = false;
	} 
	if (mid == NULL){//no MID? kill
		ret[1] = false; ret[0] = false; return ret;
	}
	else {
		mid->setMIDseq(true);
	}

	//2nd: check if basic DNA signatures are valid

	//3rd: check if all quality criteria are met, BC is true etc.
	//if (FilterUse == -1){
	//ret = curFil->check_pairs(p1, p2, mid, ret, changePHead);
	//}	else { //mutithread, to avoid race conditions etc use separate filter
	//	ret = subFilter[FilterUse]->check_pairs(p1, p2, mid, ret, changePHead);
	//}
	//if (ret[0]){ p1->prepareWrite(fastQoutVer); }
	//if (ret[1]){	p2->prepareWrite(fastQoutVer);}

	return ret;
}

bool OutputStreamer::checkFastqHeadVersion(shared_ptr<DNA> d,bool disable) {
	if (!b_checkedHeaderChange) {
		return b_changeFQheadVer;
	}
	b_checkedHeaderChange = true;
	if (pairedSeq == 1) { b_changeFQheadVer = false; return false; }
	drpMTX.lock();
	b_changeFQheadVer = false;
	int fastQheadVer = 0;
	string head = d->getId();
	int shouldVer = MFil->FQheadV();
	if (head.find("/1") != string::npos || head.find("/2") != string::npos){
		fastQheadVer = 1;
	}
	if (head.find(" 1:") != string::npos|| head.find(" 2:") != string::npos){
		fastQheadVer = 2;
	}
	if (shouldVer!=0 && shouldVer != fastQheadVer){
		b_changeFQheadVer=true;
	}
	bool ret = b_changeFQheadVer;
	if (disable){
		b_changeFQheadVer = false;
	}
	drpMTX.unlock();
	return ret;
}

void OutputStreamer::closeOutFilesDemulti() {
//    cout << "FILES TO t DEMULTIPLEX OUT " << demultiSinglFiles.size() << endl;
    //cout << "domeulti: " << bDoDemultiplexIntoFiles << endl;

	cdbg("OutputStreamer::closeOutFilesDemulti");
	if (!bDoDemultiplexIntoFiles) { return; }
	for (size_t i = 0; i < demultiSinglFiles.size(); i++) {
	    //cout << " do fancy demultiplexing for " << i << endl;
		if (demultiSinglFiles[i][0] != nullptr) {
		    //cout << "yaas" << endl;
		    delete demultiSinglFiles[i][0];
			demultiSinglFiles[i][0] = nullptr;
		}
		if (demultiSinglFiles[i][1] != nullptr)
		{
		    delete demultiSinglFiles[i][1];
			demultiSinglFiles[i][1] = nullptr;
		}
	}
	for (size_t i = 0; i < demultiMergeFiles.size(); i++) {
		if (demultiMergeFiles[i] != nullptr) {
			delete demultiMergeFiles[i];
			demultiMergeFiles[i] = nullptr;
		}
	}
    //cout << " closeOutFilesDemulti" << endl;
    return;
}

/*
void OutputStreamer::writeAllStoredDNA(){
#ifdef DEBUG
	printStorage();
	cerr << "Writting stored DNA" << DNApairStreams[0].size() <<" " << DNApairStreams[1].size()<<endl;
#endif

#ifdef _THREADED
	if (writeThreadStatus>0){
		if (writeThreadStatus>1){wrThread.join();}
		writeThreadStatus++;
		wrThread = thread(&OutputStreamer::writeAllStoredDNA2t,this);
	} else {
		writeAllStoredDNA2();
	}
#else
		writeAllStoredDNA2();
#endif
}
void OutputStreamer::writeAllStoredDNA2(){
	mem_used=false;
	//
	//this is handled by ofbufstream now..
	
	if (b_writeGreenQual) {

		//
		if (DNAsP1.size() != DNAsP2.size()) {
			for (unsigned int i = 0; i < DNAsP1.size(); i++) {
				writeAndDel(DNAsP1[i], 0);
				ReadsWritten++;
			}
			for (unsigned int i = 0; i < DNAsP2.size(); i++) {
				writeAndDel(DNAsP2[i], 1);
			}
		} else {
			for (unsigned int i = 0; i < DNAsP1.size(); i++) {
				writeAndDel(DNAsP1[i], 0);
				writeAndDel(DNAsP2[i], 1);
				ReadsWritten++;
			}
		}
		for (unsigned int i = 0; i < DNAsS1.size(); i++) {
			writeAndDel(DNAsS1[i], 2);
		}

		for (unsigned int i = 0; i < DNAsS2.size(); i++) {
			writeAndDel(DNAsS2[i], 3);
		}
		DNAsP1.resize(0); DNAsP2.resize(0);
		DNAsS1.resize(0); DNAsS2.resize(0);
	} 
	if (b_writeYellowQual) {

		//Xtra file
		if (DNAsP1_alt.size() != DNAsP2_alt.size()) {
			for (unsigned int i = 0; i < DNAsP1_alt.size(); i++) {
				writeAndDel(DNAsP1_alt[i], 0);
				ReadsWritten++;
			}
			for (unsigned int i = 0; i < DNAsP2_alt.size(); i++) {
				writeAndDel(DNAsP2_alt[i], 1);
			}
		} else {
			for (unsigned int i = 0; i < DNAsP1_alt.size(); i++) {
				writeAndDel(DNAsP1_alt[i], 0);
				writeAndDel(DNAsP2_alt[i], 1);
				ReadsWritten++;
			}
		}
		for (unsigned int i = 0; i < DNAsS1_alt.size(); i++) {
			writeAndDel(DNAsS1_alt[i], 2);
		}

		for (unsigned int i = 0; i < DNAsS2_alt.size(); i++) {
			writeAndDel(DNAsS2_alt[i], 3);
		}
		DNAsP1_alt.resize(0); DNAsP2_alt.resize(0);
		DNAsS1_alt.resize(0); DNAsS2_alt.resize(0);
	} 
	//just to be on the safe side..
	delAllDNAvectors();
#ifdef DEBUG
	cerr <<  " .. Finished" << endl;
#endif
}
*/


#ifdef _THREADED
void OutputStreamer::writeAllStoredDNA2t(){
    std::lock_guard<std::mutex> guard(mutex);
	write2File=true;mem_used=false;
	if (DNAsP1.size() != DNAsP2.size()){
		for (unsigned int i=0; i<DNAsP1.size();i++){
			writeAndDel(DNAsP1[i],0);
		}
		for (unsigned int i=0; i<DNAsP2.size();i++){
			writeAndDel(DNAsP2[i],1);
		}
	}else {
		for (unsigned int i=0; i<DNAsP1.size();i++){
			writeAndDel(DNAsP1[i],0);
			writeAndDel(DNAsP2[i],1);
		}
	}
	for (unsigned int i=0; i<DNAsS1.size();i++){
		writeAndDel(DNAsS1[i],2);
	}

	for (unsigned int i=0; i<DNAsS2.size();i++){
		writeAndDel(DNAsS2[i],3);
	}
	DNAsP1.resize(0);DNAsP2.resize(0);
	DNAsS1.resize(0);DNAsS2.resize(0);
}
#endif


void OutputStreamer::write2Demulti(shared_ptr<DNA> d1, shared_ptr<DNA> d2, int BCoffset,
	int curThread) {
	if (!this->Demulti2Fls()) {
		return;
	}
	
	bool demultFini = demultiBPperSR > 0 && BPwrittenInSR > demultiBPperSR;
	bool demultMrgFini = demultiBPperSR > 0 && BPwrittenInSRmerg > demultiBPperSR;
	if (demultFini && demultMrgFini) {
		return;
	}
	int idx = d1->getBarcodeNumber() - BCoffset; //correct for BC offset as well..

	if (idx < 0) {
		return;
	}
	bool green1(d1->isGreenQual());
	bool green2(false);
	if (d2 != nullptr) {
		green2 = d2->isGreenQual();
	}
	bool mergeWr(false);

	//merging attempts/prep, requires paired reads
	if (this->isPEseq() ==2 && (green1 || green2) && b_merge_pairs_demulti_ && d1->merge_seed_pos_ > 0) {
		shared_ptr<DNA> dna_merged = merger[curThread]->merge(d1, d2);
		if (dna_merged) {
			// write out merged DNA
			dna_merged->prepareWrite(fastQoutVer);
			//dmltMTX.lock();
			*(demultiMergeFiles[idx]) << dna_merged->writeFastQ();
			BPwrittenInSRmerg += d1->length();
			//dmltMTX.unlock();
			//mergeWr = true;
		}
	}
	if (mergeWr) {
		return;
	}

	//check for paired reads
	if ((onlyCompletePairsDemulti && !green1) ||
		(onlyCompletePairsDemulti && !green2)) {
		return;
	}
	/*
	string x = d2->getShortId();
	if (x.find("M04428:252:000000000-BM3M5:1:2106:11089:1656") != string::npos) {
		bool x = true; // DEBUG
	}
	/**/

	if (green1) {
		//cout << "trouble" << endl;
		d1->prepareWrite(fastQoutVer);
		//dmltMTX.lock();
		*(demultiSinglFiles[idx][0]) << d1->writeFastQ( false);
		BPwrittenInSR += d1->length();
		//dmltMTX.unlock();
	}
	if (green2) {
		//cout << "trouble" << endl;
		d2->prepareWrite(fastQoutVer);
		//dmltMTX.lock();
		*(demultiSinglFiles[idx][1]) << d2->writeFastQ( false);
		//dmltMTX.unlock();
	}


}

void OutputStreamer::write2Demulti(shared_ptr<DNA> d, int p, int BCoffset) {//this->getBCoffset()
	if (!this->Demulti2Fls()) {
		return;
	}
	if (demultiBPperSR > 0 && BPwrittenInSR > demultiBPperSR) {
		return;
	}
	int idx = d->getBarcodeNumber() - BCoffset; //correct for BC offset as well..

	if (idx < 0 || !d->isGreenQual()) {
		return;
	}

	d->prepareWrite(fastQoutVer);
	*(demultiSinglFiles[idx][p]) << d->writeFastQ(false);
	BPwrittenInSR += d->length();

}


void OutputStreamer::generateDemultiOutFiles(string path, Filters* fil, std::ios_base::openmode writeStatus) {
	//fill in demultiSinglFiles vector
	vector<ofbufstream*> empVec(2, NULL);
	//vector<string> empVec2(2, "");

	bool doMC = Nthrds > 1;

	struct stat info;
	if (stat(path.c_str(), &info) != 0) {
		cerr << "Output path for demultiplexed files does not exist, please create this directory:\n" << path << endl;
		exit(833);
	}
	else if (info.st_mode & S_IFDIR) { // S_ISDIR() doesn't exist on my windows 
		cerr << "Writing demultiplexed files to: " << path << endl;// printf("%s is a directory\n", path);
	}
	else {
		cerr << path << " is no directory\n"; exit(834);
	}
	path += "/";

	bDoDemultiplexIntoFiles = true;
	bool openOstreams = true; uint ostrCnt(0);
	size_t bufS = 30000;
	
	demultiSinglFiles.resize(fil->SampleID.size(), empVec);
	demultiMergeFiles.resize(fil->SampleID.size(), nullptr);
	//demultiSinglFilesF.resize(fil->SampleID.size(), empVec2);
	for (size_t i = 0; i < fil->SampleID.size(); i++) {
		//actually needs to know if paired files..
		//if (ostrCnt > maxFileStreams) {openOstreams = false;}
		if (pairedSeq == 1 || pairedSeq == -1) {
			string nfile = path + fil->SampleID[i] + ".fq";
			if (openOstreams) { demultiSinglFiles[i][0] = new ofbufstream(nfile.c_str(), writeStatus, doMC, bufS); }
			//demultiSinglFilesF[i][0] = nfile;
			ostrCnt++;
		}
		else {
			string nfile = path + fil->SampleID[i] + ".1.fq";
			if (openOstreams) { demultiSinglFiles[i][0] = new ofbufstream(nfile.c_str(), writeStatus, doMC,bufS*0.8); }
			//demultiSinglFilesF[i][0] = nfile;
			nfile = path + fil->SampleID[i] + ".2.fq";
			if (openOstreams) { demultiSinglFiles[i][1] = new ofbufstream(nfile.c_str(), writeStatus, doMC,bufS*1.2); }
			//demultiSinglFilesF[i][1] = nfile;
			nfile = path + fil->SampleID[i] + ".merg.fq";
			if (openOstreams) { demultiMergeFiles[i] = new ofbufstream(nfile.c_str(), writeStatus, doMC, bufS); }
			
			ostrCnt += 2;
		}
	}

}
void OutputStreamer::incrementOutputFile(){
    //TODO has to be made threadsafe
	MFil->incrementFileIncrementor();
	this->closeOutStreams(true);
	//this is definetely a new file

    //TODO needs to be threadsafe
	this->openOutStreams(locCmdArgs,MFil->getFileIncrementor(),ios_base::out);
	ReadsWritten = 0;
}


// CHECK THIS METHOD AFTER IMPLEMENTING MULTIPLE STATS OBJECTs
//
//void Filters::addToStatistics(shared_ptr<DNA> d, Statistics &statistics) {
//    //here should be the only place to count Barcodes!
//    statistics.total2++;
//    if (d->isGreenQual() || d->isYellowQual()) {
////        this->DNAstatLQ(d, easyPair, d->isYellowQual());
//        statistics.totalSuccess++;
//    } else {
//        statistics.totalRejected++;
//    }
//
//    //some general stats that always apply:
//    if (d->QualCtrl.PrimerFail) {
//        statistics.PrimerFail++;
//    }
//    if (d->QualCtrl.PrimerRevFail) {
//        statistics.PrimerRevFail++;
//    }
//    if (d->QualCtrl.minLqualTrim) {
//        statistics.minLqualTrim++;
//    }
//    if (d->QualCtrl.TagFail) {
//        statistics.TagFail++;
//    }
//    if (d->QualCtrl.fail_correct_BC) {
//        statistics.fail_correct_BC++;
//    }
//    if (d->QualCtrl.suc_correct_BC) {
//        statistics.suc_correct_BC++;
//    }
//    if (d->QualCtrl.RevPrimFound) {
//        statistics.RevPrimFound++;
//    }
//    if (d->QualCtrl.QWinTrimmed || d->QualCtrl.AccErrTrimmed) {
//        statistics.Trimmed++;
//    }
//    if (d->getTA_cut()) {
//        statistics.adapterRem++;
//    }
//
////exit(0);
//
//    if (d->isGreenQual() || d->isYellowQual()) {
////        countBCdetected(d->getBarcodeNumber(), easyPair, false);
//        //and register as success
//    } else {
//        if (d->getBarcodeDetected()) {
//            //DNA is no longer useful
////            failedStats2(d, easyPair);
//        }
//        //delete d;
//        if (d->QualCtrl.AvgQual) {
//            statistics.AvgQual++;
//        }
//        if (d->QualCtrl.minL) {
//            statistics.minL++;
//        }
//        if (d->QualCtrl.maxL) {
//            statistics.maxL++;
//        }
//        if (d->QualCtrl.HomoNT) {
//            statistics.HomoNT++;
//        }
//        if (d->QualCtrl.MaxAmb) {
//            statistics.MaxAmb++;
//        }
//        if (d->QualCtrl.BinomialErr) {
//            statistics.BinomialErr++;
//        }
//        if (d->QualCtrl.QualWin) {
//            statistics.QualWin++;
//        }
//    }
//    if (d->isDereplicated()) {
//        if (d->getBarcodeDetected() && !d->isGreenQual() && !d->isYellowQual()) {
//            this->statAddDerepBadSeq(d->getBarcodeNumber());
//        }
//    }
//}

void Filters::addDNAtoCStats(shared_ptr<DNA> d,int Pair) {
	//here should be the only place to count Barcodes!
	int easyPair = Pair < 3 ? Pair - 1 : Pair - 3;
	
	csMTX[easyPair]->lock();
	collectStatistics[easyPair]->total2++;


	if (d->isGreenQual() || d->isYellowQual()) {
		this->DNAstatLQ(d, easyPair, d->isYellowQual());
		collectStatistics[easyPair]->totalSuccess++;
	} else {
		collectStatistics[easyPair]->totalRejected++;
	}

	//some general stats that always apply:
	if (d->QualCtrl.PrimerFail) {
		collectStatistics[easyPair]->PrimerFail++;
	}
	if (d->QualCtrl.PrimerRevFail) {
		collectStatistics[easyPair]->PrimerRevFail++;
	}
	if (d->QualCtrl.minLqualTrim) {
		collectStatistics[easyPair]->minLqualTrim++;
	}
	if (d->QualCtrl.TagFail) {
		collectStatistics[easyPair]->TagFail++;
	}
	if (d->QualCtrl.fail_correct_BC) {
		collectStatistics[easyPair]->fail_correct_BC++;
	}
	if (d->QualCtrl.suc_correct_BC) {
		collectStatistics[easyPair]->suc_correct_BC++;
	}
	if (d->QualCtrl.RevPrimFound) {
		collectStatistics[easyPair]->RevPrimFound++;
	}
	if (d->QualCtrl.QWinTrimmed || d->QualCtrl.AccErrTrimmed) {
		collectStatistics[easyPair]->Trimmed++;
	}
	if (d->getTA_cut()) {
		collectStatistics[easyPair]->adapterRem++;
	}

//exit(0);
	
	if (d->isGreenQual() || d->isYellowQual()) {
		countBCdetected(d->getBarcodeNumber(), easyPair, false);
		//and register as success
	} else {
		if (d->getBarcodeDetected()) {
			//DNA is no longer useful
			failedStats2(d, easyPair);
		}
		//delete d; 
		if (d->QualCtrl.AvgQual) {
			collectStatistics[easyPair]->AvgQual++;
		}
		if (d->QualCtrl.minL) {
			collectStatistics[easyPair]->minL++;
		}
		if (d->QualCtrl.maxL) {
			collectStatistics[easyPair]->maxL++;
		}
		if (d->QualCtrl.HomoNT) {
			collectStatistics[easyPair]->HomoNT++;
		}
		if (d->QualCtrl.MaxAmb) {
			collectStatistics[easyPair]->MaxAmb++;
		}
		if (d->QualCtrl.BinomialErr) {
			collectStatistics[easyPair]->BinomialErr++;
		}
		if (d->QualCtrl.QualWin) {
			collectStatistics[easyPair]->QualWin++;
		}
	}
	if (d->isDereplicated()) {
		if (d->getBarcodeDetected() && !d->isGreenQual() && !d->isYellowQual()) {
			this->statAddDerepBadSeq(d->getBarcodeNumber());
		}
	}
	csMTX[easyPair]->unlock();
}



void Filters::addDNAtoCStatsMT(shared_ptr<DNA> d, int pair, int thread_id) {
    //here should be the only place to count Barcodes!
    int easyPair = pair < 3 ? pair - 1 : pair - 3;

    auto &statistics = statistics_[thread_id].main_read_stats_;

    statistics[easyPair]->total2++;


    if (d->isGreenQual() || d->isYellowQual()) {
        this->DNAstatLQ(d, easyPair, d->isYellowQual());
        statistics[easyPair]->totalSuccess++;
    } else {
        statistics[easyPair]->totalRejected++;
    }

    //some general stats that always apply:
    if (d->QualCtrl.PrimerFail) {
        statistics[easyPair]->PrimerFail++;
    }
    if (d->QualCtrl.PrimerRevFail) {
        statistics[easyPair]->PrimerRevFail++;
    }
    if (d->QualCtrl.minLqualTrim) {
        statistics[easyPair]->minLqualTrim++;
    }
    if (d->QualCtrl.TagFail) {
        statistics[easyPair]->TagFail++;
    }
    if (d->QualCtrl.fail_correct_BC) {
        statistics[easyPair]->fail_correct_BC++;
    }
    if (d->QualCtrl.suc_correct_BC) {
        statistics[easyPair]->suc_correct_BC++;
    }
    if (d->QualCtrl.RevPrimFound) {
        statistics[easyPair]->RevPrimFound++;
    }
    if (d->QualCtrl.QWinTrimmed || d->QualCtrl.AccErrTrimmed) {
        statistics[easyPair]->Trimmed++;
    }
    if (d->getTA_cut()) {
        statistics[easyPair]->adapterRem++;
    }

//exit(0);

    if (d->isGreenQual() || d->isYellowQual()) {
        countBCdetected(d->getBarcodeNumber(), easyPair, false);
        //and register as success
    } else {
        if (d->getBarcodeDetected()) {
            //DNA is no longer useful
            failedStats2(d, easyPair);
        }
        //delete d;
        if (d->QualCtrl.AvgQual) {
            statistics[easyPair]->AvgQual++;
        }
        if (d->QualCtrl.minL) {
            statistics[easyPair]->minL++;
        }
        if (d->QualCtrl.maxL) {
            statistics[easyPair]->maxL++;
        }
        if (d->QualCtrl.HomoNT) {
            statistics[easyPair]->HomoNT++;
        }
        if (d->QualCtrl.MaxAmb) {
            statistics[easyPair]->MaxAmb++;
        }
        if (d->QualCtrl.BinomialErr) {
            statistics[easyPair]->BinomialErr++;
        }
        if (d->QualCtrl.QualWin) {
            statistics[easyPair]->QualWin++;
        }
    }
    if (d->isDereplicated()) {
        if (d->getBarcodeDetected() && !d->isGreenQual() && !d->isYellowQual()) {
            this->statAddDerepBadSeq(d->getBarcodeNumber());
        }
    }
}


bool  OutputStreamer::saveForWrite(shared_ptr<DNA> d,int Pair, int thr) {
	//second most important part: collect stats on DNA passing through here (should be all read)
	//most important part: save DNA to be written later (or discard)
	if (d == NULL || stopAll) {
		return !stopAll;
	}
	int Cstream = 0;
	bool writen(false);
	
	//threadsafe
	Filters* curFil = this->getFilters(thr);

	curFil->addDNAtoCStats(d, Pair);

	if (suppressOutWrite == 3) {
		return !stopAll;
	}

	if (curFil->doFilterAtAll()) {
		if ((b_writeGreenQual && d->isGreenQual() ) || 
			(b_writeYellowQual && d->isYellowQual()) ) {
			d->prepareWrite(fastQoutVer);
			if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
				d->changeHeadPver(curFil->FQheadV());
			}
			writen = true;
			
			if (d->isYellowQual()) {
				Cstream = 1;
			}
		}
		else {
			return !stopAll;
		}
	}
	//sqfqostrMTX.lock();
	if (BWriteFastQ) {//write in fastq format
		*fqFile[Cstream][Pair-1] << d->writeFastQ();
	} else {//write in fasta (and maybe qual) format
	    *sFile[Cstream][Pair-1] << d->writeSeq(b_oneLinerFasta);
		if (BWriteQual) {
			*qFile[Cstream][Pair-1] << d->writeQual( b_oneLinerFasta);
		}
	}
	if (writen) { ReadsWritten++; }
	if (maxReadsPerOFile > 0 && ReadsWritten + DNAinMem >= maxReadsPerOFile) {
		//cerr << "ReadsWritten " << ReadsWritten << " DNAinMem " << DNAinMem << endl;
		DNAinMem = 0;
		//TODO multithread output file
		incrementOutputFile();
	}
	//sqfqostrMTX.unlock();


	return !stopAll;
/*
	if( d->isGreenQual()){//DNA is of good qual_ and should be written out //green
		if (b_writeGreenQual){ 

			//lock OutputStreamer
#ifdef _THREADED
			//Joachim: still needed? ofbusfstream should handle the mutex now..
			//std::lock_guard<std::mutex> guard(mutex);
#endif
			//dereplicate & create copy of DNA?

			//mem_used = true;
			if (Pair == 1){
				DNAsP1.push_back(d);
			}
			else if (Pair == 2){
				DNAsP2.push_back(d);
			}
			else if (Pair == 3){
				DNAsS1.push_back(d);
				curFil->collectStatistics[0]->singleton++;
			}
			else if (Pair == 4){
				DNAsS2.push_back(d);
				curFil->collectStatistics[1]->singleton++;
			}
			DNAinMem++;
		}
		
	} else if (d->isYellowQual()){//yelllow

		if (b_writeYellowQual){
			d->prepareWrite(fastQoutVer);
			mem_used = true;
			if (Pair == 1){//yellow P1
				DNAsP1_alt.push_back(d);
			}
			else if (Pair == 2){//yellow P2
				DNAsP2_alt.push_back(d);
			}
			else if (Pair == 3){ //yellow P1
				curFil->statAddition[0]->singleton++;
				DNAsS1_alt.push_back(d);
			}
			else if (Pair == 4){//yellow P2
				curFil->statAddition[1]->singleton++;
				DNAsS2_alt.push_back(d);
			}
			DNAinMem++;
		}

	} 
	//automatic mechanism to write to File, once enough DNA is in memory
	if (write2File && DNAinMem > DNA_MAX_IN_MEM){
        //TODO multithread writeAllStoredDNA
		writeAllStoredDNA();
		DNAinMem=0;
	}
	if (maxRdsOut > 0 && ReadsWritten + DNAinMem >= maxRdsOut) {
		writeAllStoredDNA();
		stopAll = true;
		
	}
	*/
}

bool OutputStreamer::saveForWrite_merge(shared_ptr<DNA> d, shared_ptr<DNA> d2,
		string newHeader,int curThread, bool elseWriteD1) {
	// CHECK WHERE TO IMPLEMENT BEST
// MERGE DNA1 AND DNA2
	shared_ptr<DNA> dna_merged = nullptr;
	if (d->merge_seed_pos_ > 0) {
		dna_merged = merger[curThread]->merge(d, d2);
	}
	if (dna_merged) {
		// write out merged DNA
		if (newHeader != "") {
			dna_merged->setNewID(newHeader);
		}
		dna_merged->prepareWrite(fastQoutVer);
		*(of_merged_fq[0]) << dna_merged->writeFastQ();
		return true;
	}
	if (!elseWriteD1) {//not passed? not my problem..
		return false;
	}
	//actually it is still my problem...
	d->prepareWrite(fastQoutVer);
	*(of_merged_fq[0]) << d->writeFastQ();

	return false;
}


/* Joachim: still needed?? Should be handled by ofbufstream?
bool  OutputStreamer::saveForWriteMT(shared_ptr<DNA> dna, int thread, int pair) {
    //second most important part: collect stats on DNA passing through here (should be all read)
    //most important part: save DNA to be written later (or discard)
    if (dna == NULL || stopAll) {
        return !stopAll;
    }

    MFil->addDNAtoCStatsMT(dna, pair, thread);

    if( dna->isGreenQual() ){//DNA is of good qual and should be written out

        if (b_writeGreenQual){
            dna->prepareWrite(fastQoutVer);
            //lock OutputStreamer
#ifdef _THREADED
            std::lock_guard<std::mutex> guard(mutex);
#endif
            //dereplicate & create copy of DNA?

            mem_used = true;
            if (pair == 1){
                //TODO multithread Vector DNAsP1
                DNAsP1.push_back(dna);
            }
            else if (pair == 2){
                //TODO multithread DNAsP2
                DNAsP2.push_back(dna);
            }
            else if (pair == 3){
                //TODO multithread DNAsP3
                DNAsS1.push_back(dna);
                MFil->statistics_[thread].passed_quality_stats_[0].singleton++;
            }
            else if (pair == 4){
                //TODO multithread DNAsP4
                DNAsS2.push_back(dna);
                MFil->statistics_[thread].passed_quality_stats_[1].singleton++;
            }
            //TODO multithread DNAinMem counter
            DNAinMem++; // Not thread safe
        }

    } else if (dna->isYellowQual()){

        if (b_writeYellowQual){
            dna->prepareWrite(fastQoutVer);
            mem_used = true;
            if (pair == 1){
                DNAsP1_alt.push_back(dna);
            }
            else if (pair == 2){
                DNAsP2_alt.push_back(dna);
            }
            else if (pair == 3){
                MFil->statistics_[thread].mid_quality_stats_[0].singleton++;
                DNAsS1_alt.push_back(dna);
            }
            else if (pair == 4){
                MFil->statistics_[thread].mid_quality_stats_[1].singleton++;
                DNAsS2_alt.push_back(dna);
            }
            DNAinMem++; // Not thread safe
        }

    }
    //automatic mechanism to write to File, once enough DNA is in memory
    if (write2File && DNAinMem > DNA_MAX_IN_MEM){
        writeAllStoredDNA();
        DNAinMem=0;
    }
    if (maxReadsPerOFile>0 && ReadsWritten+DNAinMem >= maxReadsPerOFile){
        //cerr << "ReadsWritten " << ReadsWritten << " DNAinMem " << DNAinMem << endl;
        DNAinMem=0;
        incrementOutputFile();
    }
    if (maxRdsOut > 0 && ReadsWritten + DNAinMem >= maxRdsOut) {
        writeAllStoredDNA();
        stopAll = true;

    }
    return !stopAll;
}
*/

/*no longer needed in new design using ofbufstream
void OutputStreamer::writeAndDel(shared_ptr<DNA> d,int Pair) {
	//ofstream tmpS, tmpQ, tmpFQ;
//	int PairC = Pair;
//	if (Pair > 1) {
//		PairC = Pair - 2;
//	}

    //DONE ersetze durch ofbufstream
	
	int ClsVec = -1;
	if (d != NULL) {
		if (MFil->doFilterAtAll()) {
			if (d->isGreenQual()) {
				ClsVec = 0;
			} else if (d->isYellowQual()) {
				ClsVec = 1;
			}
		}
		else {
			ClsVec = 0;
		}
	}
	if (ClsVec >= 0) {
		if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
			d->changeHeadPver(MFil->FQheadV());
		}
		if (BWriteFastQ) {
			if (fqFile[ClsVec][Pair ] == NULL ) {
				//cerr << "_";
				openOFstreamFQ(fqFileStr[ClsVec][Pair], wrMode, ClsVec, Pair , "Appending");
			}
			//cerr << "X";
			d->writeFastQ(*(fqFile[ClsVec][Pair]));
		} else {
			if (sFile[ClsVec][Pair]==NULL ) {
				openOFstreamFNA(sFileStr[ClsVec][Pair], wrMode, ClsVec, Pair, "Appending");
			}

			d->writeSeq(*sFile[ClsVec][Pair ], b_oneLinerFasta);
			if (BWriteQual) { 
				if (qFile[ClsVec][Pair]==NULL) {
					openOFstreamQL(qFileStr[ClsVec][Pair], wrMode, ClsVec, Pair, "Appending");
				}
				d->writeQual(*qFile[ClsVec][Pair], b_oneLinerFasta);
			}
		}
	}

//	delete d;

}
*/
void OutputStreamer::writeSelectiveStream(shared_ptr<DNA> d, int Pair,int FS) {
	//ofstream tmpS, tmpQ, tmpFQ;
	if (d == 0) {
		return;
	}
	d->prepareWrite(fastQoutVer);
	int PairC = Pair;
	if (Pair > 1) {
		PairC = Pair - 2;
	}
	if  (d->isGreenQual()) {
		MFil->DNAstatLQ(d, PairC, false);
	} else if ( d->isYellowQual()) {
		MFil->DNAstatLQ(d, PairC, true);
	}
	if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
		d->changeHeadPver(MFil->FQheadV());
	}
	if (BWriteFastQ) {
/*		if (fqFile[FS][Pair]==NULL||!*fqFile[FS][Pair]) {
			openOFstreamFQ(fqFileStr[FS][Pair], ios_base::app, FS, Pair, "Appending");
		}*/
		*(fqFile[FS][Pair]) << d->writeFastQ();
		delete fqFile[FS][Pair]; fqFile[FS][Pair] = NULL;
	} else {
		//if (sFile[FS][Pair]==NULL||!*sFile[FS][Pair]) {			openOFstreamFNA(sFileStr[FS][Pair], ios_base::app, FS, Pair, "Appending");		}
		*sFile[FS][Pair]<< d->writeSeq(b_oneLinerFasta);
		if (BWriteQual) { 
			//if (qFile[FS][Pair]==NULL||!*qFile[FS][Pair]) {				openOFstreamQL(sFileStr[FS][Pair], ios_base::app, FS, Pair, "Appending");			}
			*qFile[FS][Pair]<<d->writeQual( b_oneLinerFasta);
		}
	}


//	delete d;

}
void OutputStreamer::writeNonBCReads(shared_ptr<DNA> d, shared_ptr<DNA> d2) {
	if (fqNoBCFile.size() == 2) {
		if (!d->getBarcodeDetected() && !d2->getBarcodeDetected()) {
			if ((!d->getBarcodeDetected() && d2->getBarcodeDetected()) ||
				(d->getBarcodeDetected() && !d2->getBarcodeDetected())) {
				cerr << "Barcode only set in 1 reads.. something wrong!\n";
			}
			//nobcostrmMTX.lock();
			*(fqNoBCFile[0]) << d->writeFastQ();
			*(fqNoBCFile[1]) << d2->writeFastQ();
			//nobcostrmMTX.unlock();
		}
	}
}

void OutputStreamer::setSubfilters(int num) {
	if (num<1){return;}
	subFilter.resize(num,NULL);
	for (uint i=0;i<subFilter.size();i++){
		subFilter[i] = new Filters (MFil,MFil->currentBCnumber(),true);
    }
}



void OutputStreamer::mergeSubFilters() {
	vector<int> idx (MFil->Barcode.size(),0);
	for (uint i=0;i<idx.size();i++){idx[i]=i;}
	for (uint i=0; i<subFilter.size();i++){
		MFil->addStats(subFilter[i],idx);
	}
}
/*
void OutputStreamer::mergeSubFiltersMT() {
    vector<int> idx (MFil->Barcode.size(),0);

    for (uint i=0;i<idx.size();i++){idx[i]=i;}
    for (uint i=0; i<subFilter.size();i++){
//        MFil->addStats(subFilter[i],idx);
        StatisticsMultithreaded::mergeStatisticsInto(MFil->collectStatistics, 
			MFil->statAddition, subFilter[i]->statistics_, idx);
    }
}
*/

void OutputStreamer::attachDereplicator(shared_ptr<Dereplicate> de) {
	if (de != NULL) {
        dereplicator = de; b_doDereplicate = true;
		dereplicator->setPaired(pairedSeq > 1);
		//insert code here to fix BC offset in filter & add to derep info on sample names from the current filter
		int curBCOffset = dereplicator->getHighestBCoffset();

		cdbg(  "BARCODE INCREMENT: " + itos(curBCOffset )+"\n");
		
		MFil->setBCoffset(curBCOffset);
		for (size_t i = 0; i < subFilter.size(); i++) {
			subFilter[i]->setBCoffset(curBCOffset);
		}
		dereplicator->BCnamesAdding(MFil);

		if (dereplicator->mergeDereRead()) {
			b_merge_pairs_ = true;
		}

	}
}

/*void OutputStreamer::dereplicateDNA(shared_ptr<DNA> tdn) {

	if (b_doDereplicate && tdn->isGreenQual()) {
		cntDerep++;
		dereplicator->addDNA(tdn);
	}
}*/
void OutputStreamer::dereplicateDNA(shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2) {
	if (false || !b_doDereplicate) {
		return;
	}
	dereplicator->addDNA(tdn, tdn2);	
	if (tdn->isDereplicated()) {
		//drpMTX.lock();
		cntDerep++;
		//drpMTX.unlock();
	}
}

void OutputStreamer::resetOutFilesAndFilter(){
#ifdef DEBUG
	cerr << "resetOutFilesAndFilter" << endl;
#endif
	//reset count stats
	MFil->resetStats();
	//reset all stored DNA
	delAllDNAvectors();
	dereplicator.reset();
	//close streams -- why? no reason, since nothing has been writen
	//this->closeOutStreams();
	totalFileStrms = 0;
}

void OutputStreamer::closeOutStreams(bool wr){
	write2File = true;
	for (size_t i = 0; i < sFile.size(); i++) {
		for (size_t j = 0; j < sFile[i].size(); j++) {
			if (sFile[i][j] != nullptr) {
				if (!wr) { sFile[i][j]->emptyStream();}
				delete sFile[i][j];
				sFile[i][j] = nullptr;
			}
		}
	}
	for (size_t i = 0; i < qFile.size(); i++) {
		for (size_t j = 0; j < qFile[i].size(); j++) {
			if (qFile[i][j] != nullptr) {
				if (!wr) { qFile[i][j]->emptyStream(); }
				delete qFile[i][j];
				qFile[i][j] = nullptr;
			}
		}
	}
	for (size_t i = 0; i < fqFile.size(); i++) {
		for (size_t j = 0; j < fqFile[i].size(); j++) {
			if (fqFile[i][j] != nullptr) {
				if (!wr) { fqFile[i][j]->emptyStream(); }
				delete fqFile[i][j];
				fqFile[i][j] = nullptr;
			}
		}
	}
	for (size_t i = 0; i < of_merged_fq.size(); i++) {
		if (of_merged_fq[i] != nullptr) {
			if (!wr) { of_merged_fq[i]->emptyStream(); }
			delete of_merged_fq[i];
			of_merged_fq[i] = nullptr;
		}
	}
	/*  completely useless, replaced with ofbufstream
	if (wr){
		this->writeAllStoredDNA();
	}
	
#ifdef DEBUG
	cerr << "closing output streams";
#endif
	for (int i = 0; i < (int)fqFile.size(); i++) {
		for (int j = 0; j < (int)fqFile[i].size(); j++) {
			if (fqFileStr[i][j] != "T") {
				delete fqFile[i][j]; fqFile[i][j] = NULL;
			}
		}
	}

	for (int i = 0; i < (int)sFile.size(); i++) {
		for (int j = 0; j < (int) sFile[i].size(); j++) {
			if (sFileStr[i][j] != "T") {
				delete sFile[i][j]; sFile[i][j] = NULL;
			}
		}
	}
	for (int i = 0; i < (int)qFile.size(); i++) {
		for (int j = 0; j < (int)qFile[i].size(); j++) {
			if (qFileStr[i][j] != "T") {
				delete qFile[i][j]; qFile[i][j] = NULL;
			}
		}
	}
	//if(qFile){qFile.close();}if(sFile){sFile.close();}if(fqFile){fqFile.close(); }
	//other housekeeping tasks
	if (pool) {
	    // If pool exists it means the program runs in multithreaded mode
	    this->mergeSubFiltersMT();
	} else {
        this->mergeSubFilters();
    }
	*/
	MFil->setWrittenReads(ReadsWritten);
#ifdef DEBUG
	cerr << ".. closed\n";
#endif

}
/*void OutputStreamer::resetOutStreams(){
	if(qFile){qFile.seekp(qFilePos);}if(sFile){sFile.seekp(sFilePos);}if(fqFile){fqFile.seekp(fqFilePos); }
	if(qFile2){qFile2.seekp(qFile2Pos);}if(sFile2){sFile2.seekp(sFile2Pos);}if(fqFile2){fqFile2.seekp(fqFile2Pos); }
	if(qFileS){qFileS.seekp(qFileSPos);}if(sFileS){sFileS.seekp(sFileSPos);}if(fqFileS){fqFileS.seekp(fqFileSPos); }
	if(qFileS2){qFileS2.seekp(qFileS2Pos);}if(sFileS2){sFileS2.seekp(sFileS2Pos);}if(fqFileS2){fqFileS2.seekp(fqFileS2Pos); }
}*/
void OutputStreamer::openOFstream(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep, int wh) {
	switch (wh) {
	case 1:
		openOFstreamFNA(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	case 0:
		openOFstreamFQ(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	case 2:
		openOFstreamQL(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	default:
		cerr << "Wrong wh specified"; exit(1002);
	}
}
void OutputStreamer::openNoBCoutstrean(const string inS) {
	vector<string> tfnaout = splitByCommas(inS);
	bool doMC = Nthrds > 1;

	fqNoBCFile.resize(2, NULL);
	fqNoBCFile[0] = new ostr(tfnaout[0], wrMode,doMC);
	fqNoBCFile[1] = new ostr(tfnaout[1], wrMode,doMC);
}

void OutputStreamer::openOFstreamFQ(const string opOF, std::ios_base::openmode wrMode, 
	int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)fqFile.size()) {
		vector<ostr*> nullVec(4, NULL);
		fqFile.resize(p1+1, nullVec);
	}
	bool doMC = Nthrds > 1;

	//if ((int)fqFileStr.size() - 1 <= p1) {		fqFileStr.resize(p1 + 1, vector<string>(4, ""));	}
	//fqFileStr[p1][p2] = opOF;
	//if (onlyPrep) { return; }
	//if (p1 == 1 && !b_writeYellowQual ){ return; }//p1==1: mid passed suppressOutWrite >= 2
	//if (p1 == 0 && !b_writeGreenQual){ return; }//suppressOutWrite == 1
	fqFile[p1][p2] = new ostr(opOF, wrMode,doMC);
	if (!onlyPrep) {
		fqFile[p1][p2]->activate();
	}
	if (!*fqFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " fastq output file " << opOF << endl << p1 << " " << p2 << " " << totalFileStrms << endl;
		exit(4);
	}
}
void OutputStreamer::openOFstreamFQ_mrg(const string opOF, std::ios_base::openmode wrMode, 
	int p1, string errMsg, bool onlyPrep) {
	//also set up potential stream for merge out
	if (p1 + 1 >= (int)of_merged_fq.size()) {
		of_merged_fq.resize(p1 + 1, nullptr);
	}
	bool doMC = Nthrds > 1;

	of_merged_fq[p1] = new ostr(opOF, wrMode,doMC);
	if (!*of_merged_fq[p1]) {
		cerr << "Could not open " << errMsg << " fastq merged output file " << opOF << endl << p1 << " " << totalFileStrms << endl;
		exit(4);
	}
}
void OutputStreamer::openOFstreamFNA(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)sFile.size()) {
		vector<ostr*> nullVec(4, NULL);
		sFile.resize(p1+1, nullVec);
	}
	/*if ((int)sFileStr.size() - 1 <= p1) {
		sFileStr.resize(p1 + 1, vector<string>(4, ""));
	}
	sFileStr[p1][p2] = opOF;*/
	bool doMC = Nthrds > 1;

	sFile[p1][p2] = new ostr(opOF, wrMode,doMC);
	if (!*sFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " fasta output file " << opOF << endl << p1 << " " << p2 << " " << totalFileStrms << endl;
		exit(4);
	}
//	if (onlyPrep ) { return; }
	//if (p1 == 1 && !b_writeYellowQual){ return; }//p1==1: mid passed suppressOutWrite >= 2
	//if (p1 == 0 && !b_writeGreenQual){ return; }//suppressOutWrite == 1

}
void OutputStreamer::openOFstreamQL(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)qFile.size()) {
		vector<ostr*> nullVec(4, NULL);
		qFile.resize(p1+1, nullVec);
	}
	bool doMC = Nthrds > 1;

	//if (onlyPrep) { return; }
	//if (p1 == 1 && !b_writeYellowQual){ return; }//p1==1: mid passed suppressOutWrite >= 2
	//if (p1 == 0 && !b_writeGreenQual){ return; }//suppressOutWrite == 1
	qFile[p1][p2] = new ostr(opOF, wrMode,doMC);
	if (!*qFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " quality output file " << opOF << endl;
		exit(4);
	}
	/*if (opOF == "T") {
		qFile[p1][p2] = &std::cout;
	}
	else if (isGZfile(opOF)) {
#ifdef _gzipread
//		qFile[p1][p2] = new ogzstream(opOF.c_str(), wrMode);
        // experimental
        qFile[p1][p2] = new zstr::ofstream(opOF.c_str(), wrMode);
#else
		cerr << "gzip outpout not supported in your sdm build\n" << opOF; exit(50);
#endif
	} else {
		qFile[p1][p2] = new ofstr(opOF, wrMode);
	}
	if (!*qFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " quality output file " << opOF << endl;
		exit(4);
	}*/
}

void OutputStreamer::openSeveralOutstreams(OptContainer& cmdArgs, shared_ptr<ReadSubset> RDS, std::ios_base::openmode wrMode) {
#ifdef DEBUG
	cerr << " openining multiple out streams" << endl;
#endif
	string path = "", fileEnd(".fna");
	vector<string> outFile = RDS->getOFiles();
	bool openStrms = true; int omode(1);
	vector<vector<ostr*>>& tmp = sFile;
	if (cmdArgs.find("-o_fastq") != cmdArgs.end() && cmdArgs["-o_fastq"] != "" && cmdArgs["-o_fna"] == "") { //write fastq
		BWriteFastQ = true;
		path = cmdArgs["-o_fastq"];
		tmp = fqFile; omode = 0;
		fileEnd = ".fastq";
	} else  if (cmdArgs["-o_fna"] != "") {
		BWriteFastQ = false;
		path = cmdArgs["-o_fna"];
	} else {
		cerr << "Could not find valid sequence output file path. Exiting\n";
		exit(55);
	}
	if (RDS->multiFile()) {
		b_multiOutStream = true;
	}
	int i(0);
	for (i = 0; i < (int)outFile.size(); i++) {
		if (totalFileStrms >= MAX_FILE_STREAMS && openStrms) {
			cerr << "Too many output file streams\nSwitching to dynamical file appending\n";
			openStrms = false;
			//this->closeOutStreams();
		}
		string baseFile = path + removeFileEnding(outFile[i]);
		if (pairedSeq > 1) {
			if (i < 4) {
				cerr << "Outfile " << i << ": "<<baseFile + ".1" << fileEnd << "\n";
			}
			openOFstream(baseFile + ".1" + fileEnd, wrMode, i, 0, "paired 1st " + itos(i), !openStrms, omode);
			//2nd pair_
			openOFstream(baseFile + ".2" + fileEnd, wrMode, i, 1, "paired 2nd " + itos(i), !openStrms, omode);
			//1st singleton
			//if (openStrms) { openOFstream(baseFile + ".1" + fileEnd + SingletonFileDescr, wrMode, i, 2, "Singleton 1 " + itos(i), omode); }
			//2nd singleton
			//if (openStrms) { openOFstream(baseFile + ".2" + fileEnd + SingletonFileDescr, wrMode, i, 3, "Singleton 2 " + itos(i), omode); }
			totalFileStrms += 2;
		} else {
			openOFstream(baseFile + fileEnd, wrMode, i, 0, "main " + itos(i), !openStrms, omode);
			totalFileStrms++;
		}

	}
	//pipe remaining reads into this file
	if (cmdArgs["-excludeFile"] != "") {
		//set a flag
		RDS->setRemainingFilepipe(i);
		if (pairedSeq == 1) {
			openOFstream(cmdArgs["-excludeFile"], wrMode, i, 0, "excludeFile file ", !openStrms, omode);
		} else {
			string baseFile = path + removeFileEnding(cmdArgs["-excludeFile"]);
			if (i < 4) {
				cerr << "out file " << i << baseFile + ".1" << fileEnd << "\n";
			}
			openOFstream(baseFile + ".1" + fileEnd, wrMode, i, 0, "paired 1st " + itos(i), !openStrms, omode);
			//2nd pair_
			openOFstream(baseFile + ".2" + fileEnd, wrMode, i, 1, "paired 2nd " + itos(i), !openStrms, omode);
			//1st singleton
			//if (openStrms) { openOFstream(baseFile + ".1" + fileEnd + SingletonFileDescr, wrMode, i, 2, "Singleton 1 " + itos(i), omode); }
			//2nd singleton
			//if (openStrms) { openOFstream(baseFile + ".2" + fileEnd + SingletonFileDescr, wrMode, i, 3, "Singleton 2 " + itos(i), omode); }
			totalFileStrms += 2;
		}
	}
}

void OutputStreamer::openOutStreams(OptContainer& cmdArgs,int fileIt,
	std::ios_base::openmode wrMode_i, 
	string fileExt, int forceFmt){
	this->setwriteMode(wrMode_i);
	if ( suppressOutWrite == 3 || (cmdArgs["-o_fastq"] == "" && cmdArgs["-o_fna"] == "" && cmdArgs["-o_qual"] == "") ){
		suppressOutWrite = 3; b_writeGreenQual = false; 
		b_writeYellowQual = false; return;
	}
	if (forceFmt != -1){
		if (forceFmt == 1 && (cmdArgs.find("-o_fna") == cmdArgs.end() || cmdArgs["-o_fna"] == "") ){//force fna, no qual, required for seed ref fastas
			cmdArgs["-o_fna"] = cmdArgs["-o_fastq"];
			cmdArgs["-o_fastq"] = "";
		}
	}
	

	if (cmdArgs.find("-o_fastq")  != cmdArgs.end() && cmdArgs["-o_fastq"] != "" && cmdArgs["-o_fna"] == ""){ //write fastq
		this->setFastQWrite(true);
		if (pairedSeq>1){ //open second pair_ + singleton
#ifdef DEBUG
			cerr << " paired fastq out " << endl;
#endif
			vector<string> tfnaout = splitByCommas(cmdArgs["-o_fastq"]);
			if (tfnaout.size()!=2){
				cerr<<"Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = "<<cmdArgs["-o_fastq"]<<endl; exit(57);
			}
			//1st pair (or main) file out
			leadingOutf = applyFileIT(tfnaout[0] + fileExt, fileIt);
			openOFstreamFQ(leadingOutf.c_str(), wrMode, 0, 0, "paired 1st");
			//2nd pair
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			//1st singleton
			openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 2, "Singleton 1", true);
			//2nd singleton
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 3, "Singleton 2", true);
			
			//merged reads
			int prefL = (int)len_common_prefix_base(tfnaout[0].c_str(), tfnaout[1].c_str());
			string prefix = tfnaout[0].substr(0,prefL);
			openOFstreamFQ_mrg(applyFileIT(prefix + "merg.fq" + fileExt, fileIt).c_str(), wrMode, 0, "merged");
			
			//additional file
			if (cmdArgs.find("-o_fastq2")  != cmdArgs.end() && cmdArgs["-o_fastq2"].length()>1){ //write fastq
				tfnaout = splitByCommas(cmdArgs["-o_fastq2"]);
				openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st");
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd");
				openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt , fileIt, SingletonFileDescr).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt , fileIt, SingletonFileDescr).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
			} else { b_writeYellowQual = false; }
		} else {
#ifdef DEBUG
			cerr << " single fastq out " << endl;
#endif
			openOFstreamFQ(applyFileIT(cmdArgs["-o_fastq"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "the main");
			leadingOutf = applyFileIT(cmdArgs["-o_fastq"] + fileExt, fileIt);
			//additional file
			if (cmdArgs.find("-o_fastq2") != cmdArgs.end() && cmdArgs["-o_fastq2"].length()>1) { //write fastq
				openOFstreamFQ(applyFileIT(cmdArgs["-o_fastq2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "the additional");
			}
		}

		return;
	}
	BWriteFastQ=false;
	if (pairedSeq==1){
#ifdef DEBUG
		cerr << " fasta singleton output " << endl;
#endif
		openOFstreamFNA(applyFileIT(cmdArgs["-o_fna"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
		leadingOutf = applyFileIT(cmdArgs["-o_fna"] + fileExt, fileIt);
		if (cmdArgs["-o_qual"] != ""){
			openOFstreamQL(applyFileIT(cmdArgs["-o_qual"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
			this->setQualWrite(true);
		} else {
			this->setQualWrite(false);
		}
		//additional file (secondary filter)
		if (cmdArgs.find("-o_fna2") != cmdArgs.end() && cmdArgs["-o_fna2"].length()>1) { //add file
			openOFstreamFNA(applyFileIT(cmdArgs["-o_fna2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
			if (cmdArgs.find("-o_qual2") != cmdArgs.end() && cmdArgs["-o_qual2"].length() > 1) { //add file
				openOFstreamQL(applyFileIT(cmdArgs["-o_qual2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
				this->setQualWrite(true);
			}
		} else { b_writeYellowQual = false; }

	} else {
#ifdef DEBUG
		cerr << " fasta paired output " << endl;
#endif

		vector<string> tfnaout = splitByCommas(cmdArgs["-o_fna"]);
		if (tfnaout.size()!=2){
			cerr<<"Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = "<<cmdArgs["-o_fna"]<<endl; exit(57);
		}
		
		openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 0, 0, "paired 1st");
		leadingOutf = applyFileIT(tfnaout[0] + fileExt, fileIt);

		openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
		openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 2, "Singleton 1", true);
		openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 3, "Singleton 2", true);
		/*setFilePos(sFile,sFilePos);
		setFilePos(sFile2,sFile2Pos);
		setFilePos(sFileS,sFileSPos);
		setFilePos(sFileS2,sFileS2Pos);*/
		if (cmdArgs.find("-o_fna2")  != cmdArgs.end() && cmdArgs["-o_fna2"].length()>1){ //write fastq
			tfnaout = splitByCommas(cmdArgs["-o_fna2"]);

			openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st", true);
			openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd", true);
			openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
			openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
		}
		if (cmdArgs["-o_qual"] != ""){
			this->setQualWrite(true);
			vector<string> tqout = splitByComma(cmdArgs["-o_qual"],true);
			openOFstreamQL(applyFileIT(tqout[0] + fileExt, fileIt).c_str(), wrMode, 0, 0, "paired 1st");
			openOFstreamQL(applyFileIT(tqout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			openOFstreamQL(applyFileIT(tqout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 2, "Singleton 1", true);
			openOFstreamQL(applyFileIT(tqout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 3, "Singleton 2", true);
			if (cmdArgs["-o_qual2"] != "") {
				vector<string> tqout = splitByComma(cmdArgs["-o_qual2"], true);
				openOFstreamQL(applyFileIT(tqout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st", true);
				openOFstreamQL(applyFileIT(tqout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd", true);
				openOFstreamQL(applyFileIT(tqout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
				openOFstreamQL(applyFileIT(tqout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
			} 
		} else {
			this->setQualWrite(false);
		}

	}
}

/*
shared_ptr<DNA> OutputStreamer::mergeDNA(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2, ReadMerger &merger) {
    //bool mergeable;

//    if (!mergeable) {
//        return dna1;
//    } else {
//        shared_ptr<DNA>
//    }
	return NULL;
}
*/

void OutputStreamer::findSeedForMerge(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2, int thrPos) {
	bool didMerge(false);
	if (dna1->length() == 0 || dna2->length() == 0) {
		total_read_preMerge_++;  return; 
	}
    if (merger[thrPos]->findSeed(dna1->getSequence(), dna2->getSequence())) {
		didMerge = true;
        dna1->merge_seed_pos_ = (int) merger[thrPos]->result.seed.pos1;
        dna1->merge_offset_ = merger[thrPos]->result.offset1;
        dna2->merge_seed_pos_ = (int) merger[thrPos]->result.seed.pos2;
        dna2->merge_offset_ = merger[thrPos]->result.offset2;
        dna2->reversed_merge_ = merger[thrPos]->result.seed.is2reversed;
    }
	total_read_preMerge_++; 
	if (didMerge) {
		merged_counter_++;
	}
	//mergStatMTX.lock(); 	mergStatMTX.unlock(); 
	return;
}

//*******************************************
//*        DNAuniqSet OBJECT
//*******************************************


void DNAuniqSet::setBest() {
	int bestCnt = 0;
	int bestPos = -1;
	//shared_ptr<DNAunique> lastBest;
	for (auto dd : DNUs) {
		if (dd.second->totalSum() > bestCnt) {
			bestCnt = dd.second->totalSum();
			bestPos = dd.first;
			bestDNU = dd.second;
		}
	}
	//completely unbiased selection of whatever has the highest counts.. could select non-merge before merge
	if (bestPos != -1 && DNUs.size() > 1) {
		//include + - 1?
		auto xx = DNUs.find((bestPos - 1));
		if (xx != DNUs.end()) {bestDNU->transferOccurence(xx->second);}
		xx = DNUs.find((bestPos + 1));
		if (xx != DNUs.end()) {bestDNU->transferOccurence(xx->second);}
	}
	bestSet = true;
}

//*******************************************
//*        DEREPLICATE OBJECT
//*******************************************
Dereplicate::Dereplicate(OptContainer& cmdArgs, Filters* mf):
        barcode_number_to_sample_id_(0), b_usearch_fmt(true), b_singleLine(true), b_pairedInput(false),
        minCopies(1,0), minCopiesStr("0"), //default minCopies accepts every derep
		totSize(0), tmpCnt(0), curBCoffset(0), b_derep_as_fasta_(true), b_derepPerSR(false),
		b_wroteMapHD(false), b_merge_pairs_derep_(false),merger(nullptr),
		mapF(""), outHQf(""), outHQf_p2(""), outRest(""),
		mainFilter(mf)
{
	outfile = cmdArgs["-o_dereplicate"];

	string baseOF = outfile.substr(0, outfile.find_last_of('.'));
	mapF = baseOF + ".map";
	outHQf = baseOF + ".1.hq.fq";
	outHQf_p2 = baseOF + ".2.hq.fq";
	outRest = outfile + ".rest";
	//clean up eventually existing files
	remove(mapF.c_str());	remove(outHQf.c_str());
	remove(outHQf_p2.c_str());	remove(outRest.c_str());


	if (cmdArgs.find("-dere_size_fmt") != cmdArgs.end() && cmdArgs["-dere_size_fmt"] == "1") {
		b_usearch_fmt = false;
	}
	if (cmdArgs["-derepPerSR"] == "1") {
		b_derepPerSR = true;
	}

	if (cmdArgs.find("-min_derep_copies") != cmdArgs.end()) {
		minCopies[0] = -1;//in this case reset to -1 the first entry..
		minCopiesStr = cmdArgs["-min_derep_copies"];
		string x = cmdArgs["-min_derep_copies"];

		vector<string> xs = splitByCommas(x, ',');
		for (size_t i = 0; i < xs.size(); i++) {
			vector<string> tmp = splitByComma(xs[i], false, ':');
			if (tmp.size() > 2) {
				cerr << "Derep string " << xs[i] << " has to be two integers seqparated by \":\"\n"; exit(623);
			}
			int pos = 1; int cnt = -1;
			if (tmp.size() == 1) {//this is the 1 sample position, if not set
				cnt = atoi(tmp[0].c_str());
			}else{
				pos = atoi(tmp[1].c_str());//number of samples
				cnt = atoi(tmp[0].c_str());//16s copy numbers
			}
			if (pos <= 0) { cerr << "wrong derplicate position (<=0):" << xs[i] << endl; exit(312); }
			if (pos > 3000) { cerr << "too large derplicate position (>3000):" << xs[i] << endl; exit(313); }
			if (minCopies.size() < (size_t)pos) {
				minCopies.resize(pos, -1);
			}
			minCopies[pos-1] = cnt;
		}
		minCopiesSiz = minCopies.size();
		//minCopies = atoi(cmdArgs["-min_derep_copies"].c_str());

	}
	// set up output format of dereplication:
	if (cmdArgs.find("-derep_format") != cmdArgs.end()) {
		auto res = cmdArgs["-derep_format"];
		if (strcmp(res.c_str(), "fa") == 0) {
			b_derep_as_fasta_ = true;
		}
		else if (strcmp(res.c_str(), "fq") == 0) {
			b_derep_as_fasta_ = false;
		}
		else {
			std::cerr << "Invalid input for option -derep_format. If specified, -derep_format can take either \"fa\" for derep_as_fasta_ and \"fq\" for fastq. ";
			exit(77);
		}
	}

	// Flag for merging paired end fsatq reads
    if (cmdArgs.find("-merge_pairs_derep") != cmdArgs.end() && 
		cmdArgs["-merge_pairs_derep"] == "1") {
		b_merge_pairs_derep_ = true;
    }
}
/*bool Dereplicate::addDNA(shared_ptr<DNA> d) {
	//1st build hash of DNA
	string seq = d->getSeqPseudo();
	int BCN = d->getBarcodeNumber();
	if (BCN >= 0) { tmpCnt++; }
	HashDNAIT spotted = Tracker.find(seq);
	if (spotted != Tracker.end()) {// found something
		int hpos ( spotted->second );
		Dnas[hpos]->incrementSampleCounter(BCN);
		if (betterPreSeed(d, NULL, Dnas[hpos])) {
			//takeOver old DNA
			DNAunique *du = new DNAunique(d, -1);
			du->saveMem();
			du->setBestSeedLength(Dnas[hpos]->getBestSeedLength());
			du->transferOccurence(Dnas[hpos]);
			delete Dnas[hpos];
			Dnas[hpos] = du;
		}
		return false;
	} else {
		//create entry
		Tracker[seq] = (int)Dnas.size();
		DNAunique *tmp = new DNAunique(d, BCN);
		tmp->saveMem();
		Dnas.push_back(tmp);
		return true;
	}
	return true;
}*/

bool Dereplicate::addDNA(shared_ptr<DNA> dna, shared_ptr<DNA> dna2) {
	//1st build hash of DNA
	if (!dna->getBarcodeDetected()) {
		return false;
	}
	// Get copy of sequence (might have already been modified)
	//string seq = dna->getSeqPseudo();

	int sample_id = dna->getBarcodeNumber();
	bool pass = dna->isGreenQual();
	//deactivate this for now..
	int MrgPos1 = dna->merge_seed_pos_;
	//int MrgPos1 = -1;
	shared_ptr<DNA> dna_merged = nullptr;
	string srchSeq("");

	if (dna2 != nullptr && dna->merge_seed_pos_ != -1) {
		MrgPos1 = dna2->merge_seed_pos_;
		dna_merged = merger->merge(dna, dna2);
	}

	if (dna_merged){
		srchSeq = dna_merged->getSeqPseudo().substr(0, dna->length());
	}
	else {
		srchSeq = dna->getSeqPseudo();
	}

    // Lock because were accessing the base_map

    // See if DNA object is present
    // auto dna_unique = Tracker.find(new_dna_unique);
	bool new_insert(true);
	map<int, shared_ptr<DNAunique>>::iterator dna_unique;
    drpMTX.lock_shared();//lock for hash
	HashDNA::iterator dna_unique1 = Tracker.find(srchSeq);
	//HashDNA::iterator TrackerEnd = Tracker.end();
	if (dna_unique1 != Tracker.end()) {// found something
		new_insert = false;
		dna_unique1->second.lockMTX.lock();
		dna_unique = dna_unique1->second.find(MrgPos1);
		//truly dereplicated? at least overlap should fit..
		if (dna_unique == dna_unique1->second.end()) {
			//new_insert = true;
			dna_unique1->second.addNewDNAuniq(dna, dna2, MrgPos1, sample_id);
		} else { // compare to existing DNA
			dna_unique->second->matchedDNA(dna, dna2, sample_id, b_derep_as_fasta_);
		} 
		dna_unique1->second.lockMTX.unlock();//just to be on safe side, lock entire section
	}
	drpMTX.unlock_shared();//lock for hash

	if (new_insert && pass) {
		drpMTX.lock(); 
		Tracker[srchSeq].addNewDNAuniq(dna, dna2, MrgPos1, sample_id);
		drpMTX.unlock();
        // Create new dna_unique object
		//cdbg("set new DNAderep ");
        /*
		dna->setMidQual(false); dna->setDereplicated();
		shared_ptr<DNAunique> new_dna_unique = make_shared<DNAunique>(dna, sample_id);
        new_dna_unique->saveMem();
        if (dna2 != nullptr) {new_dna_unique->attachPair(make_shared<DNAunique>(dna2, sample_id));} 
		drpMTX.lock();
		Tracker[srchSeq][MrgPos1] = new_dna_unique;
		drpMTX.unlock();
		*/
    }
	//drpMTX.unlock();
    return new_insert;//didn't do a thing..
}

void Dereplicate::BCnamesAdding(Filters* fil) {
	vector<string> refSID = fil->SampleID;
	for (size_t i = 0; i < refSID.size(); i++) {
		barcode_number_to_sample_id_.push_back(refSID[i]);
	}
	//this will in the end be used to synchronize the different barcodes between samples
	curBCoffset = (int)barcode_number_to_sample_id_.size();
	//cerr << curBCoffset << " Added BC list\n";
}
void Dereplicate::reset() {
//	for (size_t i = 0; i < Dnas.size(); i++) { delete Dnas[i]; }
//	Dnas.resize(0);
	totSize = 0; tmpCnt = 0;
	//barcode_number_to_sample_id_.resize(0);
	Tracker.clear();
}
bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) { 
	return l->totalSum() > r->totalSum(); 
}

void Dereplicate::finishMap() {
	//at this point we can onl be sure that barcode_number_to_sample_id_ is finished
	//hence now it the point to add this to map
	totSize = 0; tmpCnt = 0;
	Tracker.clear();
	std::ifstream inputFile(mapF.c_str(),ios::in);
	std::ofstream outputFile((mapF+"t").c_str(), ios::out);

	outputFile << "#SMPLS";
	bool bCombiSmpl = mainFilter->combineSamples();
	vector<int> smplId2comb(0, 0);
	unordered_map<string, int> & combiMapCollectGrp = mainFilter->combiMapCollectGrp;
	if (!bCombiSmpl) {
		for (size_t i = 0; i < barcode_number_to_sample_id_.size(); i++) {
			outputFile << "\t" + itos((int)i) + ":" + barcode_number_to_sample_id_[i];
		}
	}
	else {
		smplId2comb = mainFilter->combiSmplConvergeVec(barcode_number_to_sample_id_);
		if (barcode_number_to_sample_id_.size() != smplId2comb.size()) {
			cerr << "FATAL: barcode_number_to_sample_id_ != smplId2comb\n"; exit(234);
		}
		for (auto IT = combiMapCollectGrp.begin(); IT != combiMapCollectGrp.end(); IT++) {
			outputFile << "\t" << itos(IT->second) + ":" + IT->first;
		}
	}
	outputFile << "\n";
	b_wroteMapHD = true;

	outputFile << inputFile.rdbuf();

	inputFile.close();
	outputFile.close();

	std::remove(mapF.c_str());
	std::rename((mapF + "t").c_str(), mapF.c_str());
}
string Dereplicate::writeDereplDNA(Filters* mf, string SRblock) {
	ofstream of, omaps, of2, ofRest, of2p2, of_merged;
	cerr << "Evaluating and writing dereplicated reads..\n";
	int fastqVer = mf->getuserReqFastqOutVer();
	
	//set the correct file names for dereplication output files
	string baseOF = outfile.substr(0, outfile.find_last_of('.'));
	string baseOF2 = baseOF;	
	if (b_derepPerSR && SRblock != "") {//if writing out per SRblock, needs different filename for map and 
		baseOF2 += "." + SRblock;
	}	
	string fileSuff = outfile.substr(outfile.find_last_of('.'));
	string outfile2 = baseOF2 + fileSuff;
    string outfile_merged = baseOF2 + ".merg" + fileSuff;

	// OPEN OUTPUT FILE STREAMS
	omaps.open(mapF.c_str(), ios::app);//| ios::binary
	//actual dereplicated reads as required by clustering algo (uparse, cdhit, dada2,...)
	of.open(outfile2.c_str(), ios::out);
	//reads that did not pass dereplication min number
	ofRest.open(outRest.c_str(), ios::app);
	//copy of high quality read pairs -> needed for Seed extension step
	of2.open(outHQf.c_str(), ios::app);
    if (b_merge_pairs_derep_) {
        of_merged.open(outfile_merged.c_str(), ios::out);
    }
	if (b_pairedInput) {
		of2p2.open(outHQf_p2, ios::app);
	}



	//sample specific derep filter
	bool dereplicate_sample_specific(false);
	vector<int> dereplicate_sample_specific_indices = mf->getDrerepSampleSpecifity();
	if (!dereplicate_sample_specific_indices.empty()) { dereplicate_sample_specific = true; }

	//combiner relevant vars
	bool bCombiSmpl = mf->combineSamples();
	vector<int> smplId2comb(0,0);
	unordered_map<string, int> & combiMapCollectGrp = mf->combiMapCollectGrp;
	if (bCombiSmpl) {
		smplId2comb = mf->combiSmplConvergeVec(barcode_number_to_sample_id_);
	}
	//vector<string> refSID = mf->SampleID;
	if (false && !b_wroteMapHD) {
		omaps << "#SMPLS";
		if (!bCombiSmpl) {
			for (size_t i = 0; i < barcode_number_to_sample_id_.size(); i++) {
				omaps << "\t" + itos((int)i) + ":" + barcode_number_to_sample_id_[i];
			}
		}
		else {
			if (barcode_number_to_sample_id_.size() != smplId2comb.size()) {
				cerr << "FATAL: barcode_number_to_sample_id_ != smplId2comb\n"; exit(234);
			}
			for (auto IT = combiMapCollectGrp.begin(); IT != combiMapCollectGrp.end(); IT++) {
				omaps << "\t" << itos(IT->second) + ":" + IT->first;
			}
		}
		omaps << "\n";
		b_wroteMapHD = true;
	}
	
	//convert to vector, that can than be written out
	vector<shared_ptr<DNAunique>> dereplicated_dnas (Tracker.size());
	int count = 0;
	HashDNA::iterator dd;
	for (dd = Tracker.begin(); dd != Tracker.end();dd++) {
		//super simple, double check later for correctness TODO
		size_t xsi = dd->second.size();
		if (xsi > 1) {
			int y = 0;
		}
        dereplicated_dnas[count] = dd->second.best();
		count++;
	}
//	bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) {	return l->getCount() < r->getCount();}
	sort(dereplicated_dnas.begin(), dereplicated_dnas.end() , DNAuPointerCompare);
	totSize = 0; size_t passed_hits(0);
	//bool thrHit = false;

	//sanity check
	vector<int> counts_per_sample(barcode_number_to_sample_id_.size(), 0);
	int total_count(0);
	ofstream* derepNowOut;
	//print unique DNAs
	for (size_t i = 0; i < dereplicated_dnas.size(); i++) {
		shared_ptr<DNAunique> dna = dereplicated_dnas[i];;
        total_count ++;
		dna->Count2Head(b_usearch_fmt);
		shared_ptr<DNA> dna_merged = nullptr;
		if ((dereplicate_sample_specific && 
			dna->pass_deprep_smplSpc(dereplicate_sample_specific_indices)) ||
            pass_deprep_conditions(dna) ) {
			//we do have a real derep that needs to be clustered..
			
			passed_hits++;
			totSize += dna->totalSum();

			//only do merge here, because these are known good dereps already
			if (b_merge_pairs_derep_ && dna->merge_seed_pos_ >= 0) {
				dna_merged = merger->merge(dna, dna->getPair());
			}
			derepNowOut = &of;
        } else {
			//dna->writeSeq(ofRest, b_singleLine);
			derepNowOut = &ofRest;
		}

		if (b_derep_as_fasta_)
			if (dna_merged) {//either or mechanic
				dna_merged->writeSeq(of_merged, b_singleLine);
			} else {
				dna->writeSeq((*derepNowOut), b_singleLine);
			}
		else {
			if (dna_merged) {
				dna_merged->prepareWrite(fastqVer);
				dna_merged->writeFastQ(of_merged);
			} else {
				dna->prepareDerepQualities(fastqVer);
				dna->writeDerepFastQ((*derepNowOut));
			}

		}


		//map is written no matter what
		dna->writeMap(omaps, dna->getId(), counts_per_sample, smplId2comb);

		//write out full seq + fastq
		dna->resetTruncation();
		dna->prepareWrite(fastqVer);
		dna->writeFastQ(of2);

		if (b_pairedInput) {
			shared_ptr<DNAunique> oD = dna->getPair();
			if (oD != nullptr) {
				oD->resetTruncation();
				oD->prepareWrite(fastqVer);
				oD->writeFastQ(of2p2);
			} else {
				shared_ptr<DNA> tmp = make_shared< DNA>("", dna->getId());
				tmp->writeFastQEmpty(of2p2);

			}
		}
	}
//	if (tmpCnt != totSize) {
//		cerr << "Counting failed\n" << tmpCnt << " " << totSize<<endl;
//	}
	of.close(); omaps.close();  ofRest.close();
	float avgSize = (float)totSize / (float)(passed_hits);
	string report = "";
	report += "Dereplication:\nAccepted " + intwithcommas((int)passed_hits) + " unique sequences ( "
              + itos(total_count) + " counts, " + minCopiesStr;
	if (dereplicate_sample_specific) { report += " & sample specific restrictions"; }
	report += " )";
	if (passed_hits > 0) {
		report += "; average size in this set is " + ftos(avgSize) + ".\nUniques with insufficient abundance : " + intwithcommas(int(Tracker.size() - passed_hits)) + " not passing derep conditions\n";
	}
	//cerr << tmpCnt << endl;
	cerr << "\n" << report << endl << endl;
	if (this->b_derepPerSR) {
		return report;
	}

	int uneqCnts(0);
	//check counts_per_sample vector
	vector<string> SampleID = mf->SampleID;
	for (unsigned int i = 0; i<SampleID.size(); i++) {
		size_t j(0); bool detected = false;
		for (; j < barcode_number_to_sample_id_.size(); j++) {
			if (SampleID[i] == barcode_number_to_sample_id_[j]) { detected = true;  break; }
		}
		if (detected && mf->collectStatistics[0]->BarcodeDetected[i] != counts_per_sample[j] && mf->collectStatistics[0]->BarcodeDetected[i] != -1) {
			//cerr << "ERROR: Unequal counts for " << SampleID[i] << ": " << mf->collectStatistics.BarcodeDetected[i] << " vs " << counts_per_sample[j] << endl;
			uneqCnts++;
			
		} else if (!detected) {
			cerr << "Could not detect sample_id_ " << barcode_number_to_sample_id_[j] << " in ref set (check that mapping file is correctly formatted?).\n"; exit(87);
		}
	}
#ifdef DEBUG
	cerr << "Derep Fin" << endl;
#endif
	if (b_pairedInput) { of2.close(); }
	if (uneqCnts>0) {
		cerr << "Unequal counts in " << uneqCnts << " cases. \n";
		//exit(66);
	}
	return report;
}


bool Dereplicate::pass_deprep_conditions(shared_ptr<DNAunique> d) {
	vector<int> x = d->getDerepMapSort(minCopiesSiz);
	int cumSum(0);
	for (size_t i = 0; i < x.size(); i++) {
		cumSum += x[i];
		//if (minCopies[i] == -1) { continue; }
		size_t yy(minCopiesSiz);
		if (yy > i) { yy = i+1; }//if checking for 1 smpl copy, don't need to check for 2 samlpe allowed copy number..
		for (size_t j = 0; j < yy; j++) {
			if (minCopies[j] == -1) { continue; }

			//if (i > minCopiesSiz) { break; }
			if (cumSum >= minCopies[j]) { return true; }
		}
	}
	return false;
}

void writeLog(string& logf) {
	ofstream of;
	of.open(logf.c_str());
	of.close();
}

string additionalFileName(const string& in){
	if (in.find(",") == std::string::npos){
		return additionalFileName2(in);
	}
	else {
		vector<string> vi = splitByCommas(in);
		string out = additionalFileName2(vi[0]);
		for (uint i = 1; i < vi.size(); i++){
			out += "," + additionalFileName2(vi[i]);
		}
		return out;
	}
}
string additionalFileName2(const string& in){
	size_t point = in.find_last_of('.');
	if (point == (size_t)-1){ return in + ".add"; }
	return in.substr(0, point) + ".add." + in.substr(point + 1);
}

string subfile(string x, string y) {
	if (y == "") { return x; }
	size_t pos = x.find_last_of(".");
	string s = "";
	if (pos != string::npos) {
		s= x.substr(0, pos) + "." + y + x.substr(pos);
	}
	else {
		s = x + "." + y;
	}
	return s;
}


//*******************************************
//*        FILTERS OBJECT
//*******************************************


Filters::Filters(OptContainer& cmdArgs1) :
        PrimerL(0), PrimerR(0), PrimerL_RC(0), PrimerR_RC(0), PrimerIdx(0),
        Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
        HeadSmplID(0),
        hetPrimer(2,vector<string>(0)),
        collectStatistics(2), csMTX(2,NULL) , statAddition(2),
        FastaF(0), QualF(0), FastqF(0), MIDfqF(0),
        derepMinNum(0),
	SequencingRun(0),

        lMD(NULL), 
        tAdapter(""), tAdapterLength(0),
        removeAdapter(false), bDoMultiplexing(true), bDoBarcode(true),
        bDoBarcode2(false), bDoBarcode2Rd1(false),
        bDoHeadSmplID(false), bBarcodeSameSize(false),
        bOneFileSample(false), curBCnumber(-1), BCoffset(0),
        bAdditionalOutput(false), b2ndRDBcPrimCk(false),
        bRevRdCk(false), bChkRdPrs(true),
        min_l(0), alt_min_l(0), min_l_p(-1.f), alt_min_l_p(-1.f),
        maxReadLength(0), norm2fiveNTs(false),
        max_l(10000), min_q(0.f), alt_min_q(0.f),
        BcutPrimer(true), alt_BcutPrimer(true), bPrimerR(false),
        bRequireRevPrim(false), alt_bRequireRevPrim(false),
        bRequireFwdPrim(false), alt_bRequireFwdPrim(false), BcutTag(true),
        bCompletePairs(false), bShortAmplicons(false),
        minBCLength1_(0), minBCLength2_(0), maxBCLength1_(0), maxBCLength2_(0), minPrimerLength_(0), maxHomonucleotide(0),
        PrimerErrs(0), alt_PrimerErrs(0), barcodeErrors_(0),
        MaxAmb(-1), alt_MaxAmb(-1),
        FQWwidth(0), EWwidth(0),
        RevPrimSeedL(5),
        b_BinFilBothPairs(false),
        BinFilErr(2.5), BinFilP(-1.f),
        alt_FQWthr(0), alt_EWthr(0),
        PEheaderVerWr(0), TrimStartNTs(0), TruncSeq(-1),
        userReqFastqVer(0), userReqFastqOutVer(33), maxAccumQP(-1),
        alt_maxAccumQP(-1),
        pairedSeq(-1),
        revConstellationN(0),
        BCdFWDREV(2),
        Xreads(-1),
        restartSet(false), b_optiClusterSeq(false),
        b_subselectionReads(false), b_doQualFilter(true),
        b_doFilter(true),
        bDoDereplicate(false), bDoCombiSamples(false),
        maxReadsPerOFile(0), ReadsWritten(0), OFileIncre(0),
		demultiBPperSR(0),
        barcodeLengths1_(0), barcodeLengths2_(0),
        cmdArgs(cmdArgs1)
		{
	csMTX[0] = new mutex();	csMTX[1] = new mutex();
	//csMTX[0].unlock(); 
	//csMTX[1].unlock();

    collectStatistics[0] = make_shared<collectstats>(); 
	collectStatistics[1] = make_shared<collectstats>();

	statAddition[0] = make_shared<collectstats>(); 
	statAddition[1] = make_shared<collectstats>();
	
	bool alt_bRequireRevPrimSet=false;

	string optF ("");
	if (cmdArgs.find("-options") != cmdArgs.end()) {
		optF = (cmdArgs)["-options"];
	}

	iniSpacer = (cmdArgs)["-sample_sep"];
	//***************************************
	//default options
	int maxAmb(0),PrimerErrs(1),TagErrs(0);
	float minQual(25);
	float minL(250.f);
	int maxL(1000);
	int QualWinWidth = 50;
	float QualWinThr = 0;
	int EndWinWidth = 15;
	float EndWinThr = 20;
	int maxHomoNT(12);
	bool keepTag(false),keepPrimer(false);
	bool addModConf = false;

	//set up some basic objects
	if ( (cmdArgs).find("-paired") != cmdArgs.end() ) {
		pairedSeq = atoi((cmdArgs)["-paired"].c_str()); //fakeEssentials();
		if ( pairedSeq<1 || pairedSeq>3 ) { cerr << "Argument \"-paired\" supplied with unknown parameter. Aborting.\n"; exit(28); }
		if ( (cmdArgs)["-onlyPair"] == "1" || (cmdArgs)["-onlyPair"] == "2" ) {
			pairedSeq = 1;
		}
	}
	if (cmdArgs.find("-normRdsToFiveNTs") != cmdArgs.end()) {
		norm2fiveNTs = true;
		cerr << "Warning: normRdsToFiveNTs is not implemented!\n";
	}
	//delimit output file size to X reads
	if (cmdArgs.find("-maxReadsPerOutput") != cmdArgs.end()) {
		maxReadsPerOFile = atoi(cmdArgs["-maxReadsPerOutput"].c_str());
	}
	if (cmdArgs.find("-DemultiBPperSR") != cmdArgs.end()) {
		stringstream ss(cmdArgs["-DemultiBPperSR"]);
		double d = 0;
		ss >> d;
		demultiBPperSR =(uint) d;
	}
	//important for fastq format
	if ( cmdArgs.find("-i_qual_offset") != cmdArgs.end() ) {
		if ( (cmdArgs)["-i_qual_offset"] == "auto" ) {
			userReqFastqVer = 0;
		} else {
			userReqFastqVer = atoi((cmdArgs)["-i_qual_offset"].c_str());
		}
	}
	//cerr<<cmdArgs["-o_qual_offset"]<<endl;
	userReqFastqOutVer = atoi((cmdArgs)["-o_qual_offset"].c_str());
	//statistic tracker
		//do new SEED sequence selection?
	if ( cmdArgs.find("-optimalRead2Cluster") != cmdArgs.end() ) {
		b_optiClusterSeq = true;
	}
	//do selection of specific reads?
	if ( (cmdArgs)["-specificReads"] != "" ) {
		b_subselectionReads = true;
	}
	if (cmdArgs.find("-binomialFilterBothPairs") != cmdArgs.end() && (cmdArgs)["-binomialFilterBothPairs"] == "1") {
		b_BinFilBothPairs = true;
	}
	//***************************************
	//read options
	ifstream opt;
	opt.open(optF.c_str(),ios::in);
	if ( !opt || optF=="") {
		cerr << "NO filtering will be done on your reads (just rewriting / log files created)." << endl;
		b_doFilter = false;
		return;
	}
	string line;
	while (getline(opt,line,'\n')){
		
		if (line.length()<=1 || line.substr(0,1)=="#"){
			continue;
		}
		
		bool addMod = false;
		if (line.substr(0,1)=="*"){
			addMod= true;
			line = line.substr(1);
		}
		string segs;
		string segs2;
		stringstream ss;
		ss << line;
		getline(ss,segs,'\t');
		getline(ss,segs2,'\t');

		if ((cmdArgs)["-XfirstReads"] != "") {
			Xreads = atoi((cmdArgs)["-XfirstReads"].c_str());
		}


		if (strcmp(segs.c_str(),"minSeqLength") == 0){
			if (addMod){ 
				 float tmp = (float)atof(segs2.c_str());
				 if (tmp>1.f) {
					 alt_min_l = (int)tmp;
				 } else {
					 alt_min_l = -1;
					 alt_min_l_p = tmp;
				 }
				if (alt_min_l != minL){	addModConf = true;	}
			} else {
				minL = (float)atof(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"maxSeqLength") == 0){
			maxL = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"minAvgQuality") == 0){
			if (addMod){ 
				alt_min_q = (float) atof(segs2.c_str());
				if (alt_min_q!= minQual){addModConf = true;}
			} else {
				minQual = (float) atof(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"maxAmbiguousNT") == 0){
			if (addMod){ 
				alt_MaxAmb = atoi(segs2.c_str());
				if (MaxAmb!= alt_MaxAmb){addModConf = true;}
			} else {
				maxAmb = atoi(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"QualWindowThreshhold") == 0){
			if (addMod){ 
				alt_FQWthr = (float) atof(segs2.c_str());
				if (alt_FQWthr!= QualWinThr){addModConf = true;}
			} else {
				QualWinThr = (float) atof(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "QualWindowWidth") == 0){
			QualWinWidth = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "BinErrorModelMaxExpError") == 0){
			BinFilErr = (float) atof(segs2.c_str());
			if (BinFilErr < 0){
				cerr << "BinErrorModelMaxExpError was set to <0. Set to 0 instead.\n";
				BinFilErr = 0;
			}
		}
		else if (strcmp(segs.c_str(), "BinErrorModelAlpha") == 0){
			BinFilP = (float)atof(segs2.c_str());
			if (BinFilP != -1.f && (BinFilP<0.f || BinFilP>1.f)){
				cerr << "BinErrorModelAlpha has to be between 0 and 1 (or -1 to deactivate).\nAborting..\n";
				exit(542);
			}

		} else if (strcmp(segs.c_str(),"TrimWindowWidth") == 0){
			EndWinWidth = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"TrimWindowThreshhold") == 0){
			if (addMod){ 
				alt_EWthr = (float) atof(segs2.c_str());
				if (alt_EWthr != EndWinThr){addModConf = true;}
			} else {
				EndWinThr = (float) atof(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"maxBarcodeErrs") == 0){
			TagErrs = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"maxPrimerErrs") == 0){
			if (addMod){ 
				alt_PrimerErrs = atoi(segs2.c_str());
				if (alt_PrimerErrs!=PrimerErrs){ addModConf = true;}
			} else {
				PrimerErrs = atoi(segs2.c_str());
			}
		} else if (strcmp(segs.c_str(),"keepBarcodeSeq") == 0){
			atoi(segs2.c_str())==0 ? keepTag=false : keepTag=true;
		} else if (strcmp(segs.c_str(),"keepPrimerSeq") == 0){
			if (addMod){ 
				atoi(segs2.c_str())==0 ? alt_BcutPrimer=false : alt_BcutPrimer=true;
				if(alt_BcutPrimer!=keepPrimer) {addModConf = true;}
			} else {
				atoi(segs2.c_str())==0 ? keepPrimer=false : keepPrimer=true;
			}
		} else if (strcmp(segs.c_str(),"maxHomonucleotide") == 0){
			maxHomoNT = atoi(segs2.c_str());
		} else if (strcmp(segs.c_str(),"maxAccumulatedError") == 0){
			if (addMod){
				alt_maxAccumQP = double(atof(segs2.c_str()));
				if(alt_maxAccumQP!=maxAccumQP) {addModConf = true;}
			} else{
				maxAccumQP = double(atof(segs2.c_str()));
			}
		} else if (strcmp(segs.c_str(),"TechnicalAdapter") == 0){
			tAdapter = segs2.c_str();
			transform(tAdapter.begin(), tAdapter.end(),tAdapter.begin(), ::toupper);
			tAdapterLength = (int) tAdapter.length();
            removeAdapter = true;
		}else if (segs == "PEheaderPairFmt"){
			PEheaderVerWr = atoi(segs2.c_str());
		} else if (segs == "TrimStartNTs"){
			TrimStartNTs = atoi(segs2.c_str());
		} else if (segs == "fastqVersion"){
			if (segs2 == "auto") {
				userReqFastqVer = 0;
			} else {
				userReqFastqVer = FastqVerMod(atoi(segs2.c_str()));
			}
		} else if (segs == "RejectSeqWithoutRevPrim"){
			if (addMod){ 
				alt_bRequireRevPrimSet=true;
				if (segs2=="T"){alt_bRequireRevPrim=true;
				} else {alt_bRequireRevPrim=false;}				
				if(alt_bRequireRevPrim!=bRequireRevPrim){addModConf = true;}
			} else {
				if (segs2=="T"){bRequireRevPrim=true;
				} else {bRequireRevPrim=false;}
			}
		} else if (segs == "RejectSeqWithoutFwdPrim"){
			if (addMod){ 
				alt_bRequireFwdPrim=true;
				if (segs2=="T"){alt_bRequireFwdPrim=true;
				} else {alt_bRequireFwdPrim=false;}				
				if(alt_bRequireFwdPrim!=bRequireFwdPrim){addModConf = true;}
			} else {
				if (segs2=="T"){bRequireFwdPrim=true;
				} else {bRequireFwdPrim=false;}
			}
		} else if (segs == "TruncateSequenceLength") {
			TruncSeq = atoi(segs2.c_str());
			if (TruncSeq != -1 && TruncSeq < (int)minL) { minL = (float)TruncSeq; }
		} else if (segs == "AmpliconShortPE") {
			if (segs2 == "T") {
				bShortAmplicons = true;
			} else { bShortAmplicons = false; }
		} else if (segs == "CheckForMixedPairs") {
			if (segs2 == "T") {
				b2ndRDBcPrimCk = true;
			} else { b2ndRDBcPrimCk = false; }
		} else if ( segs == "CheckForReversedSeqs" ) {
			if ( segs2 == "T" ) {
				bRevRdCk = true;
			} else { bRevRdCk = false; }
		}
		else if (segs == "SyncReadPairs") {
			if (segs2 == "T") {
				bChkRdPrs = true;
			}
			else { bChkRdPrs = false; }
		}
		

	}
	
	//report some non-std options
	if (bShortAmplicons){
		cerr << "Checking for reverse primers on 1st read.\n";
	}
	if (b2ndRDBcPrimCk){
		cerr << "Checking for switched pairs.\n";
	}

	opt.close();
	//set in filter object
	this->setSeqLength(minL,maxL);
	this->setPrimerErrs(PrimerErrs);
	this->setTagErrs(TagErrs);
	this->removePrimer(!keepPrimer);
	this->removeTag(!keepTag);
	this->setMaxAmb(maxAmb);
	this->setAvgMinQual(minQual);
	this->setFloatingQWin(QualWinWidth,QualWinThr);
	this->setFloatingEWin(EndWinWidth,EndWinThr);
	this->setMaxHomo(maxHomoNT);
	//alternative options (mid qual filtering)
	if (addModConf){
		if (!alt_bRequireRevPrimSet){alt_bRequireRevPrim = bRequireRevPrim;}
		if (cmdArgs.find("-o_fna")  != cmdArgs.end() && (cmdArgs)["-o_fna"].length()>1){
			if (cmdArgs.find("-o_fna2")  == cmdArgs.end()){
				(cmdArgs)["-o_fna2"] = additionalFileName((cmdArgs)["-o_fna"]);
				//cmdArgs["-o_fna2"] = cmdArgs["-o_fna"].substr(0,cmdArgs["-o_fna"].length()-4)+".add.fna";
			}
		} else if (cmdArgs.find("-o_fastq")  != cmdArgs.end() && (cmdArgs)["-o_fastq"].length()>1){
			if (cmdArgs.find("-o_fastq2")  == cmdArgs.end()){
				(cmdArgs)["-o_fastq2"] = additionalFileName((cmdArgs)["-o_fastq"]);
			}
		}
		bAdditionalOutput=true;
	}
}



Filters::Filters(Filters* of, int BCnumber, bool takeAll, size_t threads) 
	:
	PrimerL(0,""), PrimerR(0,""),
        PrimerL_RC(0,""), PrimerR_RC(0,""),
        PrimerIdx(of->PrimerIdx),
        Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
	    HeadSmplID(0), hetPrimer(2, vector<string>(0)),
		collectStatistics(2,nullptr), csMTX(2, NULL), statAddition(2,nullptr),
        FastaF(of->FastaF), QualF(of->QualF), FastqF(of->FastqF),
        MIDfqF(of->MIDfqF), derepMinNum(of->derepMinNum),
        lMD(nullptr), 
        tAdapter(of->tAdapter), tAdapterLength(of->tAdapterLength),
        removeAdapter(of->removeAdapter), bDoMultiplexing(of->bDoMultiplexing),
        bDoBarcode(of->bDoBarcode), bDoBarcode2(of->bDoBarcode2), bDoBarcode2Rd1(of->bDoBarcode2Rd1),
		bDoHeadSmplID(of->bDoHeadSmplID),
        bBarcodeSameSize(of->bBarcodeSameSize),
        bOneFileSample(of->bOneFileSample), curBCnumber(BCnumber), BCoffset(of->BCoffset),
        bAdditionalOutput(of->bAdditionalOutput), b2ndRDBcPrimCk(of->b2ndRDBcPrimCk),
        bRevRdCk(of->bRevRdCk), bChkRdPrs(of->bChkRdPrs),
        min_l(of->min_l), alt_min_l(of->alt_min_l), min_l_p(of->min_l_p), alt_min_l_p(of->alt_min_l_p),
        maxReadLength(0), norm2fiveNTs(of->norm2fiveNTs),
        max_l(of->max_l), min_q(of->min_q), alt_min_q(of->alt_min_q),
        BcutPrimer(of->BcutPrimer), alt_BcutPrimer(of->alt_BcutPrimer),
        bPrimerR(of->bPrimerR),
        bRequireRevPrim(of->bRequireRevPrim), alt_bRequireRevPrim(of->alt_bRequireRevPrim),
        bRequireFwdPrim(of->bRequireFwdPrim), alt_bRequireFwdPrim(of->alt_bRequireFwdPrim),
        BcutTag(of->BcutTag),
        bCompletePairs(of->bCompletePairs), bShortAmplicons(of->bShortAmplicons),
        minBCLength1_(of->minBCLength1_), minBCLength2_(of->minBCLength2_), maxBCLength1_(of->maxBCLength1_), maxBCLength2_(of->maxBCLength2_), minPrimerLength_(of->minPrimerLength_), maxHomonucleotide(of->maxHomonucleotide),
        PrimerErrs(of->PrimerErrs), alt_PrimerErrs(of->alt_PrimerErrs), barcodeErrors_(of->barcodeErrors_),
        MaxAmb(of->MaxAmb), alt_MaxAmb(of->alt_MaxAmb),
        FQWwidth(of->FQWwidth), EWwidth(of->EWwidth),
        RevPrimSeedL(of->RevPrimSeedL),
        b_BinFilBothPairs(of->b_BinFilBothPairs),
        BinFilErr(of->BinFilErr), BinFilP(of->BinFilP),
        FQWthr(of->FQWthr), EWthr(of->EWthr),
        alt_FQWthr(of->alt_FQWthr), alt_EWthr(of->alt_EWthr),
        PEheaderVerWr(of->PEheaderVerWr), TrimStartNTs(of->TrimStartNTs),
        TruncSeq(of->TruncSeq),
        iniSpacer(of->iniSpacer), userReqFastqVer(of->userReqFastqVer),
        userReqFastqOutVer(of->userReqFastqOutVer), maxAccumQP(of->maxAccumQP),
        alt_maxAccumQP(of->alt_maxAccumQP),
        //BChit, BCrevhit initialize to 0 - new set, new luck
        pairedSeq(of->pairedSeq),
        revConstellationN(0),
        BCdFWDREV(of->BCdFWDREV),
        restartSet(false),
        b_optiClusterSeq(of->b_optiClusterSeq), b_subselectionReads(of->b_subselectionReads),
        b_doQualFilter(of->b_doQualFilter),
        b_doFilter(of->b_doFilter),
        bDoDereplicate(of->bDoDereplicate),
        bDoCombiSamples(of->bDoCombiSamples),
        maxReadsPerOFile(of->maxReadsPerOFile),
		demultiBPperSR(of->demultiBPperSR),
        //ReadsWritten(of->ReadsWritten), OFileIncre(of->OFileIncre),
        barcodeLengths1_(0), barcodeLengths2_(0), SequencingRun(0)
{
	cdbg("New Filter object from copy\n");
	ReadsWritten = of->writtenReads();
	OFileIncre = of->getFileIncrementor();
    BCdFWDREV[0].fix(); BCdFWDREV[1].fix();
	//collectStatistics.resize(2); statAddition.resize(2);
	csMTX[0] = new mutex(); csMTX[1] = new mutex();

	collectStatistics[0] = make_shared<collectstats>(); 
	collectStatistics[1] = make_shared<collectstats>();
    statAddition[0] = make_shared<collectstats>(); statAddition[1] = make_shared<collectstats>(); ;
	cdbg("New Filter::resize collectStatistics done\n");
	if (takeAll) {
        this->allResize((uint)of->PrimerIdx.size());
        PrimerIdxRev = of->PrimerIdxRev;
        PrimerIdx = of->PrimerIdx;
        Barcode = of->Barcode;
        Barcode2 = of->Barcode2;
        SampleID = of->SampleID;
        SampleID_Combi = of->SampleID_Combi;
        HeadSmplID = of->HeadSmplID;
        PrimerL = of->PrimerL;
        PrimerR = of->PrimerR;
        PrimerL_RC = of->PrimerL_RC;
        PrimerR_RC = of->PrimerR_RC;
        hetPrimer = of->hetPrimer;
        lMD = of->lMD;
        barcodeLengths1_ = of->barcodeLengths1_;
        barcodeLengths2_ = of->barcodeLengths2_;
		SequencingRun = of->SequencingRun;
		SequencingRun2id = of->SequencingRun2id;
BarcodePreStats();

	}

}

/*
//copy of filter, but leave stat vals empty
Filters::Filters(Filters* of, int BCnumber, bool takeAll, size_t threads) :
		PrimerL(0,""), PrimerR(0,""),
		PrimerL_RC(0,""), PrimerR_RC(0,""),
		PrimerIdx(of->PrimerIdx),
		Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
		HeadSmplID(0), hetPrimer(2, vector<string>(0)),
		collectStatistics(2), statAddition(2), statistics_(threads),
		FastaF(of->FastaF), QualF(of->QualF), FastqF(of->FastqF),
		MIDfqF(of->MIDfqF), derepMinNum(of->derepMinNum),
		lMD(NULL),
		tAdapter(of->tAdapter), tAdapterLength(of->tAdapterLength),
		removeAdapter(of->removeAdapter), bDoMultiplexing(of->bDoMultiplexing),
		bDoBarcode(of->bDoBarcode), bDoBarcode2(of->bDoBarcode2), bDoBarcode2Rd1(of->bDoBarcode2Rd1),
		bDoHeadSmplID(of->bDoHeadSmplID),
		bBarcodeSameSize(of->bBarcodeSameSize),
		bOneFileSample(of->bOneFileSample), curBCnumber(BCnumber), BCoffset(0),
		bAdditionalOutput(of->bAdditionalOutput), b2ndRDBcPrimCk(of->b2ndRDBcPrimCk),
		bRevRdCk(of->bRevRdCk), bChkRdPrs(of->bChkRdPrs),
		min_l(of->min_l), alt_min_l(of->alt_min_l), min_l_p(of->min_l_p), alt_min_l_p(of->alt_min_l_p),
		maxReadLength(0), norm2fiveNTs(of->norm2fiveNTs),
		max_l(of->max_l), min_q(of->min_q), alt_min_q(of->alt_min_q),
		BcutPrimer(of->BcutPrimer), alt_BcutPrimer(of->alt_BcutPrimer),
		bPrimerR(of->bPrimerR),
		bRequireRevPrim(of->bRequireRevPrim), alt_bRequireRevPrim(of->alt_bRequireRevPrim),
		bRequireFwdPrim(of->bRequireFwdPrim), alt_bRequireFwdPrim(of->alt_bRequireFwdPrim),
		BcutTag(of->BcutTag),
		bCompletePairs(of->bCompletePairs), bShortAmplicons(of->bShortAmplicons),
		minBCLength1_(of->minBCLength1_), minBCLength2_(of->minBCLength2_), maxBCLength1_(of->maxBCLength1_), maxBCLength2_(of->maxBCLength2_), minPrimerLength_(of->minPrimerLength_), maxHomonucleotide(of->maxHomonucleotide),
		PrimerErrs(of->PrimerErrs), alt_PrimerErrs(of->alt_PrimerErrs), barcodeErrors_(of->barcodeErrors_),
		MaxAmb(of->MaxAmb), alt_MaxAmb(of->alt_MaxAmb),
		FQWwidth(of->FQWwidth), EWwidth(of->EWwidth),
		RevPrimSeedL(of->RevPrimSeedL),
		b_BinFilBothPairs(of->b_BinFilBothPairs),
		BinFilErr(of->BinFilErr), BinFilP(of->BinFilP),
		FQWthr(of->FQWthr), EWthr(of->EWthr),
		alt_FQWthr(of->alt_FQWthr), alt_EWthr(of->alt_EWthr),
		PEheaderVerWr(of->PEheaderVerWr), TrimStartNTs(of->TrimStartNTs),
		TruncSeq(of->TruncSeq),
		iniSpacer(of->iniSpacer), userReqFastqVer(of->userReqFastqVer),
		userReqFastqOutVer(of->userReqFastqOutVer), maxAccumQP(of->maxAccumQP),
		alt_maxAccumQP(of->alt_maxAccumQP),
		//BChit, BCrevhit initialize to 0 - new set, new luck
		pairedSeq(of->pairedSeq),
		revConstellationN(0),
		BCdFWDREV(of->BCdFWDREV),
		restartSet(false),
		b_optiClusterSeq(of->b_optiClusterSeq), b_subselectionReads(of->b_subselectionReads),
		b_doQualFilter(of->b_doQualFilter),
		b_doFilter(of->b_doFilter),
		bDoDereplicate(of->bDoDereplicate),
		bDoCombiSamples(of->bDoCombiSamples),
		maxReadsPerOFile(of->maxReadsPerOFile),
		ReadsWritten(of->ReadsWritten), OFileIncre(of->OFileIncre),
		barcodeLengths1_(0), barcodeLengths2_(0), SequencingRun(0)
		{
	BCdFWDREV[0].fix(); BCdFWDREV[1].fix();
	collectStatistics[0] = collectstats(); collectStatistics[1] = collectstats();
	collectStatistics[0]->ini_repStat(); collectStatistics[1]->ini_repStat();

	if (bAdditionalOutput) {
		statAddition[0] = collectstats(); statAddition[1] = collectstats();
		statAddition[0]->ini_repStat(); statAddition[1]->ini_repStat();

	}
	if (takeAll) {
		this->allResize((uint)of->PrimerIdx.size());
		PrimerIdxRev = of->PrimerIdxRev;
		PrimerIdx = of->PrimerIdx;
		Barcode = of->Barcode;
		Barcode2 = of->Barcode2;
		SampleID = of->SampleID;
		SampleID_Combi = of->SampleID_Combi;
		HeadSmplID = of->HeadSmplID;
		PrimerL = of->PrimerL;
		PrimerR = of->PrimerR;
		PrimerL_RC = of->PrimerL_RC;
		PrimerR_RC = of->PrimerR_RC;
		hetPrimer = of->hetPrimer;
		lMD = of->lMD;
		barcodeLengths1_ = of->barcodeLengths1_;
		barcodeLengths2_ = of->barcodeLengths2_;
		SequencingRun = of->SequencingRun;
		SequencingRun2id = of->SequencingRun2id;
		BarcodePreStats();

	}

}
*/

Filters::~Filters() {
	cdbg("Deleting filter .. ");
	for (size_t i = 0; i < csMTX.size(); i++){
		delete csMTX[i];
	}

//	for (size_t i = 0; i < 2; i++) { delete PostFilt[i]; delete RepStatAddition[i]; }
//	delete PreFiltP1; delete PreFiltP2;
	cdbg("Done\n");
}


Filters* Filters::filterPerBCgroup(const vector<int> idxi) {
	cdbg("filterPerBCgroup::start : "+ itos(idxi[0])+"\n");

	// get filter from main filter object passing an index for mapping?!
//	shared_ptr<Filters> filter = make_shared<Filters>(shared_from_this(), idxi[0]);
	Filters* filter = new Filters(this, idxi[0]);
	cdbg("filterPerBCgroup::star2t\n");

	// number of mapping file lines associated with that unique fastx
	unsigned int tarSize = (unsigned int)idxi.size();
	filter->allResize(tarSize);

	int tarID = -1;
	bool isDoubleBarcoded = this->doubleBarcodes();
	cdbg("filterPerBCgroup::Go over BCs\n");

	// iterate over every occurence of unique fa
	for (unsigned int j = 0; j < tarSize; j++) { //fill in filter
		// Get id_ of file in tar
		tarID = idxi[j];
		if (this->PrimerIdx[tarID] > -1) {
			filter->addPrimerL(this->PrimerL[this->PrimerIdx[tarID]], j);
		}
		if (this->doReversePrimers() && this->PrimerIdxRev[tarID] > -1) {
			filter->addPrimerR(this->PrimerR[this->PrimerIdxRev[tarID]], j);
		}
		filter->Barcode[j] = this->Barcode[tarID];
		if (isDoubleBarcoded) {
			filter->Barcode2[j] = this->Barcode2[tarID];
		}
		filter->SampleID[j] = this->SampleID[tarID];
		filter->SampleID_Combi[j] = this->SampleID_Combi[tarID];
		filter->HeadSmplID[j] = this->HeadSmplID[tarID];
	}
	cdbg("filterPerBCgroup::check 4 doubles\n");
	//sanity check no double barcodes..
	filter->checkDoubleBarcode();
	return filter;
}


//service function to ini what needs to be done
UClinks * Filters::ini_SeedsReadsDerep(UClinks *ucl, shared_ptr<ReadSubset>& RDSset, 
	shared_ptr<Dereplicate>& Dere) {
	if (this->doOptimalClusterSeq()) {
		ucl = new UClinks(cmdArgs);
		if (cmdArgs.find("-mergedPairs") != cmdArgs.end() && cmdArgs["-mergedPairs"] == "1") {
			ucl->pairedSeqsMerged();
			this->setFloatingEWin(0, 0.f);
		}
		else {
			this->setFloatingEWin(10, 25);
		}
		//are fallback fasta sequences available?
		if (cmdArgs["-OTU_fallback"] != "") {
			shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true,
				this->getuserReqFastqVer(), cmdArgs["-ignore_IO_errors"], cmdArgs["-pairedRD_HD_out"],1);
			FALL->setupFna(cmdArgs["-OTU_fallback"]);
			ucl->setupDefSeeds(FALL, SampleID);
		}
		ReadMerger* merg = new ReadMerger(); //create special object for these functions
		ucl->attachMerger(merg);
	}
	else if (this->doSubselReads()) {
		//this will select a list of reads and distribute these into multiple files
		RDSset = make_shared<ReadSubset>(cmdArgs["-specificReads"], "");
	}
	else if (this->doDereplicate()) {
		Dere = make_shared<Dereplicate>(cmdArgs, this);
		ReadMerger* merg = new ReadMerger(); //create special object for these functions
		Dere->attachMerger(merg);
	}
	return ucl;
}



//simulates that in mapping file links to sequence file was given.
bool Filters::setcmdArgsFiles(){

	if (FastqF.size()==0 && QualF.size()==0 && FastaF.size() > 0){
		//fasta entry but no qual entries

		string path="";
		if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() >= 1){
			path=cmdArgs["-i_path"] + string("/");
		}

		QualF.resize(FastaF.size());
		for (unsigned int i=0; i< FastaF.size(); i++){
			string newQ = FastaF[i];
			int pos = (int) newQ.find_last_of(".");
			newQ = newQ.substr(0,pos);
			newQ += string(".qual_");
			fstream fin;
			string fullQ = path + newQ;
			fin.open(fullQ.c_str(),ios::in);
			if( fin.is_open() )	{
				cerr<<"Using quality file: "<<fullQ <<endl;
			} else if (cmdArgs.find("-number")!=  cmdArgs.end() && cmdArgs["-number"] =="T"){
				;
			}else {
				cerr<<"You did not supply a quality file for"<<path+ FastaF[i]<<". \nPlease give the path to your quality file as command line argument:\n  -i_qual <PathToQualityFile>\n";
				newQ = "";
				//fin.close();return false;
			}
			fin.close();
			QualF[i] = newQ;
		}
	}

	int fileSiz = (int)Barcode.size();
	//instead of max:
	if (!bDoMultiplexing){fileSiz=1;}

	if (FastaF.size()==0 && FastqF.size()==0){
		//set up fasta/fastq vector specific to corresponding BC (that should be in this file)
		if (cmdArgs.find("-i_fastq")  == cmdArgs.end()){
			FastaF.resize(fileSiz);
			QualF.resize(fileSiz);
			for (unsigned int i=0; i< FastaF.size(); i++){
				FastaF[i] = cmdArgs["-i_fna"];
				QualF[i] = cmdArgs["-i_qual"];
			}
		} else {//fastq input
			vector<string> fqTmp (1,cmdArgs["-i_fastq"]);
			if (cmdArgs["-i_fastq"].find(";") != string::npos) {//";" denotes several files
				if (fileSiz == 1) {//no BC, 
					fqTmp = splitByCommas(cmdArgs["-i_fastq"], ';');
					this->allResize((uint) fqTmp.size());
					fileSiz = (int) fqTmp.size();
					cerr << "Detected " << fileSiz << " input files (pairs)." << endl;
					FastqF = fqTmp;
				} else {
					cerr << "Fastq string contains symbol \";\". Not allowed in input string"; exit(32);
				}
			} else {
				FastqF.resize(fileSiz, cmdArgs["-i_fastq"]);
			}
		}
	}


	if (MIDfqF.size() == 0)
		if (cmdArgs.find("-i_MID_fastq") != cmdArgs.end()) {
		MIDfqF.resize(fileSiz, "");
		for (unsigned int i = 0; i < MIDfqF.size(); i++) {
			MIDfqF[i] = cmdArgs["-i_MID_fastq"];
		}
	}
		

	if (cmdArgs["-o_dereplicate"] != "") {
		//check if file could exist
		ofstream temp;
		temp.open(cmdArgs["-o_dereplicate"].c_str(), ios::out);
		if (!temp) { cerr << "Could not open outstream to dereplicated sequences:\n" << cmdArgs["- o_dereplicate"] << endl; exit(78); }
		temp.close();
		bDoDereplicate = true;
	}

	return true;
}


//only does BC 1
void Filters::reverseTS_all_BC(){
//	for (int i=0; i<Barcode.size();i++){
//		reverseTS(Barcode[i]);
//	}
	Barcode = revBarcode;
	revBarcode.resize(0);
	barcodes1_.clear();
	for (uint i = 0; i < Barcode.size(); i++){
        barcodes1_[Barcode[i]] = i;
        barcodeLengths1_[i] = (int) Barcode[i].length();
	}

}
void Filters::reverseTS_all_BC2() {
	//	for (int i=0; i<Barcode.size();i++){
	//		reverseTS(Barcode[i]);
	//	}
	Barcode2 = revBarcode2;
	revBarcode2.resize(0);
	barcodes2_.clear();
	for ( uint i = 0; i < Barcode2.size(); i++ ) {
        barcodes2_[Barcode2[i]] = i;
        barcodeLengths2_[i] = (int)Barcode2[i].length();
	}

}
bool Filters::swapReverseDNApairs(vector< shared_ptr<DNA>>& tdn){
	if (!this->checkRevRd()) {
		return false;
	}
	int tagIdx = 0;//just try

	if (BcutPrimer) {
		//1test if fwd read has primer 1
		if (cutPrimer(tdn[0], PrimerIdx[tagIdx], false, 0) ||
			//test if read2 has primer2
			cutPrimerRev(tdn[1], PrimerIdxRev[tagIdx], false)) {
			return false;
		}


		if (cutPrimerRev(tdn[0], PrimerIdxRev[tagIdx], false) ||
			cutPrimer(tdn[1], PrimerIdxRev[tagIdx], false, 0)) {
			//swap out
			shared_ptr<DNA> x = tdn[1];
			tdn[0] = tdn[1];
			tdn[1] = x;
			return true;
		}
		//reversed?
		tdn[0]->reverse_transcribe();
		tdn[1]->reverse_transcribe();
		if (cutPrimer(tdn[0], PrimerIdx[tagIdx], false, 0) ||
			cutPrimerRev(tdn[1], PrimerIdxRev[tagIdx], false)) {
			return true;
		}

		//no? back to normal..
		tdn[0]->reverse_transcribe();
		tdn[1]->reverse_transcribe();
		return false;
	}
	tagIdx = -2;
	string presentBC = ""; int c_err = 0; int chkRev1=false;
	tagIdx = findTag(tdn[0], presentBC, c_err, true, chkRev1);


	/*
	if (true && checkReversedRead && (tagIdx2 < 0 && tagIdx < 0)) {
		tdn[0]->reverse_transcribe(); tdn[1]->reverse_transcribe();
		Pr1 = curFil->findPrimer(tdn[0], 0, false, 0);
		Pr2 = curFil->findPrimer(tdn[1], 0, false, 0);
		tagIdx = curFil->findTag(tdn[0], presentBC, c_err, true, chkRev1);
		tagIdx2 = curFil->findTag(tdn[1], presentBC, c_err, true, chkRev2);
		revT = true;
	}
	//this is all about barcodes..
	if (checkReversedRead && tdn[0] != NULL && tagIdx < 0) {
		if (!MIDuse) { tagIdx = -2; }
		//		curFil->sTotalMinus(0);
		tdn[0]->reverse_transcribe();
		MD->analyzeDNA(tdn[0], -1, 0, tagIdx, curThread);
		ch1 = tdn[0]->isGreenQual();
		isReversed = ch1;
		if (!isReversed) {//reset
			tdn[0]->reverse_transcribe();
		}
	}
	*/


	return false;
}

void Filters::preFilterSeqStat(shared_ptr<DNA> d, int pair) {
	if (d == NULL)
		return;
	int easyPair = 1;
	if (pair <= 0) {
		easyPair = 0;
	}
	//csMTX[easyPair]->lock();
	collectStatistics[easyPair]->PreFilt->addDNAStats(d);
	updateMaxSeqL(d->length());
	//csMTX[easyPair]->unlock();
}

void Filters::preFilterSeqStatMT(shared_ptr<DNA> d, int pair, uint thread) {
    if (d == NULL)
        return;


    if (pair <= 0) {
        statistics_[thread].main_read_stats_[0]->PreFilt->addDNAStats(d);
    } else if (pair == 1) {
        statistics_[thread].main_read_stats_[1]->PreFilt->addDNAStats(d);
    }
    updateMaxSeqL(d->length());
}

std::mutex updateMaxSeqMutex;
void Filters::updateMaxSeqL(int x) {
    {
        std::lock_guard<std::mutex> updateMaxSeqLLck(updateMaxSeqMutex);
        if (x < maxReadLength) { return; }
        maxReadLength = x;
        if (min_l_p != -1.f) {
            min_l = (int) ((float) maxReadLength * min_l_p);
        }
    }
}
void Filters::setSeqLength(float minL, int maxL) {
	if (minL>1.f) {
		min_l = (int)minL;
		min_l_p = -1.f;
	} else {
		min_l = -1;
		min_l_p = minL;
	}
	max_l = maxL;
}


//ever_best is the best %id_ that was ever observed for this cluster match
bool Filters::betterSeed(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNA> ref, shared_ptr<DNA> ref2, float ever_best, 
	uint bestL, int usePair, bool checkBC) {
	float d1pid(d1->getTempFloat()), refpid(ref->getTempFloat());
	int TagIdx(0);
	if (checkBC) {
		TagIdx = -2;
	}
	//0.2% difference is still ok, but within 0.5% of the best found seed (prevent detoriating sequence match)
	//float blen = (float)ref->length() + (float)d1->length();
	uint curL = d1->length();
	if (d2 != NULL) {		curL += d2->length();	}
	if (float(curL) / float(bestL) < BestLengthRatio) { return false; }
	if (d1pid<refpid - 0.4f || d1pid < ever_best - 1){ return false; }

	//*** DNA1
	//needs to quality filter first
	if (!check(d1, true, usePair, TagIdx)) {
		return false;
	}
	//at least 90% length of "good" hit
	if (d1->length() / ref->length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//1 added to qual, in case no qual DNA is used
	float thScore = (1+d1->getAvgQual())*(d1pid ) * log((float)d1->length() );
	float rScore = (1+ref->getAvgQual())*(refpid ) * log((float)ref->length() );
	if (thScore > rScore){
		//also check for stable lowest score
		if (d1->minQual() > ref->minQual() - MinQualDiff && (d2 == NULL || ref2 == NULL)) { return true; }
	}
	if (d2 == NULL || ref2 == NULL) {
		return false;
	}
	//*** DNA2
	//second pair_ likely to be of worse qual_, but only direct comparison relevant here

	check(d2, true, 1, TagIdx);
	//at least 90% length of "good" hit
	if (d2->length() / ref2->length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//weigh with average id_ to OTU seed
	thScore +=(1+ d2->getAvgQual()) * log((float)d2->length()) * 97;
	rScore += (1+ref2->getAvgQual()) * log((float)ref2->length()) * 97;
	if (thScore > rScore) {
		return true;
	}

	return false;
}


bool Filters::check(shared_ptr<DNA> d, bool doSeeding, int pairPre,
	int &tagIdx) {
    
    //corrects for -1 (undefined pair_) to set to 0
	int pair = max(0, pairPre);
	unsigned int hindrance = 0;

	if (check_length(d->length())) {
		d->QualCtrl.minL = true;
		return false;
	}
	//remove technical adapter
	if (pair == -1 && removeAdapter)
	    remove_adapter(d);

	//BC already detected (e.g. MID)?
	if (tagIdx == -2){
		tagIdx = detectCutBC(d, pair == 0); //barcode 2nd part
	}
	if ((bDoBarcode2 || bDoBarcode) && tagIdx < 0) {
		 d->QualCtrl.TagFail = true;
		 d->failed();
		 return false;
	}
	
	//fwd primer only on pair_ 0
	if (BcutPrimer){
		if (pair != 1 && !cutPrimer(d, PrimerIdx[tagIdx], false,pair) && bRequireFwdPrim) {//0 or -1
		    d->failed();
		    return false;
		}
		else if (bShortAmplicons) {
		    //pair_ == 1, check for fwd primer in pair_ 2 (rev-compl)
			cutPrimer(d, PrimerIdx[tagIdx], true,pair);
		}
	}


	if (check_length(d->length())){
		d->QualCtrl.minL = true;
		return false;
	}
	if (max_l!=0 && d->length()-hindrance > max_l){
		d->QualCtrl.maxL = true;
		return false;
	}

	//rev primer is the first that needs to be looked for
	//makes it slower, as higher chance for low qual_ and this routine is costly.. however more important to get good lock on rev primer
	if ((pair != 0 || bShortAmplicons) && bPrimerR) {
		//removal of reverse primer
		bool revCheck = pair == -1 || pair == 0;//1:false for RC, else always a reverse check
		cutPrimerRev(d, PrimerIdxRev[tagIdx], revCheck);
		if (d->getRevPrimCut()) {
			//check length
			if (check_lengthXtra(d, hindrance)) {
				d->QualCtrl.minL = true; //sMinLength(pair_);
				d->failed();
				return false;
			}
		} else  {//stats, but only for 2nd pair_
			if (pair == 1 && bRequireRevPrim) {//failed to find reverse primer
				return false;
			}
		}
	}


	if (doSeeding){
		//cut off low qual, hard limits
		d->qualWinPos(EWwidth, EWthr);
		return true;
	}
	
	
	//if seq needs to be cut, than here
	if (TruncSeq>0){
		d->cutSeqPseudo(TruncSeq);
		if ( check_length(d->length()) ){
			d->QualCtrl.minL = true; //sMinLength(pair_);
			return false;
		}
	}
	if (b_doQualFilter) {
		//second cut off low qual
		d->qualWinPos(EWwidth, EWthr);
		//cut off accumulation error larger than maxAccumQP
		d->qualAccumTrim(maxAccumQP);
		//if (check_length(d->length(),hindrance) ){sMinLength();	return false;}
		if (check_length(d->length())) {
			d->QualCtrl.minLqualTrim = true; return false;//sMinQTrim(pair_);
		}
		int rea = 2;
		if ((min_q > 0 || FQWthr > 0) && d->qualWinfloat(FQWwidth, FQWthr, rea) < min_q) {
			d->QualCtrl.AvgQual = true;//sAvgQual(pair_);
			return false;
		}
		if (rea == 1) {
			d->QualCtrl.QualWin = true; //sQualWin(pair_);
			return false;
		}
		//binomial filter here
		if (b_BinFilBothPairs || pair != 1 ){
			float ExpErr = d->binomialFilter((int)BinFilErr, BinFilP);
			if (ExpErr > BinFilErr){
				d->QualCtrl.BinomialErr = true;
				//sBinomError(pair_, ExpErr);
				d->failed(); return false;
			}
		}
	}

	if (MaxAmb!=-1 && d->numACGT() > MaxAmb){
		d->QualCtrl.MaxAmb = true;//sMaxAmbig(pair_);
		return false;
	}
	if (maxHomonucleotide!=0 && !d->HomoNTRuns(maxHomonucleotide)){
		d->QualCtrl.HomoNT = true;// sHomoNT(pair_);
		return false;
	}

	//adapter removed, quality filtering done. If no base_map is provided, that is all that is needed
	if (!bDoMultiplexing){
		if (TrimStartNTs>0){
			if (d->length()-TrimStartNTs > max_l){//length check
				d->QualCtrl.maxL = true; //sMaxLength(pair_);
				return false;
			}
			//remove start NTs
			d->cutSeq(0,TrimStartNTs);

		}
		d->setPassed(true);
		return true;
	}
	d->setPassed(true);

	//keep control over passed / not as close as possible to source
	return true;
}
//DNA qual_ check, and some extra parameters
//should be safe to call from different threads
bool Filters::checkYellowAndGreen(shared_ptr<DNA> d, int pairPre, int &tagIdx) {

	unsigned int hindrance = 0;
	int pair = max(0, pairPre);//corrects for -1 (undefined pair_) to set to 0
	
	//remove technical adapter
	if (pairPre == -1 && removeAdapter)
	    remove_adapter(d);
	
	//set in outer routines that check mid (-1) or needs to be checked here(-2)
	if (tagIdx == -2){
		if ( bDoBarcode2 && pair == 1 ) {
			tagIdx = detectCutBC(d, false); //barcode 2nd part
		} else if(bDoBarcode && pair == 0) {
			tagIdx = detectCutBC(d,true); //barcode
		}
		else {
			tagIdx = 0;
		}
	}
	if ((bDoBarcode || bDoBarcode2) && tagIdx < 0) {
	    d->QualCtrl.TagFail = true;
	    return false;
	}
	//bShortAmplicons checks for reverse primer on 1st read
	if (BcutPrimer) {
		if (pair != 1  ) {//0 or -1
			cutPrimer(d, PrimerIdx[tagIdx], false, pair);
			if (bShortAmplicons) {//also check other end of primer..
				cutPrimerRev(d, PrimerIdxRev[tagIdx], true);
			}
			if (!d->getFwdPrimCut() && bRequireFwdPrim) {
				d->failed(); return false;
			}
		}
		else if (pair != 0) {//pair_ == 1, check for fwd primer in pair_ 2 (rev-compl)
			bool revCheck = pair == -1 || pair == 0;//1:false for RC, else always a reverse check
			cutPrimerRev(d, PrimerIdxRev[tagIdx], revCheck);
			if (bShortAmplicons) {//also check other end of primer..
				cutPrimer(d, PrimerIdx[tagIdx], true, pair);
			}
			if (!d->getRevPrimCut()) {
				if (pair != 0 && bRequireRevPrim) {//failed to find reverse primer
					if (alt_bRequireRevPrim) {
						d->failed(); return false;
					}
					else {
						d->setMidQual(true);
					}
				}
			}
		}
	}
	if (check_lengthXtra(d)){
		d->failed(); return false;
	}

	//if seq needs to be cut, than here
	if (TruncSeq>0){
		d->cutSeqPseudo(TruncSeq);
		if ( check_lengthXtra(d) ){
			d->failed(); return false;
		}
	}

	if (b_doQualFilter) {
		//second cut off low qual_
		d->qualWinPos(EWwidth, EWthr);	// { qualWinTrim = true; }
		//cut off accumulation error larger than maxAccumQP
		if (maxAccumQP != -1.0) {
			int cP = d->qualAccumulate(maxAccumQP);
			if (check_lengthXtra(d, 0, cP)) {
				d->isYellowQual();  d->QualCtrl.minLqualTrim = true;//sMinQTrim(pair_);
				cP = d->qualAccumulate(alt_maxAccumQP);
				if (check_lengthXtra(d, 0, cP)) {//check if passes alt
					d->failed(); return false;
				}
			} else {
				d->qualAccumTrim(maxAccumQP);//) {  AccErrTrim = true; }
			}
		}


		int rea (2), rea2 (2);
		if ((min_q > 0 || FQWthr > 0) && d->qualWinfloat(FQWwidth, FQWthr, rea) < min_q) {
			d->QualCtrl.AvgQual = true; //sAvgQual(pair_);
			if ((alt_min_q > 0 || alt_FQWthr > 0) && d->qualWinfloat(FQWwidth, alt_FQWthr, rea2) < alt_min_q) {
				d->QualCtrl.AvgQual= true; //statAddition.AvgQual++;
				d->failed(); return false;
			} else {
				d->QualCtrl.AvgQual = false;
				d->setMidQual(true);
			}
		}
		if (rea == 1) {
			d->QualCtrl.QualWin = true; //sQualWin(pair_);
			if (rea2 == 1) {
				d->failed(); return false;
			}
		}
		if (b_BinFilBothPairs || pair == 0){
			float ExpErr = d->binomialFilter((int)BinFilErr, BinFilP);
			if (ExpErr > BinFilErr){
				d->QualCtrl.BinomialErr = true;
				d->failed(); return false;
			}
		}
	}
	int ambNTs = d->numACGT();
	if (MaxAmb!=-1 && ambNTs > MaxAmb){
		d->QualCtrl.MaxAmb = true;
		if (alt_MaxAmb!=-1 && ambNTs>= alt_MaxAmb){
			d->QualCtrl.MaxAmb = true; //statAddition.MaxAmb++;
			d->failed(); return false;
		}
	}
	if (maxHomonucleotide!=0 && !d->HomoNTRuns(maxHomonucleotide)){
		d->QualCtrl.HomoNT = true;//sHomoNT(pair_);
		d->failed(); return false;
	}

	//adapter removed, quality filtering done. If no map is provided, that is all that is needed
	if (!bDoMultiplexing){
		if (TrimStartNTs>0){
			if (d->length()-TrimStartNTs > max_l){//length check
				d->QualCtrl.maxL = true; //sMaxLength(pair_);
				d->failed();
				return false;
			}
			//remove start NTs
			d->cutSeq(0,TrimStartNTs);

		}
		d->setPassed(true);
		return true;
	}

	if (!d->isYellowQual()) {
		d->setPassed(true);
	}

	return true;
}
void Filters::noMapMode(){
	string noMapTxt = "sdm run in No Map Mode.";
	if (cmdArgs.find("-paired")  != cmdArgs.end() && (cmdArgs["-paired"]=="2" || cmdArgs["-paired"]=="2")){
		pairedSeq = 2; //fakeEssentials();
		noMapTxt += " Using paired end sequencing files.";
	}

	BcutPrimer = false; bDoBarcode = false; bDoBarcode2 = false; bDoBarcode2Rd1 = false;
    removeAdapter=false;bDoMultiplexing=false;
	bDoHeadSmplID=false;
	fakeEssentials();
    minBCLength1_ = 0; minBCLength2_ = 0; maxBCLength1_ = 0; maxBCLength2_ = 0; minPrimerLength_ = 0;
	cerr<<noMapTxt<<endl;
}
void Filters::fakeEssentials(){
	//create fake entries
	PrimerIdx.push_back(0);Barcode.push_back("NA");
	barcodeLengths1_.push_back(0);
	barcodeLengths2_.push_back(0);
	SequencingRun.push_back("");
	SequencingRun2id[""] = vector<int>(1, 0);
	PrimerL.push_back(""); PrimerL_RC.push_back(""); SampleID.push_back("NA"); SampleID_Combi.push_back("NA");
	HeadSmplID.push_back("");bDoHeadSmplID=false;
	collectStatistics[0]->BarcodeDetected.push_back(-1);
	collectStatistics[1]->BarcodeDetected.push_back(-1);
	collectStatistics[0]->BarcodeDetectedFail.push_back(-1);
	collectStatistics[1]->BarcodeDetectedFail.push_back(-1);
	
}
void Filters::allResize(unsigned int x){
	//cerr<<"resize "<<x<<endl;
	PrimerIdx.resize(x,0);
	PrimerIdxRev.resize(x,0);
	Barcode.resize(x, "");
	Barcode2.resize(x, "");
	barcodeLengths1_.resize(x, 0);
	barcodeLengths2_.resize(x, 0);
	SampleID.resize(x, "");
	SampleID_Combi.resize(x, "");
	SequencingRun.resize(x, "");
	
	collectStatistics[0]->BarcodeDetected.resize(x, 0);
	collectStatistics[1]->BarcodeDetected.resize(x, 0);
	collectStatistics[0]->BarcodeDetectedFail.resize(x, 0);
	collectStatistics[1]->BarcodeDetectedFail.resize(x, 0);

	statAddition[0]->BarcodeDetected.resize(x, 0);
	statAddition[0]->BarcodeDetectedFail.resize(x, 0);
	statAddition[1]->BarcodeDetected.resize(x, 0);
	statAddition[1]->BarcodeDetectedFail.resize(x, 0);
	HeadSmplID.resize(x, "");
	vector<ofbufstream*> emptVec(2, NULL);
	vector<string> emptVec2(2, "");
}

bool Filters::remove_adapter(shared_ptr<DNA> d){ //technical adapter
	//allows for 0 errors, no shifts
	if (d->getTA_cut()) {
		return true;
	}
	const string& se = d->getSequence();
	for (unsigned int i=0;i<tAdapterLength;i++){
		if (se[i]!=tAdapter[i] ){
			return false;
		}
	}
	d->cutSeq(0,tAdapterLength);
	d->setTA_cut(true);
	return true;
}
//only identifies based on dual BCding
void Filters::dblBCeval(int& tagIdx, int& tagIdx2, string presentBC, shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2) {
	//bool BCfail = false;// , BCfail2 = false;
	if ( tagIdx < 0 || tagIdx2 < 0 || !tdn->getBarcodeDetected() || !tdn2->getBarcodeDetected()) {
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != NULL) { 
			tdn->setPassed(false); /*BCfail = true; */
			tdn->setBCnumber(tagIdx, BCoffset); tdn->setMidQual(false);
		} 
		if (tdn2 != NULL) { tdn2->setPassed(false); tdn2->setMidQual(false); tdn2->setBCnumber(tagIdx2, BCoffset);}
		
		collectStatistics[0]->dblTagFail++;
		return;
	}
	string BC1 = Barcode[tagIdx];
	string BC2 = Barcode2[tagIdx2];
	bool hit(false);
	//this routine finds two matching barcodes (as several combinations are possible)
	for ( uint i = 0; i < Barcode.size(); i++ ) {
		if ( Barcode[i] == BC1 && Barcode2[i] == BC2 ) {
			tagIdx = i; tagIdx2 = i; hit = true; break;
		}
	}

	if ( !hit ) {
		//no BC, useless
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != NULL) { tdn->setPassed(false); tdn->setMidQual(false); tdn->setBCnumber(tagIdx, BCoffset);	}
		if (tdn2 != NULL) { tdn2->setPassed(false); tdn2->setMidQual(false); tdn2->setBCnumber(tagIdx2, BCoffset);	}
		return;
	}
	presentBC = BC1 + "|" + BC2;
	//add new BC info to DNA
	//also reset BC in DNA
	if ( tdn != NULL ) {
		BCintoHead(tagIdx, tdn, presentBC, -1, false, true);
		//already done in BCintoHead
	}
	if ( tdn2 != NULL ) {
		BCintoHead(tagIdx2, tdn2, presentBC, -1, true, true);
	}
}

//cuts & identifies - version is just for mid sequences
int Filters::detectCutBC(shared_ptr<DNA> d, string&presentBC, int& c_err, bool isPair1) {
	if (bDoHeadSmplID) {
		for (unsigned int i = 0; i<HeadSmplID.size(); i++) {
			size_t pos = d->getOldId().find(HeadSmplID[i]);
			if (pos != string::npos) {
				SampleIntoHead(i, d, pos);
				return i;
			}
		}
		return -1;
	}
	/*BCdecide & locBCD(BCdFWD);
	if ( !isPair1 ) {
	locBCD = BCdREV;
	}*/
	int start(-1), stop(-1);
	int idx(-1);
	int scanRegion = 4; //dna region to scan for Tag sequence_
	if (!d->getTA_cut() && isPair1) {//no technical adapter found / given by user: scan wider region for barcode
		scanRegion = 14; //arbitary value
	}
	if (d->isMIDseq()) {
		if (d->length() < minBCLength1_) { return -1; }
		scanRegion =  d->length() - minBCLength1_ + 1;
	}

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
	if (!BCdFWDREV[!isPair1].b_BCdirFix) {
		if (start == -1) {//check reverse transcription
						  //d->reverseTranscribe();
			scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
			if (start != -1) {
				BCdFWDREV[!isPair1].BCrevhit++;
			}
		}
		else {
			BCdFWDREV[!isPair1].BChit++;
		}
		//check if BC direction can be fixed
		if (BCdFWDREV[!isPair1].BCrevhit + BCdFWDREV[!isPair1].BChit > DNA_MAX_IN_MEM) {
			if (!eval_reversingBC(isPair1)) { return -1; }
		}
	}
	if (start != -1) {
		if (BcutTag && !d->isMIDseq()) {
			//remove tag from DNA
			d->cutSeq(start, stop);
			d->setBarcodeCut();
		}
	}
	else {
		idx = -1;
	}
	return idx;
}

int Filters::findTag(shared_ptr<DNA> d, string&presentBC, int& c_err, 
			bool isPair1, int& revChecks) {
    
    //cout << "FIND TAG DO HEAD " << bDoHeadSmplID << endl;
	if (bDoHeadSmplID) {
		for (unsigned int i = 0; i<HeadSmplID.size(); i++) {
			size_t pos = d->getOldId().find(HeadSmplID[i]);
			if (pos != string::npos) {
				SampleIntoHead(i, d, pos);
				return i;
			}
		}
		return -1;
	}
	/*BCdecide & locBCD(BCdFWD);
	if ( !isPair1 ) {
	locBCD = BCdREV;
	}*/
	int start(-1), stop(-1);
	int idx(-1);
	int scanRegion = 16; //dna region to scan for Tag sequence_

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
	return -1;
	
	if (!BCdFWDREV[!isPair1].b_BCdirFix) {
		if (start == -1) {//check reverse transcription
						  //d->reverseTranscribe();
			scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
			if (start != -1) {
				BCdFWDREV[!isPair1].BCrevhit++;
			}
		}
		else {
			BCdFWDREV[!isPair1].BChit++;
		}
		//check if BC direction can be fixed
		if (BCdFWDREV[!isPair1].BCrevhit + BCdFWDREV[!isPair1].BChit > DNA_MAX_IN_MEM) {
			if (!eval_reversingBC(isPair1)) { return -1; }
		}
	}
	else if (idx < 0 && revChecks > 0) {
		scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
		if (idx >= 0) {
			revChecks = 0;
		}
	}
	if (start == -1) {
		idx = -1;
	}
	return idx;
}

int Filters::detectCutBC(shared_ptr<DNA> d, bool isPair1) {
	//seq too short for BC
	if (d->length() < minBCLength1_ ) {
	    return -1;
	}
	//already detected barcode
	if (d->getBarcodeCut()) {
		return d->getBarcodeNumber() - BCoffset;
	}
	if ((isPair1 && !bDoBarcode) || (!isPair1 && !bDoBarcode2)) {
		d->setBCnumber(0, BCoffset);
		return BCoffset; //not failed, just not requested
	}

	//ok, really start looking for BC in seq
	int idx(-1);
	if (bDoHeadSmplID){

	    // must be locked for multithreading
        unsigned int i = 0;
        for (; i < HeadSmplID.size(); i++) {
            size_t pos = d->getOldId().find(HeadSmplID[i]);
            if (pos != string::npos) {
                if (!bDoBarcode2) {
                    SampleIntoHead(i, d, pos);
                }//this has to be done AFTER two BCs are read (on a higher lvl)
                else {
                    d->setBCnumber(i, BCoffset);
                }
                idx = i;
                break;
            }
        }
        return idx;
	} else if (bOneFileSample) {
        // if theres only one sample per file then there is no need for demultiplexing
		d->setBCnumber(0,this->currentBCnumber());

		// needs to be locked when multithreading
		BCintoHead(0, d, "FileName", isPair1, false);

        return 0;
	}
	
	int start(-1);
	int stop(-1);
	string presentBC;
	int c_err(0);
	int scanRegion=4; //dna region to scan for Tag sequence_
	
	if (!d->getTA_cut() && isPair1){//no technical adapter found / given by user: scan wider region for barcode
		scanRegion = 22; //arbitary value
	}
	if (d->isMIDseq() || d->length() < scanRegion){
		scanRegion = d->length() - minBCLength1_ + 1;
	}

	// needs to be locked when multithreaded
	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);

    
	if ( !BCdFWDREV[!isPair1].b_BCdirFix ) {
		if (start == -1){//check reverse transcription
			//d->reverseTranscribe();
			scanBC_rev(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);
            //cout << "scanBCrev: " << start << "," << stop << "," << idx << "," << c_err << "," << presentBC << endl;
			if (start!=-1){
				BCdFWDREV[!isPair1].BCrevhit++;
			}
		} else {
			BCdFWDREV[!isPair1].BChit++;
		}
		//check if BC direction can be fixed
		if ( BCdFWDREV[!isPair1].BCrevhit + BCdFWDREV[!isPair1].BChit > 5000 ) {
			eval_reversingBC(isPair1);//){return -1;}
		}
	}
	if (start < 0) {
		return (-1);
	}

	//check also for reverse BC on same read??? (PacBio)
	int startX(-1), stopX(-1);
	string presentBCX(""); int c_errX(0);
	int scanRegionX = 4;
	int idxX = -2;
	if (bDoBarcode2Rd1) {
		scanBC_rev(d, startX, stopX, idxX, c_errX, scanRegionX, presentBCX, isPair1);
        //cout << "scanBCrev: " << start << "," << stop << "," << idx << "," << c_err << "," << presentBC << endl;
	}
	

	d->setBCnumber(idx, BCoffset);

	if (BcutTag && !d->isMIDseq()) {
		//remove tag from DNA
		d->cutSeq(start, stop);
		d->setBarcodeCut();

		// needs to be locked when multithreaded
		BCintoHead(idx, d, presentBC, c_err, isPair1);
        
    }
	return idx;
}

void Filters::BCintoHead(int idx, shared_ptr<DNA> d,const string presentBC,
						 const int c_err, bool pair1, bool atEnd){
	vector<string> * locBC = (!pair1) ? &Barcode2 : &Barcode;
	string on (d->getId());
	string spacee (" ");
	//keep only until first space
	//remove">"
	if (!atEnd){
		on = on.substr(0,on.find_first_of(' ',0));
	}
	string nID ( SampleID[idx]);
	if (bDoCombiSamples){
		nID = SampleID_Combi[idx];
	}
	nID += iniSpacer +
		on + spacee + string("orig_bc=") + (*locBC)[idx];
	if (c_err > 0 || atEnd) {//atEnd: dbl barcode
		//convert c_err to s_c_err;
		//string s_c_err; stringstream conve;
		//conve << c_err;
		//s_c_err = conve.str();
		nID += spacee + string("new_bc=") + presentBC +
			spacee + string("bc_diffs=") + itos(c_err);
	}
	d->setNewID(nID);
	d->setBCnumber(idx, BCoffset);
}

void Filters::SampleIntoHead(const int idx, shared_ptr<DNA> d, const size_t pos){
	string on = d->getId(),spacee=" ";
	size_t pos2 = on.find_first_of(" ",pos);
	string on2 = on.substr(0,pos)+on.substr(pos2+1);
	string on3 = on.substr(pos,pos2);

	string nID(SampleID[idx]);
	if (bDoCombiSamples){
		nID = SampleID_Combi[idx];
	}

	nID +=  iniSpacer +
		on2 + spacee + string("orig_hdPart=")+on3;
	d->setNewID(nID);
	d->setBCnumber(idx, BCoffset);
}
void Filters::countBCdetected(int BC, int Pair, bool MidQ) {
    
	if (!bDoMultiplexing) { return; }
	if (BC - BCoffset < 0) {
		exit(132);
	}
	if (Pair < 0) { Pair = 0; }
	if (!MidQ) {
		collectStatistics[Pair]->BarcodeDetected[BC - BCoffset]++;
	}
	else {
		statAddition[Pair]->BarcodeDetected[BC - BCoffset]++;
	}
}
bool Filters::eval_reversingBC(bool fwd){
	if ( !fwd && !bDoBarcode2 ) {
		return true;
	}
	/*BCdecide lbcd(BCdFWD);
	if ( !fwd ) {
		lbcd = BCdREV;
	}*/
	if ( BCdFWDREV[!fwd].b_BCdirFix ) { return true; }
	BCdFWDREV[!fwd].b_BCdirFix = true; lMD->setBCfixed(true, fwd);
	if ( BCdFWDREV[!fwd].BCrevhit> BCdFWDREV[!fwd].BChit * 8 ) {//use reversed_ BC ..
		BCdFWDREV[!fwd].reversedBCs = true;
		if ( fwd ) {
			reverseTS_all_BC();
		} else {
			reverseTS_all_BC2();
		}
		if ( BCdFWDREV[!fwd].BChit > 0 ) {
			restartSet=true;
			return false;
		}
	} else if ( BCdFWDREV[!fwd].BCrevhit>0 ) {
		restartSet=true;
		return false;
	}
	return true;
}

void Filters::scanBC_rev(shared_ptr<DNA> d,int& start,int& stop,int& idx,int c_err, 
					 int scanRegion,string & presentBC,
					 bool fwdStrand) {
	vector<string> emptyV(0), emptyV2(0);

	vector<string>& localBarcodesRev(emptyV2);
	vector<string>&  locBC(emptyV);
	if ( !fwdStrand ) {
        localBarcodesRev = revBarcode2;
		locBC = Barcode2;
	} else {
        localBarcodesRev = revBarcode;
		locBC = Barcode;
	}
	int BCs = (int) localBarcodesRev.size();
	
	if (BCs==0){
        localBarcodesRev = locBC;
		BCs = (int)localBarcodesRev.size();
		for (int i=0; i< BCs; i++){
			reverseTS(localBarcodesRev[i]);
		}
		if ( !fwdStrand ) {//copy over BC
			revBarcode2 = localBarcodesRev;
		} else {
			revBarcode = localBarcodesRev;
		}
	}
	//check each possible BC for a match
	if (barcodeErrors_ == 0){
		for (idx=0; idx< BCs; idx++){
			start = d->matchSeq_tot(localBarcodesRev[idx], 0, scanRegion, c_err);
			if (start!=-1){
				presentBC = localBarcodesRev[idx];
				stop = start+ (int)localBarcodesRev[idx].length();
				break;
			}
		}
	} else {
		vector<int> stars(0),idxses(0);
		bool zeroErr = false;
		//this version tries all BC's and if there are more than one possible match, will reject all matches
		for (idx=0; idx< BCs; idx++){
			start = d->matchSeq_tot(localBarcodesRev[idx], barcodeErrors_, scanRegion, c_err);
			if (start!=-1){
				if (c_err==0){
					stop = start+ (int) localBarcodesRev[idx].length();
					presentBC = localBarcodesRev[idx];
					zeroErr = true;
					break;
				}
				stars.push_back(start);
				idxses.push_back(idx);
			}
		}
		if (!zeroErr && stars.size()>0){
			//int pair_ = (int)!fwdStrand;//d->getReadMatePos();
			if (stars.size() > 1){//too many matches, thus true seq can't be found
				//currently only have only one BC, could be changed in future
				//sTagNotCorrected(pair_);
				d->QualCtrl.fail_correct_BC = true;
				idx=-1; start = -1;
				return;
			}
			d->QualCtrl.suc_correct_BC = true;
			//sTagCorrected(pair_);// collectStatistics.suc_correct_BC++;
			start = stars[0];
			idx = idxses[0];
			stop = start+(int)localBarcodesRev[idx].length();
			presentBC = d->getSubSeq(start,stop);
		}

	}
}

void Filters::scanForBarcode(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err, int scanRegion, string& barcode, bool fwdStrand) {
//    BarcodeMap* barcodes;
//    if (fwdStrand) {
//        barcodes = &barcodesFwdStrand_;
//    } else {
//        barcodes = &barcodesRevStrand_;
//    }
}

void Filters::scanBC(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err,
	int scanRegion, string& presentBC, bool fwdStrand) {
	if (d->length() < minBCLength1_) { return; }
    bool leaveFunction = false;


//    cout << "scanBC variables: " << endl;
//    cout << "start: " << start << endl;
//    cout << "stop: " << stop << endl;
//    cout << "idx: " << idx << endl;
//    cout << "scanRegion: " << scanRegion << endl;
//    cout << "presentBC: " << presentBC << endl;
        
        //check each possible BC for a match
        //TODO: check for BC using suffix tree
        //vector<string> emptyV(0);
        
        //BarcodeMap &localBarcodes(emptyBarcodes);
        //BarcodeMap& localBarcodes(emptyBarcodes);
        //BarcodeMap* localBarcodes = nullptr;
        BarcodeMap* localBarcodes = &emptyBarcodes;
        vector<int> *localBarcodeLengths;


        // was locked with an omp pragma
        unsigned int maxBCLength;
        unsigned int minBCLength;
        if (fwdStrand) {
            localBarcodes = &barcodes1_;
            localBarcodeLengths = &barcodeLengths1_;
            maxBCLength = maxBCLength1_;
            minBCLength = minBCLength1_;
        } else {
            localBarcodes = &barcodes2_;
            localBarcodeLengths = &barcodeLengths2_;
            maxBCLength = maxBCLength2_;
            minBCLength = minBCLength2_;
        }
    {
//    cout << "minBCLength: " << minBCLength << endl;
//    cout << "maxBCLength: " << maxBCLength << endl;
//    cout << "barcodeErrs: " << barcodeErrors_ << endl;
//
        //this for loop needs to be threadsafe



        
        //bool found = false;

        // was locked with pragma for multithreading
        for (start = 0; start < scanRegion; start++) {
            const string test = d->getSubSeq(start, maxBCLength);
            //cout << "start: " << start << " maxBCLength: " << maxBCLength << " test: " << test << endl;
            //cout << test << endl;
            auto barcodeIterator = localBarcodes->find(test);

            // found barcode
            if (barcodeIterator != localBarcodes->end()) {
                // set index if found
                idx = (*barcodeIterator).second;
                stop = start + (int) (*localBarcodeLengths)[idx];
                presentBC = test; // (*locBC)[idx];
                //found = true;
                //break;
                return;
            }
        }

        
        
        start = -1;
        if (barcodeErrors_ != 0 || minBCLength != maxBCLength) {
            vector<int> stars(0), indices(0);
            bool zeroErr = false;
            
            //this version tries all BC's and if there are more than one possible match, will reject all matches
            for (auto jx = localBarcodes->begin(); jx != localBarcodes->end(); jx++) {
                start = d->matchSeq_tot((*jx).first, barcodeErrors_, scanRegion, c_err);
                if (start != -1) {
                    idx = (*jx).second;
                    if (c_err == 0) {
                        stop = start + (int) (*localBarcodeLengths)[idx];
                        presentBC = d->getSubSeq(start, maxBCLength); // (*locBC)[idx];
                        zeroErr = true;
                        break;
                    }
                    stars.push_back(start);
                    indices.push_back(idx);
                }
            }
            if (!zeroErr && stars.size() > 0) {
                if (stars.size() > 1) {//too many matches, thus true seq can't be found
                    //sTagNotCorrected(pair_);
                    d->QualCtrl.fail_correct_BC = true;
                    idx = -1;
                    start = -1;
                    return;
                    //leaveFunction = true;
                }
                if (!leaveFunction) {
                    d->QualCtrl.suc_correct_BC = true;
                    //sTagCorrected(pair_);// collectStatistics.suc_correct_BC++;
                    
                    start = stars[0];
                    idx = indices[0];
                    stop = start + (int) (*localBarcodeLengths)[idx];
                    presentBC = d->getSubSeq(start, stop);
                }
            }
        
            
        }
    }
    
    return;
}
//cuts primers, tags
bool Filters::cutPrimer(shared_ptr<DNA> d,int primerID,bool RC,int pair){
	//only adapted to singular BC
	if (PrimerL.size()==0 || PrimerL[0].length()==0){return true;}
	if (d->getFwdPrimCut()) {
		return true;
	}
	int start(-1) ,stop(-1);
	int tolerance(30), startSearch(0);
	if (!d->getBarcodeCut() && maxBCLength1_ > 0) { tolerance = maxBCLength1_ + 4;
	} else { tolerance = 22; }//in this case nothing is known about 5' end

	if (!BcutTag){
		//Tag was not cut out of sequence_, take this into account
		startSearch = minBCLength1_ - 2;
		tolerance += (maxBCLength1_ - minBCLength1_) + 4;
	}
	if (!RC) {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, tolerance, startSearch);
		stop = start + (int)PrimerL[primerID].length();
	} else {
		int QS = d->length();int limit = max(QS >> 1, QS - 150); stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
	}
	if (start == -1){//failed to match primer
		d->QualCtrl.PrimerFail = true;
		//sPrimerFail(pair_);// max(0, (int)d->getReadMatePos()));
		if (alt_PrimerErrs!= 0 && PrimerErrs < alt_PrimerErrs){
			if (!RC) {
				start = d->matchSeq(PrimerL[primerID], alt_PrimerErrs, tolerance, startSearch);
				stop = start + (int)PrimerL[primerID].length();
			} else {
				int QS = d->length(); int limit = max(QS >> 1, QS - 150); stop = QS;
				start = d->matchSeqRev(PrimerL_RC[primerID], alt_PrimerErrs, limit, startSearch);
			}
		}
		if (start== -1){
			//statAddition.PrimerFail++;
			return false;
		}else if (pair!=1){ //2nd read shouldnt be affected by fwd primer (but still checked in short read mode)
			d->setMidQual(true);
		}
	}
	//remove everything before/after primer cut
	if (!BcutPrimer){
		if ( !RC ) { d->cutSeq(0, start); }
		else { d->cutSeq(stop,-1); }
		return true;
	}
	
	//remove the primer, if confimed before
	if ( !RC ) { d->cutSeq(0, stop); }
	else { d->cutSeq(start, -1); }
	if (!RC) {
		d->setFwdPrimCut();
	}
	return true;
}
bool Filters::findPrimer(shared_ptr<DNA> d, int primerID, bool RC, int pair) {
	//only adapted to singular BC
	if (PrimerL[0].length() == 0) { return true; }
	int start(-1);// , 
	//int stop(-1);
	int tolerance(22), startSearch(0);
	if (!d->getBarcodeCut() && maxBCLength1_ > 0) {
		tolerance = maxBCLength1_ + 4;
		startSearch = minBCLength1_ - 4;
	}
	else { tolerance = 16; }//in this case nothing is known about 5' end
	if (!RC) {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, tolerance, startSearch);
		//stop = start + (int)PrimerL[primerID].length();
	}
	else {
		int QS = d->length(); int limit = max(QS >> 1, QS - 150); //stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
	}
	if (start == -1) {//failed to match primer
		return false;
	}
	return true;
}
bool Filters::cutPrimerRev(shared_ptr<DNA> d,int primerID,bool RC){
	//const string& se = d->getSequence();
	if (d->QualCtrl.PrimerRevFail) {
		return false;
	}
	if (d->getRevPrimCut()){
		return true;
	}

	int start(-1) ,stop(d->length());
	int QS = d->length(); 
	int limit=max(QS>>1,QS-150);
	//int limit = QS>>1;

	if (!RC) {
		start = d->matchSeq(PrimerR[primerID] , PrimerErrs, 15,0);
		stop = start + (int)PrimerR[primerID].length();
	} else {
		start = d->matchSeqRev(PrimerR_RC[primerID] , PrimerErrs, limit);
	}
	


	if (start == -1){//failed to match primer
		d->QualCtrl.PrimerRevFail = true;
		return false;
	} 

	if ( !BcutPrimer ) { //found it, but no cut
		if ( !RC ) { d->cutSeq(0, start); }
		else { d->cutSeq(stop,-1); }
		return true; 
	}

	//remove the primer, if confimed before
	if ( !RC ) {
		d->cutSeq(0, stop);//start  everything in front has to be removed
	} else {
		d->cutSeq(start, -1); // everything in the end has to be removed
	}
	//string neSe = se.substr(0,start) + se.substr(stop);
	if (!RC) {
		d->setRevPrimCut();
	}

	return true;
}
bool Filters::readMap(){//core routine to read map info
	if (cmdArgs.find("-map")  == cmdArgs.end()){
		this->noMapMode(  );
		return true;
	}
	string MapF = cmdArgs["-map"];

	string path = ""; bool pathMode = false;
	if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() >= 1){
		path=cmdArgs["-i_path"] + string("/");
		pathMode = true;//check later if mapping file contains fasta/fastq
	}
    
    minBCLength1_ = 100000; minBCLength2_ = 1000000; maxBCLength1_ = 0; maxBCLength2_ = 0; minPrimerLength_ = 100000;
	string line;
	ifstream in(MapF.c_str());
	if (!in){
		cerr<<"Could not find "<<MapF<<" mapping file. Exiting.\n";exit(2);
	}
	int ini_ColPerRow(0),cnt(0),skips(0);

	//check MAP format
	//while(getline(in,line,'\n')) {
	while(!safeGetline(in,line).eof()) {
		if(line.substr(0,1) == "#" || line.length()==0){skips++;continue;}
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss,segments,'\t')) {
			ColsPerRow++;
		}
		if (segments == "") { ColsPerRow++; }

		if (cnt==0){
			ini_ColPerRow = ColsPerRow;
		} else {
			if (ColsPerRow != ini_ColPerRow){
				cerr<<"Number of columns on line "<<cnt+skips<<" is "<<ColsPerRow<<". Expected "<<ini_ColPerRow<<" columns.\n";
				return false;
			}
		}
		cnt++;
	}
	if (ini_ColPerRow==0){
		cerr<<"Mapping File exists, but appears to be badly formated (0 columns detected). Exiting\n";exit(2);
	}
	if (cnt==0){
		cerr<<"Mapping File exists, but appears to be badly formated (0 lines detected). Exiting\n";exit(2);
	}
	

	PrimerIdx.resize(cnt,-1); PrimerIdxRev.resize(cnt,-1);
	Barcode.resize(cnt, "");
	Barcode2.resize(cnt, "");
	SampleID.resize(cnt, "");
	SampleID_Combi.resize(cnt, "");
	barcodeLengths1_.resize(cnt, 0);
	barcodeLengths2_.resize(cnt, 0);
	SequencingRun.resize(cnt, "");

	collectStatistics[0]->BarcodeDetected.resize(cnt, 0);
	collectStatistics[1]->BarcodeDetected.resize(cnt, 0);
	collectStatistics[0]->BarcodeDetectedFail.resize(cnt, 0);
	collectStatistics[1]->BarcodeDetectedFail.resize(cnt, 0);
	HeadSmplID.resize(cnt, "");
	statAddition[0]->BarcodeDetected.resize(cnt, 0);
	statAddition[0]->BarcodeDetectedFail.resize(cnt, 0);
	statAddition[1]->BarcodeDetected.resize(cnt, 0);
	statAddition[1]->BarcodeDetectedFail.resize(cnt, 0);
	hetPrimer[0].resize(cnt, ""); hetPrimer[1].resize(cnt, "");
	in.clear();
	in.seekg(0, ios::beg);   

	//extract MAP content
	cnt=0;
	bool hasQualityColumn=false;
	vector<string> terms(16);terms[0]="SampleID";
	terms[1] = "BarcodeSequence"; terms[2]="LinkerPrimerSequence";
	terms[3] = "ReversePrimer"; terms[4] = "fastqFile";
	terms[5] = "fnaFile"; terms[6] = "qualFile";
	terms[7] = "SampleIDinHead"; terms[8] = "MIDfqFile";
	terms[9] = "CombineSamples"; terms[10] = "ForwardPrimer";
	terms[11] = "Barcode2ndPair";
	terms[12] = "HetSpacerFwd";	terms[13] = "HetSpacerRev";
	terms[14] = "derepMin"; terms[15] = "SequencingRun";
	bool hetOneSide = false;

	vector<int> termIdx(terms.size(),-1);

	while(!safeGetline(in,line).eof()) {
//	while(getline(in,line,'\n')) {
		if(cnt!=0 && line.substr(0,1) == "#"){continue;}
		if (line.length()<10){continue;}
		if (cnt==0){
			line = line.substr(1);
		}

		string segments;
		stringstream ss;
		ss << line;
		int tbcnt=0;

		//cmdArgs["-i_MID_fastq"]
		while (getline(ss,segments,'\t')) {
			trim(segments);
			if (cnt==0){ //search for header
				//Primer, BarcodeSequence, LinkerPrimerSequence
				//PrLCol(-1), PrRCol(-1), BCCol(-1), SIDCol(-1);
				for (unsigned int i=0; i<terms.size();i++){
					if (segments == terms[i]){
						if (i==6){hasQualityColumn=true;}
						if (i == 12){ if (hetOneSide){ doHetPrimerExplicit = true; } else{ hetOneSide = true; } }
						if (i == 13){ if (hetOneSide){ doHetPrimerExplicit = true; } else{ hetOneSide = true; } }
						if (termIdx[i] != -1){
							cerr<<"Header contains ambiguous entries: "<<segments<<"\n";
							exit(9);
						}
						termIdx[i] = tbcnt;
					}
				}
			} else { //read data into entries
				for (uint k=0; k<terms.size();k++){
					if (termIdx[k]==tbcnt){
						extractMap(k, cnt - 1, tbcnt, segments, hasQualityColumn);
					}
				}
			}
			tbcnt++;
		}
		cnt++;
	}
	//check some prerequisites 
	if (HeadSmplID.size() ==0 && Barcode.size() ==0){
		cerr<<"Could not find in mapping file either (1) valid Barcodes or (2) valid SampleID's in Sequence header. Exiting\n"; 
		exit(2);
	}
	if (pathMode) {
		if (termIdx[5] == -1 && termIdx[4] == -1) {
			cerr << "The defined input through a directory (-i_path) requries either \"fnaFile\" or \"fastqFile\" columns in the mapping file. \nAborting..\n";
			exit(55);
		}
	}

	decideHeadBC();
	
	//check for duplicate barcodes, but only if no list provided
	if(termIdx[4]!=-1 &&termIdx[5]!=-1 &&termIdx[6]!=-1){
		checkDoubleBarcode();
	}
	if (termIdx[9] != -1){ bDoCombiSamples = true; }
	checDoubleSampleID();

	//eleminate required checks for primers, if there is simply no primer given
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0 || PrimerL[0].substr(0, 1) == " ") {
		BcutPrimer = false;
		bRequireFwdPrim = false;
		alt_bRequireFwdPrim = false;
	}
	if (PrimerR.size() == 0 || PrimerR[0].length() == 0 || PrimerR[0].substr(0, 1) == " ") {
		alt_bRequireRevPrim = false;
		bRequireRevPrim = false;
	}
	if (pathMode) {
		this->sanityCheckFilesSRs();
	} else {
		this->removeSRs();
	}

	this->BarcodePreStats();

	return true;
}
void Filters::decideHeadBC(){
	bDoMultiplexing=true;
	if (HeadSmplID[0].length()>0 && Barcode[0].length()==0){
		bDoBarcode = true; bDoHeadSmplID = true; minPrimerLength_ = 0; return;
	} else if ( HeadSmplID[0].length() == 0 && Barcode[0].length() > 0 && Barcode2[0].length() > 0 ) {
		bDoBarcode = true; bDoHeadSmplID = false; bDoBarcode2 = true;  return;
	} else if ( HeadSmplID[0].length() == 0 && Barcode[0].length() > 0 ) {
		bDoBarcode = true; bDoHeadSmplID = false; bDoBarcode2 = false;  return;
	} else if ( HeadSmplID[0].length() == 0 && Barcode[0].length() == 0 ) {
		//simply check if each filename is different
		if (FastqF.size() > 0 ){
			for (uint i = 0; i<FastqF.size(); i++){
				for (uint j = i+1; j<FastqF.size(); j++){
					if (FastqF[i] == FastqF[j]){
						cerr << "File names " << i << " and " << j << " are equal - no identification by filename possible.\n   Aborting..\n"; exit(55);
					}
				}
			}
			bOneFileSample = true; bDoBarcode = true; bDoHeadSmplID = false;
			return;
		}
		if (FastaF.size() > 0){
			for (uint i = 0; i < FastaF.size(); i++){
				for (uint j = i + 1; j < FastaF.size(); j++){
					if (FastaF[i] == FastaF[j]){
						cerr << "File names " << i << " and " << j << " are equal - no identification by filename possible.\n   Aborting..\n"; exit(55);
					}
				}
			}

			bOneFileSample = true; bDoBarcode = true; bDoHeadSmplID = false;
			return;
		}
	}

	bDoMultiplexing = false;
	cerr << "No Barcode and no id_ in header defined.. aborting\n";
	exit(53);

}


void Filters::checkDoubleBarcode(){
	if (!bDoBarcode || bDoHeadSmplID || Barcode.size()==0){ return; }
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	if ( bDoBarcode2 ) {
		if ( Barcode.size() != Barcode2.size() ) {
			cerr << "Unequal Barcode vector sizes in dual barcoding controls. Exiting.."; exit(45);
		}
		for ( unsigned int i = 0; i < Barcode.size(); i++ ) {
			for ( unsigned int j = i + 1; j < Barcode.size(); j++ ) {
				if ( strcmp(Barcode[i].c_str(), Barcode[j].c_str()) == 0 && strcmp(Barcode2[i].c_str(), Barcode2[j].c_str()) == 0 ) {
					empty[0] = i; empty[1] = j ; doubles.push_back(empty);
				}
			}
		}
		if (doubles.size() > 0){
			for (uint x = 0; x < doubles.size(); x++){
				int i = doubles[x][0]; int j = doubles[x][1];
				cerr << "Duplicate dual Barcode detected: Barcode1 " << i + 1 << " (" << Barcode[i] << ") and " << j + 1 << " (" << Barcode[j] << ")  as well as Barcode1 " << i + 1 << " (" << Barcode2[i] << ") and " << j + 1 << " (" << Barcode2[j] << ") are equal.\n";
			}
			exit(8);
		}
	}
	else {
		for ( unsigned int i = 0; i < Barcode.size(); i++ ) {
			for ( unsigned int j = i + 1; j < Barcode.size(); j++ ) {
				if ( strcmp(Barcode[i].c_str(), Barcode[j].c_str()) == 0 ) {
					empty[0] = i; empty[1] = j; doubles.push_back(empty);
				}
			}
		}
		if (doubles.size() > 0){
			for (uint x = 0; x < doubles.size(); x++){
				int i = doubles[x][0]; int j = doubles[x][1];
				cerr << "Duplicate Barcode detected: Barcode " << i + 1 << " (" << Barcode[i] << ") and " << j + 1 << " (" << Barcode[j] << ") are equal.\n";
			}
			exit(8);
		}

	}
}
void Filters::checkDoubleSampleIDHead(){
	if (!bDoHeadSmplID){return;}
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	for (unsigned int i=0; i<HeadSmplID.size();i++){
		for (unsigned int j=i+1;j<HeadSmplID.size();j++){
			if(strcmp(HeadSmplID[i].c_str(),HeadSmplID[j].c_str()) == 0){
				empty[0] = i; empty[1] = j; doubles.push_back(empty);
			}
		}
	}
	if (doubles.size() > 0){
		for (uint x = 0; x < doubles.size(); x++){
			int i = doubles[x][0]; int j = doubles[x][1];
			cerr << "Duplicate Header2split detected: pattern " << i + 1 << " and " << j + 1 << " are equal.\n";
		}
		exit(8);
	}

}

void Filters::checDoubleSampleID(){
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	for (unsigned int i = 0; i<SampleID.size(); i++){
		for (unsigned int j=i+1;j<SampleID.size();j++){
			if(strcmp(SampleID[i].c_str(),SampleID[j].c_str()) == 0){
				empty[0] = i; empty[1] = j; doubles.push_back(empty);
			}
		}
	}
	if (doubles.size() > 0){
		for (uint x = 0; x < doubles.size(); x++){
			int i = doubles[x][0]; int j = doubles[x][1];
			cerr << "Duplicate SampleID detected: SampleID " << i + 1 << " and " << j + 1 << " are equal.\n";
		}
		exit(8);
	}
	if (!bDoCombiSamples){
		return;
	}
	bDoCombiSamples = false;
	if (SampleID_Combi.size() <= 1){
		return;
	}
	
	string prevCSID = SampleID_Combi[0];
	for (unsigned int i = 1; i < SampleID_Combi.size(); i++){
		if (SampleID_Combi[i] != prevCSID){
			bDoCombiSamples = true; 
		}
		if (SampleID_Combi[i] == ""){
			SampleID_Combi[i] = SampleID[i];
		}
	}
}

void Filters::removeSRs(void) {
	SequencingRun2id.clear();
	SequencingRun2id[""] = vector(0,0);
	auto XX = SequencingRun2id.find("");
	for (int k = 0; k < (int)SequencingRun.size(); k++) {
		SequencingRun[k] = "";
		XX->second.push_back(k);
		
	}
}
void Filters::sanityCheckFilesSRs(void) {
	//just checks that SampleRun and input files+BCs do not collide
	//so same SampleRunshould be in the same file(s)
	unordered_map<string, int> files;
	int curBlk = -1;
	for (auto YY : SequencingRun2id) {
		vector<int> ids = YY.second;
		curBlk++;
		for (auto J : ids) {
			string FQ("");
			if (FastqF.size() == 0) {
				FQ = FastaF[J];
			} else {
				FQ = FastqF[J];
			}
			auto XX = files.find(FQ);
			if (XX == files.end()) {
				files[FQ] = curBlk;
			} else if (XX->second != curBlk) {
				cerr << "Sample " << SampleID[J] << "("<< FQ <<") has SequencingRun " << SequencingRun[J] << " already found for other files. This is not allowed.\n";
				exit(323);
			}
		}
	}
}

void Filters::BarcodePreStats(){
    minBCLength1_=100000;maxBCLength1_=0;
	for (unsigned int i=0; i<Barcode.size();i++ ){
		if (Barcode[i].length() < minBCLength1_) minBCLength1_= (unsigned int) Barcode[i].length();
		if (Barcode[i].length() > maxBCLength1_) maxBCLength1_= (unsigned int) Barcode[i].length();
		//initialize Barcodes
		barcodes1_[Barcode[i]] = i;
        barcodeLengths1_[i] = (int) Barcode[i].length();
	}
	if (minBCLength1_ == maxBCLength1_){
		bBarcodeSameSize = true;
	}
    minBCLength2_ = 100000; maxBCLength2_ = 0;
	for ( unsigned int i = 0; i < Barcode2.size(); i++ ) {
		//create index
		barcodes2_[Barcode2[i]] = i;
		if (Barcode2[i].length() < minBCLength2_ ) minBCLength2_ = (unsigned int)Barcode2[i].length();
		if (Barcode2[i].length() > maxBCLength2_ ) maxBCLength2_ = (unsigned int)Barcode2[i].length();
        barcodeLengths2_[i] = (int)Barcode2[i].length();
	}
	//fix empty last column specifically for derepMinNum
	if (derepMinNum.size() > 0) {
		derepMinNum.resize(Barcode.size(), -1);
	}

}
void Filters::resetStats(){
	for (size_t i = 0; i < 2; i++) {
		 collectStatistics[i]->reset();
		 statAddition[i].reset();
	}

}

void Filters::failedStats2(shared_ptr<DNA> d,int pair){
	int pa = max(pair, 0);
	if (bDoMultiplexing){
		int idx = d->getBarcodeNumber() - BCoffset;
		if ( bOneFileSample ) {
			collectStatistics[pa]->BarcodeDetectedFail[0]++;
		} else if (idx >= 0) {

			
#ifdef DEBUG
			if (pa < 0 || pa>1) { cerr << "Pair in failedStats2 set to:" << pa << endl; }
			if (idx >= (int)collectStatistics[pa]->BarcodeDetectedFail.size()) {
				cerr << "idx in failedStats2 too big:" << idx << endl;
			}
#endif // DEBUG
				collectStatistics[pa]->BarcodeDetectedFail[idx]++;
		}
	}

}
void Filters::prepStats() {
	float remSeqs = float(collectStatistics[0]->total - collectStatistics[0]->totalRejected);
	for (size_t i = 0; i < 2; i++) {
		collectStatistics[i]->PostFilt->calcSummaryStats(remSeqs, min_l, min_q);
		if (bAdditionalOutput) {
			remSeqs = float(statAddition[0]->total - statAddition[0]->totalRejected);
			statAddition[0]->PostFilt->calcSummaryStats(remSeqs, min_l, min_q);//yellow

		}
		collectStatistics[i]->PreFilt->calcSummaryStats(1, min_l, min_q);
	}
	
}


void Filters::addPrimerL(string segments, int cnt){
	int used = -1;
	for (unsigned int i=0; i<PrimerL.size();i++){
		if (segments==PrimerL[i]){used=i;}
	}
	if (used == -1){
		PrimerL.push_back(segments);
		PrimerL_RC.push_back(reverseTS2(segments));
		PrimerIdx[cnt] = (int)PrimerL.size() - 1;
		if (segments.length() < minPrimerLength_) minPrimerLength_= (unsigned int) segments.length();
	} else {
		PrimerIdx[cnt] = used;
	}
}
void Filters::addPrimerR(string segments, int cnt){
	bPrimerR=true;
	int used = -1;
	for (unsigned int i=0; i<PrimerR.size();i++){
		if (segments==PrimerR[i]){used=i;}
	}
	if (used == -1){
		PrimerR.push_back(segments);
		PrimerR_RC.push_back(reverseTS2(segments));
		PrimerIdxRev[cnt] = (int)PrimerR.size() - 1;
	} else {
		PrimerIdxRev[cnt] = used;
	}
}


void Filters::extractMap(int k, int cnt, int tbcnt, string & segments,
	bool hasQualityColumn){
	//terms[4] = "fastqFile";terms[5] = "fnaFile"; terms[6] = "qualFile";
	switch(k)
	{
	case 4: // fastq file
		trim(segments);
		FastqF.push_back(segments);
		break;
	case 5: // fna file
		trim(segments);
		FastaF.push_back(segments);
		if (!hasQualityColumn  ){//code to create artificial quality file name
				string newQ = segments;
				size_t pos = newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual_");
				QualF.push_back(newQ);
		}
		break;
	case 6: // qual_ file
		trim(segments);
		QualF.push_back(segments);
		break;
	case 2: //left primer 		
	case 10:
		trim(segments);
		transform(segments.begin(), segments.end(),segments.begin(), ::toupper);
		this->addPrimerL(segments,cnt);
		break;
	case 3: // right primer
		trim(segments);
		transform(segments.begin(), segments.end(),segments.begin(), ::toupper);
		this->addPrimerR(segments,cnt);
		break;
	
	case 0: //id_
		trim(segments);
		SampleID[cnt] = segments;
		break;

	case 1: //Barcode
		trim(segments);
		transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
		Barcode[cnt] = segments;
		break;
	case 11: //Barcode rev
		trim(segments);
		transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
		Barcode2[cnt] = segments;
		break;
	case 7: //sample id_ in head
		trim(segments);
		HeadSmplID[cnt] = segments;
		break;
	case 8://mid xtra fq
		trim(segments);
		MIDfqF.push_back(segments);
		break;
	case 9://combine samples
		trim(segments);
		SampleID_Combi[cnt] = segments;
		break;
	case 14://demultiplex num
		if (segments.length() == 0) {
			derepMinNum.push_back(-1);
		}else if (!is_digits(segments)){
			cerr << "Wrong map entry \"" << segments << "\". For header derepMin only number can be used.\n"; exit(313);
		} else {
			int nint = atoi(segments.c_str());
			derepMinNum.push_back( nint);
		}
		break;
	case 12://het primer fw
		if (!doHetPrimerExplicit){break;}
		hetPrimer[0][cnt] = segments;
		break;
	case 13://het primer rv
		if (!doHetPrimerExplicit){ break; }
		hetPrimer[1][cnt] = segments;
		break;
	case 15: //SequencingRun
		trim(segments);
		SequencingRun[cnt] = segments;
		auto XX = SequencingRun2id.find(segments);
		if (XX != SequencingRun2id.end()) {
			XX->second.push_back(cnt);
		} else {
			SequencingRun2id[segments] = vector<int>(1, cnt);
		}
		break;
	}

	if (k==6 || k==5){//a qual_ pushback was "" (was empty); takeOver
		if (QualF.size() == FastaF.size() && QualF.back()==""){
				string newQ = FastaF.back();
				size_t pos =  newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual_");
				QualF.back() = newQ;
		}
	}

}

void Filters::printHisto(ostream& give,int which, int set){
	bool p2stat = pairedSeq > 1 ;

	if (set == 0) {
		vector<uint> colStats(  collectStatistics[0]->PostFilt->get_rstat_Vmed(which));
		vector<size_t> ra(collectStatistics[0]->PostFilt->getVrange(which) );

		if (which == 1) {
			give << "qual_\tFilterObs" << endl;
		} else {
			give << "Length\tFilterObs" << endl;
		}
		for (size_t i = ra[0]; i < ra[1] ; i++) {
				give << i << "\t" << colStats[i] << endl;
		}
	} else if (set == 1) {
		vector<size_t> ra(2,0), tra;
		vector<uint> stat;
		vector<bool> skips(6, false);
		vector<vector<uint>> matHist;
		if (which == 1) {	give << "#qual_\t";
		} else { give << "#Length\t"; }
		if (p2stat && b_doFilter) { give << "FilteredP1\tFilteredP2\t"; } 
		else if ( b_doFilter){ give << "Filtered\t"; skips[1] = true; }
		else { skips[1] = true; skips[0] = true; }
		
		if ( bAdditionalOutput && b_doFilter ) {
			if (p2stat) { give << "AddFilterP1\tAddFilterP2\t"; }
			else { give << "AddFilter\t"; skips[3] = true; }
		} else {
			skips[2] = true; skips[3] = true;
		}
		if (p2stat) { give << "RawReadsP1\tRawReadsP2"; }
		else { give << "RawReads"; skips[5] = true; }
		give << endl;
		ra = collectStatistics[0]->PostFilt->getVrange(which);
		if (p2stat) { tra = collectStatistics[1]->PostFilt->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		if (!skips[2]) {
			tra = statAddition[0]->PostFilt->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
			if (p2stat) { tra = statAddition[0]->PostFilt->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		} 
		tra = collectStatistics[0]->PreFilt->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
		if (p2stat) { tra = collectStatistics[1]->PreFilt->getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		vector<uint> empt(ra[1], 0);
		matHist = vector<vector<uint>>(6, empt);
		for (size_t kk = 0; kk < 6; kk++) {
			if (skips[kk]) { continue; }
			switch (kk) {
			case 0:	stat = collectStatistics[0]->PostFilt->get_rstat_Vmed(which); break;
			case 1:	stat = collectStatistics[1]->PostFilt->get_rstat_Vmed(which); break;
			case 2:	stat = statAddition[0]->PostFilt->get_rstat_Vmed(which); break;
			case 3:	stat = statAddition[1]->PostFilt->get_rstat_Vmed(which); break;
			case 4:	stat = collectStatistics[0]->PreFilt->get_rstat_Vmed(which); break;
			case 5:	stat = collectStatistics[1]->PreFilt->get_rstat_Vmed(which); break;
			}
			for (size_t i = 0; i < stat.size(); i++) {
				if (i>=ra[1]) {break;}
				matHist[kk][i] = stat[i];
			}
		}
		for (size_t i = ra[0]; i < ra[1] ; i++) {
			give << i ;
			for (size_t kk = 0; kk < matHist.size(); kk++) {
				if (skips[kk]) { continue; }
				give << "\t" << matHist[kk][i];
			}
			give << endl;
		}

	}
}

vector<int> Filters::combiSmplConvergeVec(const vector<string>& inNames){
	vector<int> retV(inNames.size(), -1);
	unordered_map<string, int> smpl2combi;
	unordered_map<string, int>::iterator s2cIT;
	int cntGrps(-1);
	for (size_t i = 0; i < SampleID_Combi.size(); i++){
		s2cIT = combiMapCollectGrp.find(SampleID_Combi[i]);
		if (s2cIT == combiMapCollectGrp.end()){
			cntGrps++;
			combiMapCollectGrp[SampleID_Combi[i]] = cntGrps;
		}
		smpl2combi[SampleID[i]] = combiMapCollectGrp[SampleID_Combi[i]];
	}
	for (size_t i = 0; i < inNames.size(); i++){
		s2cIT = smpl2combi.find(inNames[i]);
		if (s2cIT == smpl2combi.end()){
			cerr << "Can't find SampleID " << inNames[i] << " in reference sample_id_ Names\n"; exit(113);
		}
		retV[i] = s2cIT->second;
	}
	return retV;
}


string Filters::shortStats( const string & file) {
	shared_ptr<collectstats> cst = collectStatistics[0];
	string ret("");
	if (file != ""){
		ret+= file + "\n";
	}
	if (pairedSeq > 1) {
		ret+= "Pair 1: ";
	}
	char buffer[50];
	float tmp = (100.f*float(cst->total - cst->totalRejected) / (float)cst->total);
	ostringstream os;
	os << tmp << "% of " << cst->total << " reads accepted (" << (100.f * float(cst->total - cst->Trimmed) / (float)cst->total) << "% end-trimmed)\n";
//	sprintf(buffer, "%.3f%% of %d", tmp, cst->total); ret += buffer;
//	sprintf(buffer," reads accepted (%.3f%% end-trimmed)\n", (100.f* float(cst->total - cst->Trimmed) / (float)cst->total)); ret += buffer;
	
	if (pairedSeq > 1) {
		shared_ptr<collectstats> cst = collectStatistics[1]; 
		os<<"Pair 2: "<< (100.f * float(cst->total - cst->totalRejected) / (float)cst->total) << "% of " << cst->total;
		os << " reads accepted ("<< (100.f * float(cst->total - cst->Trimmed) / (float)cst->total) <<"% end - trimmed)\n";
//		sprintf(buffer,"Pair 2: %.3f%% of %d", (100.f*float(cst->total - cst->totalRejected) / (float)cst->total), cst->total); ret += buffer;
//		sprintf(buffer," reads accepted (%.3f%% end - trimmed)\n", (100.f* float(cst->total - cst->Trimmed) / (float)cst->total)); ret += buffer;
	}
	ret = os.str();
	return ret;
}
void Filters::printGC(ostream& os,int Npair) {
	os << "Subset\t\tOccurence\t\tAvg.Quality\n";
	os << "\tA\tT\tG\tC\tA\tT\tG\tC\n";
	os << "R1 pre-filter";
	collectStatistics[0]->PreFilt->printGCstats(os);
	if ( Npair > 1 ) {
		os << "R2 pre-filter";
		collectStatistics[1]->PreFilt->printGCstats(os);
	}
	if ( !b_doFilter ) {return;}
	os << "R1 filtered";
	collectStatistics[0]->PostFilt->printGCstats(os);
	if ( Npair > 1 ) {
		os << "R2 filtered";
		collectStatistics[1]->PostFilt->printGCstats(os);
	}
}

void Filters::printStats(ostream& give, string file, string outf, bool greenQualStats) {
	//TODO switch min_l to min_l_add
	shared_ptr<collectstats> cst = collectStatistics[0];
	shared_ptr<collectstats> cst2 = collectStatistics[1];
	if (cst->total != cst->total2) {
		cerr << "Unequal read numbers recorded " << cst->total << "," << cst->total2 << endl;
	}
	if (!greenQualStats) {
		cst = statAddition[0];
	}
	bool p2stat = pairedSeq > 1 && greenQualStats;
	give << "sdm " << sdm_version << " " << sdm_status << endl;
	if (file.length()>0){
		give<<"Input File:  "<<file<<endl;
	}
	if (outf == "-") {
		give << "Output File: stdout\n";
	}else if (outf.length() > 0) {
		give<<"Output File: "<<outf<<endl;
	}
	if ( !b_doFilter ) {
		give << "No valid Filter file provided; no filtering done on files\n";
		return;
	}
	if (!greenQualStats) {
		give << "Statistics of reads that passed the mid qual filter\n";
	} else {
		give << "Statistics of high quality reads\n";
	}
	float remSeqs = float (cst->total-cst->totalRejected);
	give << endl;
	if (!greenQualStats){
		give << "Reads not High qual_: " << intwithcommas((int)cst->totalRejected);
	} else {
		give << "Reads processed: " << intwithcommas((int)cst->total);
		if (p2stat) {
			give << "; " << intwithcommas((int)cst2->total) << " (pair 1;pair 2)";
		}
	}
	give << endl;
	//int numAccept = (int)(cst->total - cst->totalRejected);
	int numAccept = (int)(cst->totalSuccess );
	if (!greenQualStats){
		give << "Rejected:" << intwithcommas((int)(collectStatistics[0]->totalRejected - numAccept)) << endl << "Accepted (Mid+High qual): " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst->Trimmed) << " were end-trimmed";
	} else {
		if (p2stat) {
			give << "Rejected: " << intwithcommas((int)cst->totalRejected) << "; " << intwithcommas((int)cst2->totalRejected) << endl << "Accepted (Mid+High qual): " << intwithcommas((int)numAccept) << "; " << intwithcommas((int)cst2->totalSuccess) << " (" << intwithcommas((int)cst->Trimmed) << "; " << intwithcommas((int)cst2->Trimmed) << " were end-trimmed";
		}
		else {
			give << "Rejected: " << intwithcommas((int)cst->totalRejected) << endl << "Accepted (Mid+High qual): " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst->Trimmed) << " were end-trimmed";
		}
	}


	if ( false && bPrimerR ) { //confusing collectStatistics
		give << ", with rev. primer: " << intwithcommas((int)cst->RevPrimFound); if ( p2stat ) { give << "; " << intwithcommas((int)cst2->RevPrimFound); }
	}
	give<<")"<<endl;

	if (pairedSeq>1) {
		give <<"Singletons among these: " << intwithcommas((int)cst->singleton) << "; " << intwithcommas((int)cst2->singleton) << endl;
	}
	give << "Bad Reads recovered with dereplication: " << intwithcommas((int)cst->DerepAddBadSeq) << endl;

	if ( bShortAmplicons ) {
		give << "Short amplicon mode.\n";
	}

	if ( checkBC2ndRd() ) {
		give << "Looked for switched read pairs (" << intwithcommas(revConstellationN) << " detected)" << endl;
	}
	if (greenQualStats) {
		collectStatistics[0]->PostFilt->printStats2(give, remSeqs,0);
		collectStatistics[1]->PostFilt->printStats2(give, remSeqs,1);
	} else {
		statAddition[0]->PostFilt->printStats2(give,remSeqs,0);
	}

	give << "Trimmed due to:\n";
	//EWwidth, EWthr  no stat for this so far
	float dval = (float)EWthr;	if (!greenQualStats) { dval = (float)alt_EWthr; }
	int Xval = EWwidth; 
	if (EWthr > 0) {
		give << "  > " << EWthr << " avg qual_ in " << Xval << " bp windows : " << spaceX(10 - digitsFlt(dval)) << intwithcommas(cst->QWinTrimmed);
			if (p2stat) { give << "; " << intwithcommas((int)cst2->QWinTrimmed); } give << endl;
	}
	dval = (float)maxAccumQP;	if (!greenQualStats) { dval = (float)alt_maxAccumQP; }
	if (maxAccumQP>0.0) {
		give << "  > (" << dval << ") acc. errors, trimmed seqs : " << spaceX(8 - digitsFlt(dval)) << intwithcommas((int)cst->AccErrTrimmed);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->AccErrTrimmed); } give << endl;
	}

	give << "Rejected due to:\n";
	float val = (float)min_l;
	if (val == -1.f) {val = min_l_p;}
	if (!greenQualStats){ val = (float)alt_min_l; }

	give << "  < min sequence_ length (" << val << ")  : " << spaceX(18 - digitsFlt(val)) << intwithcommas((int)cst->minL);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->minL); } give << endl;
	if (cst->minLqualTrim>0){//this is failed because seq was too short after trimming
		give << "       -after Quality trimming : " << spaceX(10) << intwithcommas((int)cst->minLqualTrim);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->minLqualTrim); } give << endl;
	}
	float valf = min_q;	if (!greenQualStats){ valf = alt_min_q; }
	give << "  < avg Quality (" << valf << ")  : " << spaceX(21 - digitsInt((int)min_q)) << intwithcommas((int)cst->AvgQual);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->AvgQual); } give << endl;
	give << "  < window (" << FQWwidth << " nt) avg. Quality (" << FQWthr << ")  : " << spaceX(5 - digitsInt(FQWwidth)) << intwithcommas((int)cst->QualWin);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->QualWin); } give << endl;
	give << "  > max sequence_ length (" << max_l << ")  : " << spaceX(18 - digitsInt(max_l)) << intwithcommas((int)cst->maxL);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->maxL); } give << endl;
	give << "  > (" << maxHomonucleotide << ") homo-nt run  : " << spaceX(21 - digitsInt(maxHomonucleotide)) << intwithcommas((int)cst->HomoNT);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->HomoNT); } give << endl;
	int val2 = MaxAmb;	if (!greenQualStats){ val2 = alt_MaxAmb; }
	give << "  > (" << val2 << ") amb. Bases  : " << spaceX(22 - digitsInt(val2)) << intwithcommas((int)cst->MaxAmb);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->MaxAmb); } give << endl;
	if (BinFilP >= 0.f){
		give << "  > (" << BinFilErr << ") binomial est. errors : " << spaceX(13 - digitsFlt(BinFilErr)) << intwithcommas((int)cst->BinomialErr);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->BinomialErr); } give << endl;
	}
	if ((removeAdapter && tAdapter != "") || (bDoMultiplexing || cst->PrimerFail > 0) || ((!greenQualStats && alt_bRequireFwdPrim) || bRequireFwdPrim)
        || bPrimerR) {
		give << "Specific sequence searches:\n";
	}
	if (removeAdapter && tAdapter != "") {
		give << "  -removed adapter (" << tAdapter << ")  : " << spaceX(18 - (uint)tAdapter.length()) << intwithcommas((int)cst->adapterRem);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->adapterRem); } give << endl;
	}
	if ( (bDoMultiplexing || cst->PrimerFail>0) || ((!greenQualStats && alt_bRequireFwdPrim) || bRequireFwdPrim) ){
		give << "  -With fwd Primer remaining (<= " << PrimerErrs << " mismatches";
		if ((!greenQualStats && alt_bRequireFwdPrim) || bRequireFwdPrim){
			give << ", required) : ";
			give << spaceX(1 - digitsInt(PrimerErrs));
		}
		else  {
			give <<") : "<< spaceX(11 - digitsInt(PrimerErrs));
		}
		give << intwithcommas((int)cst->PrimerFail);
		if ( p2stat ) { give << "; " << intwithcommas((int)cst2->PrimerFail) << endl; }
		else { give << endl; }
	} 
	if (bPrimerR){
		give<<"  -With rev Primer remaining (<= "<<PrimerErrs <<" mismatches";
		if ((!greenQualStats && alt_bRequireRevPrim) || bRequireRevPrim){ give << ", required) : "; 
			give << spaceX(1 - digitsInt(PrimerErrs));
		}	else  {
			give << ") : "<<spaceX(11 - digitsInt(PrimerErrs));
		}
		give << intwithcommas((int)cst->PrimerRevFail);
		if ( p2stat ) { give << "; " << intwithcommas((int)cst2->PrimerRevFail) << endl; }
		else { give << endl; }
	}
	if (bDoMultiplexing){
		if (bDoBarcode){
			give << "  -Barcode unidentified (max " << barcodeErrors_ << " errors) : " << spaceX(19 - digitsInt(barcodeErrors_)) << intwithcommas((int)cst->TagFail);
			if (p2stat && (cst2->TagFail > 0 || doubleBarcodes())) { give << "; " << intwithcommas((int)cst2->TagFail); give << " (" << intwithcommas((int)cst->dblTagFail) << " pairs failed)"; }
			give << endl;

			if (barcodeErrors_ > 0){
				give << "    -corrected barcodes: " << spaceX(18) << intwithcommas((int)cst->suc_correct_BC);
				if (p2stat){give << "; " << intwithcommas((int)cst2->suc_correct_BC);	}
				give << endl;
				//<< ", failed to correct barcode: " << spaceX(5 - digitsInt(FQWwidth)) << intwithcommas((int)cst->fail_correct_BC) << endl;
			}
			
			if ( bDoBarcode2 ) {
				give << "    -used dual index barcodes";
				if ( BCdFWDREV[0].reversedBCs || BCdFWDREV[1].reversedBCs ) {
					give << " (reversed_ ";
					if ( BCdFWDREV[1].reversedBCs && BCdFWDREV[0].reversedBCs ) {
						give << " fwd & rev";
					} else	if ( BCdFWDREV[0].reversedBCs ) {
						give << " fwd";
					} else	if ( BCdFWDREV[1].reversedBCs ) {
						give << " rev";
					}
					give << " BCs)" << endl;
				}
				
			} else if ( BCdFWDREV[0].reversedBCs ) {
				give << "    -reversed_ all barcodes" << endl;
			}
			give << endl << "SampleID";
			if (bDoCombiSamples){
				give << "\tSampleGroup";
			}
			give << "\tBarcode";
			if ( bDoBarcode2 ) {give << "\tBarcode2";}
			give << "\tInstances\n";
			for (unsigned int i =0; i<Barcode.size();i++){
				give << SampleID[i] << "\t";
				if (bDoCombiSamples){ give << SampleID_Combi[i] << "\t"; }
				give << Barcode[i];
				if ( bDoBarcode2 ) {
					give << "\t"<<Barcode2[i];
				}
				give << "\t" << intwithcommas((int)cst->BarcodeDetected[i]) << endl;
			}
		} else if (bDoHeadSmplID){
			give << "  -Failed to assign sequences to header tag : " << intwithcommas((int)barcodeErrors_ ) << endl;
			give << endl << "SampleID\t";
			if (bDoCombiSamples){	give << "\tSampleGroup";	}
			give << "\tSampleID\tInstances\n"; 
			for (unsigned int i =0; i<Barcode.size();i++){
				give << SampleID[i] << "\t";
				if (bDoCombiSamples){ give << SampleID_Combi[i] << "\t"; }
				give << HeadSmplID[i] << "\t" << intwithcommas((int)cst->BarcodeDetected[i]) << endl;
			}
		}
	}
}

//statistics for each single sample 
void Filters::SmplSpecStats(ostream & give){
	shared_ptr<collectstats> cst = collectStatistics[0];
	shared_ptr<collectstats> cst2 = collectStatistics[1];
	bool p2stat = pairedSeq > 1 ;
	give << std::setprecision(3);
	give  << "SampleID";
	if ( bDoCombiSamples ) {
		give << "\tSampleGroup";
	}
	if ( bDoHeadSmplID ) {
		give << "\tSampleID";

	} else {
		if ( bDoBarcode2 ) { give << "\tBarcode\tBarcode2"; }
		else { give << "\tBarcode"; }
	}
	if ( p2stat ) {
		give << "\tRead1Accepted\tRead1Filtered\tRead1PassedFrac\tRead2Accepted\tRead2Filtered\tRead2PassedFrac\n";
	} else {
		give << "\tReadsAccepted\tReadsFailed\tPassed%\n";
	}
	for ( unsigned int i = 0; i<Barcode.size(); i++ ) {
		give << SampleID[i] << "\t";
		if ( bDoCombiSamples ) { give << SampleID_Combi[i] << "\t"; }
		if ( p2stat ) {
			give << Barcode[i] << "\t";
			if ( bDoBarcode2 ) { give << Barcode2[i] << "\t"; }
			give <<	cst->BarcodeDetected[i] << "\t" << cst->BarcodeDetectedFail[i] << "\t";
			float totSum = (float(cst->BarcodeDetected[i]) + float(cst->BarcodeDetectedFail[i]));
			if ( totSum > 0 ) {
				give << float(cst->BarcodeDetected[i]) / totSum << "\t";
			} else {give << "NA\t";	}

			give << cst2->BarcodeDetected[i] << "\t" << cst2->BarcodeDetectedFail[i] <<"\t" ;
			totSum = (float(cst2->BarcodeDetected[i]) + float(cst2->BarcodeDetectedFail[i]));
			if ( totSum > 0 ) {
				give << float(cst2->BarcodeDetected[i]) / totSum << "";
			} else { give << "NA"; }
			give<<endl;
		} else {
			give << Barcode[i] << "\t";
			if ( bDoBarcode2 ) { give << Barcode2[i] << "\t"; }
			give << (int)cst->BarcodeDetected[i]
				<< "\t" << cst->BarcodeDetectedFail[i] << "\t";
				
			float totSum = (float(cst->BarcodeDetected[i]) + float(cst->BarcodeDetectedFail[i]));
			if ( totSum > 0 ) {
				give << float(cst->BarcodeDetected[i]) / totSum << "";
			} else { give << "NA"; }
				give<< endl;
		}
	}



}

void ReportStats::calcSummaryStats(float remSeqs, unsigned int min_l, float min_q){
	if (remSeqs == 0){ return; }
	if (bMedianCalcs){
		rstat_Smed = (int) calc_median(rstat_VSmed,0.5f);
		rstat_Qmed = (int) calc_median(rstat_VQmed,0.5f);
		USQS=0.f;
	}
	RSQS = ( ( (float(rstat_NTs)/remSeqs)/(float)min_l ) +
			( (float(rstat_qualSum)/remSeqs) / min_q) ) / 2.f;
}

void Filters::addStats(Filters* fil, vector<int>& idx){
	for (size_t i = 0; i < 2; i++) {
		collectStatistics[i]->addStats(fil->collectStatistics[i], idx); 
		if (bAdditionalOutput) {
			statAddition[i]->addStats(fil->statAddition[i], idx);
		}
	}
	maxReadsPerOFile = fil->maxReadsPerOFile;
	ReadsWritten = fil->writtenReads();//the idea here is to have a number of reads in CURRENT file, not total reads
	OFileIncre = fil->getFileIncrementor();
	revConstellationN += fil->revConstellationN;
}




void collectstats::addStats(shared_ptr<collectstats> cs, vector<int>& idx){
	if (BarcodeDetected.size() > (uint)10000){
		cerr<<"Unrealistic number of barcodes (>10000) in addStats\n"; exit(79);}
	int BCS = (int)BarcodeDetected.size();
	for (unsigned int i=0;i<idx.size();i++){
		if (idx[i] >= BCS){ return; }
		//assert(idx[i] < BCS);
		BarcodeDetected[idx[i]] += cs->BarcodeDetected[i];
		BarcodeDetectedFail[idx[i]] += cs->BarcodeDetectedFail[i];
	}
	maxL += cs->maxL;	PrimerFail += cs->PrimerFail ;
	AvgQual += cs->AvgQual; HomoNT += cs->HomoNT;
	PrimerRevFail += cs->PrimerRevFail;
	minL += cs->minL ; minLqualTrim+= cs->minLqualTrim; TagFail += cs->TagFail;
	MaxAmb += cs->MaxAmb ; QualWin += cs->QualWin;
	Trimmed += cs->Trimmed ; AccErrTrimmed+= cs->AccErrTrimmed; total += cs->total;
	QWinTrimmed += cs->QWinTrimmed;
	totalRejected += cs->totalRejected;
	fail_correct_BC += cs->fail_correct_BC; suc_correct_BC += cs->suc_correct_BC ;
	failedDNAread += cs->failedDNAread; adapterRem += cs->adapterRem ;
	RevPrimFound += cs->RevPrimFound;
	singleton += cs->singleton;
	BinomialErr += cs->BinomialErr;
	dblTagFail += cs->dblTagFail;
	DerepAddBadSeq += cs->DerepAddBadSeq;
	total2 += cs->total2; totalSuccess += cs->totalSuccess;
	PostFilt->addRepStats(cs->PostFilt);
	PreFilt->addRepStats(cs->PreFilt);
}

void collectstats::reset(){
	singleton=0;
	size_t BCsiz = BarcodeDetected.size();
	for (unsigned int i=0; i<BCsiz;i++){
		BarcodeDetected[i]=0;
		BarcodeDetectedFail[i] = 0;
	}
	maxL=0; PrimerFail=0;AvgQual=0; HomoNT=0;
	PrimerRevFail=0;
	minL=0; TagFail=0; MaxAmb=0; QualWin=0;
	Trimmed=0; total=0; totalRejected=0;
	fail_correct_BC=0; suc_correct_BC=0;failedDNAread=0;
	adapterRem = 0; RevPrimFound = 0; DerepAddBadSeq = 0;
	total2 = 0; totalSuccess = 0;

	//reportStats
	PostFilt->reset(); PreFilt->reset();
}


inline void ReportStats::addDNAStats(shared_ptr<DNA> d){
	stats_mutex.lock();
	//pretty fast
	addMeanStats(d->length(),(int) d->getAvgQual(), (float)d->getAccumError());
	//NT specific quality scores
	
//    d->ntSpecQualScores(QperNT, NTcounts); // Thread safe

    // Test this function
    addNtSpecQualScores(d);
	
	//more memory intensive
	if (bMedianCalcs){
		//quali
		uint avq = (uint) (d->getAvgQual() + 0.5f);
        addMedian2Histo(avq, rstat_VQmed); // Thread safe (?)
        addMedian2Histo(d->length(), rstat_VSmed); // Thread safe
	}
	stats_mutex.unlock();
}


inline void ReportStats::mergeStats(data_MT &data) {
    rstat_NTs += data.total_nts;
    rstat_totReads += data.total_reads;
    rstat_qualSum += data.qual_sum;
    rstat_accumError += (float) data.accum_error;

    if (data.per_base_quality_sum.size() > QperNT.size())
        QperNT.resize(data.per_base_quality_sum.size(), 0);
    for (size_t i = 0; i < data.per_base_quality_sum.size(); i++) {
        QperNT[i] += data.per_base_quality_sum[i];
    }
    if (data.nucleotide_counter.size() > NTcounts.size())
        NTcounts.resize(data.nucleotide_counter.size(), 0);
    for (size_t i = 0; i < data.nucleotide_counter.size(); i++) {
        NTcounts[i] += data.nucleotide_counter[i];
    }
}

inline void ReportStats::addDNAStatsMT(shared_ptr<DNA> d, data_MT *data){
    data->total_nts += d->length();
    data->qual_sum += uint(d->getAvgQual()+0.5f);
    data->accum_error += d->getAccumError();
    ++data->total_reads;

    //NT specific quality scores
    d->ntSpecQualScores(data->per_base_quality_sum, data->nucleotide_counter);


    // Not yet with separate arrays
    //more memory intensive
    if (bMedianCalcs){
        uint avq = (uint) (d->getAvgQual() + 0.5f);
        addMedian2Histo(avq, rstat_VQmed);
        addMedian2Histo(d->length(), rstat_VSmed);
    }
}



void ReportStats::reset() {
	rstat_totReads = 0; rstat_NTs = 0; rstat_qualSum=0;
	rstat_Qmed = 0; rstat_Smed = 0; 
	RSQS = 0.f; USQS = 0.f; rstat_accumError = 0.f;
	QperNT.resize(1000,0); NTcounts.resize(1000,0);
	std::fill(QperNT.begin(), QperNT.end(), 0);
	std::fill(NTcounts.begin(), NTcounts.end(), 0);
	rstat_VQmed.resize(0); rstat_VSmed.resize(0);
}
unsigned int ReportStats::lowest(const vector<uint>& in){
	for (int i=0;i<(int)in.size();i++){
		if (in[i]>0){return i;}
	}
	return(0);
}
unsigned int ReportStats::highest(const vector<uint>& in){
	if (in.size()==0){return 0;}
	for (int i=(int)in.size()-1;i>=0;i--){
		if (in[i]>0){return i;}
	}
	return 0;
}
void ReportStats::printGCstats(ostream& give) {
	//NT_POS['A'] = 0; NT_POS['T'] = 1; NT_POS['G'] = 2; NT_POS['C'] = 3;	NT_POS['N'] = 4;
	vector<string> NTs(6, "X"); NTs[0] = "A"; NTs[1] = "T";
	NTs[2] = "G"; NTs[3] = "C"; NTs[4] = "N";
	for ( uint i = 0; i < 4; i++ ) {
		give << "\t" << NTcounts[i];
	}
	//give << endl;
	for ( uint i = 0; i < 4; i++ ) {
		give << "\t" << float(QperNT[i]) / float(NTcounts[i]);
	}
	give << endl;
}
void ReportStats::printStats2(ostream& give, float remSeqs,int pair){
	if ( pair == 1 ) {
		return;//deactivate for now
	}
	if ( pair == 0 ) {
		if ( bMedianCalcs ) {
			unsigned int minS = lowest(rstat_VSmed);
			unsigned int maxS = highest(rstat_VSmed);
			unsigned int minQ = lowest(rstat_VQmed);
			unsigned int maxQ = highest(rstat_VQmed);
			give << "Min/Avg/Max stats Pair 1";// -RSQS : "<<RSQS;
			if ( remSeqs == 0 ) {
				give << "\n     - sequence Length : " << "0/0/0"
					<< "\n     - Quality :   " << "0/0/0";
			} else {
				give << "\n     - sequence Length : " << minS << "/" << float(rstat_NTs) / (float)rstat_totReads << "/" << maxS
					<< "\n     - Quality :   " << minQ << "/" << float(rstat_qualSum) / (float)rstat_totReads << "/" << maxQ;
			}
		} else {
			give << "Average Stats - RSQS : " << RSQS;
			if ( remSeqs == 0 ) {
				give << "\n     - sequence Length : " << "0/0/0"
					<< "\n     - Quality :   " << "0/0/0";
			} else {
				give << "\n     - sequence Length : " << float(rstat_NTs) / (float)rstat_totReads
					<< "\n     - Quality :   " << float(rstat_qualSum) / (float)rstat_totReads;
			}
		}
	} else {
		give << "Pair 2 stats";
	}
	
		if (bMedianCalcs){
			//give << "Median Stats Pair 1";// -USQS : " << USQS;
			give << "\n     - Median sequence Length : " << rstat_Smed << ", Quality : " << rstat_Qmed;// << "\n";
		}
		give << "\n     - Accum. Error " << (rstat_accumError / (float)rstat_totReads) << "\n";
}

//add the stats from a different ReportStats object
void ReportStats::addRepStats(const ReportStats* RepoStat){
	//report stats:
	rstat_NTs += RepoStat->rstat_NTs; rstat_totReads += RepoStat->rstat_totReads;
	rstat_qualSum += RepoStat->rstat_qualSum;
	rstat_accumError += RepoStat->rstat_accumError;
	for ( uint i = 0; i < 6; i++ ) {
		QperNT[i] += RepoStat->QperNT[i];
		NTcounts[i] += RepoStat->NTcounts[i];
	}
	if (bMedianCalcs){
		//vectors for median calcs
		if (rstat_VQmed.size() < RepoStat->rstat_VQmed.size()){
			rstat_VQmed.resize(RepoStat->rstat_VQmed.size(),0);
			assert(rstat_VQmed.size() < 10000);
		}
		for (unsigned int i=0; i<RepoStat->rstat_VQmed.size(); i++){
			rstat_VQmed[i] += RepoStat->rstat_VQmed[i];
		}
		if (rstat_VSmed.size() < RepoStat->rstat_VSmed.size()){
			rstat_VSmed.resize(RepoStat->rstat_VSmed.size(),0);
			//times change.. wrong assert here
			//assert(rstat_VSmed.size() < 10000);
		}
		for (unsigned int i=0; i<RepoStat->rstat_VSmed.size(); i++){
			rstat_VSmed[i] += RepoStat->rstat_VSmed[i];
		}
	}
}
//calculate median value from data stored as histogram-vector
// for median use perc = 0.5f
float ReportStats::calc_median(vector<unsigned int>& in, float perc){
	unsigned int sum = 0;
	for (unsigned int i=0; i<in.size(); i++){
		sum += in[i];
	}
	unsigned int threshold = (unsigned int) (((float)sum) * perc);
	sum = 0;
	for (unsigned int i=0; i<in.size(); i++){
		sum += in[i];
		if (sum >= threshold){
			return (float) i;
		}
	}

	return 0.f;
}
void ReportStats::add_median2histo(vector<unsigned int>& in, vector<unsigned int>& histo)
{
	unsigned int max = *max_element(in.begin(),in.end());
	if (max> histo.size()){
		if (max > 10000){cerr<<"max bigger 10000.\n"; exit(77);}
		histo.resize(max,0);
	}
	for (unsigned int i=0; i<histo.size(); i++){
		histo[ in[i] ] ++;
	}
}
vector<size_t> ReportStats::getVrange(int which) {
	if (which == 1) {
		return medVrange(rstat_VQmed);
	} else {
		return medVrange(rstat_VSmed);
	}
}
vector<size_t> ReportStats::medVrange(const vector<uint> x) {
	vector<size_t> ret(2, 0); ret[1] = 0; bool empty = true;
	for (size_t i = 0; i < x.size(); i++) {
		if (x[i]>0) {
			if (empty) { ret[0] = i; empty = false; }
			ret[1] = i;
		}
	}
	ret[1]++;
	return ret;
}

void ReportStats::addMedian2Histo(unsigned int in, vector<unsigned int>& histo)
{
    {
        if (in >= histo.size()) {
            histo.resize(in + 3, 0);
            assert(in < 1e6);
        }
        histo[in] += 1;
    }
}


////////////  UCLINKS  ///////////////
UClinks::~UClinks(){
	ucf.close();
	mapdere.close();
	if (merger != nullptr) {
		delete merger;
		merger = nullptr;
	}
/*
	for (uint i=0; i< bestDNA.size();i++){
		delete bestDNA[i];
	}
//	for (std::map<string, shared_ptr<DNA>>::iterator iterator = unusedID.begin(); iterator != unusedID.end(); iterator++) {
//		delete (*iterator).second;
//	}
	for (uint i = 0; i<oldDNA.size(); i++) {
		if (oldDNA[i] != NULL) { delete oldDNA[i]; }
	}
	for (uint i = 0; i<oldDNA2.size(); i++) {
		if (oldDNA2[i] != NULL) { delete oldDNA2[i]; }
	}
	*/
}
UClinks::UClinks( OptContainer& cmdArgs):
	CurSetPair(-1),//maxOldDNAvec(20000),
	DNAunusedPos(0), derepMapFile(""),
	bestDNA(0, NULL), oriKey(0), bestPID(0), bestLEN(0),
	clusCnt(0), uclines(0),
	SEP(""), 
	UCread(false), pairsMerge(false), MAPread(false),
	b_derepAvailable(false),
	UPARSE8up(false), UPARSE9up(false), UPARSE11up(false),
	UpUcFnd(false),
	otuTerm("OTU"), otuOUTterm("OTU_"),
	RefDBmode(false), RefDBotuStart(-1),
	SeedsAreWritten(false),
	OTUmat(0), unregistered_samples(false),
	doChimeraCnt(false), OTUnumFixed(true),
	b_merge_pairs_optiSeed_(false),merger(nullptr),
	totalDerepCnt(0)

{
	//read in UC file and assign clusters
	SEP = cmdArgs["-sample_sep"];
	string str = cmdArgs["-optimalRead2Cluster"];
	

	ucf.open(str.c_str(),ios::in);
	if (!ucf){
		UCread = true;
		cerr<<"Could not open uc file\n"<< str<<endl; exit(46);
	}
	if (cmdArgs.find("-derep_map") != cmdArgs.end()) {
		derepMapFile = cmdArgs["-derep_map"];
	}
	if (cmdArgs.find("-count_chimeras") != cmdArgs.end() &&
		cmdArgs["-count_chimeras"] == "T") {
		doChimeraCnt = true;
	}
	if (cmdArgs.find("-merge_pairs_seed") != cmdArgs.end() &&
		cmdArgs["-merge_pairs_seed"] == "1") {
		b_merge_pairs_optiSeed_ = true;
	}
	if (cmdArgs["-uparseVer"] != "") {
		if (cmdArgs["-uparseVer"] == "N11") {//UNOISE v11
			UPARSE8up = true;
			otuTerm = "Zot";
			otuOUTterm = "Zotu";
		}
		else {
			int upVer = atoi(cmdArgs["-uparseVer"].c_str());
			if (upVer >= 8 && upVer < 9) {
				UPARSE8up = true; UpUcFnd = true;
			}
			else if (upVer >= 9 && upVer < 11) {
				UPARSE8up = true; UPARSE9up = true; UpUcFnd = true;
			}
			else if (upVer >= 11) {
				UPARSE8up = true; UPARSE9up = true; UPARSE11up = true; UpUcFnd = true;
				otuTerm = "otu";
			}
		}
	}
}


void UClinks::readDerepInfo(const string dereM) {
	mapdere.open(dereM.c_str(), ios::in);
	if (!mapdere) {
		cerr << "Can't open " << dereM << ". \nAborting\n"; exit(54);
	}
	int cnt(0);
	totalDerepCnt = 0;
	b_derepAvailable = true;
	string line("");
	//only read header
	while (!safeGetline(mapdere, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (line.length() < 3) { continue; }
		string segments;
		stringstream ss;
		ss << line;
		int tbcnt = -1;

		//cmdArgs["-i_MID_fastq"]
		while (getline(ss, segments, '\t')) {
			tbcnt++;
			if (cnt == 0) { //search for header
				if (tbcnt == 0) {
					if (segments != "#SMPLS") { cerr << "First line in dereplicate map has to start with #SMPLS, corrupted file.\n" << dereM << endl; exit(98); }
					continue;
				}
				vector<string> spl = header_string_split(segments, ":");
				int ccnt = stoi(spl[0]);
				SmplIDs[spl[1]] = ccnt;				
			} else {//actual derep info per sorted line
				cerr << "wrong mapping reading\n"; exit(64);
			}
		}
		if (cnt == 0) {
			OTUmat.resize(SmplIDs.size(), vector<matrixUnit>(clusCnt + 1, 0));
		}
		break;
	}
	
}

//read in dereplicated info from derep.map and derep.hq.fq (in IS) -> they are in the same order
int UClinks::oneDerepLine(shared_ptr<DNAunique> d) {
	
	if (MAPread){
		return 0;
	}
	string line("");
	int cnt;
	while (!safeGetline(mapdere, line).eof()) {
		//if (line.length() < 3) { continue; }
		if (line[0] == '#') {continue;}
		string segments;
		//stringstream ss;
		//string head("");
		//ss << line;
		int tbcnt = -1;
		long curCnt(0),curCnt2(0);
		size_t strpos(line.find_first_of('\t')), lastpos(0);
		while (strpos != string::npos){
		//while (getline(ss, segments, '\t')) {
			segments = line.substr(lastpos, strpos - lastpos);
			lastpos = strpos+1;
			strpos = line.find_first_of('\t', lastpos);
			tbcnt++;
			if (tbcnt == 0) {
				if (!d->sameHead(segments)) {
					cerr << segments << " is not " << d->getId() << endl;
					//this is actually a fatal error and currently not covered
					//derep.map and fq [which one??] need to have reads in the same order..
					exit(738);
				}
				size_t idx = segments.find(";size=");
				if (idx != string::npos) {
					curCnt = atoi(segments.substr(idx + 6, segments.find(";", idx + 5) - (idx + 6)).c_str());
					//curCnt
				}
				continue;
			}
			size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl+1));
            d->incrementSampleCounterBy(stoi(segments.substr(0, spl)), occur);
			curCnt2 += occur;
			//string tmp = segments.substr(0, spl);
		}
		//last round
		segments = line.substr(lastpos);
		size_t spl = segments.find(":"); long occur = stoi(segments.substr(spl + 1));
        d->incrementSampleCounterBy(stoi(segments.substr(0, spl)), occur);
		curCnt2 += occur;

		if (curCnt2 != curCnt) {
			cerr << "ERROR: Mapping file unique abundance reconstruction: failed to find correct count (" << curCnt2 << " vs " << curCnt<<"):\n" << line << endl;
			//exit(67); //DEBUG
		}
		totalDerepCnt += curCnt2;
		cnt = curCnt2;
		break;
	}

	if (mapdere.eof()) {// end of file
		mapdere.close();
		MAPread = true;
	}
	return cnt;
	//
}

//finishes reading the map, loading the sequences as DNAunique into mem (no DNA, mind)
void UClinks::finishMAPfile(){
	if (MAPread){
		return;
	}
	string line("");
	//bool hardAdd = false;

	int finishMapCnt = 0;

	//go through ucl file line by line
	while (!safeGetline(mapdere, line).eof()) {	
		if (line[0] == '#') { continue; }
		string segments;
		int tbcnt = -1;
		int curCnt(0), curCnt2(0);
		size_t strpos(line.find_first_of('\t')), lastpos(strpos + 1);
		segments = line.substr(0, strpos);
		//1st setup empty DNA for each remaining line
		shared_ptr<DNAunique> d = make_shared<DNAunique>("", segments);
		d->seal();
		size_t idx = segments.find(";size=");
		if (idx != string::npos) {
			curCnt = atoi(segments.substr(idx + 6, segments.find(";", idx + 5) - (idx + 6)).c_str());
		}
		strpos = line.find_first_of('\t', lastpos);

		while (strpos != string::npos){
			segments = line.substr(lastpos, strpos - lastpos);
			lastpos = strpos + 1;
			strpos = line.find_first_of('\t', lastpos);
			tbcnt++;
			size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl + 1));
            d->incrementSampleCounterBy(stoi(segments.substr(0, spl)), occur);
			curCnt2 += occur;
			//string tmp = segments.substr(0, spl);
		}
		//last round
		segments = line.substr(lastpos);
		size_t spl = segments.find(":"); int occur = stoi(segments.substr(spl + 1));
        d->incrementSampleCounterBy(stoi(segments.substr(0, spl)), occur);
		curCnt2 += occur;

		if (curCnt2 != curCnt) {
			cerr << "ERROR: Mapping file unique abundance reconstruction: failed to find correct count (" << curCnt2 << " vs " << curCnt << "):\n" << line << endl;
			exit(67);//DEBUG
		}
		finishMapCnt += curCnt;

		//add d to oldDNA
		string curID = d->getId();
		curID = curID.substr(0, curID.find_first_of(' '));
		unusedID[curID] = DNAunusedPos;
		oldDNA[DNAunusedPos] = d;
		DNAunusedPos++;
			
	}

	if (mapdere.eof()) {// end of file
		mapdere.close();
		MAPread = true;
	}

	cerr << "Found " << finishMapCnt << " counts derep map\n";

}

//functions selects optimal read to represent OTU
void UClinks::findSeq2UCinstruction(shared_ptr<InputStreamer> IS, bool readFQ, 
	Filters* fil){
	if (UCread){return;}

	if (IS->hasMIDseqs()){
		CurSetPair = 0;
	} else {
		CurSetPair = -1;
	}
#ifdef DEBUG
	cerr << "UC seed extension" << endl;
#endif
	int UCcnt = 0;
	DNAidmapsIT fndDNAinOld;
	std::list<string>::iterator listIT;
	shared_ptr<DNAunique> match(NULL); shared_ptr<DNA>  match2(NULL);
	bool cont(true);//,cont2(true);
	string segs("");
	string segs2;
	float perID;
	vector<int> curCLID(0,0);
	//bool sync(false); // syncing of 2 read pairs; not implemented for this function yet

	while ( getUCFlineInfo(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		//int subcnt = 0;
		if ( curCLID.size() == 0 ) { continue; }
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}
		//goes through UC file
		while(cont){
			shared_ptr<DNA> tmpDNA = IS->getDNA(0);
			if (tmpDNA == NULL) { cont = false; break; }//signal that at end of file
			match.reset( new DNAunique(tmpDNA, -1));
			//delete tmpDNA;
			//assummes in original implementation, that we can get derep.map lines with the same 
			//ordering as fq derep (which works normally, just not for dada2 mode)
			UCcnt+= oneDerepLine(match);
			match2 = IS->getDNA(1);
			string curID = match->getId();
			curID = curID.substr(0,curID.find_first_of(' '));
			//check if dnaTemp1 a) matches id_ b) is better
			if (curID != segs){
				//block to store unused DNA & find this id_ in this block
				unusedID[curID] = DNAunusedPos;
				oldDNA[DNAunusedPos] = match;	oldDNA2[DNAunusedPos] = match2;
				DNAunusedPos++;
			} else {
				//assign % identity score to DNA object
#ifdef DEBUG
				cerr << "UC Hit";
#endif
				match->setTempFloat(perID);
				besterDNA(curCLID, match, match2, fil);
				curCLID.resize(0);
				break;
			}
			
			//if (tdn!=NULL && ch1 != tdn->isGreenQual()){cerr<<"isGreenQual is != ch1! Aborting..\n";exit(12);}
		}
	}
#ifdef DEBUG
	cerr << "UC seed initialized" << endl;
#endif
	
	
	//cerr << "\nRead " << totalDerepCnt << " read counts from derep.map\n";
	cerr << "Associated " << totalDerepCnt << ", "<< OTUmatSum() <<" counts to primary clusters\n";

}

void UClinks::finishUCfile(Filters* fil, string addUC, bool bSmplHd){
	if (UCread){
		addUCdo(addUC,bSmplHd);
		UCread = true;
		return;
	}
	string segs;
	string segs2;
	float perID(0.f);
	vector<int> curCLID(0);

	while ( getUCFlineInfo(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}
		cerr << segs << " ";
	}
	
	addUCdo(addUC, bSmplHd);
	UCread = true;
}
void UClinks::addUCdo(string addUC,bool SmplHd) {
	if (addUC == "") {
		return;
	}
	string segs;
	string segs2;
	float perID;
	int cntsAddUC=0;
	vector<int> curCLID(0);
	ucf.open(addUC.c_str(), ios::in);
	if (!ucf) {
		UCread = true;
		cerr << "Could not find additional uc file: " << addUC << endl;
	} else {
		std::cerr << "Reading " << addUC << endl;
	}
	UCread = false;
	if (SmplHd){
		while (getUCFlineInfo(segs, segs2, perID, curCLID, SmplHd)) { curCLID.resize(0); }
		cntsAddUC++;
	}else {//complicated..
		int cnt(0);
		while (getUCFlineInfo(segs, segs2, perID, curCLID, SmplHd)) { 
			//find segs in remaining dereps
			cnt++; //if (cnt < 61403) { continue; }
			if (curCLID.size() == 0) { continue; }
			if (uclInOldDNA_simple(segs, curCLID, cntsAddUC)) {
				curCLID.resize(0);
				continue;
			}

			curCLID.resize(0); 
		}
	}

	cout << "Added " << cntsAddUC << " counts from uc add file\n";

}

//used for specific situation where only empty DNA string was read with derep info attached
bool UClinks::uclInOldDNA_simple(const string& segs,const vector<int>& curCLID, int& countsAdd) {
	//check if DNA is in group of unmatched DNA's ?
	if (segs == ""){ return false; }//empty ucl line
	DNAidmapsIT unusedIT = unusedID.find(segs);
	if (unusedIT != unusedID.end()){// found something
		int mID = (*unusedIT).second;
		//if (mID >= (int)oldDNA.size()) { cerr << "SEC MID too high\n"; exit(55); }
		if (oldDNA[mID] == NULL){ return false; }
		//give sequence a chance to be selected
		matrixUnit matchSiz = (matrixUnit)curCLID.size();
		for (int k = 0; k < matchSiz; k++) {
			add2OTUmat(oldDNA[mID], curCLID[k], matchSiz);
			countsAdd += matchSiz;
		}
		//was matched once to an OTU seed. Even if several later matches, doesn't mater - delete

		unusedID.erase(unusedIT);
//		delete oldDNA[mID]; if (oldDNA2[mID] != NULL){ delete oldDNA2[mID]; }

		oldDNA.erase(mID); oldDNA2.erase(mID);// = NULL;oldDNA2[mID] = NULL;
		return true;
	}
	return false;
}
bool UClinks::uclInOldDNA(const string& segs,const vector<int>& curCLID, float perID,
	Filters* fil) {
	//check if DNA is in group of unmatched DNA's ?
	if (segs == ""){ return false; }//empty ucl line
	DNAidmapsIT unusedIT = unusedID.find(segs);
	if (unusedIT != unusedID.end()){// found something
		//shared_ptr<DNA> fndDNA = (*unusedIT).second;
		//fndDNA->setTempFloat(perID);
		//besterDNA(curCLID, fndDNA, fil);
		//unusedID.erase(unusedIT);

		int mID = (*unusedIT).second;
		//give sequence a chance to be selected
		oldDNA[mID]->setTempFloat(perID);
		besterDNA(curCLID, oldDNA[mID], oldDNA2[mID], fil);
		//remove all trace
		unusedID.erase(unusedIT);
		oldDNA.erase(mID); oldDNA2.erase(mID);
		//oldDNA[mID] = NULL;		oldDNA2[mID] = NULL;
		return true;
	}
	return false;
}

bool UClinks::getUCFlineInfo(string& segs, string& segs2,float& perID, 
	vector<int>& curCLID,  bool addFromHDstring) {
	//reads UC file line by line
	//can also be used to delineate UC's


	if (UCread){return false;}
	//close all file streams
	if (ucf.eof()){
		UCread=true;
		ucf.close();
		return false;
	}
	string line; 
	std::unordered_map<string, int>::iterator itCL;
	while (getline(ucf, line, '\n')) {
	//	cerr<<line<<endl;
		uclines++;
		if (line.length() <= 1){ continue; }	
		if (!UpUcFnd){
			if (line.substr(1, 1) == "\t"){ UPARSE8up = false; 
			}else {	UPARSE8up = true;	}

			if (UPARSE8up){//could also be uparse9

				stringstream ss2;
				ss2 << line;   int tabCnt = 1;
				while (getline(ss2, segs, '\t')) {
					tabCnt++;
				}
				if (tabCnt >= 3) {
					UPARSE9up = true;
					cerr << "Switching to Uparse 10+ style map file.\n";
				}
			}
			UpUcFnd = true;
		}
		
		stringstream ss;
		ss << line;
		bool chimera = false;
		vector<string>tarsV; //saves hits to OTUs
		//2 ways to get to a) hit info b) query & otu
		if (!UPARSE8up){ //uparse 7
			if ( (line.substr(0, 1) != "H")) {
				continue;
			}
			for (uint i = 0; i < 4; i++){//jump to pos X
				getline(ss, segs, '\t');
			}
			perID = (float)atof(segs.c_str());
			for (uint i = 0; i < 5; i++){//jump to pos X
				getline(ss, segs, '\t');
			}
			getline(ss, segs2, '\t');
		} else if (!UPARSE9up){ // uparse 8
			//query first entry
			string tmp;
			getline(ss, segs, '\t');//0
			getline(ss, tmp, '\t');//1
			//should be "match"
			if ( tmp == "chimera") {
				if ( !doChimeraCnt ) {continue;}
				chimera = true; }// segs = ""; continue;}
			if (tmp == "otu"){
				segs2 = segs;
				perID = 100.f;
			} else {
				getline(ss, tmp, '\t');//2
				perID = (float)atof(tmp.c_str());
				//indicator if hit
				//OTU last entry
				getline(ss, segs2, '\t');//3
				getline(ss, segs2, '\t');//4
			}
		} else { //UP9, uparse 10, uparse 11, dada2 fake .uc
				 //query first entry
			string tmp;
			getline(ss, segs, '\t');//0
			getline(ss, tmp, '\t');//1
								   //should be "match"
			if (tmp == "chimera" ) {
				continue;
			} else if ( tmp == "noisy_chimera" ||  tmp == "good_chimera") { //tmp == "perfect_chimera" ||
				if (!doChimeraCnt) { continue; }
				chimera = true;
			}else if (tmp == "perfect_chimera") {
				removeSizeStr(segs);
				perfectChims.insert(segs);
				continue;
			}// segs = ""; continue;}
			if (tmp.substr(0,3) == otuTerm ){
				segs2 = segs;
				perID = 100.f;
			} 
			else {//match or perfect_chimera case
				getline(ss, tmp, '\t');//2
				//dqt=1;top=GZV0ATA01ANJXZ;size=14;(99.6%);
				size_t p1(tmp.find("top=")+4);//4
				size_t p2(tmp.find(";(", p1)+2);
				size_t p3(tmp.find("%);", p2));
				segs2 = tmp.substr(p1, p2 - p1-1);
				//string xx = tmp.substr(p2, p3 - p2);
				perID = (float)atof(tmp.substr(p2,p3-p2).c_str());
				if (false && chimera) {//just use up9 top hit
//					p1(tmp.find(";top=") + 5,p2);
//					p2(tmp.find(";(", p1) + 2);

				}
			}
		}
		//remove spaces
		segs = segs.substr(0,segs.find_first_of(' '));
		//also remove sample identifier in string
		string smplID = "";
		removeSampleID(segs, SEP, smplID);

		if ( chimera && UPARSE8up) {
			tarsV = splitByComma(segs2, false, '+');
			for ( uint kk = 0; kk < tarsV.size(); kk++ ) {
				tarsV[kk] = tarsV[kk].substr(0,tarsV[kk].find_last_of("("));
			}
		} else if (!chimera){
			tarsV.push_back(segs2);
		}
		matrixUnit splChim = (matrixUnit)tarsV.size();
		curCLID.resize(tarsV.size());

		for ( uint kk = 0; kk < tarsV.size(); kk++ ) {
			//goes over all chimeric hits
			string oriClKey (tarsV[kk]);
			segs2 = oriClKey;
			//remove ;size= argument
			removeSizeStr(segs2);
			removeSampleID(segs2, SEP);
			itCL = seq2CI.find(segs2);

			if ( itCL == seq2CI.end() ) {//new OTU?
				if (OTUnumFixed) {
					//cerr << "XX\n";
					if (perfectChims.find(segs2) == perfectChims.end()) {
						cerr << "Unkown OTU entry found:" << segs2 << endl;
					}
				}	else {
					//not found in known clusters.. create entry
					bestDNA.push_back(NULL);
					bestDNA2.push_back(NULL);
					oriKey.push_back(oriClKey);
					bestPID.push_back(0.f);
					bestLEN.push_back(0);
					clusCnt = (int)bestDNA.size() - 1;
					curCLID[kk] = clusCnt;
					seq2CI[segs2] = clusCnt;
				}
			} else {//cluster exists
				curCLID[kk] = ((*itCL).second);
			}
			bestDNA[curCLID[kk]]->totalSum();
#ifdef matrix_sum
			if ( addFromHDstring ) {
				//cerr << segs << "\t" << segs2 << endl;
				add2OTUmat(smplID, (*itCL).second, splChim);
			} 
#endif
		}
		return true;
	}
	return true;
}
void UClinks::writeOTUmatrix(string outf) {
	cerr << "Writing OTU matrix to " << outf << endl;
	this->setOTUnms();
	ofstream MA;
	MA.open(outf);
	//first write all sample IDs
	std::unordered_map<string, int>::iterator OTUid;
	MA << "OTU";
	matrixUnit totalCnts = 0;
	
	vector<string> smpls(SmplIDs.size(),"");
	unordered_map<string, int>::iterator smplNit;
	for (smplNit = SmplIDs.begin(); smplNit != SmplIDs.end(); smplNit++) {
		smpls[smplNit->second] = smplNit->first;
	}
	
	for (size_t i = 0; i < smpls.size(); i++) {
		MA<<"\t"<<smpls[i];
	}
	//write OTU values
	for (OTUid = seq2CI.begin(); OTUid != seq2CI.end(); OTUid++) {
		MA << endl << oriKey[OTUid->second];
		int rowI = OTUid->second;
		for (size_t i = 0; i < smpls.size(); i++) {
			matrixUnit val = OTUmat[ SmplIDs[smpls[i]] ][rowI];
			totalCnts += val;
			if ( isnan(val) ) {
				MA << "\t0";
			} else {
				MA << "\t"<<round(val);
			}
		}
	}
	MA << endl;
	MA.close();
	cerr << "Recruited " << (int)totalCnts << " reads in OTU matrix\n";
}
void UClinks::setOTUnms(){
	int cnt = 0;
	std::unordered_map<string, int>::iterator OTUid;
	for (OTUid = seq2CI.begin(); OTUid != seq2CI.end(); OTUid++) {
		string newlySetID = otuOUTterm + itos(cnt);
		/*if (newlySetID == "OTU_7711"){
			int x = 0;	cerr << "7711 id_ " << OTUid->second; if (bestDNA[OTUid->second] == NULL){ cerr << "no DNA\n"; }
			else { cerr << "YES\n"; }
		}*/
		oriKey[OTUid->second] = newlySetID;
		cnt++;
	}

}
void UClinks::add2OTUmat(const string& smplID, int curCLID, matrixUnit spl) {

	std::unordered_map<string, int>::iterator smplNit = SmplIDs.find(smplID);
	if (smplNit == SmplIDs.end()) {//create entry & expand matrix
		SmplIDs[smplID] = (int) OTUmat.size();
		OTUmat.push_back(vector<matrixUnit>(clusCnt + 1, (matrixUnit)0));
		smplNit = SmplIDs.find(smplID);
		cerr << "New sample_id_ id_ in uc file detected, that is not present in map: " << smplID<< endl;
		unregistered_samples = true;
	}
	//easy, now add in
	if ( spl <1 ) { spl = 1; }
	OTUmat[(*smplNit).second][curCLID] += (matrixUnit)1 / spl;

}
void UClinks::add2OTUmat(shared_ptr<DNAunique> d, int curCLID, matrixUnit rep) {
	if (d == NULL) { cerr << " add2OTUmat::d is NULL\n"; exit(85); }
	read_occ map = d->getDerepMap();
	if ( rep <1 ) { rep = 1; }
	for (auto iterator = map.begin(); iterator != map.end(); iterator++) {
		//assert(OTUmat.size() <= iterator->first);
		OTUmat[iterator->first][curCLID] += (matrixUnit)iterator->second / rep;
	}
}

void UClinks::setupDefSeeds(shared_ptr<InputStreamer> FA, const vector<string>& smpls) {
	bool contRead = true; //bool sync(false);
	while (contRead) {
		shared_ptr<DNA> tmpDNA = FA->getDNA( 0);
		if (tmpDNA == NULL) { contRead = false; break; }
		shared_ptr<DNAunique> tmp = make_shared<DNAunique>(tmpDNA, -1);
//		delete tmpDNA;
		//second pair_
		//shared_ptr<DNA> tmp2 = FA->getDNA(contRead, 1);
		string segs2 = tmp->getId();
		string oriClKey = segs2;
		//remove ;size= argument
		removeSizeStr(segs2);
		//also remove sample identifier in string
		removeSampleID(segs2, SEP);

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			bestDNA.push_back(tmp);
			bestDNA2.push_back(NULL);
			oriKey.push_back(oriClKey);
			bestPID.push_back(100.f);
			bestLEN.push_back(tmp->length());
			clusCnt = (int)bestDNA.size() - 1;
			seq2CI[segs2] = clusCnt;
		} else {//cluster exists
			cerr << "Found double id_ " << segs2 << endl;//
			if (false) {
				cerr << "Aborting.." << endl;
				exit(74);
			}
		}
	}
	//std::map<string, int>::iterator smplNit;
	if (derepMapFile != "") {
		readDerepInfo(derepMapFile);
	} else {
//		const vector<string>& smpls = fil->SampleID;
		for (uint i = 0; i < smpls.size(); i++) {
			SmplIDs[smpls[i]] = (int)OTUmat.size();
			if (i != OTUmat.size()) { cerr << "Err in setupDefSeeds\n"; exit(453); }
			OTUmat.push_back(vector<matrixUnit>(clusCnt + 1, (matrixUnit)0));
		}
	}
}

void UClinks::addDefSeeds(shared_ptr<InputStreamer> FA, Filters* fil) {
	bool contRead = true;
	int addCnt = 0; //bool sync(false);
	while (contRead) {
		shared_ptr<DNA> tmpDNA = FA->getDNA( 0);
		if (tmpDNA == NULL) { contRead = false; break; }
		shared_ptr<DNAunique> tmp = make_shared<DNAunique>(tmpDNA, -1);
		//delete tmpDNA;
		//second pair_
		//shared_ptr<DNA> tmp2 = FA->getDNA(contRead, 1);
		string segs2 = tmp->getId();
		string oriClKey = segs2;
		//remove ;size= argument
		removeSizeStr(segs2);
		removeSampleID(segs2, SEP);

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			bestDNA.push_back(tmp);
			bestDNA2.push_back(NULL);
			oriKey.push_back(oriClKey);
			bestPID.push_back(100.f);
			bestLEN.push_back(tmp->length());
			clusCnt = (int)bestDNA.size() - 1;
			seq2CI[segs2] = clusCnt;
		}
		else {//cluster exists
			cerr << "Found double id_ " << segs2 << endl << "Aborting.." << endl;
			exit(74);
		}
		addCnt++;
	}
	int siz = (int)bestDNA.size();
	//resize complete matrix to take these up
	
	for (uint i = 0; i < OTUmat.size(); i++){
		OTUmat[i].resize(siz, (matrixUnit)0);
	}
}

void UClinks::besterDNA(const vector<int>& curCLIDpre, shared_ptr<DNAunique> tdn1, 
	shared_ptr<DNA> tdn2, Filters* fil) {
	bool checkBC = true; 
	int TagIdx(-2);
	if (tdn2 != NULL) {//fix for pairs assuming midSeqs
		checkBC = false; //TagIdx = 0;
	}
	matrixUnit matchSiz = (matrixUnit)curCLIDpre.size();
	if ( matchSiz>1 ) {
		for ( int k = 0; k < matchSiz; k++ ) {
			add2OTUmat(tdn1, curCLIDpre[k], matchSiz);
		}
		return;
	}
	int curCLID = (int) curCLIDpre[0];
	add2OTUmat(tdn1, curCLID,1);
	if (SeedsAreWritten){
		//delete dnaTemp1; if (dnaTemp2 != NULL){ delete dnaTemp2; }
		return;
	}
	//general routine, matching tdn found
	if (bestDNA[curCLID] == NULL) { //just fill with current DNA
		if (tdn2 == NULL) {
			if ( fil->doReversePrimers() && !fil->check(tdn1, true, CurSetPair, TagIdx) ) {
				return;//delete dnaTemp1;
			}
			bestDNA[curCLID] = tdn1;
			bestPID[curCLID] = tdn1->getTempFloat();
			bestLEN[curCLID] = tdn1->length();
		} else {
			if (fil->doReversePrimers() && !fil->check(tdn1, true, CurSetPair, TagIdx)) {
				 return;//delete dnaTemp2; delete dnaTemp1;
			}
			bestDNA[curCLID] = tdn1;
			bestPID[curCLID] = tdn1->getTempFloat();
			bestLEN[curCLID] = tdn1->length() + tdn2->length();
			fil->check(tdn2, true, 1, TagIdx);
			bestDNA2[curCLID] = tdn2;
		}
	}//already a candidate sequence? check who is better..
	else if (
		fil->betterSeed(tdn1, tdn2, bestDNA[curCLID], bestDNA2[curCLID], bestPID[curCLID], bestLEN[curCLID], CurSetPair, checkBC)
		){
//		delete bestDNA[curCLID];
		bestDNA[curCLID] = tdn1; 
		if (tdn2 != NULL) {
//			if (bestDNA2[curCLID] != NULL) {delete bestDNA2[curCLID];	}
			bestDNA2[curCLID] = tdn2;
		}
		if (bestPID[curCLID] < tdn1->getTempFloat()){
			bestPID[curCLID] = tdn1->getTempFloat();
		}
		uint curL = tdn1->length();
		if (tdn2 != NULL) { curL += tdn2->length(); }
		if (bestLEN[curCLID] < curL){
			bestLEN[curCLID] = curL;
		}
	} else {
//		delete dnaTemp1;delete dnaTemp2;
	}
}

void UClinks::removeSampleID(string& w, const string &SEP) {
	size_t pos = w.find(SEP);
	if (pos != std::string::npos) {
		w = w.substr(pos + SEP.length());
	}
}
void UClinks::removeSampleID(string& w, const string &SEP, string & SMplID) {
	size_t pos = w.find(SEP);
	if (pos != std::string::npos) {
		SMplID = w.substr(0,pos);
		w = w.substr(pos + SEP.length());
	}
}
void UClinks::removeSizeStr(string& w) {
	size_t idx = w.find(";size=");
	/*#ifdef matrix_sum   //not required, sanity check that is not really working out with usearch mappings
	int OTUsize = atoi( segs2.substr(idx+6).substr(0,-1).c_str() );
	#endif*/
	w = w.substr(0, idx);
}

void UClinks::writeNewSeeds(shared_ptr<OutputStreamer> MD, Filters* fil, 
			bool refSeeds, bool printLnk) {
	if (!RefDBmode && refSeeds){ return; }
	//ofstream O(outf.c_str());
	int paired = fil->isPaired();
	//MD->printStorage();
	ofstream links;
	uint st(0); uint to((uint)oriKey.size());
	//in case of refDB mode, need to start from point where refDBs were added in..
	if (refSeeds && RefDBmode && RefDBotuStart >= 0){
		cerr << "Writing ref DB sequences..";
		st = RefDBotuStart;
		//and also write a link between OTU_name and refSeq
		if (printLnk){
			string refLinkF = MD->leadOutFile() + ".lnks";
			links.open(refLinkF.c_str(), ios::out);
		}
	} else if (!refSeeds && RefDBmode && RefDBotuStart >= 0){
		cerr << "Writing new OTU seeds..";
		to = RefDBotuStart;
	}
	shared_ptr<DNA> d;
	
	for (uint i = st; i < to; i++) {
		//if (i == 6628){			int x = 0;		}
		//if check DNA was done before (remove rev primer), do it now
		if (bestDNA[i] == NULL) {
			//this should not happen on a regular lotus run
			MD->closeOutStreams();
			cerr << "No seed sequence found for DNA " << oriKey[i] << ". Aborting..\n"; // " << bestDNA[i]->getOldId()<<"(
			exit(54);
		}
		//fil->check(bestDNA[i],true);
		string newH = oriKey[i];
		d = bestDNA[i];
		if (printLnk){
			string oriH = d->getShortId(); removeSizeStr(oriH);
			links << newH << "\t" << oriH << endl;
		}
		//takeOver id_ to allow for later linkup with cluster
		if (paired == 2) {
			//first check for readmerging
			d->setPassed(true);
			d->setNewID(newH + ".1");

			if (b_merge_pairs_optiSeed_ && bestDNA2[i] != NULL) {
				MD->findSeedForMerge(d, bestDNA2[i],0);
				bool didMerge(false);
				didMerge = MD->saveForWrite_merge(d, bestDNA2[i], newH , 0, true);
				if (!didMerge) {//not merged? we want to add this DNA nonetheless to output
					//better to do this in saveForWrite already
				}
			}

			if (bestDNA2[i] != NULL) {
			    //sorted by pair1,2, quality yellow green, singleton (pair 1 2 3 4) etc
				MD->saveForWrite(d, 1,-1);
				bestDNA2[i]->setPassed(true);
				bestDNA2[i]->setNewID(newH + ".2");
				MD->saveForWrite(bestDNA2[i], 2,-1);
			} else {
				MD->saveForWrite(d, 3,-1);
			}
		} else {
			d->setNewID(newH);
			d->setPassed(true);
			MD->saveForWrite(d, 1,-1);
		}
	}
	MD->closeOutStreams();
	SeedsAreWritten = true;
	if (printLnk){ links.close(); }
	cerr << "Done" << endl;
}
void UClinks::printStats(ostream& os){
	if (oriKey.size() == 0) {
		os << "No OTU's to re-seed\n";
		return;
	}
	uint numCls = 0;
	float avgQ(0.f), minQ(10000.f), maxQ(0.f),
		avgA(0.f), minA(10000.f), maxA(0.f),
		avgS(0.f), minS(10000.f), maxS(0.f);
	uint avgL = 0, minL(10000), maxL(0);
	
	//uint div(0);// = (uint)oriKey.size();
	vector<uint> lengths;
	vector<float> quals,accums,sims;
	//oriKey.size() - not if there are only refs
	uint to((uint)oriKey.size());
	if (RefDBmode && RefDBotuStart >= 0){
		to = RefDBotuStart;
	}

	for (uint i = 0; i < to; i++){
		if (bestDNA[i] == NULL){ continue; }
		float curQ = bestDNA[i]->getAvgQual();
		if (bestDNA2[i] != NULL) {
			curQ += bestDNA2[i]->getAvgQual(); curQ /= 2.f;
		}
		uint curL = (uint) bestDNA[i]->length();
		if (bestDNA2[i] != NULL) {
			curL += bestDNA2[i]->length(); 
		}

		if (curL < minL){ minL = curL; }
		if (curL > maxL){ maxL = curL; }
		avgL += curL;
		lengths.push_back(curL);
		
		if (curQ < 1) {//no new Seed found, default seed
			continue;
		}
		float sc = bestDNA[i]->getTempFloat();

		if (curQ < minQ){ minQ = curQ; }
		if (curQ > maxQ){ maxQ = curQ; }
		avgQ += curQ;

		float curA = (float)bestDNA[i]->getAccumError();
		if (bestDNA2[i] != NULL) {
			curA += (float) bestDNA2[i]->getAccumError(); 
		}

		quals.push_back(curQ);
		accums.push_back(curA); 
		if (curA < minA){ minA = curA; }
		if (curA > maxA){ maxA = curA; }
		avgA += curA;	

		sims.push_back(sc);
		if (sc < minS){ minS = sc; }
		if (sc > maxS){ maxS = sc; }
		avgS += sc;
		numCls++; 
	}
	os << "Found " << numCls << " seeds of " << oriKey.size() << " OTU's in " << uclines << " mappings." << endl;

	//sort vectors
	std::sort(lengths.begin(), lengths.end());
	std::sort(quals.begin(), quals.end());
	std::sort(accums.begin(), accums.end());
	std::sort(sims.begin(), sims.end());
	os << "Stats of Seed sequences (0th/10th/50th/90th/100th) percentile:\n";
	if (lengths.size() > 0) { os << "\n     - sequence_ Length :   " << minL << "/" << calc_median2(lengths, 0.1f) << "/" << calc_median2(lengths, 0.5f) << "/" << calc_median2(lengths, 0.9f) << "/" << maxL; }
	if (quals.size() > 0) { os << "\n     - Quality :      " << minQ << "/" << calc_median2(quals, 0.1f) << "/" << calc_median2(quals, 0.5f) << "/" << calc_median2(quals, 0.9f) << "/" << maxQ; }
	if (accums.size() > 0) { os << "\n     - Accum. Error : " << minA << "/" << calc_median2(accums, 0.1f) << "/" << calc_median2(accums, 0.5f) << "/" << calc_median2(accums, 0.9f) << "/" << maxA;	}
	if (sims.size() > 0) {	os << "\n     - Sim2Consensus: " << minS << "/" << calc_median2(sims, 0.1f) << "/" << calc_median2(sims, 0.5f) << "/" << calc_median2(sims, 0.9f) << "/" << maxS;		}
	os << endl;
}
