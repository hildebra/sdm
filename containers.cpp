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

void ReadSubset::findMatches(shared_ptr<InputStreamer> IS, OutputStreamer* MD, bool mocatFix) {
	vector<shared_ptr<DNA>> match; //shared_ptr<DNA> match2(NULL);
	bool cont(true), cont2(true);
	int idx(0);
	int pairs = IS->pairNum();
	unordered_map <string, int>::iterator SEEK;
	bool b_doHD = newHD.size() > 0;
	//bool sync(false);//meaningless placeholder
	while (cont) {
		match = IS->getDNApairs();
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





OutputStreamer::OutputStreamer(Filters* fil, OptContainer* cmdArgs,
        std::ios_base::openmode writeStatus, 
		shared_ptr<ReadSubset> RDSset,
		int numThreads, string fileExt, int forceFmt) :
        MFil(nullptr), subFilter(0), //DNAsP1(0), DNAsP2(0), DNAsS1(0), DNAsS2(0),
        //DNAsNoHead(0), DNAsP1_alt(0), DNAsP2_alt(0), DNAsS1_alt(0), DNAsS2_alt(0),
		//locBufGreen(2,""),
		suppressOutWrite(0), write2File(true),// mem_used(false),
        DNAinMem(0), writeThreadStatus(0),
        fastQver(33),
        fastQoutVer(33), BWriteQual(false),
        BWriteFastQ(false), b_multiOutStream(false), pairedSeq(-1),
		b_changeFQheadVer(false),
		b_checkedHeaderChange(false),
        b_oneLinerFasta(false), b_doDereplicate(false), b_writeGreenQual(true), b_writeYellowQual(true),
        maxReadsPerOFile(0), 
		demultiBPperSR(0), 
		BPwrittenInSR(0), BPwrittenInSRmerg(0),
		ReadsWritten(0),
        maxRdsOut(-1), stopAll(false),
        leadingOutf(""), locCmdArgs(cmdArgs), dereplicator(nullptr), cntDerep(0), wrMode(ios::out),
        sFile(0), qFile(0), fqFile(0),
		sPairFile(0), qPairFile(0), fqPairFile(0),
		of_merged_fq(0),
        //sFileStr(0), qFileStr(0), fqFileStr(0), fqNoBCFile(0), 
		totalFileStrms(0),doTIO(true),
        bDoDemultiplexIntoFiles(false), demultiSinglFiles(0), //demultiSinglFilesF(0),
		demultiMergeFiles(0),
        onlyCompletePairsDemulti(false),
		b_merge_pairs_demulti_(false), b_merge_pairs_filter_(false),
		b_merge_pairs_(false), mergers(0,nullptr),
		Nthrds(numThreads),_benchmark(nullptr)
{
	MFil = fil;
	
	fastQver = fil->getuserReqFastqVer();
	fastQoutVer = fil->getuserReqFastqOutVer();

	maxReadsPerOFile = fil->maxReadsOutput();
	demultiBPperSR = fil->getDemultiBPperSR();
	ReadsWritten = fil->writtenReads();

 	pairedSeq = MFil->isPaired();
	if (cmdArgs->find("-suppressOutput") != cmdArgs->end()) {
		suppressOutWrite = atoi( (*cmdArgs)["-suppressOutput"].c_str() );
	}
	if (cmdArgs->find("-pairedDemulti") != cmdArgs->end() && (*cmdArgs)["-pairedDemulti"] == "1") { //only demultiplex proper pairs
		onlyCompletePairsDemulti = true;
	}
	if (cmdArgs->find("-oneLineFastaFormat") != cmdArgs->end() && (*cmdArgs)["-oneLineFastaFormat"] == "1") {
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
	if ((*cmdArgs)["-o_fastq_noBC"] != "") {
		openNoBCoutstrean((*cmdArgs)["-o_fastq_noBC"]);
	}
	if (cmdArgs->find("-o_demultiplex") != cmdArgs->end()) {//demulti: always write mode out
		bool gzDemulti = false;
		if (cmdArgs->find("-o_demultiplex_gz") != cmdArgs->end()) {
			if ((*cmdArgs)["-o_demultiplex_gz"] == "1") {
				gzDemulti = true;
				//cerr << "Gzip demulti output\n";
			}
			else {
				//cerr << "\n\nNo Gzip demulti output\n";
			}
		}
		generateDemultiOutFiles((*cmdArgs)["-o_demultiplex"],fil, ios::out, gzDemulti);
	}
	if (cmdArgs->find("-merge_pairs_filter") != cmdArgs->end()
		&& (*cmdArgs)["-merge_pairs_filter"] == "1") {
		b_merge_pairs_filter_ = true;
	}
	if (cmdArgs->find("-merge_pairs_demulti") != cmdArgs->end() 
		&& (*cmdArgs)["-merge_pairs_demulti"] == "1") {
		b_merge_pairs_demulti_ = true;
	}
	if (b_merge_pairs_demulti_ || b_merge_pairs_filter_) {
		b_merge_pairs_ = true;
	}
	setSubfilters(Nthrds);

	if (this->pairedSeq != 2) {
		onlyCompletePairsDemulti = false;
	}
	if ((*cmdArgs)["-threadIO"] == "1") { doTIO = true; }
	if ((*cmdArgs)["-threadIO"] == "0") { doTIO = false; }


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
		if (subFilter[i] != nullptr) { delete subFilter[i]; }
	}
	for (size_t x = 0; x < mergers.size(); x++) {
		if (mergers[x] != nullptr) { delete mergers[x]; mergers[x] = nullptr; }
	}
	dmltMTX.unlock();
	cdbg("Subfilters deleted, streams closed\n" );
	//for (size_t x = 0; x < merger.size(); x++) { delete merger[x]; } merger.clear(); 
	//delete optim;
}
void OutputStreamer::delAllDNAvectors() {
#ifdef DEBUG
	cerr << "cleaning MD..";
#endif
	//locBufGreen[0].clear(); locBufGreen[1].clear();
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
	//s/q/fqPairFile  destructor..
	for (size_t i = 0; i < fqPairFile.size(); i++) {
		if (fqPairFile[i] != nullptr) {
			delete fqPairFile[i];
		}
	}
	for (size_t i = 0; i < qPairFile.size(); i++) {
		if (qPairFile[i] != nullptr) {
			delete qPairFile[i];
		}
	}
	for (size_t i = 0; i < sPairFile.size(); i++) {
		if (sPairFile[i] != nullptr) {
			delete sPairFile[i];
		}
	}

}

void OutputStreamer::analyzeDNA(shared_ptr<DNA> d, int FilterUse, int pair, int& idx, int thr) {
	if (!d){
		return ;
	}
	if (d->isFailed()) { return; }
	Filters * curFil = this->getFilters(thr);
	bool isP1 = max(0, pair) == 0;

	if (isP1 && curFil->getcut5PR1()) {
		d->cutSeq(0,curFil->getcut5PR1());

	} else if (!isP1 && curFil->getcut5PR2()) {
		d->cutSeq(0, curFil->getcut5PR2());
	}

	if ( !curFil->doFilterAtAll() ) {
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

	//pricipally safe to call from thread
    curFil->checkYellowAndGreen(d, pair, idx, false);

	//bit redundant.. but better safe than sorry
	if (!curFil->secondaryOutput() && d->isYellowQual() && ! d->isGreenQual()) {
		d->failed();
	}
	if (idx == -1 ) {
		d->setBarcodeDetected(false); 
	}
	return;

	/*
	if (curFil->secondaryOutput() ){
		//make this primary routine..

	} else if (FilterUse==-1){
		curFil->check(d, false, pair, idx);
	} else {
		cerr << "Invalid control path in analyzeDNA\n"; exit(55);
		//Filters* f = subFilter[FilterUse];
		//subFilter[FilterUse].check(d, false, pair_, idx);
        //f->check(d, false, pair, idx);
    }
	*/

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
	
	bool demultFini = demultiBPperSR > 0 && (uint)BPwrittenInSR > demultiBPperSR;
	bool demultMrgFini = demultiBPperSR > 0 && (uint) BPwrittenInSRmerg > demultiBPperSR;
	if (demultFini && demultMrgFini) {
		return;
	}
	int idx = d1->getBCnumber() - BCoffset; //correct for BC offset as well..

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
	if (this->isPEseq() == 2 && 
		(green1 || green2) && b_merge_pairs_demulti_ && d1->merge_seed_pos_ > 0) {
		shared_ptr<DNA> dna_merged = mergers[curThread]->merge(d1, d2);
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
	if ((onlyCompletePairsDemulti && (!green1 || !green2) )) {
		return;
	}
	/*
	string x = d2->getShortId();
	if (x.find("M04428:252:000000000-BM3M5:1:2106:11089:1656") != string::npos) {
		bool x = true; // DEBUG
	}
	/**/
	string dn1(""), dn2("");
	if (green1) {
		//cout << "trouble" << endl;
		d1->prepareWrite(fastQoutVer);
		//dmltMTX.lock();
		BPwrittenInSR += d1->length();
		//dmltMTX.unlock();
		d1->writeFastQ(dn1,false);
	}
	if (green2) {
		//cout << "trouble" << endl;
		d2->prepareWrite(fastQoutVer);
		//dmltMTX.lock();
		//dmltMTX.unlock();
		d2->writeFastQ(dn2,false);
	}
	//cerr << dn2 << "WEW";
	//move them closer together..
	dmltMTX.lock();
	*(demultiSinglFiles[idx][0]) << dn1;
	//cdbg("write2Demulti4"+itos(idx) + "X "+itos(demultiSinglFiles[idx].size())+"Y"+ dn2);
	*(demultiSinglFiles[idx][1]) << dn2;
	dmltMTX.unlock();


}

void OutputStreamer::write2Demulti(shared_ptr<DNA> d, int p, int BCoffset) {//this->getBCoffset()
	if (!this->Demulti2Fls()) {
		return;
	}
	if (demultiBPperSR > 0 && (uint)BPwrittenInSR > demultiBPperSR) {
		return;
	}
	int idx = d->getBCnumber() - BCoffset; //correct for BC offset as well..

	if (idx < 0 || !d->isGreenQual()) {
		return;
	}

	d->prepareWrite(fastQoutVer);
	*(demultiSinglFiles[idx][p]) << d->writeFastQ(false);
	BPwrittenInSR += d->length();

}


void OutputStreamer::generateDemultiOutFiles(string path, Filters* fil, 
	std::ios_base::openmode writeStatus, bool demulti2gz) {
	//fill in demultiSinglFiles vector
	vector<ofbufstream*> empVec(2, NULL);
	//vector<string> empVec2(2, "");

	bool doMC =  Nthrds > 1;
	string gzSuffiz = "";
	string gzReport = "";
	if (demulti2gz) {
		gzSuffiz = ".gz";
		gzReport = "gzipped ";
	}

	bool displayInfo(false);
	if (fil->getReadsWritten() == 0) {
		displayInfo = true;
	}

	struct stat info;
	if (stat(path.c_str(), &info) != 0) {
		cerr << "Output path for demultiplexed files does not exist, please create this directory:\n" << path << endl;
		exit(833);
	}
	else if (info.st_mode & S_IFDIR) { // S_ISDIR() doesn't exist on my windows 
		if (displayInfo) {
			cerr << "\nWriting " << gzReport << "demultiplexed files to : " << path << endl;// printf(" % s is a directory\n", path);
		}
	}
	else {
		cerr << path << " is no directory\n"; exit(834);
	}
	path += "/";

	bDoDemultiplexIntoFiles = true;
	bool openOstreams = true; uint ostrCnt(0);
	size_t bufS = OUTPUT_BUFFER_SIZE;
	//save as .gz? -> no for dada2..
	


	demultiSinglFiles.resize(fil->SampleID.size(), empVec);
	demultiMergeFiles.resize(fil->SampleID.size(), nullptr);
	//demultiSinglFilesF.resize(fil->SampleID.size(), empVec2);
	for (size_t i = 0; i < fil->SampleID.size(); i++) {
		//actually needs to know if paired files..
		//if (ostrCnt > maxFileStreams) {openOstreams = false;}
		if (pairedSeq == 1 || pairedSeq == -1) {
			string nfile = path + fil->SampleID[i] + ".fq" + gzSuffiz;
			if (openOstreams) { 
				if (demultiSinglFiles[i][0] != nullptr) { delete demultiSinglFiles[i][0]; }
				demultiSinglFiles[i][0] = DBG_NEW ofbufstream(nfile.c_str(), writeStatus, doMC, (size_t)bufS); 
			}
			//demultiSinglFilesF[i][0] = nfile;
			ostrCnt++;
		}
		else {
			string nfile = path + fil->SampleID[i] + ".1.fq" + gzSuffiz;
			if (openOstreams) { 
				if (demultiSinglFiles[i][0] != nullptr) { delete demultiSinglFiles[i][0]; }
				demultiSinglFiles[i][0] = DBG_NEW ofbufstream(nfile.c_str(), writeStatus, doMC, size_t(bufS * 0.8));
			}
			//demultiSinglFilesF[i][0] = nfile;
			nfile = path + fil->SampleID[i] + ".2.fq" + gzSuffiz;
			if (openOstreams) { 
				if (demultiSinglFiles[i][1] != nullptr) { delete demultiSinglFiles[i][1]; }
				demultiSinglFiles[i][1] = DBG_NEW ofbufstream(nfile.c_str(), writeStatus, doMC, size_t(bufS * 1.2)); }
			//demultiSinglFilesF[i][1] = nfile;
			nfile = path + fil->SampleID[i] + ".merg.fq" + gzSuffiz;
			if (openOstreams) { 
				if (demultiMergeFiles[i] != nullptr) { delete demultiMergeFiles[i]; }
				demultiMergeFiles[i] = DBG_NEW ofbufstream(nfile.c_str(), writeStatus, doMC, (size_t)bufS);
			}
			
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
//    if (d->QualCtrl.PrimerFwdFail) {
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
////        countBCdetected(d->getBCnumber(), easyPair, false);
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
//            this->statAddDerepBadSeq(d->getBCnumber());
//        }
//    }
//}

void Filters::addDNAtoCStats(shared_ptr<DNA> d,int Pair) {
	//here should be the only place to count Barcodes!
	int easyPair = Pair < 3 ? Pair - 1 : Pair - 3;
	
	//csMTX[easyPair]->lock();
	collectStatistics[easyPair]->total2++;


	if (d->isGreenQual() || d->isYellowQual()) {
		this->DNAstatLQ(d, easyPair, d->isYellowQual());
		collectStatistics[easyPair]->totalSuccess++;
		if (d->isYellowQual()) {
			collectStatistics[easyPair]->totalMid++;
		}
	} else {
		collectStatistics[easyPair]->totalRejected++;
	}

	//some general stats that always apply:
	if (d->QualCtrl.PrimerFwdFail) {
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
	
	if (d->isGreenQual() ){//|| d->isYellowQual()) {
		countBCdetected(d->getBCnumber(), easyPair, false);
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
		if (d->QualCtrl.HomoNTtrimmed) {
			collectStatistics[easyPair]->HomoNTtrimmed++;
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
			this->statAddDerepBadSeq(d->getBCnumber());
		}
	}
	//csMTX[easyPair]->unlock();
}


/*
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
    if (d->QualCtrl.PrimerFwdFail) {
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
        countBCdetected(d->getBCnumber(), easyPair, false);
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
            this->statAddDerepBadSeq(d->getBCnumber());
        }
    }
}

*/

bool  OutputStreamer::saveForWrite(shared_ptr<DNA> d,int Pair, int thr,int& Cstream, bool write) {
	//second most important part: collect stats on DNA passing through here (should be all read)
	//most important part: save DNA to be written later (or discard)
	if (d == NULL || stopAll) {
		Cstream = 100;
		return !stopAll;
	}
	//int Cstream = 0;
	//bool writen(false);
	
	//threadsafe
	Filters* curFil = this->getFilters(thr);

	curFil->addDNAtoCStats(d, Pair);

	if (suppressOutWrite == 3) {
		Cstream = 100; //set to impossibly high value
		return !stopAll;
	}

	if (curFil->doFilterAtAll()) {
		if ((b_writeGreenQual && d->isGreenQual() ) || 
			(b_writeYellowQual && d->isYellowQual()) ) {
			if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
				d->changeHeadPver(curFil->FQheadV());
			}
			if (d->isYellowQual()) {
				Cstream = 1;
			} else { Cstream = 0; }
		} else {
			Cstream = 100;
			return !stopAll;
		}
	} else {
		Cstream = 0;
	}
	d->prepareWrite(fastQoutVer);

	/*if (curFil->passedReads(25000)) {
		if (pairedSeq > 1) {
			cerr << "25k read pairs/" << _benchmark->now_interval_secs(2) << "s\r";
		} else {
			cerr << "25k reads/" << _benchmark->now_interval_secs(2) << "s\r";
		}
		//cerr.flush();
	}*/

	//sqfqostrMTX.lock();
	if (write) {
			writeForWrite(d, Pair, Cstream, nullptr, -1, 100);
	}
	return !stopAll;

	
}

void OutputStreamer::writeForWrite(shared_ptr<DNA> d1, int Pair1, int Cstream1,
	shared_ptr<DNA> d2, int Pair2, int Cstream2) {
	
	//Cstream 100 means not to write DNA
	if (Cstream1 >= 100 && Cstream2 >= 100) {
		//both reads not written..
		return;
	}
	if (maxRdsOut > 0 &&  ReadsWritten >= maxRdsOut) {
		return;
	}
	 
	
	//incrementing output file number
	if (maxReadsPerOFile > 0 && ReadsWritten + DNAinMem >= maxReadsPerOFile) {
		//cerr << "ReadsWritten " << ReadsWritten << " DNAinMem " << DNAinMem << endl;
		DNAinMem = 0;
		sqfqostrMTX.lock();
		incrementOutputFile();
		sqfqostrMTX.unlock();
	}
	bool fq1cool = Cstream1 < (int)fqFile.size();
	bool fq2cool = Cstream2 < (int)fqFile.size();
	if (BWriteFastQ && !fq1cool && !fq2cool) {
		return;
	}
	//string str1(""); if (fq1cool) { str1 = d1->writeFastQ();}
	//string str2(""); if (fq2cool) { str2 = d2->writeFastQ(); }
	
	/*if (Cstream1 == 0 && Pair1 < 2 && BWriteFastQ) {
		sqfqostrMTX.lock();//lock to ensure read pairs written together
		locBufGreen[0] += str1; locBufGreen[1] += str2;
		if (locBufGreen[0].length() > 100000) {
			*fqFile[Cstream1][0] << locBufGreen[0];
			*fqFile[Cstream1][1] << locBufGreen[1];
			locBufGreen[0].clear(); locBufGreen[1].clear();
		}
		sqfqostrMTX.unlock();
		return;
	}
	*/
	//case of yellow qual or unpaired reads..
	//write out read pair 1
	bool doMutex(Pair1 < 3 && pairedSeq > 1);
	if (doMutex) {
		if (BWriteFastQ) {
			/*fqPairFile[Cstream1]->write(d1->writeFastQ(), 0);
			fqPairFile[Cstream1]->write(d2->writeFastQ(), 1);*/
			fqPairFile[Cstream1]->write2(d1->writeFastQ(),d2->writeFastQ());
			ReadsWritten++; 
			return;
		}
	}

	if (doMutex) {
		sqfqostrMTX.lock();
	}
	if (BWriteFastQ) {//write in fastq format
		if (fq1cool) {
			*fqFile[Cstream1][Pair1 - 1] << d1->writeFastQ();
			ReadsWritten++;
		}
	} else {//Cstream1 in fasta (and maybe qual) format
		if (Cstream1 < (int) sFile.size()) {
			ReadsWritten++; 
			*sFile[Cstream1][Pair1 - 1] << d1->writeSeq(b_oneLinerFasta);
			if (BWriteQual) {
				*qFile[Cstream1][Pair1 - 1] << d1->writeQual(b_oneLinerFasta);
			}
		}
	}

	//write out read pair 2
	if (BWriteFastQ) {//write in fastq format
		if (fq2cool) {
			*fqFile[Cstream2][Pair2 - 1] << d2->writeFastQ();
			ReadsWritten++;
		}
	} else {//Cstream1 in fasta (and maybe qual) format
		if (Cstream2 < (int) sFile.size()) {
			*sFile[Cstream2][Pair2 - 1] << d2->writeSeq(b_oneLinerFasta);
			if (BWriteQual) {
				*qFile[Cstream2][Pair2 - 1] << d2->writeQual(b_oneLinerFasta);
			}
			ReadsWritten++;
		}
	}
	if (doMutex) {
		sqfqostrMTX.unlock();
	}
}

bool OutputStreamer::saveForWrite_merge(shared_ptr<DNAunique> d, 
		string newHeader,int curThread, bool elseWriteD1) {
	// CHECK WHERE TO IMPLEMENT BEST
// MERGE DNA1 AND DNA2
	shared_ptr<DNA> d2 = d->getPair();
	shared_ptr<DNAunique> dna_merged = d->getMerge();
	bool didMerge(false);
	
	if (!dna_merged && d2 != nullptr) {
		findSeedForMerge(d, d2, 0);
		if (d->merge_seed_pos_ > 0) {
			shared_ptr<DNA> tdM = mergers[curThread]->merge(d, d2);
			if (tdM != nullptr) {
				dna_merged = make_shared<DNAunique>(tdM, -1);
			}
		}
		d->attachMerge(dna_merged);
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
		*(fqFile[FS][Pair]) << d->writeFastQ();
	} else {
		*sFile[FS][Pair]<< d->writeSeq(b_oneLinerFasta);
		if (BWriteQual) { 
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
			string d1s,d2s;
			d->writeFastQ(d1s);
			d2->writeFastQ(d2s);
			dmltMTX.lock();
			*(fqNoBCFile[0]) << d1s;
			*(fqNoBCFile[1]) << d2s;
			dmltMTX.unlock();
			//nobcostrmMTX.unlock();
		}
	}
}

void OutputStreamer::setSubfilters(int num) {
	if (num<1){return;}
	subFilter.resize(num,NULL);
	for (uint i=0;i<subFilter.size();i++){
		if (subFilter[i] != nullptr) { delete subFilter[i]; }
		subFilter[i] = DBG_NEW Filters(MFil, MFil->currentBCnumber(), true);
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
				if (wr) { sFile[i][j]->emptyStream();}
				delete sFile[i][j];
				sFile[i][j] = nullptr;
			}
		}
	}
	for (size_t i = 0; i < qFile.size(); i++) {
		for (size_t j = 0; j < qFile[i].size(); j++) {
			if (qFile[i][j] != nullptr) {
				if (wr) { qFile[i][j]->emptyStream(); }
				delete qFile[i][j];
				qFile[i][j] = nullptr;
			}
		}
	}
	for (size_t i = 0; i < fqFile.size(); i++) {
		for (size_t j = 0; j < fqFile[i].size(); j++) {
			if (fqFile[i][j] != nullptr) {
				if (wr) { fqFile[i][j]->emptyStream(); }
				delete fqFile[i][j];
				fqFile[i][j] = nullptr;
			}
		}
	}
	for (size_t i = 0; i < of_merged_fq.size(); i++) {
		if (of_merged_fq[i] != nullptr) {
			if (wr) { of_merged_fq[i]->emptyStream(); }
			delete of_merged_fq[i];
			of_merged_fq[i] = nullptr;
		}
	}
	for (size_t i = 0; i < fqPairFile.size(); i++) {
		if (fqPairFile[i] != nullptr) {
			if (wr) { fqPairFile[i]->emptyStreams(true); }
			delete fqPairFile[i]; fqPairFile[i] = nullptr;
		}
	}
	for (size_t i = 0; i < sPairFile.size(); i++) {
		if (sPairFile[i] != nullptr) {
			if (wr) { sPairFile[i]->emptyStreams(true); }
			delete sPairFile[i]; sPairFile[i] = nullptr;
		}
	}
	for (size_t i = 0; i < qPairFile.size(); i++) {
		if (qPairFile[i] != nullptr) {
			if (wr) { qPairFile[i]->emptyStreams(true); }
			delete qPairFile[i]; qPairFile[i] = nullptr;
		}
	}



	if (MFil != nullptr) {
		MFil->setWrittenReads(ReadsWritten);
	}
#ifdef DEBUG
	cerr << ".. closed\n";
#endif

}

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
	//bool doMC = Nthrds > 1;

	fqNoBCFile.resize(2, NULL);
	if (fqNoBCFile[0] != nullptr) { delete fqNoBCFile[0]; }
	fqNoBCFile[0] = DBG_NEW ostr(tfnaout[0], wrMode, doTIO);
	if (fqNoBCFile[1] != nullptr) { delete fqNoBCFile[1]; }
	fqNoBCFile[1] = DBG_NEW ostr(tfnaout[1], wrMode, doTIO);
}

void OutputStreamer::openOFstreamFQpair(const string opOF, const string opOF2, std::ios_base::openmode wrMode,
	int p1,  string errMsg, bool onlyPrep) {
	//bool doMC = Nthrds >= 1;
	assert(p1 < 2);
	while (fqPairFile.size() <= p1) {
		//fqPairFile.resize(p1);
		fqPairFile.push_back(new dualOfBufStream());
	}
	if (!fqPairFile[p1]->open(opOF, wrMode, 0, doTIO,250000)) {
		cerr << "Could not open " << errMsg << " pair 1fastq output file " << opOF << endl << p1 <<  " " << totalFileStrms << endl;
		exit(5);
	}	
	if (!fqPairFile[p1]->open(opOF2, wrMode, 1, doTIO,350000)) {
		cerr << "Could not open " << errMsg << " pair 2 fastq output file " << opOF << endl << p1 << " " << totalFileStrms << endl;
		exit(5);
	}	
	fqPairFile[p1]->activate();
}

void OutputStreamer::openOFstreamFQ(const string opOF, std::ios_base::openmode wrMode, 
	int p1, int p2, string errMsg, bool onlyPrep) {
	//p1: 0:green, 1:yellow
	//p2: 0:pair1, 1:pair2, 2:single1, 3:single2
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1+1 >= (int)fqFile.size()) {
		vector<ostr*> nullVec(4, NULL);
		fqFile.resize(p1+1, nullVec);
	}
	//bool doMC = Nthrds > 1;

	if (fqFile[p1][p2] != nullptr) { delete fqFile[p1][p2]; }
	fqFile[p1][p2] = DBG_NEW ostr(opOF, wrMode, doTIO);
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
	//bool doMC = Nthrds > 1;
	if (of_merged_fq[p1] != nullptr) { delete of_merged_fq[p1]; }
	of_merged_fq[p1] = DBG_NEW ostr(opOF, wrMode, doTIO);
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
	//bool doMC = Nthrds > 1;

	if (sFile[p1][p2] != nullptr) { delete sFile[p1][p2]; }
	sFile[p1][p2] = DBG_NEW ostr(opOF, wrMode, doTIO);
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
	//bool doMC = Nthrds > 1;

	//if (onlyPrep) { return; }
	//if (p1 == 1 && !b_writeYellowQual){ return; }//p1==1: mid passed suppressOutWrite >= 2
	//if (p1 == 0 && !b_writeGreenQual){ return; }//suppressOutWrite == 1
	if (qFile[p1][p2] != nullptr) { delete qFile[p1][p2]; }
	qFile[p1][p2] = DBG_NEW ostr(opOF, wrMode, doTIO);
	if (!*qFile[p1][p2]) {
		cerr << "Could not open " << errMsg << " quality output file " << opOF << endl;
		exit(4);
	}

}

void OutputStreamer::openSeveralOutstreams(OptContainer* cmdArgs, shared_ptr<ReadSubset> RDS, std::ios_base::openmode wrMode) {
#ifdef DEBUG
	cerr << " openining multiple out streams" << endl;
#endif
	string path = "", fileEnd(".fna");
	vector<string> outFile = RDS->getOFiles();
	bool openStrms = true; int omode(1);
	if (cmdArgs->find("-o_fastq") != cmdArgs->end() && (*cmdArgs)["-o_fastq"] != "" && (*cmdArgs)["-o_fna"] == "") { //write fastq
		BWriteFastQ = true;
		path = (*cmdArgs)["-o_fastq"];
		omode = 0;
		fileEnd = ".fastq";
	} else  if ((*cmdArgs)["-o_fna"] != "") {
		BWriteFastQ = false;
		path = (*cmdArgs)["-o_fna"];
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
	if ((*cmdArgs)["-excludeFile"] != "") {
		//set a flag
		RDS->setRemainingFilepipe(i);
		if (pairedSeq == 1) {
			openOFstream((*cmdArgs)["-excludeFile"], wrMode, i, 0, "excludeFile file ", !openStrms, omode);
		} else {
			string baseFile = path + removeFileEnding((*cmdArgs)["-excludeFile"]);
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

void OutputStreamer::openOutStreams(OptContainer* cmdArgs,int fileIt,
	std::ios_base::openmode wrMode_i, 
	string fileExt, int forceFmt){
	this->setwriteMode(wrMode_i);
	if ( suppressOutWrite == 3 || ((*cmdArgs)["-o_fastq"] == "" && (*cmdArgs)["-o_fna"] == "" && (*cmdArgs)["-o_qual"] == "") ){
		suppressOutWrite = 3; b_writeGreenQual = false; 
		b_writeYellowQual = false; return;
	}
	if (forceFmt != -1){
		if (forceFmt == 1 && (cmdArgs->find("-o_fna") == cmdArgs->end() || (*cmdArgs)["-o_fna"] == "") ){//force fna, no qual, required for seed ref fastas
			(*cmdArgs)["-o_fna"] = (*cmdArgs)["-o_fastq"];
			(*cmdArgs)["-o_fastq"] = "";
		}
	}
	

	if (cmdArgs->find("-o_fastq")  != cmdArgs->end() && (*cmdArgs)["-o_fastq"] != "" && (*cmdArgs)["-o_fna"] == ""){ //write fastq
		this->setFastQWrite(true);
		if (pairedSeq>1){ //open second pair_ + singleton
#ifdef DEBUG
			cerr << " paired fastq out " << endl;
#endif
			vector<string> tfnaout = splitByCommas((*cmdArgs)["-o_fastq"]);
			if (tfnaout.size()!=2){
				cerr<<"Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = "<<(*cmdArgs)["-o_fastq"]<<endl; exit(57);
			}
			leadingOutf = applyFileIT(tfnaout[0] + fileExt, fileIt);
			//now save pairs in dedicated object:
			openOFstreamFQpair(leadingOutf.c_str(), applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, "paired out");

			/*
			//1st pair (or main) file out
			openOFstreamFQ(leadingOutf.c_str(), wrMode, 0, 0, "paired 1st");
			//2nd pair
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			//1st singleton
			*/
			openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 2, "Singleton 1", true);
			//2nd singleton
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 3, "Singleton 2", true);
			
			//merged reads
			int prefL = (int)len_common_prefix_base(tfnaout[0].c_str(), tfnaout[1].c_str());
			string prefix = tfnaout[0].substr(0,prefL);
			openOFstreamFQ_mrg(applyFileIT(prefix + "merg.fq" + fileExt, fileIt).c_str(), wrMode, 0, "merged");
			
			//additional file
			if (cmdArgs->find("-o_fastq2")  != cmdArgs->end() && (*cmdArgs)["-o_fastq2"].length()>1){ //write fastq
				tfnaout = splitByCommas((*cmdArgs)["-o_fastq2"]);
				openOFstreamFQpair(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, "additional paired out");
				/*openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st");
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd");
				*/
				openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt , fileIt, SingletonFileDescr).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt , fileIt, SingletonFileDescr).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
			} else { b_writeYellowQual = false; }
		} else {
#ifdef DEBUG
			cerr << " single fastq out " << endl;
#endif
			openOFstreamFQ(applyFileIT((*cmdArgs)["-o_fastq"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "the main");
			leadingOutf = applyFileIT((*cmdArgs)["-o_fastq"] + fileExt, fileIt);
			//additional file
			if (cmdArgs->find("-o_fastq2") != cmdArgs->end() && (*cmdArgs)["-o_fastq2"].length()>1) { //write fastq
				openOFstreamFQ(applyFileIT((*cmdArgs)["-o_fastq2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "the additional");
			}
			else { b_writeYellowQual = false; }
		}

		return;
	}
	BWriteFastQ=false;
	if (pairedSeq==1){
#ifdef DEBUG
		cerr << " fasta singleton output " << endl;
#endif
		openOFstreamFNA(applyFileIT((*cmdArgs)["-o_fna"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
		leadingOutf = applyFileIT((*cmdArgs)["-o_fna"] + fileExt, fileIt);
		if ((*cmdArgs)["-o_qual"] != ""){
			openOFstreamQL(applyFileIT((*cmdArgs)["-o_qual"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
			this->setQualWrite(true);
		} else {
			this->setQualWrite(false);
		}
		//additional file (secondary filter)
		if (cmdArgs->find("-o_fna2") != cmdArgs->end() && (*cmdArgs)["-o_fna2"].length()>1) { //add file
			openOFstreamFNA(applyFileIT((*cmdArgs)["-o_fna2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
			if (cmdArgs->find("-o_qual2") != cmdArgs->end() && (*cmdArgs)["-o_qual2"].length() > 1) { //add file
				openOFstreamQL(applyFileIT((*cmdArgs)["-o_qual2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
				this->setQualWrite(true);
			}
		} else { b_writeYellowQual = false; }

	} else {
#ifdef DEBUG
		cerr << " fasta paired output " << endl;
#endif

		vector<string> tfnaout = splitByCommas((*cmdArgs)["-o_fna"]);
		if (tfnaout.size()!=2){
			cerr<<"Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = "<<(*cmdArgs)["-o_fna"]<<endl; exit(57);
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
		if (cmdArgs->find("-o_fna2")  != cmdArgs->end() && (*cmdArgs)["-o_fna2"].length()>1){ //write fastq
			tfnaout = splitByCommas((*cmdArgs)["-o_fna2"]);

			openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st", true);
			openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd", true);
			openOFstreamFNA(applyFileIT(tfnaout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
			openOFstreamFNA(applyFileIT(tfnaout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
		}
		if ((*cmdArgs)["-o_qual"] != ""){
			this->setQualWrite(true);
			vector<string> tqout = splitByComma((*cmdArgs)["-o_qual"],true);
			openOFstreamQL(applyFileIT(tqout[0] + fileExt, fileIt).c_str(), wrMode, 0, 0, "paired 1st");
			openOFstreamQL(applyFileIT(tqout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			openOFstreamQL(applyFileIT(tqout[0] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 2, "Singleton 1", true);
			openOFstreamQL(applyFileIT(tqout[1] + fileExt + SingletonFileDescr, fileIt).c_str(), wrMode, 0, 3, "Singleton 2", true);
			if ((*cmdArgs)["-o_qual2"] != "") {
				vector<string> tqout = splitByComma((*cmdArgs)["-o_qual2"], true);
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
	if (dna1->length() == 0 || dna2->length() == 0) {
		total_read_preMerge_++;  return; 
	}
	bool didMerge = mergers[thrPos]->findSeedForMerge(dna1, dna2);
    /*if (merger[thrPos]->findSeed(dna1->getSequence(), dna2->getSequence())) {
		didMerge = true;
        dna1->merge_seed_pos_ = (int) merger[thrPos]->result.seed.pos1;
        dna1->merge_offset_ = merger[thrPos]->result.offset1;
        dna2->merge_seed_pos_ = (int) merger[thrPos]->result.seed.pos2;
        dna2->merge_offset_ = merger[thrPos]->result.offset2;
        dna2->reversed_merge_ = merger[thrPos]->result.seed.is2reversed;
    }
	*/
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


void DNAuniqSet::setBest(bool addCnts) {
	int bestCnt = 0;
	int bestPos = -1;
	if (DNUs.size() == 1) {
		bestSet = true;
		bestDNU = DNUs.begin()->second;
		return;
	}
	//shared_ptr<DNAunique> lastBest;
	for (auto dd : DNUs) {
		if (dd.second == nullptr) {
			continue;
		}
		bool nhM = dd.second->getMerge() != nullptr;
		//weigh by whether any has a merge
		float modN = 1.f; float modB = 1.f; float ratMLs(1.f);
		if (nhM && bestHasMerge) {//compare merge length
			ratMLs = (float)dd.second->getMerge()->length() / (float)bestDNU->getMerge()->length();
		}
		else {
			if (!nhM) { modN = 0.8f; }
			if (!bestHasMerge) { modB = 0.8f; }
		}
		float ratCns = ((float)dd.second->totalSum() * modN )/  ( (float)bestCnt * modB );
		ratCns *= ratMLs;
		if (ratCns > 1  ) {
			bestCnt = dd.second->totalSum();
			bestPos = dd.first;
			bestDNU = dd.second;
		}
	}
	//completely unbiased selection of whatever has the highest counts.. could select non-merge before merge
	if (bestPos != -1 && DNUs.size() > 1 && !cntsAdded2best) {

		if (addCnts) {
			for (auto dd : DNUs) {
				if (bestPos != dd.first) {
					bestDNU->transferOccurence(dd.second);
				}
			}
		} else {
			//include + - 1?
			auto xx = DNUs.find((bestPos - 1));
			if (xx != DNUs.end()) { bestDNU->transferOccurence(xx->second); }
			xx = DNUs.find((bestPos + 1));
			if (xx != DNUs.end()) { bestDNU->transferOccurence(xx->second); }
		}
		cntsAdded2best = true;
	}
	bestSet = true;
}

//*******************************************
//*        DEREPLICATE OBJECT
//*******************************************
Dereplicate::Dereplicate(OptContainer* cmdArgs, Filters* mf):
        barcode_number_to_sample_id_(0), b_usearch_fmt(true), b_singleLine(true), b_pairedInput(false),
        minCopies(1,0), minCopiesStr("0"), //default minCopies accepts every derep
		tmpCnt(0), curBCoffset(0), b_derep_as_fasta_(true), b_derepPerSR(false),
		b_wroteMapHD(false), b_merge_pairs_derep_(false),merger(nullptr),
		mapF(""), outHQf(""), outHQf_p2(""), outRest(""),
		mainFilter(mf)
{
	outfile = (*cmdArgs)["-o_dereplicate"];

	string baseOF = outfile.substr(0, outfile.find_last_of('.'));
	mapF = baseOF + ".map";
	outHQf = baseOF + ".1.hq.fq";
	outHQf_p2 = baseOF + ".2.hq.fq";
	outRest = outfile + ".rest";
	//clean up eventually existing files
	remove(mapF.c_str());	remove(outHQf.c_str());
	remove(outHQf_p2.c_str());	remove(outRest.c_str());


	if (cmdArgs->find("-dere_size_fmt") != cmdArgs->end() && (*cmdArgs)["-dere_size_fmt"] == "1") {
		b_usearch_fmt = false;
	}
	if ((*cmdArgs)["-derepPerSR"] == "1") {
		b_derepPerSR = true;
	}

	if (cmdArgs->find("-min_derep_copies") != cmdArgs->end()) {
		minCopies[0] = -1;//in this case reset to -1 the first entry..
		minCopiesStr = (*cmdArgs)["-min_derep_copies"];
		string x = (*cmdArgs)["-min_derep_copies"];

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
		//minCopies = atoi((*cmdArgs)["-min_derep_copies"].c_str());

	}
	// set up output format of dereplication:
	if (cmdArgs->find("-derep_format") != cmdArgs->end()) {
		auto res = (*cmdArgs)["-derep_format"];
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
    if (cmdArgs->find("-merge_pairs_derep") != cmdArgs->end() && 
		(*cmdArgs)["-merge_pairs_derep"] == "1") {
		b_merge_pairs_derep_ = true;
    }
}


bool Dereplicate::addDNA(shared_ptr<DNA> dna, shared_ptr<DNA> dna2) {
	//1st build hash of DNA
	if (!dna->getBarcodeDetected()) {
		return false;
	}
	// Get copy of sequence (might have already been modified)
	//string seq = dna->getSeqPseudo();

	int sample_id = dna->getBCnumber();
	bool pass = dna->isGreenQual();
	//deactivate this for now..

	//completely deactive the merge position searches..reactivate later
	//int MrgPos1 = -1;// dna->merge_seed_pos_;
	int MrgPos1 = dna->merge_seed_pos_;
	//int MrgPos1 = -1;
	shared_ptr<DNA> dna_merged = nullptr;
	string srchSeq("");

	bool searchWithMerg = true;

	if (searchWithMerg && dna2 != nullptr && dna->merge_seed_pos_ != -1) {
		//MrgPos1 = dna2->merge_seed_pos_;
		dna_merged = merger->merge(dna, dna2);
	}

	if (dna_merged){
		srchSeq = dna_merged->getSeqPseudo().substr(0, dna->length());
		dna_merged = nullptr;
		//merge can be a lot shorter, potentially leading to problems searching this seq
/*		if (false && srchSeq.length() != dna->length()) {
			//int y = 1;
			srchSeq = dna->getSeqPseudo();
		}*/
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
			dna_unique1->second.addNewDNAuniq(dna, dna2, dna_merged, MrgPos1, sample_id);
		} else { // compare to existing DNA
			dna_unique->second->matchedDNA(dna, dna2, dna_merged, sample_id, b_derep_as_fasta_);
		} 
		dna_unique1->second.lockMTX.unlock();//just to be on safe side, lock entire section
	}
	drpMTX.unlock_shared();//lock for hash

	if (new_insert && pass) {
		drpMTX.lock(); 
		Tracker[srchSeq].addNewDNAuniq(dna, dna2, dna_merged, MrgPos1, sample_id);
		drpMTX.unlock();
        // Create new dna_unique object
		//cdbg("set new DNAderep ");
        /*
		dna->setYellowQual(false); dna->setDereplicated();
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
	//passedSize = 0; 
	tmpCnt = 0;
	//barcode_number_to_sample_id_.resize(0);
	Tracker.clear();
}
bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) { 
	return l->totalSum() > r->totalSum(); 
}

void Dereplicate::finishMap() {
	//at this point we can onl be sure that barcode_number_to_sample_id_ is finished
	//hence now it the point to add this to map
	//passedSize = 0; 
	tmpCnt = 0;
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
	int x = std::rename((mapF + "t").c_str(), mapF.c_str());

	return;
}
string Dereplicate::writeDereplDNA(Filters* mf, string SRblock) {
	ofstream of, omaps, of2, ofRest, of2p2, of_merged;
	cerr << "\nEvaluating and writing dereplicated reads..\n";
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
        dereplicated_dnas[count] = dd->second.best(true);
		count++;
	}
//	bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r) {	return l->getCount() < r->getCount();}
	sort(dereplicated_dnas.begin(), dereplicated_dnas.end() , DNAuPointerCompare);
	size_t passedSize = 0; size_t notPassedSize = 0; size_t passed_hits(0);
	//bool thrHit = false;

	//sanity check
	vector<int> counts_per_sample(barcode_number_to_sample_id_.size(), 0);
	int total_count(0);
	int passed_count(0);
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
			passedSize += dna->totalSum();

			//only do merge here, because these are known good dereps already
			if (b_merge_pairs_derep_ && dna->merge_seed_pos_ >= 0) {
				dna_merged = merger->merge(dna, dna->getPair());
			}
			derepNowOut = &of;
        } else {
			notPassedSize += dna->totalSum();
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
//	if (tmpCnt != passedSize) {
//		cerr << "Counting failed\n" << tmpCnt << " " << passedSize<<endl;
//	}
	of.close(); omaps.close();  ofRest.close();
	float avgSize = (float)passedSize / (float)(passed_hits);
	string report = "";
	string N_notPassed = intwithcommas(int(dereplicated_dnas.size() - passed_hits));
	string N_total = intwithcommas(int(Tracker.size() ));
	string N_passed = intwithcommas((int)passed_hits);
	report += "Dereplication: " + N_passed +
		" unique sequences (avg size " + ftos(avgSize, 2) + "; "+ intwithcommas(passedSize) + " counts)\n";
	
		if (passed_hits > 0) {
		report += N_notPassed + "/" + N_total + " not passing derep conditions (" + intwithcommas(notPassedSize) + " counts; "+ minCopiesStr;
		if (dereplicate_sample_specific) { report += " & sample specific restrictions"; }
	report += ")";
	}

		
		
	cerr  << report << endl << endl;
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


Filters::Filters(OptContainer* cmdArgs1) :
        PrimerL(0), PrimerR(0), PrimerL_RC(0), PrimerR_RC(0), PrimerIdx(0),
        Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
        HeadSmplID(0),
        hetPrimer(2,vector<string>(0)),
        collectStatistics(2), statAddition(2),
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
		BextensivePrimerChecks(false),
        bRequireRevPrim(false), alt_bRequireRevPrim(false),
        bRequireFwdPrim(false), alt_bRequireFwdPrim(false), BcutTag(true),
        bCompletePairs(false), bShortAmplicons(false),
        minBCLength1_(0), minBCLength2_(0), maxBCLength1_(0), maxBCLength2_(0), minPrimerLength_(0), maxHomonucleotide(0), trimHomonucleotide(0),
		cut5PR1(0), cut5PR2(0),
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
        //revConstellationN(0),
        BCdFWDREV(2),
		firstXreadsW(-1), firstXreadsR(-1),
        restartSet(false), b_optiClusterSeq(false),
        b_subselectionReads(false), b_doQualFilter(true),
        b_doFilter(true),
        bDoDereplicate(false), bDoSeedExtension(false), bDoCombiSamples(false),
        maxReadsPerOFile(0), ReadsWritten(0), OFileIncre(0),
		demultiBPperSR(0),
        barcodeLengths1_(0), barcodeLengths2_(0),
		illuPEfwd(""), illuPErev(""), illuSEuni(""), illuSEidx(""), 
		Bcheck4illuAdapts(false), doGoldAxe(false),
		GoldAxeMinAmpli(-1), GoldAxeMaxAmpli(-1),
        cmdArgs(cmdArgs1), passed_interval_reads(0)
		{
	//csMTX[0].unlock(); 
	//csMTX[1].unlock();

	//set up objects to collect statistics on run
    collectStatistics[0] = make_shared<collectstats>(); 
	collectStatistics[1] = make_shared<collectstats>();

	statAddition[0] = make_shared<collectstats>(); 
	statAddition[1] = make_shared<collectstats>();
	
	GAstatistics = make_shared<GAstats>();

	mergeStats = make_shared<MEstats>();
	


	bool alt_bRequireRevPrimSet=false;

	string optF ("");
	if (cmdArgs->find("-options") != cmdArgs->end()) {
		optF = (*cmdArgs)["-options"];
	}

	iniSpacer = (*cmdArgs)["-sample_sep"];
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
	int maxHomoNT(12); int trimHomoNT(12);
	bool keepTag(false),keepPrimer(false);
	bool addModConf = false;

	//set up some basic objects
	if ( (*cmdArgs).find("-paired") != cmdArgs->end() ) {
		pairedSeq = atoi((*cmdArgs)["-paired"].c_str()); //fakeEssentials();
		if ( pairedSeq<1 || pairedSeq>3 ) { cerr << "Argument \"-paired\" supplied with unknown parameter. Aborting.\n"; exit(28); }
		if ( (*cmdArgs)["-onlyPair"] == "1" || (*cmdArgs)["-onlyPair"] == "2" ) {
			pairedSeq = 1;
		}
	}
	if (cmdArgs->find("-normRdsToFiveNTs") != cmdArgs->end()) {
		norm2fiveNTs = true;
		cerr << "Warning: normRdsToFiveNTs is not implemented!\n";
	}
	if ((*cmdArgs)["-logLvsQ"].c_str() != "") {
		collectStatistics[0]->setbLvsQlogsPreFilt(true);
	}
	if ((*cmdArgs)["-GoldenAxe"] == "1") { this->setGoldAxe(true, stoi((*cmdArgs)["-GoldenAxeMaxAmpli"]), stoi((*cmdArgs)["-GoldenAxeMinAmpli"])); }


	//delimit output file size to X reads
	if (cmdArgs->find("-maxReadsPerOutput") != cmdArgs->end()) {
		maxReadsPerOFile = atoi((*cmdArgs)["-maxReadsPerOutput"].c_str());
	}
	if (cmdArgs->find("-DemultiBPperSR") != cmdArgs->end()) {
		stringstream ss((*cmdArgs)["-DemultiBPperSR"]);
		double d = 0;
		ss >> d;
		demultiBPperSR =(uint) d;
	}
	//important for fastq format
	if ( cmdArgs->find("-i_qual_offset") != cmdArgs->end() ) {
		if ( (*cmdArgs)["-i_qual_offset"] == "auto" ) {
			userReqFastqVer = 0;
		} else {
			userReqFastqVer = atoi((*cmdArgs)["-i_qual_offset"].c_str());
		}
	}

	cut5PR1= atoi((*cmdArgs)["-5PR1cut"].c_str());
	cut5PR2 = atoi((*cmdArgs)["-5PR2cut"].c_str());
	//cerr<<(*cmdArgs)["-o_qual_offset"]<<endl;
	userReqFastqOutVer = atoi((*cmdArgs)["-o_qual_offset"].c_str());
	//statistic tracker
		//do new SEED sequence selection?
	if ( cmdArgs->find("-optimalRead2Cluster") != cmdArgs->end() ) {
		b_optiClusterSeq = true;
	}
	//do selection of specific reads?
	if ( (*cmdArgs)["-specificReads"] != "" ) {
		b_subselectionReads = true;
	}
	if (cmdArgs->find("-binomialFilterBothPairs") != cmdArgs->end() && (*cmdArgs)["-binomialFilterBothPairs"] == "1") {
		b_BinFilBothPairs = true;
	}

	if ((*cmdArgs)["-illuminaClip"] == "1") {
		Bcheck4illuAdapts = true;
	}
	if ((*cmdArgs)["-XfirstReadsWritten"] != "") {
		firstXreadsW = atoi((*cmdArgs)["-XfirstReadsWritten"].c_str());
	}
	if ((*cmdArgs)["-XfirstReadsRead"] != "") {
		firstXreadsR = atoi((*cmdArgs)["-XfirstReadsRead"].c_str());
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
		}else if (strcmp(segs.c_str(),"trimHomonucleotide") == 0){
			trimHomoNT = atoi(segs2.c_str());
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
		} else if (segs == "ExtensivePrimerChecks") {
			if (segs2 == "T") {
				BextensivePrimerChecks = true;
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
		} else if (segs == "SyncReadPairs") {
			if (segs2 == "T") {
				bChkRdPrs = true;
			}
			else { bChkRdPrs = false; }
		} else if (segs == "illuminaFwd") {
			illuPEfwd = segs2;
		} else if (segs == "illuminaRev") {
			illuPErev = segs2;
		} else if (segs == "illuminaSngUni") {
			illuSEuni = segs2;
		}else if (segs == "illuminaSngIdx") {
			illuSEidx = segs2;
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
	this->setTrimHomo(trimHomoNT);

	//alternative options (mid qual filtering)
	if (addModConf){
		if (!alt_bRequireRevPrimSet){alt_bRequireRevPrim = bRequireRevPrim;}
		if (cmdArgs->find("-o_fna")  != cmdArgs->end() && (*cmdArgs)["-o_fna"].length()>1){
			if (cmdArgs->find("-o_fna2")  == cmdArgs->end()){
				(*cmdArgs)["-o_fna2"] = additionalFileName((*cmdArgs)["-o_fna"]);
				//(*cmdArgs)["-o_fna2"] = (*cmdArgs)["-o_fna"].substr(0,(*cmdArgs)["-o_fna"].length()-4)+".add.fna";
			}
		} else if (cmdArgs->find("-o_fastq")  != cmdArgs->end() && (*cmdArgs)["-o_fastq"].length()>1){
			if (cmdArgs->find("-o_fastq2")  == cmdArgs->end()){
				(*cmdArgs)["-o_fastq2"] = additionalFileName((*cmdArgs)["-o_fastq"]);
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
		collectStatistics(2,nullptr),  statAddition(2,nullptr),
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
		BextensivePrimerChecks(of->BextensivePrimerChecks),
		bRequireFwdPrim(of->bRequireFwdPrim), alt_bRequireFwdPrim(of->alt_bRequireFwdPrim),
        BcutTag(of->BcutTag),

        bCompletePairs(of->bCompletePairs), bShortAmplicons(of->bShortAmplicons),
        minBCLength1_(of->minBCLength1_), minBCLength2_(of->minBCLength2_), maxBCLength1_(of->maxBCLength1_), maxBCLength2_(of->maxBCLength2_), minPrimerLength_(of->minPrimerLength_), maxHomonucleotide(of->maxHomonucleotide), trimHomonucleotide(of->trimHomonucleotide),
		cut5PR1(of->cut5PR1), cut5PR2(of->cut5PR2),
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
        //revConstellationN(0),
        BCdFWDREV(of->BCdFWDREV),
		firstXreadsW(of->firstXreadsW), firstXreadsR(of->firstXreadsR),
		restartSet(false),
        b_optiClusterSeq(of->b_optiClusterSeq), b_subselectionReads(of->b_subselectionReads),
        b_doQualFilter(of->b_doQualFilter),
        b_doFilter(of->b_doFilter),
        bDoDereplicate(of->bDoDereplicate),
		bDoSeedExtension(of->bDoSeedExtension),
        bDoCombiSamples(of->bDoCombiSamples),
        maxReadsPerOFile(of->maxReadsPerOFile),
		demultiBPperSR(of->demultiBPperSR),
        //ReadsWritten(of->ReadsWritten), OFileIncre(of->OFileIncre),
        barcodeLengths1_(0), barcodeLengths2_(0), 
	
		illuPEfwd(of->illuPEfwd), illuPErev(of->illuPErev), illuSEuni(of->illuSEuni), illuSEidx(of->illuSEidx),
		Bcheck4illuAdapts(of->Bcheck4illuAdapts),
		doGoldAxe(of->doGoldAxe),
		GoldAxeMinAmpli(of->GoldAxeMinAmpli), GoldAxeMaxAmpli(of->GoldAxeMaxAmpli),
	
		SequencingRun(0),cmdArgs(of->cmdArgs), passed_interval_reads(0)
{
	cdbg("New Filter object from copy\n");
	ReadsWritten = of->writtenReads();
	OFileIncre = of->getFileIncrementor();
    BCdFWDREV[0].reset(); BCdFWDREV[1].reset();
	//collectStatistics.resize(2); statAddition.resize(2);
	//csMTX[0] = new mutex(); csMTX[1] = new mutex();

	collectStatistics[0] = make_shared<collectstats>(); 
	collectStatistics[1] = make_shared<collectstats>();
	collectStatistics[0]->setbLvsQlogsPreFilt (of->collectStatistics[0]->getbLvsQlogsPreFilt());
    statAddition[0] = make_shared<collectstats>(); statAddition[1] = make_shared<collectstats>(); 
	GAstatistics = make_shared<GAstats>();
	mergeStats = make_shared<MEstats>();


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


Filters::~Filters() {
	cdbg("Deleting filter .. ");
//	for (size_t i = 0; i < csMTX.size(); i++){
//		delete csMTX[i];
//	}

//	for (size_t i = 0; i < 2; i++) { delete PostFilt[i]; delete RepStatAddition[i]; }
//	delete PreFiltP1; delete PreFiltP2;
	cdbg("Done\n");
}


Filters* Filters::newFilterPerBCgroup(const vector<int> idxi) {
	
	if (idxi.size() < 1) {
		return nullptr;}
	cdbg("newFilterPerBCgroup::start : "+ itos(idxi[0])+"\n");

	// get filter from main filter object passing an index for mapping?!
//	shared_ptr<Filters> filter = make_shared<Filters>(shared_from_this(), idxi[0]);
	Filters* filter = DBG_NEW Filters(this, idxi[0]);
	cdbg("newFilterPerBCgroup::star2t\n");

	// number of mapping file lines associated with that unique fastx
	unsigned int tarSize = (unsigned int)idxi.size();
	filter->allResize(tarSize);

	int tarID = -1;
	bool isDoubleBarcoded = this->doubleBarcodes();
	cdbg("newFilterPerBCgroup::Go over BCs\n");

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
	cdbg("newFilterPerBCgroup::check 4 doubles\n");
	//sanity check no double barcodes..
	filter->checkDoubleBarcode();
	filter->singReadBC2();


	return filter;
}

void Filters::miniCheckDNA(shared_ptr<DNA> d, shared_ptr<DNA> d2) {

	int BCoffs = getBCoffset();
	bool checkSwitchedRdPairs = this->checkSwitchedRdPairs();
	bool dualBCs = this->doubleBarcodes();
	bool doBCsAtAll = this->doBarcodes();
	int pairedRd = this->isPaired();

	//needs some basic cleanups here..
	bool wasReversed;
	if (pairedRd == 2) {
		vector< shared_ptr<DNA>>tdn(0); tdn.push_back(d);
		if (d2 != nullptr) { tdn.push_back(d2); }
		else { tdn.push_back(nullptr); }
		wasReversed = this->swapReverseDNApairs(tdn);
	}
	else if (d2 == nullptr) {
		wasReversed = this->isReversedAmplicon(d);
	}

	//remove base pairs 3' or 5' ?
	if (getcut5PR1()) { d->cutSeq(0, getcut5PR1()); }
	if (getcut5PR2()) { d2->cutSeq(0, getcut5PR2()); }
	if (removeAdapter) { remove_adapter(d); }


	if (BcutPrimer || Bcheck4illuAdapts) {
		int tagIdx = 0; //for now just set to 0.. if an experiment uses > 1 primers this would need to change
		bool fwdRC = false; bool revRC = true;//if this gets changed, the read needs to be reverse complemented
		cutPrimer(d, PrimerIdx[tagIdx], fwdRC, 0, BextensivePrimerChecks);
		if (bShortAmplicons) {//also check other end of primer.. also use for PacBio amplicons
			//case for 1) long read 2) look for both primers 3) rev required
			cutPrimerRev(d, PrimerIdxRev[tagIdx], revRC, BextensivePrimerChecks);
		}
		if (d2 != nullptr) {//pair_ == 1, check for fwd primer in pair_ 2 (rev-compl)
			bool revCheck = false;// pair == -1 || pair == 0;//1:false for RC, else always a reverse check
			cutPrimerRev(d2, PrimerIdxRev[tagIdx], revCheck, false);
			if (bShortAmplicons) {//also check other end of primer..
				cutPrimer(d, PrimerIdx[tagIdx], revCheck, 1);
			}
		}
		//conditions for failing read on not finding fwd primer
		if (!d->getFwdPrimCut() && bRequireFwdPrim) {
			d->setYellowQual(true);
		}
		//conditions for failing read on not finding rev primer
		if (d2 != nullptr && bRequireRevPrim && !d->getRevPrimCut()) {//failed to find reverse primer
			d->setYellowQual(true);
		}
	}

}

//service function to ini OTU Seed extension 
UClinks * Filters::ini_SeedsReadsDerep(UClinks *ucl, shared_ptr<ReadSubset>& RDSset, 
	shared_ptr<Dereplicate>& Dere) {
	if (this->doOptimalClusterSeq()) {
		ucl = DBG_NEW UClinks(cmdArgs);
		if (cmdArgs->find("-mergedPairs") != cmdArgs->end() && (*cmdArgs)["-mergedPairs"] == "1") {
			ucl->pairedSeqsMerged();
			this->setFloatingEWin(0, 0.f);
		}
		else {
			this->setFloatingEWin(10, 25);
		}
		//are fallback fasta sequences available?
		if ((*cmdArgs)["-OTU_fallback"] != "") {
			shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true,
				this->getuserReqFastqVer(), (*cmdArgs)["-ignore_IO_errors"], (*cmdArgs)["-pairedRD_HD_out"],1);
			FALL->setupFna((*cmdArgs)["-OTU_fallback"]);
			ucl->setupDefSeeds(FALL, SampleID);
		}
		ucl->activateMerger();
	}
	else if (this->doSubselReads()) {
		//this will select a list of reads and distribute these into multiple files
		RDSset = make_shared<ReadSubset>((*cmdArgs)["-specificReads"], "");
	}
	else if (this->doDereplicate()) {
		Dere = make_shared<Dereplicate>(cmdArgs, this);
	//	ReadMerger* merg = DBG_NEW ReadMerger(); //create special object for these functions
		Dere->activateMerger();
	}
	return ucl;
}



//simulates that in mapping file links to sequence file was given.
bool Filters::setcmdArgsFiles(){

	if (FastqF.size()==0 && QualF.size()==0 && FastaF.size() > 0){
		//fasta entry but no qual entries

		string path="";
		if (cmdArgs->find("-i_path")  != cmdArgs->end() && (*cmdArgs)["-i_path"].length() >= 1){
			path=(*cmdArgs)["-i_path"] + string("/");
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
			} else if (cmdArgs->find("-number")!=  cmdArgs->end() && (*cmdArgs)["-number"] =="T"){
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
		if (cmdArgs->find("-i_fastq")  == cmdArgs->end()){
			FastaF.resize(fileSiz);
			QualF.resize(fileSiz);
			for (unsigned int i=0; i< FastaF.size(); i++){
				FastaF[i] = (*cmdArgs)["-i_fna"];
				QualF[i] = (*cmdArgs)["-i_qual"];
			}
		} else {//fastq input
			vector<string> fqTmp (1,(*cmdArgs)["-i_fastq"]);
			if ((*cmdArgs)["-i_fastq"].find(";") != string::npos) {//";" denotes several files
				if (fileSiz == 1) {//no BC, 
					fqTmp = splitByCommas((*cmdArgs)["-i_fastq"], ';');
					this->allResize((uint) fqTmp.size());
					fileSiz = (int) fqTmp.size();
					cerr << "Detected " << fileSiz << " input files (pairs)." << endl;
					FastqF = fqTmp;
				} else {
					cerr << "Fastq string contains symbol \";\". Not allowed in input string"; exit(32);
				}
			} else {
				FastqF.resize(fileSiz, (*cmdArgs)["-i_fastq"]);
			}
		}
	}


	if (MIDfqF.size() == 0)
		if (cmdArgs->find("-i_MID_fastq") != cmdArgs->end()) {
		MIDfqF.resize(fileSiz, "");
		for (unsigned int i = 0; i < MIDfqF.size(); i++) {
			MIDfqF[i] = (*cmdArgs)["-i_MID_fastq"];
		}
	}
		

	if ((*cmdArgs)["-o_dereplicate"] != "") {
		//check if file could exist
		ofstream temp;
		temp.open((*cmdArgs)["-o_dereplicate"].c_str(), ios::out);
		if (!temp) { cerr << "Could not open outstream to dereplicated sequences:\n" << (*cmdArgs)["- o_dereplicate"] << endl; exit(78); }
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

bool Filters::isReversedAmplicon( shared_ptr<DNA> tdn) {
	if ((!checkRevRd() ) ) {
		return false;
	}


	//BC should be already cut at this point..
	//int tagIdx = 0;//just try

	//method 1: just check if primer is found reversed, most basic and seems to work fine..
	//simple check if fwd rev primer is in reverse position: then reverse transcribe

	bool fwdRd1Primer = checkIfPrimerHits(tdn, 0, 0);
	bool fwdRd1PrimerRev = checkIfRevPrimerHits(tdn, 0, 0);
	if (fwdRd1Primer ) {
		return false;
	}
	if (fwdRd1PrimerRev) {
		tdn->reverse_compliment();
		collectStatistics[0]->reversedRds++;//take stats on this
		return true;
	} 

	return false;
}


vector<shared_ptr<DNA>>  Filters::GoldenAxe(vector< shared_ptr<DNA>>& tdn) {
	vector<shared_ptr<DNA>> retDNA(0);
	if (!this->isGoldAxe() || this->isPaired() != 1) {
		return retDNA;
	}
	if (tdn[0]->length() < 100) { return retDNA; }
	int idx = tdn[0]->getBCnumber();
	
	///first check if amplicon wrongly oriented..
	isReversedAmplicon(tdn[0]);

	if (idx < 0) {
		//string presentBC(""); int c_err(0);
		//idx = this->findTag(tdn[0], presentBC, c_err, true,0,true);
		idx = detectCutBC(tdn[0],true);
		if (idx < 0) {
			if (isReversedAmplicon(tdn[0])) {
				//string presentBC(""); int c_err(0);
				//idx = this->findTag(tdn[0], presentBC, c_err, true,0,true);
				idx = detectCutBC(tdn[0], true);
			}
		} 
		/*if (idx >= 0) { //this should be done within "detectCutBC"
			if (this->doubleBarcodes() ) {
				int tagIdx2(-2);
				tagIdx2 = this->findTag2(tdn[0], presentBC, c_err, false, -1);
				this->dblBCeval(idx, tagIdx2, presentBC, tdn[0], nullptr);
				if (idx != tagIdx2) {
					cerr << "GA Double BC eval unsuccesful!\n"; exit(828);
				}
			}
			tdn[0]->setBCnumber(idx, getBCoffset());
		}
		*/
		if (idx >= 0) { tdn[0]->setBCnumber(idx, getBCoffset()); }
	}
	if (idx < 0) {
		tdn[0]->failed();	//retDNA.push_back(tdn[0]);
		//this->addGAstats(tdn[0], retDNA);
		GAstatistics->addBaseGAStats(tdn[0], retDNA, 0);

		return retDNA;
	}
	shared_ptr<DNA> dn = tdn[0];

	int limitF = 0;
	int limitF2 = -1;
	int limitR = 0;
	int SearchL = 6000;

	vector<int> posF(0), posR(0);
	vector<bool> isRC(0), isProblem(0);
	

	//do fwd amplicon search
	while (1) {
		if (limitF2 > limitR) { limitF = limitF2; 
		} else {
			limitF = dn->matchSeq(PrimerL[PrimerIdx[idx]], PrimerErrs, SearchL + limitR, limitR);
		}
		if (limitF < 0 ) { break; }
		limitF2 = dn->matchSeq(PrimerL[PrimerIdx[idx]], PrimerErrs, SearchL + limitF+10, limitF+10);
		//records reverse-searched primer positions
		limitR = dn->matchSeq(PrimerR_RC[PrimerIdx[idx]], PrimerErrs, SearchL + limitF, limitF + 1);
		
		if (limitF2>0 && limitF2 < limitR) {//something went wrong.. couldn't detect correct reverse primer?
			limitF2 = -1;
			isProblem.push_back(true);
		} else {
			isProblem.push_back(false);
		}
		
		if ( limitR < 0) { break; }
		posF.push_back(limitF);	posR.push_back(limitR); isRC.push_back(false);
	}
	//do rev amplicon search
/*	limitF = 0; limitR = 0;
	while (1) {
		limitF = dn->matchSeq(PrimerL_RC[PrimerIdx[idx]], PrimerErrs, SearchL + limitR, limitR);
		//records reverse-searched primer positions
		limitR = dn->matchSeq(PrimerR[PrimerIdx[idx]], PrimerErrs, SearchL + limitF, limitF + 1);
		if (limitF < 0 || limitR < 0) { break; }
		posF.push_back(limitF);	posR.push_back(limitR); isRC.push_back(false);
	}
*/

	//create new DNA objects from each subset..
	int missedGAs(0);
	for (size_t i = 0; i < posF.size(); i++) {
		retDNA.push_back(
			dn->getDNAsubseq(posF[i], posR[i] + PrimerR_RC[PrimerIdx[idx]].length(),
				dn->getId() + "_" + itos(i))
			);
		BCintoHead(idx, retDNA.back(), "", -1, false, true);

		if (i > 0) {
			int disGAs = posF[i] - posR[i - 1] - (int)PrimerR[PrimerIdx[idx]].length();
			if (disGAs > 12) {
				//cerr << ("disGA:"+itos(disGA));
				//std::cout << "disGA" << " ";
				missedGAs++;
			}
		} else if (posF[i] > 100) {
			missedGAs++;
		}

		if (isProblem[i]) {
			dn->setYellowQual(true);
		}
	}
	//int X = 0;

	//collect some stats
	if (retDNA.size() == 0) {//failure to find any GA sequences..
		tdn[0]->failed();	//retDNA.push_back(tdn[0]);
	}
	else if ((GoldAxeMinAmpli != -1 && retDNA.size() < GoldAxeMinAmpli)
				||
			(GoldAxeMaxAmpli != -1 && retDNA.size() > GoldAxeMaxAmpli ) ) {
		tdn[0]->failed();
		retDNA.resize(0);
	}
	//this->addGAstats(dn, retDNA);
	GAstatistics->addBaseGAStats(dn, retDNA, missedGAs);
	//GAstatistics->addMissedGAs(missedGAs);

	return retDNA;

}

bool Filters::swapReverseDNApairs(vector< shared_ptr<DNA>>& tdn){
	if ((!checkRevRd() && !checkSwitchedRdPairs()) || tdn[1] == nullptr) {
		return false;
	}
	int tagIdx = 0;//just try

	//method 1: just check if primer is found reversed, most basic and seems to work fine..
	//simple check if fwd rev primer is in reverse position: then reverse transcribe
	
	bool fwdRd1Primer = checkIfPrimerHits(tdn[0], 0, 0);
	if (fwdRd1Primer) {
		return false;
	}
	if ( checkIfRevPrimerHits(tdn[0], 0, 0)) {
		tdn[0]->reverse_compliment();
		collectStatistics[0]->reversedRds++;
		if (tdn[1] != nullptr) { 
			tdn[1]->reverse_compliment(); 
		//take stats on this
			collectStatistics[1]->reversedRds++;
			return true;
		}
	} else if (checkSwitchedRdPairs() && isPaired() == 2 && tdn[1] != nullptr) {
		//more complex: check if second pair has reversed primer: whole pair swap
		if (checkIfRevPrimerHits(tdn[1], 0, 0)) {//switched pairs && reversed
			tdn[0]->reverse_compliment();
			tdn[1]->reverse_compliment(); 
			swap(tdn[1], tdn[0]);
			tdn[1]->setpairREV();		tdn[0]->setpairFWD();
			collectStatistics[0]->swappedRds++;
			collectStatistics[0]->reversedRds++;
			collectStatistics[1]->reversedRds++;
			//but redundant logging..
			tdn[1]->constellationPairRev(true);
			tdn[0]->constellationPairRev(true);
			return true;
		} else if (checkIfPrimerHits(tdn[1], 0, 0)) {//switched pairs only
			swap(tdn[1], tdn[0]);
			tdn[1]->setpairREV();		tdn[0]->setpairFWD();
			collectStatistics[0]->swappedRds++;
			return true;

		}
	}



	/*
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
		tdn[0]->reverse_compliment();
		tdn[1]->reverse_compliment();
		if (cutPrimer(tdn[0], PrimerIdx[tagIdx], false, 0) ||
			cutPrimerRev(tdn[1], PrimerIdxRev[tagIdx], false)) {
			return true;
		}

		//no? back to normal..
		tdn[0]->reverse_compliment();
		tdn[1]->reverse_compliment();
		return false;
	}
	tagIdx = -2;
	string presentBC = ""; int c_err = 0; int chkRev1=false;
	tagIdx = findTag(tdn[0], presentBC, c_err, true, chkRev1);
	*/

	/*
	if (true && checkReversedRead && (tagIdx2 < 0 && tagIdx < 0)) {
		tdn[0]->reverse_compliment(); tdn[1]->reverse_compliment();
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
		tdn[0]->reverse_compliment();
		MD->analyzeDNA(tdn[0], -1, 0, tagIdx, curThread);
		ch1 = tdn[0]->isGreenQual();
		isReversed = ch1;
		if (!isReversed) {//reset
			tdn[0]->reverse_compliment();
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
	collectStatistics[easyPair]->addPreFilt(d);// PreFilt.addDNAStats(d);
	updateMaxSeqL(d->length());
	//csMTX[easyPair]->unlock();
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
	if (max_l < 0) {
		max_l = (uint)1e9; //100 mil should be good enough for infinite length
	}
}


//ever_best is the best %id_ that was ever observed for this cluster match
bool Filters::betterSeed(shared_ptr<DNAunique> d1,
	shared_ptr<DNAunique> ref,  float ever_best,
	 int usePair, bool checkBC) {
	int TagIdx(0);
	if (checkBC) {
		TagIdx = -2;
	}
	//0.2% difference is still ok, but within 0.5% of the best found seed (prevent detoriating sequence match)
	//float blen = (float)ref->length() + (float)d1->length();

	//*** DNA1
	//needs to quality filter first
	if (!checkYellowAndGreen(d1, usePair, TagIdx, true)) {
		return false;
	}
	if (d1->getPair() != nullptr) {
		checkYellowAndGreen(d1->getPair(), 1, TagIdx, true);
	}
	/*float d1pid(d1->getTempFloat()), refpid(ref->getTempFloat());
	if (d1pid<refpid - 0.4f || d1pid < ever_best - 1){ return false; }
	*/
	//at least 90% length of "good" hit
//	if (d1->length() / ref->length() < RefLengthRatio) { return false; }

	return whoIsBetter(d1, d1->getPair(),d1->getMerge(), 
		ref, ref->getPair(), ref->getMerge(),  ever_best,true);
	
	//checks if the new DNA has a better overall quality
	//1 added to qual, in case no qual DNA is used
	/*
	float thScore = (1+d1->getAvgQual())*(d1pid ) * log((float)d1->length() );
	float rScore = (1+ref->getAvgQual())*(refpid ) * log((float)ref->length() );
	if (thScore > rScore){
		//also check for stable lowest score
		if (d1->minQual() > ref->minQual() - MinQualDiff && (d2 == NULL || ref2 == NULL)) { return true; }
	}
	if (d2 == NULL || ref2 == NULL) {
		return false;
	}
	*/
	//*** DNA2
	//second pair_ likely to be of worse qual_, but only direct comparison relevant here
	
	
	
	/*d2 irrelevant when working primarily with merged reads..
	//at least 90% length of "good" hit
	if (d2->length() / ref2->length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//weigh with average id_ to OTU seed
	thScore +=(1+ d2->getAvgQual()) * log((float)d2->length()) * 97;
	rScore += (1+ref2->getAvgQual()) * log((float)ref2->length()) * 97;
	if (thScore > rScore) {
		return true;
	}
	*/

	return false;
}


//DNA qual_ check, and some extra parameters
//should be safe to call from different threads
bool Filters::checkYellowAndGreen(shared_ptr<DNA> d, int pairPre, 
	int &tagIdx, bool doSeeding) {
	int tagIdx2(-2);
	unsigned int hindrance = 0;
	int pair = max(0, pairPre);//corrects for -1 (undefined pair_) to set to 0

	// We trim for poly-G / homonucleotide tail here, before filtering for length, and before trimming adapter sequences etc. 
	if (trimHomonucleotide != 0){
		unsigned int homoNTnewlength = d->HomoNTTrim(trimHomonucleotide);
		if (homoNTnewlength > 0) {
			d->cutSeqPseudo(homoNTnewlength);
			d->QualCtrl.HomoNTtrimmed = true;
		}
	}

	//remove technical adapter
	if (pairPre == -1 && removeAdapter) {
		remove_adapter(d);
	}
	
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
	
	if (BextensivePrimerChecks) {
		if (checkIfRevPrimerHits(d, PrimerIdx[tagIdx], 0, bShortAmplicons)) {
			d->reverse_compliment(false);
		}
	}
	//bShortAmplicons checks for reverse primer on 1st read
	if (BcutPrimer || Bcheck4illuAdapts) {
		bool fwdRC = false; bool revRC = true;//if this gets changed, the read needs to be reverse complemented
		if (pair != 1  ) {//0 or -1
			cutPrimer(d, PrimerIdx[tagIdx], fwdRC, pair, BextensivePrimerChecks);
			if (bShortAmplicons) {//also check other end of primer.. also use for PacBio amplicons
				cutPrimerRev(d, PrimerIdxRev[tagIdx], revRC, BextensivePrimerChecks);
				//case for 1) long read 2) look for both primers 3) rev required
			}
		} else if (pair != 0) {//pair_ == 1, check for fwd primer in pair_ 2 (rev-compl)
			bool revCheck = pair == -1 || pair == 0;//1:false for RC, else always a reverse check
			cutPrimerRev(d, PrimerIdxRev[tagIdx], revCheck,false);
			if (bShortAmplicons) {//also check other end of primer..
				cutPrimer(d, PrimerIdx[tagIdx], revCheck, pair);
			}
		}
		//conditions for failing read on not finding fwd primer
		if (pair != 1 && !d->getFwdPrimCut() && bRequireFwdPrim) {
			if (alt_bRequireRevPrim) {
				d->QualCtrl.PrimerFwdFail = true;
				d->failed(); return false;
			}
			d->setYellowQual(true);
		}
		//conditions for failing read on not finding rev primer
		if (bRequireRevPrim && !d->getRevPrimCut()
				&& (pair != 0 || isPaired() == 1)  //only pair2 or 1-read-PacBio will fail
			) {//failed to find reverse primer
			if (alt_bRequireRevPrim) {
				d->QualCtrl.PrimerRevFail = true;
				d->failed(); return false;
			} else {
				d->setYellowQual(true);
			}
		}

		if (fwdRC != false || revRC != true) {// matched reverse primer (unexpected).. need to rev compliment
			//d->reverse_compliment(false);
		}
	}


	if (doSeeding) {
		//cut off low qual, hard limits
		d->qualWinPos(EWwidth, EWthr);
		return true;
	}

	//if seq needs to be cut, than here
	if (TruncSeq>0){
		d->cutSeqPseudo(TruncSeq);
	}

	if (check_lengthXtra(d)) {
		d->failed(); return false;
	}

	if (b_doQualFilter) {
		//second cut off low qual_
		d->qualWinPos(EWwidth, EWthr);	// { qualWinTrim = true; }
		//cut off accumulation error larger than maxAccumQP
		if (maxAccumQP > 0.0) {
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
		float avgQ(min_q); 
		if (min_q > 0 || FQWthr > 0){
			avgQ = d->qualWinfloat(FQWwidth, FQWthr, rea);
		}
		if ( avgQ < min_q) {
			d->QualCtrl.AvgQual = true; //sAvgQual(pair_);
			float avgQalt(alt_min_q);
			if (alt_min_q > 0 || alt_FQWthr > 0) {
				avgQalt = d->qualWinfloat(FQWwidth, alt_FQWthr, rea2);
			}
			if (avgQalt < alt_min_q) {
				d->QualCtrl.AvgQual= true; //statAddition.AvgQual++;
				d->failed(); return false;
			} else {
				d->QualCtrl.AvgQual = false;
				d->setYellowQual(true);
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
	if (MaxAmb >= 0 && ambNTs > MaxAmb){
		d->QualCtrl.MaxAmb = true;
		
		if (alt_MaxAmb!=-1 && ambNTs>= alt_MaxAmb){
			d->QualCtrl.MaxAmb = true; //statAddition.MaxAmb++;
			d->failed(); return false;
		} else {
			d->setYellowQual(true);
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
	if (cmdArgs->find("-paired")  != cmdArgs->end() && ((*cmdArgs)["-paired"]=="2" || (*cmdArgs)["-paired"]=="2")){
		pairedSeq = 2; //fakeEssentials();
		noMapTxt += " Using paired end sequencing files.";
	}

	BcutPrimer = false; bDoBarcode = false; bDoBarcode2 = false; bDoBarcode2Rd1 = false;
    removeAdapter=false;bDoMultiplexing=false;
	bDoHeadSmplID=false;
    minBCLength1_ = 0; minBCLength2_ = 0; maxBCLength1_ = 0; maxBCLength2_ = 0; minPrimerLength_ = 0;
	cerr<<noMapTxt<<endl;

	///very similar in principle but easier:
	//needs to correct some parts..
	if (Bcheck4illuAdapts) {
		//BcutPrimer = true;
		fakeEssentials(false);
		PrimerIdxRev.resize(1, 0); PrimerIdx.resize(1, 0);
		if (pairedSeq) {
			string segments = illuPErev;
			trim(segments);	transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
			this->addPrimerR(segments, 0);
			segments = illuPEfwd;
			trim(segments);	transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
			this->addPrimerL(segments, 0);
		}
		else {
			string segments = illuSEuni;
			trim(segments);	transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
			this->addPrimerL(segments, 0);
		}
	}
	else {
		fakeEssentials(true);
	}


}
void Filters::fakeEssentials(bool all){
	//create fake entries
	Barcode.push_back("NA");
	barcodeLengths1_.push_back(0);
	barcodeLengths2_.push_back(0);
	SequencingRun.push_back("");
	SequencingRun2id[""] = vector<int>(1, 0);
	SampleID.push_back("NA"); SampleID_Combi.push_back("NA");
	HeadSmplID.push_back("");bDoHeadSmplID=false;
	collectStatistics[0]->BarcodeDetected.push_back(-1);
	collectStatistics[1]->BarcodeDetected.push_back(-1);
	collectStatistics[0]->BarcodeDetectedFail.push_back(-1);
	collectStatistics[1]->BarcodeDetectedFail.push_back(-1);
	if (all) {
		PrimerIdx.push_back(0); PrimerL.push_back(""); PrimerL_RC.push_back("");
	}
	
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
void Filters::dblBCeval(int& tagIdx, int& tagIdx2, string presentBC, 
	shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2) {
	//bool BCfail = false;// , BCfail2 = false;

	if ( tagIdx < 0 || tagIdx2 < 0 || !tdn->getBarcodeDetected() || 
		(tdn2 != nullptr && !tdn2->getBarcodeDetected()) ) {
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != nullptr) {
			tdn->setPassed(false); /*BCfail = true; */
			tdn->setBCnumber(tagIdx, BCoffset); tdn->setYellowQual(false);
		} 
		if (tdn2 != nullptr) { tdn2->setPassed(false); tdn2->setYellowQual(false); tdn2->setBCnumber(tagIdx2, BCoffset);}
		
		collectStatistics[0]->dblTagFail++;
		return;
	}
	string BC1 = Barcode[tagIdx];
	string BC2 = Barcode2[tagIdx2];
	bool hit(false);
	if (tagIdx == tagIdx2) {
		hit = true;
	} else {
		//this routine finds two matching barcodes (as several combinations are possible)
		for (uint i = 0; i < Barcode.size(); i++) {
			if (Barcode[i] == BC1 && Barcode2[i] == BC2) {
				tagIdx = i; tagIdx2 = i; hit = true; break;
			}
		}
	}

	if ( !hit ) {
		//no BC, useless
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != nullptr) { tdn->setPassed(false); tdn->setYellowQual(false); tdn->setBCnumber(tagIdx, BCoffset);	}
		if (tdn2 != nullptr) { tdn2->setPassed(false); tdn2->setYellowQual(false); tdn2->setBCnumber(tagIdx2, BCoffset);	}
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
/*int Filters::detectCutBC(shared_ptr<DNA> d, string& presentBC, int& c_err, bool isPair1) {
	int start = findTag(d, presentBC, c_err, isPair1,0);

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
*/

//2nd BC on same DNA sequence (from the 3' end)
//deactivated, as now implemented in findTag()
/*
int Filters::findTag2(shared_ptr<DNA> d, string& presentBC, int& c_err,
	bool isPair1, int revChecks) {
	//cerr << "findTag2\n"; exit(2316);
	int start(-1), stop(-1);
	int idx(-1);
	int scanRegion = 34; //dna region to scan for Tag sequence_

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1,false,true);

	if (revChecks) {
	}

	return -1;
}
*/
int Filters::findTag(shared_ptr<DNA> d, string&presentBC, int& c_err, 
			bool isPair1, int revChecks,bool cutBC,bool endCheck) {
    
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
	int scanRegion(4); //dna region to scan for Tag sequence_
	if (!d->getTA_cut() && isPair1) {//no technical adapter found / given by user: scan wider region for barcode
		scanRegion = 14; //arbitary value
	}
	if (d->isMIDseq()) {
		if (d->length() < minBCLength1_) { return -1; }
		scanRegion = d->length() - minBCLength1_ + 1;
	}

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1,false, endCheck);
	
	if (!BCdFWDREV[!isPair1].b_BCdirFix) {
		if (start == -1) {//check reverse transcription
						  //d->reverseTranscribe();
			scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1,true, endCheck);
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
	} else if (idx < 0 && revChecks > 0) {
		scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1,true, endCheck);
		if (idx >= 0) {
			revChecks = 0;
		}
	}
	if (start == -1) {
		idx = -1;
	}
	if (idx != -1) {
		if (cutBC&& BcutTag && !d->isMIDseq()) {
			//remove tag from DNA
			if (endCheck) {
				d->cutSeq(start, -1);
			} else {
				d->cutSeq(0, stop);//start,stop
			}
			d->setBarcodeCut();
			// needs to be locked when multithreaded
			BCintoHead(idx, d, presentBC, c_err, isPair1);
		}
	}

	return idx;
}

//somewhat redundant with findTag function..
int Filters::detectCutBC(shared_ptr<DNA> d, bool isPair1) {
	//seq too short for BC
	if (d->length() < minBCLength1_ ) {
	    return -1;
	}
	//already detected barcode
	if (d->getBarcodeCut()){// && !scndBC) {
		return d->getBCnumber() - BCoffset;
	}
	if ((isPair1 && !bDoBarcode) || (!isPair1 && !bDoBarcode2)) {
		d->setBCnumber(0, BCoffset);
		return BCoffset; //not failed, just not requested
	}

	//ok, really start looking for BC in seq
	int idx(-1);
	if (bDoHeadSmplID){
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



	
	//int start(-1);int stop(-1);
	string presentBC;int c_err(0);
	bool useBC1 = isPair1; bool useBC2 = !useBC1;
	idx = findTag(d, presentBC, c_err, useBC1, 0,true,false);

	/*
	int scanRegion=4; //dna region to scan for Tag sequence_
	
	if (!d->getTA_cut() && isPair1){//no technical adapter found / given by user: scan wider region for barcode
		scanRegion = 22; //arbitary value
	}
	if (d->isMIDseq() || (int)d->length() < scanRegion){
		scanRegion = d->length() - minBCLength1_ + 1;
	}

	// needs to be locked when multithreaded
	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1);


    
	if ( !BCdFWDREV[!isPair1].b_BCdirFix ) {
		if (start == -1){//check reverse transcription
			//d->reverseTranscribe();
			scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1,true);
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
	if (BcutTag && !d->isMIDseq()) {
		//remove tag from DNA
		d->cutSeq(0, stop);
		d->setBarcodeCut();

		// needs to be locked when multithreaded
		BCintoHead(idx, d, presentBC, c_err, isPair1);
	}
	*/

	//check also for reverse BC on same read??? (PacBio)
	if (bDoBarcode2Rd1) {
		string presentBCX(""); int c_errX(0);
		//scan for reverse barcode..
		int idxX = findTag(d, presentBCX, c_err, useBC2, 0, true, true);

		/*
		//int startX(-1), stopX(-1);
		//int scanRegionX = 44;
		//int idxX = idx;
		//scanBC(d, startX, stopX, idxX, c_errX, scanRegionX, presentBCX,false,false,true);
		scanBC_back(d, startX, stopX, idxX, c_errX, scanRegionX, presentBCX, false, true);
		if (idxX >= 0 && BcutTag) {
			d->cutSeq(startX);
		}
		*/


		dblBCeval(idx, idxX, presentBC, d, nullptr);
		c_err = -1;

		//check a second time that barcode was correctly identified, just to be double sure...
		if (idx != idxX ) {
			cerr << "(2) Unequal BC numbers:" << idx << " : " << idxX << "; in object: " << d->getBCnumber() << endl;
			cerr << "In read:" << d->getId() << endl;
			exit(835);
		}

        //cout << "scanBCrev: " << start << "," << stop << "," << idx << "," << c_err << "," << presentBC << endl;
	}
	

	d->setBCnumber(idx, BCoffset);

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
	if (presentBC != "" && (c_err > 0 || atEnd)) {//atEnd: dbl barcode
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


void Filters::reverse_all_BC() {

	if (Barcode2.size() != 0 && revBarcode2.size() == 0) {
		revBarcode2 = Barcode2;
		for (int i = 0; i < (int)Barcode2.size(); i++) {
			reverseTS(revBarcode2[i]);
			//and also prep search vector..
			revBarcodes2_[revBarcode2[i]] = i;
		}
	}
	if (Barcode.size() != 0 && revBarcode.size() == 0) {
		revBarcode = Barcode;
		for (size_t i = 0; i < Barcode.size(); i++) {
			reverseTS(revBarcode[i]);
			revBarcodes1_[revBarcode[i]] = i;
		}
	}

}

/*
void Filters::scanBC_rev(shared_ptr<DNA> d,int& start,int& stop,int& idx,int c_err, 
					 int scanRegion,string & presentBC,
					 bool fwdStrand) {
	if (d->length() < minBCLength1_) { return; }
	//vector<string> emptyV(0), emptyV2(0);
	reverse_all_BC();
	vector<string>* localBarcodesRev;
	if ( !fwdStrand ) {
        localBarcodesRev = &(revBarcode2);
	} else {
        localBarcodesRev = &(revBarcode);
	}
	int BCs = (int) localBarcodesRev->size();	

	//check each possible BC for a match
	if (barcodeErrors_ == 0){
		for (idx=0; idx< BCs; idx++){
			start = d->matchSeq_tot((*localBarcodesRev)[idx], 0, scanRegion, c_err);
			if (start!=-1){
				presentBC = (*localBarcodesRev)[idx];
				stop = start+ (int)(*localBarcodesRev)[idx].length();
				break;
			}
		}
	} else {
		vector<int> stars(0),idxses(0);
		bool zeroErr = false;
		//this version tries all BC's and if there are more than one possible match, will reject all matches
		for (idx=0; idx< BCs; idx++){
			start = d->matchSeq_tot((*localBarcodesRev)[idx], barcodeErrors_, scanRegion, c_err);
			if (start!=-1){
				if (c_err==0){
					stop = start+ (int) (*localBarcodesRev)[idx].length();
					presentBC = (*localBarcodesRev)[idx];
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
			stop = start+(int)(*localBarcodesRev)[idx].length();
			presentBC = d->getSubSeq(start,stop);
		}
	}
	if (start == -1) {
		idx = -1;
	}

}
*/
/*
void Filters::scanBC_back(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err,
	int scanRegion, string& presentBC,
	bool useBC1, bool revBC) {
	if (d->length() < minBCLength1_) { return; }
	//vector<string> emptyV(0), emptyV2(0);
	if (revBC) {
		reverse_all_BC();
	}
	vector<string>* locBC;
	if (!useBC1) {
		if (revBC) {
			locBC = &(revBarcode2);
		} else {
			locBC = &(Barcode2);
		}
	}
	else {
		if (revBC) {
			locBC = &(revBarcode);
		}
		else {
			locBC = &(Barcode);
		}
	}
	int BCs = (int)locBC->size();

	//check each possible BC for a match
	if (barcodeErrors_ == 0) {
		for (; idx < BCs; idx++) {
			start = d->matchSeqRev((*locBC)[idx], 0, scanRegion, c_err);
			if (start != -1) {
				presentBC = (*locBC)[idx];
				stop = start + (int)(*locBC)[idx].length();
				break;
			}
		}
	}
	else {
		vector<int> stars(0), idxses(0);
		bool zeroErr = false;
		//this version tries all BC's and if there are more than one possible match, will reject all matches
		for (; idx < BCs; idx++) {
			start = d->matchSeqRev((*locBC)[idx], (*locBC)[idx].size(), scanRegion, c_err);
			if (start != -1) {
				if (c_err == 0) {
					stop = start + (int)(*locBC)[idx].length();
					presentBC = (*locBC)[idx];
					zeroErr = true;
					break;
				}
				stars.push_back(start);
				idxses.push_back(idx);
			}
		}
		if (!zeroErr && stars.size() > 0) {
			//int pair_ = (int)!fwdStrand;//d->getReadMatePos();
			if (stars.size() > 1) {//too many matches, thus true seq can't be found
				//currently only have only one BC, could be changed in future
				//sTagNotCorrected(pair_);
				d->QualCtrl.fail_correct_BC = true;
				idx = -1; start = -1;
				return;
			}
			d->QualCtrl.suc_correct_BC = true;
			//sTagCorrected(pair_);// collectStatistics.suc_correct_BC++;
			start = stars[0];
			idx = idxses[0];
			stop = start + (int)(*locBC)[idx].length();
			presentBC = d->getSubSeq(start, stop);
		}

	}
	if (start < 0) {
		idx = -1;
	}
}

*/

void Filters::scanBC(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err,
	int scanRegion, string& presentBC, bool fwdStrand, bool revBC, bool endScan ) {
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
	if (revBC) {
		reverse_all_BC();
	}


    BarcodeMap* localBarcodes = &emptyBarcodes;
    vector<int> *localBarcodeLengths;


    // was locked with an omp pragma
    unsigned int maxBCLength;
    unsigned int minBCLength;
	if (fwdStrand) {
		if (revBC) {
			localBarcodes = &revBarcodes1_;
		}else {
			localBarcodes = &barcodes1_;
		}
        localBarcodeLengths = &barcodeLengths1_;
        maxBCLength = maxBCLength1_;
        minBCLength = minBCLength1_;
    } else {
		if (revBC) {
			localBarcodes = &revBarcodes2_;
		}else { localBarcodes = &barcodes2_; }
        localBarcodeLengths = &barcodeLengths2_;
        maxBCLength = maxBCLength2_;
        minBCLength = minBCLength2_;
    }
	if (d->length() < maxBCLength) {
		return ;
	}

	int seqLen = d->mem_length();
    
//    cout << "minBCLength: " << minBCLength << endl;
//    cout << "maxBCLength: " << maxBCLength << endl;
//    cout << "barcodeErrs: " << barcodeErrors_ << endl;
//
        
    //bool found = false;


    for (int start2 = 0; start2 < scanRegion; start2++) {
		start = endScan ? (seqLen - start2 - maxBCLength) : start2;
		if (start < 0) {break;}
		//cerr << start << " ";
        const string test = d->getSubSeq(start, maxBCLength);
        auto barcodeIterator = localBarcodes->find(test);

        if (barcodeIterator != localBarcodes->end()) {// set index if found
            idx = (*barcodeIterator).second;
            stop = start + (int) (*localBarcodeLengths)[idx];
            presentBC = test; // (*locBC)[idx];
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
    
	if (start == -1) {
		idx = -1;
	}

    
    return;
}

bool Filters::checkIfRevPrimerHits(shared_ptr<DNA> d, int primerID, int pair,bool twoPrimers) {
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0) { return true; }
	if (d->getFwdPrimCut()) {
		return true;
	}
	int start(-1), stop(-1);
	int tolerance(30); //, startSearch(0);
	int QS = d->length(); int limit = max(QS >> 1, QS - 150); stop = QS;
	if (pair == 1) {
		start = d->matchSeqRev(PrimerR_RC[primerID], PrimerErrs, limit);
		if (start != -1) {
			return true;
		}
	}
	if (pair!=1 || twoPrimers) {
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit);
	}

	if (start != -1) {
		return true;
	}

	return false;

}
bool Filters::passedReads(int n) {
	if (collectStatistics[0]->total - passed_interval_reads  >(uint) n) {
		passed_interval_reads = collectStatistics[0]->total;
		return true;
	}
	return false;
}

bool Filters::checkIfPrimerHits(shared_ptr<DNA> d, int primerID, int pair) {
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0) { return true; }
	if (d->getFwdPrimCut()) {
		return true;
	}
	int start(-1), stop(-1);
	int tolerance(30), startSearch(0);
	int QS = d->length(); 
	int limit = max(QS >> 1, QS - 150); stop = QS;
	if (QS < 20) {
		return false;
	}
	if (pair == 1) {
		start = d->matchSeq(PrimerR[primerID], PrimerErrs, tolerance, startSearch);
	}
	else {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, tolerance, startSearch);
	}

	if (start != -1) {
		return true;
	}

	return false;

}

//cuts primers, tags
bool Filters::cutPrimer(shared_ptr<DNA> d,int primerID,bool& RC,int pair, bool extensivePrimerCheck){
	//only adapted to singular BC
	if (PrimerL.size()==0 || PrimerL[0].length()==0){return true;}
	if (d->getFwdPrimCut()) {
		return true;
	}
	int start(-1) ,stop(-1);
	int limit(50), startSearch(0),limit2(50);
	if (Bcheck4illuAdapts) {
		limit = 40;
	} else	if (!d->getBarcodeCut() && maxBCLength1_ > 0) {
		limit = maxBCLength1_ + 4;
	} else { limit = 22; }//in this case nothing is known about 5' end
	if (extensivePrimerCheck) { 
		int QS = (int)d->length(); limit = min(QS - 1, max(((int)((float)d->length() * 0.75)), 250));
		limit2 = (d->length() - 100);;
	}

	if (!BcutTag){
		//Tag was not cut out of sequence_, take this into account
		startSearch = minBCLength1_ - 2;
		limit += (maxBCLength1_ - minBCLength1_) + 4;
	}
	if (!RC) {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, limit, startSearch);
		stop = start + (int)PrimerL[primerID].length();
	} else {
	//if (1 && start == -1){
		int QS = d->length();int limit = max(int(QS *0.75), QS - 250); stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
//		start = d->matchSeq(PrimerL_RC[primerID], PrimerErrs, tolerance, startSearch);
		if (0 && start != -1) {
			int y = 0;//debug
		}
	}
	if (start < 0){//failed to match primer
		//this should only be set if it led to failing to read
		//d->QualCtrl.PrimerFwdFail = true;
		//sPrimerFail(pair_);// max(0, (int)d->getReadMatePos()));
		if (alt_PrimerErrs!= 0 && PrimerErrs < alt_PrimerErrs){
			if (!RC) {
				start = d->matchSeq(PrimerL[primerID], alt_PrimerErrs, limit, startSearch);
				stop = start + (int)PrimerL[primerID].length();
			} else {
				int QS = d->length();  stop = QS;
				start = d->matchSeqRev(PrimerL_RC[primerID], alt_PrimerErrs, limit2, startSearch);
			}
		}
		if (start < 0){
			//statAddition.PrimerFail++;
			return false;
		}else if (pair!=1){ //2nd read shouldnt be affected by fwd primer (but still checked in short read mode)
			d->setYellowQual(true);
		}
	}
	//remove everything before/after primer cut
	if (!BcutPrimer){
		if ( !RC ) { d->cutSeq(0, start); }
		else { d->cutSeq(stop,-1); }
		return true;
	} else {
		//remove the primer, if confimed before
		if (!RC) { d->cutSeq(0, stop); }
		else { d->cutSeq(start, -1); }
	}
	d->setFwdPrimCut();
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
		int QS = d->length(); int limit = max(QS / 2, QS - 150); //stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
	}
	if (start == -1) {//failed to match primer
		return false;
	}
	return true;
}
bool Filters::cutPrimerRev(shared_ptr<DNA> d,int primerID,bool& RC,bool extraCheck){
	//const string& se = d->getSequence();
	if (d->getRevPrimCut() || PrimerR.size()==0){
		return true;
	}

	int start(-1) ,stop(d->length());
	int QS = d->length(); 
	int limit=max(QS>>1,QS-250);
	int limit2 = 50;
	if (true) {
		limit = min(QS - 1, max(((int)((float)d->length() * 0.75)), 250));
		limit2 = (d->length() - 100);
	}

	//int limit = QS>>1;

	if (RC) {
		start = d->matchSeqRev(PrimerR_RC[primerID], PrimerErrs, limit, extraCheck);
	} else {
		start = d->matchSeq(PrimerR[primerID], PrimerErrs, limit2, 0, extraCheck);
		stop = start + (int)PrimerR[primerID].length();
		if (start >= 0) {
			RC = false;
		}
	}

	if (start < 0){//failed to match primer
		//d->QualCtrl.PrimerRevFail = true;
		return false;
	} 

	if ( !BcutPrimer ) { //found it, but no cut
		if ( !RC ) { d->cutSeq(0, start); }
		else { d->cutSeq(stop,-1); }
		return true; 
	} else {
		//remove the primer, if confimed before
		if (!RC) {
			d->cutSeq(0, stop);//start  everything in front has to be removed
		}
		else {
			d->cutSeq(start, -1); // everything in the end has to be removed
		}
	}
	//string neSe = se.substr(0,start) + se.substr(stop);
	d->setRevPrimCut();

	return true;
}
bool Filters::readMap(){//core routine to read map info
	if (cmdArgs->find("-map")  == cmdArgs->end()){
		this->noMapMode();
		return true;
	}
	string MapF = (*cmdArgs)["-map"];
	if (cmdArgs->find("-optimalRead2Cluster") != cmdArgs->end()){
		bDoSeedExtension = true;
	}



	string path = ""; bool pathMode = false;
	if (cmdArgs->find("-i_path")  != cmdArgs->end() && (*cmdArgs)["-i_path"].length() >= 1){
		path=(*cmdArgs)["-i_path"] + string("/");
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

		//(*cmdArgs)["-i_MID_fastq"]
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
		int idx = d->getBCnumber() - BCoffset;
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

void Filters::addMergeStats(OutputStreamer* out) {
	mergeStats->BPwritten = (uint)out->getBPwrittenInSR();
	mergeStats->BPmergeWritte =(uint) out->getBPwrittenInSRmerg();
	mergeStats->total_read_preMerge_ = (int)out->total_read_preMerge_;
	mergeStats->merged_counter_ = (int)out->merged_counter_;



}
void Filters::prepStats() {
	float remSeqs = float(collectStatistics[0]->total - collectStatistics[0]->totalRejected);
	for (size_t i = 0; i < 2; i++) {
		collectStatistics[i]->PostFilt.calcSummaryStats(remSeqs, min_l, min_q);
		if (bAdditionalOutput) {
			remSeqs = float(statAddition[0]->total - statAddition[0]->totalRejected);
			statAddition[0]->PostFilt.calcSummaryStats(remSeqs, min_l, min_q);//yellow

		}
		collectStatistics[i]->PreFilt.calcSummaryStats(1, min_l, min_q);
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

void Filters::printLenVsQual(ostream& give) {
	collectStatistics[0]->PreFilt.printLvsQ(give);
}

void Filters::printHisto(ostream& give,int which, int set){
	bool p2stat = pairedSeq > 1 ;

	if (set == 0) {
		vector<uint> colStats(  collectStatistics[0]->PostFilt.get_rstat_Vmed(which));
		vector<size_t> ra(collectStatistics[0]->PostFilt.getVrange(which) );

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
		ra = collectStatistics[0]->PostFilt.getVrange(which);
		if (p2stat) { tra = collectStatistics[1]->PostFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		if (!skips[2]) {
			tra = statAddition[0]->PostFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
			if (p2stat) { tra = statAddition[0]->PostFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		} 
		tra = collectStatistics[0]->PreFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
		if (p2stat) { tra = collectStatistics[1]->PreFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		vector<uint> empt(ra[1], 0);
		matHist = vector<vector<uint>>(6, empt);
		for (size_t kk = 0; kk < 6; kk++) {
			if (skips[kk]) { continue; }
			switch (kk) {
			case 0:	stat = collectStatistics[0]->PostFilt.get_rstat_Vmed(which); break;
			case 1:	stat = collectStatistics[1]->PostFilt.get_rstat_Vmed(which); break;
			case 2:	stat = statAddition[0]->PostFilt.get_rstat_Vmed(which); break;
			case 3:	stat = statAddition[1]->PostFilt.get_rstat_Vmed(which); break;
			case 4:	stat = collectStatistics[0]->PreFilt.get_rstat_Vmed(which); break;
			case 5:	stat = collectStatistics[1]->PreFilt.get_rstat_Vmed(which); break;
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
void Filters::FileEssentials(filesStr& files, OptContainer* cmdArgs){//unordered_map<string, int>& UFF) {
	
	files.FastaF = this->getFastaFiles();
	files.QualF = this->getQualFiles();
	files.FastqF = this->getFastqFiles();
	files.MIDfq = this->getMIDfqFiles();



	//set up some log structures
	files.deLog = "";//dereplication main log
	files.logF = (*cmdArgs)["-log"];
	files.logFA = (*cmdArgs)["-log"].substr(0, (*cmdArgs)["-log"].length() - 3) + "add.log";


	// Set folder path
	if (cmdArgs->find("-i_path") != cmdArgs->end() && (*cmdArgs)["-i_path"].length() > 2) {
		files.path = (*cmdArgs)["-i_path"] + string("/");
	}

	// Set up b_derep_as_fasta_ or fastq way and save file vector in tar in case it is zipped
	if (files.FastaF.size() > 0) { // If b_derep_as_fasta_ vector contains elements
		files.fastXtar = files.FastaF;
		files.isFastq = false; // Set boolean Fastq to false
	}
	else { // If no files.FastaF present assume there are Fastq files
		files.fastXtar = files.FastqF;
		if (files.FastqF.size() == 0) { // no Fasta and no Fastq files -> abort
			cerr << "No FastQ or Fasta file given.\n  Aborting..\n";
			exit(12);
		}
	}



	// We dont know if it is a tar yet, but we call it tar
	//this routine is important for managing the blocks of files to be read together
	for (unsigned int i = 0; i < files.fastXtar.size(); i++) {
		bool suc = false;
		auto XX = files.uniqueFastxFiles.find(files.fastXtar[i]);
		if (XX == files.uniqueFastxFiles.end()) {//no entry for this fastq yet
			files.uniqueFastxFiles[files.fastXtar[i]] = (int)files.uniqueFastxFiles.size();
			files.idx.push_back(vector<int>(1, i));
		}
		else {//exists already..
			files.idx[XX->second].push_back(i);
		}

	}


	//files.uniqFxFls = mapToVector(files.uniqueFastxFiles);
	
	
	
	
	
	
	
	
	
	//is SeqRun covered at all?
	if (SequencingRun.size() < files.uniqueFastxFiles.size()) { 
		SequencingRun.resize(files.uniqueFastxFiles.size(), ""); 
	}
	//transfer  uniqueFastxFiles to vector with SR info
	vector<pair<string, string>>SR2File;
	for (auto uFX : files.uniqueFastxFiles) {
		int tarID = files.idx[uFX.second][0];
		string SR = this->SequencingRun[tarID];
		pair<string, string> tmp (SR, uFX.first);
		SR2File.push_back(tmp);
	}
	//sort vector
	std::sort(SR2File.begin(), SR2File.end());

	for (auto fx : SR2File) {
		pair<string, int> tmp(fx.second, files.uniqueFastxFiles[fx.second]);
		files.uniqFxFls.push_back(tmp);
	}

	if (files.uniqFxFls.size() != files.uniqueFastxFiles.size()) {
		cerr << "Wrong size files.uniqFxFls vs files.uniqueFastxFiles\nAborting..\n";
		exit(623);
	}


	//unique Fas files set up.. check for their existence
	shared_ptr<InputStreamer> testFiles =
		make_shared<InputStreamer>(!files.isFastq, this->getuserReqFastqVer(), "1", "1", 1);
	// For each unique Fa file, create to see if path etc are right
	for (auto uFX : files.uniqFxFls) {
		int tarID = files.idx[uFX.second][0]; string tmp;
		string x = testFiles->setupInput(files.path, tarID, uFX.first, files, this->setPaired(), (*cmdArgs)["-onlyPair"], tmp, true);
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
	//char buffer[50];
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
	collectStatistics[0]->PreFilt.printGCstats(os);
	if ( Npair > 1 ) {
		os << "R2 pre-filter";
		collectStatistics[1]->PreFilt.printGCstats(os);
	}
	if ( !b_doFilter ) {return;}
	os << "R1 filtered";
	collectStatistics[0]->PostFilt.printGCstats(os);
	if ( Npair > 1 ) {
		os << "R2 filtered";
		collectStatistics[1]->PostFilt.printGCstats(os);
	}
}

void Filters::printStats(ostream& give, string file, string outf, bool greenQualStats) {
	//TODO switch min_l to min_l_add
	shared_ptr<collectstats> cst = collectStatistics[0];
	shared_ptr<collectstats> cst2 = collectStatistics[1];
	if (cst2->total != cst->total) {
		//cerr << "Unequal read numbers recorded " << cst->total << "," << cst2->total << endl;
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
	/*if (!greenQualStats) {
		give << "Statistics of reads that passed the mid qual filter\n";
	} else {
		give << "Statistics of high quality reads\n";
	}*/
	float remSeqs = float (cst->total-cst->totalRejected);
	
	give << endl;

	string ReadTag = "Reads";
	if (isGoldAxe()) { ReadTag = "(Sub)Reads"; }
	
	if (!greenQualStats){
		give << ReadTag<<" not High qual_: " << intwithcommas((int)cst->totalRejected);
	} else {
		give << "Reads processed: " << intwithcommas((int)cst->total);
		if (p2stat) {
			give << "; " << intwithcommas((int)cst2->total) << " (pair 1;pair 2)";
		}
	}
	give << endl;
	if (cst->reversedRds > 0 || cst->swappedRds > 0) {
		if (cst->reversedRds > 0) {
			give << cst->reversedRds;
			if (p2stat) {				give << ";" << cst2->reversedRds;			}
			give<<  " reads reverse-translated";
			if (cst->swappedRds > 0) {
				give << ", ";
			}
		}
		if (cst->swappedRds > 0) {
			give << cst->swappedRds << " read pairs swapped";
		}
		give << endl;
	}

	//int numAccept = (int)(cst->total - cst->totalRejected);
	int numAccept = (int)(cst->totalSuccess - cst->totalMid);
	int numMid = (int)cst->totalMid;
	if (!greenQualStats){
		give << "Rejected:" << intwithcommas((int)(cst->totalRejected)) << endl;
		give << "Accepted (High qual): " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst->Trimmed) << " end-trimmed)\n";
		give << "Accepted (Mid qual): " << intwithcommas(numMid) << endl;
	} else {
		if (p2stat) {
			int numMid2 = (int)cst2->totalMid;
			int numAccept2 = int(cst2->totalSuccess - cst2->totalMid);
			give << "Rejected: " << intwithcommas((int)cst->totalRejected) << "; " << intwithcommas((int)cst2->totalRejected) << endl;
			give << "Accepted (High qual): " << intwithcommas((int)numAccept) << "; " << intwithcommas((int)numAccept2) << " (" << intwithcommas((int)cst->Trimmed) << "; " << intwithcommas((int)cst2->Trimmed) << " end-trimmed)\n";
			give << "Accepted (Mid qual): " << intwithcommas(numMid) << ";" << intwithcommas(numMid2)<< endl;
		}
		else {
			give << "Rejected: " << intwithcommas((int)cst->totalRejected) << endl;
			give << "Accepted (High qual): " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst->Trimmed) << " end-trimmed)\n";
			give << "Accepted (Mid qual): " << intwithcommas(numMid) << endl;
		}
	}


	if ( false && bPrimerR ) { //confusing collectStatistics
		give << ", with rev. primer: " << intwithcommas((int)cst->RevPrimFound); if ( p2stat ) { give << "; " << intwithcommas((int)cst2->RevPrimFound); }
	}

	if (pairedSeq>1) {
		give <<"Singletons among these: " << intwithcommas((int)cst->singleton) << "; " << intwithcommas((int)cst2->singleton) << endl;
	}
	give << "Bad Reads recovered with dereplication: " << intwithcommas((int)cst->DerepAddBadSeq) << endl;

	if ( bShortAmplicons ) {
		give << "Short amplicon mode.\n";
	}
	if ( checkSwitchedRdPairs() ) {
		//give << "Looked for switched read pairs (" << intwithcommas(revConstellationN) << " detected)" << endl;
	}
	if (greenQualStats) {
		cst->PostFilt.printStats2(give, remSeqs,0);
		collectStatistics[1]->PostFilt.printStats2(give, remSeqs,1);
	} else {
		statAddition[0]->PostFilt.printStats2(give,remSeqs,0);
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

	give << "  > (" << trimHomonucleotide << ") homo-nt trimmed  : " << spaceX(17 - digitsInt(maxHomonucleotide)) << intwithcommas((int)cst->HomoNTtrimmed);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->HomoNTtrimmed); } give << endl;

	give << "Rejected due to:\n";
	float val = (float)min_l;
	if (val == -1.f) {val = min_l_p;}
	if (!greenQualStats){ val = (float)alt_min_l; }

	give << "  < min Sequence length (" << val << ")  : " << spaceX(18 - digitsFlt(val)) << intwithcommas((int)cst->minL);
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
	give << "  > max Sequence length (" << max_l << ")  : " << spaceX(18 - digitsInt(max_l)) << intwithcommas((int)cst->maxL);
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
	if (bDoMultiplexing) {
		if (bDoBarcode) {
			give << "  -Barcode unidentified (max " << barcodeErrors_ << " errors) : " << spaceX(19 - digitsInt(barcodeErrors_)) << intwithcommas((int)cst->TagFail);
			if (p2stat && (cst2->TagFail > 0 || doubleBarcodes())) { give << "; " << intwithcommas((int)cst2->TagFail); give << " (" << intwithcommas((int)cst->dblTagFail) << " pairs failed)"; }
			give << endl;

			if (barcodeErrors_ > 0) {
				give << "    -corrected barcodes: " << spaceX(18) << intwithcommas((int)cst->suc_correct_BC);
				if (p2stat) { give << "; " << intwithcommas((int)cst2->suc_correct_BC); }
				give << endl;
				//<< ", failed to correct barcode: " << spaceX(5 - digitsInt(FQWwidth)) << intwithcommas((int)cst->fail_correct_BC) << endl;
			}

			if (bDoBarcode2) {
				give << "    -used dual index barcodes";
				if (BCdFWDREV[0].reversedBCs || BCdFWDREV[1].reversedBCs) {
					give << " (reversed_ ";
					if (BCdFWDREV[1].reversedBCs && BCdFWDREV[0].reversedBCs) {
						give << " fwd & rev";
					}
					else	if (BCdFWDREV[0].reversedBCs) {
						give << " fwd";
					}
					else	if (BCdFWDREV[1].reversedBCs) {
						give << " rev";
					}
					give << " BCs)" << endl;
				}

			}
			else if (BCdFWDREV[0].reversedBCs) {
				give << "    -reversed_ all barcodes" << endl;
			}


		}
		else if (bDoHeadSmplID) {
			give << "  -Failed to assign sequences to header tag : " << intwithcommas((int)barcodeErrors_) << endl;
		}
	}
	
	
	
	if (isGoldAxe()) {
		GAstatistics->setBCs(SampleID, Barcode, Barcode2);
		GAstatistics->printSummary(give);
		//GAstatistics->printBCtabs(give);
	}

	mergeStats->print(give);


	if (bDoMultiplexing) {

		if (bDoBarcode) {
			give << endl << "SampleID";
			if (bDoCombiSamples) {
				give << "\tSampleGroup";
			}
			give << "\tBarcode";
			if (bDoBarcode2) { give << "\tBarcode2"; }
			give << "\tInstances\n";
			for (unsigned int i = 0; i < Barcode.size(); i++) {
				give << SampleID[i] << "\t";
				if (bDoCombiSamples) { give << SampleID_Combi[i] << "\t"; }
				give << Barcode[i];
				if (bDoBarcode2) {
					give << "\t" << Barcode2[i];
				}
				give << "\t" << intwithcommas((int)cst->BarcodeDetected[i]) << endl;
			}

		} else if (bDoHeadSmplID) {
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

	GAstatistics->addStats(fil->GAstatistics);
	mergeStats->addStats(fil->mergeStats);
	//revConstellationN += fil->revConstellationN;
}




////////////  UCLINKS  ///////////////
UClinks::~UClinks(){
	ucf.close();
	mapdere.close();
	if (merger != nullptr) {
		delete merger;
		merger = nullptr;
	}
}
UClinks::UClinks( OptContainer* cmdArgs):
	CurSetPair(-1),//maxOldDNAvec(20000),
	DNAunusedPos(0), derepMapFile(""),
	bestDNA(0, NULL), oriKey(0), 
	mapLines(0),
	bestPID(0), bestLEN(0),
	clusCnt(0), uclines(0),
	SEP(""), 
	UCread(false), pairsMerge(false), MAPread(false),
	b_derepAvailable(false),
	UPARSE8up(false), UPARSE9up(false), UPARSE11up(false),
	UpUcFnd(false), cdhit(false), vsearch(false),
	ucispaf(false),
	otuTerm("OTU"), otuOUTterm("OTU_"),
	RefDBmode(false), RefDBotuStart(-1),
	SeedsAreWritten(false),
	OTUmat(0), unregistered_samples(false),
	doChimeraCnt(false), OTUnumFixed(true),
	totalDerepCnt(0),
	qCovThr(0.8f), perIDmatch(97.f),
	b_merge_pairs_optiSeed_(false),merger(nullptr)

{
	//read in UC file and assign clusters
	SEP = (*cmdArgs)["-sample_sep"];
	string str = (*cmdArgs)["-optimalRead2Cluster"];
	

	if (str.find(".paf") != string::npos) {
		ucispaf = true; UpUcFnd = true;
	}

	ucf.open(str.c_str(),ios::in);
	if (!ucf){
		UCread = true;
		cerr<<"Could not open uc file\n"<< str<<endl; exit(46);
	}
	if (cmdArgs->find("-derep_map") != cmdArgs->end()) {
		derepMapFile = (*cmdArgs)["-derep_map"];
	}
	if (cmdArgs->find("-minQueryCov") != cmdArgs->end()) {
		qCovThr = (float) atof((*cmdArgs)["-minQueryCov"].c_str());
	}
	if (cmdArgs->find("-id") != cmdArgs->end()) {
		perIDmatch =(float) atof((*cmdArgs)["-id"].c_str());
	}
	if (cmdArgs->find("-count_chimeras") != cmdArgs->end() &&
		(*cmdArgs)["-count_chimeras"] == "T") {
		doChimeraCnt = true;
	}
	if (cmdArgs->find("-merge_pairs_seed") != cmdArgs->end() &&
		(*cmdArgs)["-merge_pairs_seed"] == "1") {
		b_merge_pairs_optiSeed_ = true;
	}

	//set version of mapping/clustering
	if ((*cmdArgs)["-uparseVer"] != "") {
		if ((*cmdArgs)["-uparseVer"] == "N11") {//UNOISE v11
			UPARSE8up = true;
			otuTerm = "Zot";
			otuOUTterm = "Zotu";
		}
		else if ((*cmdArgs)["-uparseVer"] == "cdhit") {
			cdhit = true; UpUcFnd = true; UPARSE8up = true;
		} else if ((*cmdArgs)["-uparseVer"] == "vsearch") {
			UpUcFnd = true; vsearch = true;
		} else {
			int upVer = atoi((*cmdArgs)["-uparseVer"].c_str());
			if (upVer >= 8 && upVer < 9) {
				UPARSE8up = true; UpUcFnd = true;
			}
			else if (upVer >= 9 && upVer < 11) {
				UPARSE8up = true; UPARSE9up = true; UpUcFnd = true;
			}
			else if (upVer == 9966) {//cdhit code
				cdhit = true; UpUcFnd = true;
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

		//(*cmdArgs)["-i_MID_fastq"]
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
	cout << "Found " << SmplIDs.size() << " samples in derep.map\n";
}

//read in dereplicated info from derep.map and derep.hq.fq (in IS) -> they are in the same order
int UClinks::oneDerepLine(shared_ptr<DNAunique> d) {
	
	if (MAPread){
		return 0;
	}
	string line("");
	int cnt(0);
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
		totalDerepCnt += curCnt2;

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

	cerr << "Found " << totalDerepCnt<< " counts in derep.map, "<< finishMapCnt<<" counts secondary\n";

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
	string segs2("");
	float perID;
	vector<int> curCLID(0,0);
	//bool sync(false); // syncing of 2 read pairs; not implemented for this function yet


	//goes through UC file
	while ( getMAPPERline(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
		//int subcnt = 0;
		if ( curCLID.size() == 0 ) { continue; }
		if ( uclInOldDNA(segs, curCLID, perID, fil) ) {
			curCLID.resize(0);
			continue;
		}
		while(cont){//match to DNA
			shared_ptr<DNA> tmpDNA = IS->getDNA(0);
			match2 = IS->getDNA(1);
			if (tmpDNA == NULL) { cont = false; break; }//signal that at end of file
			
			//does some basic quality check, such as removing primers, reversing etc
			fil->miniCheckDNA(tmpDNA, match2);


			match = make_shared<DNAunique> (tmpDNA, -1);
			if (match2 != NULL) {
				match->attachPair(make_shared<DNAunique>(match2, -1));
				//merge read
				merger->findSeedForMerge(match, match2);
				if (match->merge_seed_pos_ > 0) {
					shared_ptr<DNA> dM = merger->merge(match, match2);
					if (dM != nullptr) {
						match->attachMerge(make_shared<DNAunique>(dM, -1));
					}
				}

			}

			//assummes in original implementation, that we can get derep.map lines with the same 
			//ordering as fq derep (which works normally, just not for dada2 mode)
			UCcnt+= oneDerepLine(match);
			string curID = match->getId();
			curID = curID.substr(0,curID.find_first_of(' '));
			//check if dnaTemp1 a) matches id_ b) is better
			if (curID != segs){
				//block to store unused DNA & find this id_ in this block
				unusedID[curID] = DNAunusedPos;
				oldDNA[DNAunusedPos] = match;	//oldDNA2[DNAunusedPos] = match2;
				DNAunusedPos++;
			} else {//find if current DNA rep needs to be replaced with a better read (pair)
				//assign % identity score to DNA object
#ifdef DEBUG
				cerr << "UC Hit";
#endif
				match->setTempFloat(perID);
				besterDNA(curCLID, match, match->getPair(), fil);
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

	while ( getMAPPERline(segs, segs2, perID, curCLID, !b_derepAvailable) ) {
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
	ucispaf = false;
	if (addUC.find(".paf") != string::npos) {
		ucispaf = true; UpUcFnd = true;
	}
	ucf.open(addUC.c_str(), ios::in);
	if (!ucf) {
		UCread = true;
		cerr << "Could not find additional uc file: " << addUC << endl;
	} else {
		std::cerr << "Reading " << addUC << endl;
	}
	UCread = false;
	if (SmplHd){
		while (getMAPPERline(segs, segs2, perID, curCLID, SmplHd)) { 
			if (curCLID.size() > 0) { cntsAddUC++; }
			curCLID.resize(0); 
		}
		
	}else {//complicated..
		int cnt(0);
		while (getMAPPERline(segs, segs2, perID, curCLID, SmplHd)) { 
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
			countsAdd += (int)matchSiz;
		}
		//was matched once to an OTU seed. Even if several later matches, doesn't mater - delete

		unusedID.erase(unusedIT);
//		delete oldDNA[mID]; if (oldDNA2[mID] != NULL){ delete oldDNA2[mID]; }

		oldDNA.erase(mID); //oldDNA2.erase(mID);// = NULL;oldDNA2[mID] = NULL;
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
		besterDNA(curCLID, oldDNA[mID], oldDNA[mID]->getPair(), fil);
		//remove all trace
		unusedID.erase(unusedIT);
		oldDNA.erase(mID); //oldDNA2.erase(mID);
		//oldDNA[mID] = NULL;		oldDNA2[mID] = NULL;
		return true;
	}
	return false;
}

void UClinks::resetMarks() {
	if (cdhit || vsearch) {
		UpUcFnd = false;
		cdhit = false; vsearch = false;
	}
}

bool UClinks::getTMPmapperLine(string& line) {
	if (mapLines.size() > 0) {
		line = mapLines.front();
		mapLines.pop_front();
	}
	else {
		if (!getline(ucf, line, '\n')) {
			return false;
		}
	}
	return true;
}

bool UClinks::getMAPPERline(string& segs, string& segs2,float& perID, 
	vector<int>& curCLID,  bool addFromHDstring) {
	//converts derep reads mapping to OTUs to counts
	//seg: read ID, segs2: OTU repr read, perID: match to repr read


	//close all file streams
	if (UCread) { resetMarks(); UpUcFnd = false; return false; }
	if (ucf.eof()){
		resetMarks();
		UCread=true;		ucf.close();		return false;
	}
	float qCov = 1.f;
	string line; 
	std::unordered_map<string, int>::iterator itCL;

	//check if any lines need to be backworked
	while (this->getTMPmapperLine(line)) {
	//	cerr<<line<<endl;
		uclines++;
		if (line.length() <= 1){ continue; }	
		if (!UpUcFnd){
			//reset
			vsearch = false;cdhit = false; UPARSE8up = false; UPARSE9up = false; UPARSE11up = false;
			if (line[0] == '>') {
				cdhit = true;
			} else if (line.substr(1, 1) == "\t"){ UPARSE8up = false; 
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
		bool chimera = false;
		vector<string>tarsV; //saves hits to OTUs

		if (cdhit) {
			if (line[0] == '>') {
				//get next line with actual mapping/clustering
				getline(ucf, line, '\n');
				repFound = false;
				segs2 = "";
			}

			//meachnism has to be slightly different for cdhit, have to save segs2 between lines
			if (!repFound  ){//&& line.back() == '*') {//report representative gene
				bool isRep = line.back() == '*';
				while (!isRep) {
					if (!isRep) {
						mapLines.push_back(line);
					}
					getline(ucf, line, '\n');
					isRep = line.back() == '*';
					if (line[0] == '>') {//should never get here
						cerr << "CD-HIT .clstr contained seq cluster with no rep (prev lines):\n" << line << endl;
						exit(163);
					}
				}
				size_t pos = line.find("nt, >");size_t pos2 = line.find("...", pos + 4);
				segs = line.substr(pos + 5, pos2 - pos - 5);
				repFound = true;
				segs2 = segs;
				perID = 100.f;
				
			} else {
				size_t pos = line.find("nt, >");size_t pos2 = line.find("...", pos + 4);
				segs = line.substr(pos + 5, pos2 - pos - 5);
				//find ID
				pos = line.find("at +/", pos2 + 3);
				string ID = line.substr(pos + 5, line.length() - pos - 6);
				perID = stof(ID);
				//cerr << line << " " << ID << endl;
			}
		} else {

			stringstream ss;
			ss << line;
			//2 ways to get to a) hit info b) query & otu
			if (ucispaf) {//minimap2 .paf file..
				getline(ss, segs, '\t');
				string query = segs;
				getline(ss, segs, '\t');
				float qLen = (float)atof(segs.c_str());
				for (uint i = 0; i < 3; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs2, '\t');//tarname
				getline(ss, segs, '\t');
				float tLen = (float) atof(segs.c_str());
				for (uint i = 0; i < 2; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs, '\t');
				float matches = (float) atof(segs.c_str());
				getline(ss, segs, '\t');
				float mapLen = (float) atof(segs.c_str());
				perID = matches / mapLen *100.f;

				qCov = mapLen / qLen;
				
				if (perID < perIDmatch || qCov < qCovThr) {
					//reject this hit
					return true;
				}
				segs = query;
			} else if (vsearch) {
				bool isSeed = line.substr(0, 1) == "S";
				if (line.substr(0, 1) != "H" && !isSeed) {
					continue;
				}
				for (uint i = 0; i < 4; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				perID = (float)atof(segs.c_str());
				for (uint i = 0; i < 5; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs2, '\t');
				if (isSeed) { segs2 = segs; perID = 100.f; }

			} else if (!UPARSE8up) { //uparse 7
				if (line.substr(0, 1) != "H" ) {
					continue;
				}
				for (uint i = 0; i < 4; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				perID = (float)atof(segs.c_str());
				for (uint i = 0; i < 5; i++) {//jump to pos X
					getline(ss, segs, '\t');
				}
				getline(ss, segs2, '\t');
			}
			else if (!UPARSE9up) { // uparse 8
			 //query first entry
				string tmp;
				getline(ss, segs, '\t');//0
				getline(ss, tmp, '\t');//1
				//should be "match"
				if (tmp == "chimera") {
					if (!doChimeraCnt) { continue; }
					chimera = true;
				}// segs = ""; continue;}
				if (tmp == "otu" || tmp == "OTU") {
					segs2 = segs;
					perID = 100.f;
				}
				else {
					getline(ss, tmp, '\t');//2
					perID = (float)atof(tmp.c_str());
					//indicator if hit
					//OTU last entry
					getline(ss, segs2, '\t');//3
					getline(ss, segs2, '\t');//4
				}

			}
			else { //UP9, uparse 10, uparse 11, dada2 fake .uc
				  //query first entry
				string tmp;
				getline(ss, segs, '\t');//0
				getline(ss, tmp, '\t');//1
									   //should be "match"
				if (tmp == "chimera") {
					continue;
				}
				else if (tmp == "noisy_chimera" || tmp == "good_chimera") { //tmp == "perfect_chimera" ||
					if (!doChimeraCnt) { continue; }
					chimera = true;
//M04428:252:000000000-BM3M5:1:1109:4503:15901;size=52;   noisy_chimera   
// dqt=20;dqm=2;div=9.0;segs=2;parents=M04428:252:000000000-BM3M5:1:2110:24763:14436(1-78/0)+M04428:252:000000000-BM3M5:1:1116:5359:16203(78-186/0);

					getline(ss, tmp, '\t');
					size_t p1(tmp.find("parents=") + 4);//4
					size_t p2(tmp.find(");", p1) + 2);
					segs2 = tmp.substr(p1, p2 - p1 - 1);
				}
				else if (tmp == "perfect_chimera") {
					removeSizeStr(segs);
					perfectChims.insert(segs);
					continue;
				}// segs = ""; continue;}
				if (tmp.substr(0, 3) == otuTerm) {
					segs2 = segs;
					perID = 100.f;
				}
				else {//match or perfect_chimera case
					getline(ss, tmp, '\t');//2
					//dqt=1;top=GZV0ATA01ANJXZ;size=14;(99.6%);
					size_t p1(tmp.find("top=") + 4);//4
					size_t p2(tmp.find(";(", p1) + 2);
					size_t p3(tmp.find("%);", p2));
					segs2 = tmp.substr(p1, p2 - p1 - 1);
					//string xx = tmp.substr(p2, p3 - p2);
					perID = (float)atof(tmp.substr(p2, p3 - p2).c_str());
					if (false && chimera) {//just use up9 top hit
	//					p1(tmp.find(";top=") + 5,p2);
	//					p2(tmp.find(";(", p1) + 2);

					}
				}
			}
		}



		//convert mapping to otu information
		//remove spaces
		segs = segs.substr(0,segs.find_first_of(' '));
		//also remove sample identifier in string
		string smplID = "";
		removeSampleID(segs, SEP, smplID);

		//finished parsing, now reformat to get matching sample

		if (chimera){
			if (UPARSE9up) {
				tarsV = splitByComma(segs2, false, '+');
				for (uint kk = 0; kk < tarsV.size(); kk++) {
					tarsV[kk] = tarsV[kk].substr(0, tarsV[kk].find_last_of("("));
				}
			}
			else if (UPARSE8up) {
				tarsV = splitByComma(segs2, false, '+');
				for (uint kk = 0; kk < tarsV.size(); kk++) {
					tarsV[kk] = tarsV[kk].substr(0, tarsV[kk].find_last_of("("));
				}
			}
		} else {
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
			removeSeqStr(segs2);
			removeCentrStr(segs2);
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
					//bestDNA2.push_back(NULL);
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
			//bestDNA[curCLID[kk]]->totalSum();
#ifdef matrix_sum
			if ( addFromHDstring && smplID != "") {
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
	cerr << "Writing OTU matrix with "<< SmplIDs.size() << " samples and " << seq2CI.size() << " OTUs to " << outf << endl;
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
		cerr << "New sample_id_ id_ in uc file detected, that is not present in map: \"" << smplID << "\"\n";// endl;
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
		tmpDNA->setAllQual(30);
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
		if (vsearch) {
			removeSeqStr(segs2);
			removeCentrStr(segs2);
		}

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			bestPID.push_back(100.f);//it's the original seq, so it's at 100%
			tmp->setTempFloat(100.f);
			bestDNA.push_back(tmp);
			//bestDNA2.push_back(NULL);
			oriKey.push_back(oriClKey);
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
		if (vsearch) {
			removeCentrStr(segs2);
			removeSeqStr(segs2);
		}

		//cluster should not exist, test
		std::unordered_map<string, int>::iterator itCL;
		itCL = seq2CI.find(segs2);

		if (itCL == seq2CI.end()) {
			//not found in known clusters.. create entry
			bestDNA.push_back(tmp);
			//bestDNA2.push_back(NULL);
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
			if ( fil->doReversePrimers() && !fil->checkYellowAndGreen(tdn1, CurSetPair, TagIdx, true) ) {
				return;//delete dnaTemp1;
			}
			bestDNA[curCLID] = tdn1;
			bestPID[curCLID] = tdn1->getTempFloat();
			bestLEN[curCLID] = tdn1->length();
		} else {
			if (fil->doReversePrimers() && !fil->checkYellowAndGreen(tdn1, CurSetPair, TagIdx, true)) {
				 return;//delete dnaTemp2; delete dnaTemp1;
			}
			bestDNA[curCLID] = tdn1;
			bestPID[curCLID] = tdn1->getTempFloat();
			bestLEN[curCLID] = tdn1->length() + tdn2->length();
			fil->checkYellowAndGreen(tdn2, 1, TagIdx, true);
			//bestDNA2[curCLID] = tdn2;
			
		}
	}//already a candidate sequence? check who is better..
	else if (
		fil->betterSeed(tdn1,  bestDNA[curCLID], 
			bestPID[curCLID], CurSetPair, checkBC) //bestLEN[curCLID], 
		){
//		delete bestDNA[curCLID];
		bestDNA[curCLID] = tdn1; 
		if (bestPID[curCLID] < tdn1->getTempFloat()){
			bestPID[curCLID] = tdn1->getTempFloat();
		}
		//uint curL = tdn1->length();
		//if (tdn2 != NULL) { curL += tdn2->length(); }
		//if (bestLEN[curCLID] < curL){			bestLEN[curCLID] = curL;		}
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
	if (idx != std::string::npos) {
		w = w.substr(0, idx);
	}
}
void UClinks::removeSeqStr(string& w) {
	size_t idx = w.find(";seqs=");
	if (idx != std::string::npos) {
		w = w.substr(0, idx);
	}
}
void UClinks::removeCentrStr(string& w) {
	size_t idx = w.find("centroid=");
	if (idx != std::string::npos) {
		w = w.substr(idx + 9);
	}
}

void UClinks::writeNewSeeds(OutputStreamer* MD, Filters* fil, 
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
	shared_ptr<DNAunique> d;


	
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
		shared_ptr<DNAunique> d2 = d->getPair();
		if (printLnk){
			string oriH = d->getShortId(); removeSizeStr(oriH);
			links << newH << "\t" << oriH << endl;
		}
		//takeOver id_ to allow for later linkup with cluster
		if (paired == 2) {
			//first check for readmerging
			d->setPassed(true);
			d->setNewID(newH + ".1");

			if (b_merge_pairs_optiSeed_ ){//&& bestDNA2[i] != NULL) {
				bool didMerge = MD->saveForWrite_merge(d, newH, 0, true);
				//didMerge = MD->saveForWrite_merge(dM, nullptr, newH, 0, true);
				if (!didMerge) {//not merged? we want to add this DNA nonetheless to output
					//better to do this in saveForWrite already
				}
			}

			int Cstream(100);
			if (d2 != nullptr) {
			    //sorted by pair1,2, quality yellow green, singleton (pair 1 2 3 4) etc
				d2->setPassed(true);
				d2->setNewID(newH + ".2");
				MD->saveForWrite(d, 1,-1, Cstream,false);
				MD->saveForWrite(d2, 2,-1, Cstream,false);
				MD->writeForWrite(d, 1, Cstream, d2, 2, Cstream);
				

			} else {
				MD->saveForWrite(d, 3,-1, Cstream,true);
			}
		} else {
			int Cstream(0);
			d->setNewID(newH);
			d->setPassed(true);
			MD->saveForWrite(d, 1,-1, Cstream,true);
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
		shared_ptr<DNAunique> d = bestDNA[i];
		if (d == nullptr) { continue; }
		shared_ptr<DNAunique> d2(nullptr);
		if (d->getMerge() != nullptr) {
			d = d->getMerge();
		} else {
			d2 = d->getPair();
		}
		float curQ = d->getAvgQual();
		if (d2 != NULL) {
			curQ += d2->getAvgQual(); curQ /= 2.f;
		}
		uint curL = (uint) d->length();
		if (d2!= NULL) {
			curL += d2->length(); 
		}

		if (curL < minL){ minL = curL; }
		if (curL > maxL){ maxL = curL; }
		avgL += curL;
		lengths.push_back(curL);
		
		if (curQ < 1) {//no new Seed found, default seed
			continue;
		}
		float sc = d->getTempFloat();

		if (curQ < minQ){ minQ = curQ; }
		if (curQ > maxQ){ maxQ = curQ; }
		avgQ += curQ;

		float curA = (float)d->getAccumError();
		if (d2 != NULL) {
			curA += (float) d2->getAccumError(); 
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
	if (lengths.size() > 0) { os << "\n     - Sequence Length :   " << calc_median2(lengths, 0.f) << "/" << calc_median2(lengths, 0.1f) << "/" << calc_median2(lengths, 0.5f) << "/" << calc_median2(lengths, 0.9f) << "/" << calc_median2(lengths, 1.f);; }
	if (quals.size() > 0) { os << "\n     - Quality :      " << calc_median2(quals, 0.f) << "/" << calc_median2(quals, 0.1f) << "/" << calc_median2(quals, 0.5f) << "/" << calc_median2(quals, 0.9f) << "/" << calc_median2(quals, 1.f);; }
	if (accums.size() > 0) { os << "\n     - Accum. Error : " << calc_median2(accums, 0.f) << "/" << calc_median2(accums, 0.1f) << "/" << calc_median2(accums, 0.5f) << "/" << calc_median2(accums, 0.9f) << "/" << calc_median2(accums, 1.f);;	}
	if (sims.size() > 0) {	os << "\n     - Sim2Consensus: " << calc_median2(sims, 0.f) << "/" << calc_median2(sims, 0.1f) << "/" << calc_median2(sims, 0.5f) << "/" << calc_median2(sims, 0.9f) << "/" << calc_median2(sims, 1.f);		}
	os << endl;
}
