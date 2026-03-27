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
#include "Common.h"
#include "containers.h"
#include "Filters.h"
#include "OutputStreamer.h"

using namespace std;





// trim and is_digits are provided by Common.cpp


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
			outputFile << "\t" + std::to_string((int)i) + ":" + barcode_number_to_sample_id_[i];
		}
	}
	else {
		smplId2comb = mainFilter->combiSmplConvergeVec(barcode_number_to_sample_id_);
		if (barcode_number_to_sample_id_.size() != smplId2comb.size()) {
			cerr << "FATAL: barcode_number_to_sample_id_ != smplId2comb\n"; exit(234);
		}
        for (auto IT = combiMapCollectGrp.begin(); IT != combiMapCollectGrp.end(); IT++) {
			outputFile << "\t" << std::to_string(IT->second) + ":" + IT->first;
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
			omaps << "\t" << std::to_string(IT->second) + ":" + IT->first;
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
	report += "Dereplication: " + N_passed + " unique sequences (avg size " + intwithcommas((int)avgSize) + "; " + intwithcommas((int)passedSize) + " counts)\n";
	
		if (passed_hits > 0) {
       report += N_notPassed + "/" + N_total + " not passing derep conditions (" + intwithcommas((int)notPassedSize) + " counts; "+ minCopiesStr;
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
		string newlySetID = otuOUTterm + std::to_string(cnt);
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
