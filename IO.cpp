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

#include <queue>
#include <mutex>
#include <thread>
#include "IO.h"
#include "Benchmark.h"
#include <future>
//#include "IOMultithreaded.h"


void read_single(OptContainer& cmdArgs, shared_ptr<OutputStreamer> MD, shared_ptr<InputStreamer> IS){
	Filters* curFil = MD->getFilters();
    curFil->singReadBC2();
    int chkDerep(0);
    bool checkReversedRead = curFil->checkRevRd();
    bool cont(true); //bool sync(false);
    
    while (cont){
        vector<shared_ptr<DNA>> tdn = IS->getDNAMC();
        
        if (tdn[0] == nullptr) {
#ifdef DEBUG
            cerr << "NULL read returned" << endl;
#endif
			cont = false;
            break;
        }


        //collect some info on general run parameters
        curFil->preFilterSeqStat(tdn[0], 0); // statistics preFilter
        curFil->sTotalPlus(0);//mutex update total counts

    
        //thread starts here
        int tagIdx(-2);
        if (true && checkReversedRead ) {
            string presentBC(""); int c_err(0);
            int chkRev(1);
            tagIdx = curFil->findTag(tdn[0], presentBC, c_err, true, chkRev);
/*			if (tagIdx < 0) { //check if on reversed_ read
				dnaTemp1->reverse_transcribe();
				tagIdx = curFil->findTag(dnaTemp1, presentBC, c_err, true);
			}*/
            
            if (chkRev==0) {//no? undo revTranscr
				tdn[0]->reverse_transcribe();
            }
        }
        tagIdx = -2;
		int curThread = -1;
        MD->analyzeDNA(tdn[0], -1, -1, tagIdx,curThread);
    
    
        //thread ends here
        //here BC has to be correctly set within DNA object
        MD->dereplicateDNA(tdn[0], nullptr);//run in extra thread?
		MD->write2Demulti(tdn[0], 0, curFil->getBCoffset());
    
        
        //first save read in mem, then write if enough reads accumulate in mem
        //extra thread for multithreading??

        // This is a quality collecting step
        if (!MD->saveForWrite(tdn[0],1, curThread)) {
            cont = false;
            break;
		}
		if (tdn[0]->isGreenQual()) {
			chkDerep++;
		}

		//if (tdn!=NULL && ch1 != tdn->isGreenQual()){cerr<<"isGreenQual is != ch1! Aborting..\n";exit(12);}
	}
	MD->closeOutStreams();
}


bool read_paired_STRready(vector< vector< string>> tmpLines,
	bool MIDuse, shared_ptr<OutputStreamer> MD, int curThread, 
	bool keepPairHd, qual_score FastqVer ) {
	vector<shared_ptr<DNA>> ret(3,nullptr);
	ret[0] = str2DNA(tmpLines[0], keepPairHd, FastqVer,0);
	ret[1] = str2DNA(tmpLines[1], keepPairHd, FastqVer,1);
	if (MIDuse) {
		ret[2] = str2DNA(tmpLines[2], keepPairHd, FastqVer,2);
	}

	return read_paired_DNAready(ret, MIDuse, MD, curThread);
}
static mutex testreadpair;
//is called from a while loop, that reads the DNA pairs
bool read_paired_DNAready(vector< shared_ptr<DNA>> tdn,
	bool MIDuse, shared_ptr<OutputStreamer> MD, int curThread) {

	if (tdn[0] == nullptr) {
	    return false;
	} //|| tdn->length()==0
	//DNA objects as they should be??
	if (MIDuse && tdn[2] == nullptr) {
		cerr << "Missing MID read pair.\n";
		exit(4);
	}
	if (tdn[1] == nullptr && tdn[0] != nullptr) {
		cerr << "Second provided file has not the same number of entries as first file.\n";
		exit(5);
	}
	//MD->checkFastqHeadVersion(tdn[0]);
	//testreadpair.lock();

	Filters* curFil = MD->getFilters(curThread);
	//register read at all with stat counter:
	curFil->sTotalPlus(0); curFil->sTotalPlus(1);
	//collect some info on general run parameters
	curFil->preFilterSeqStat(tdn[0], 0);
	curFil->preFilterSeqStat(tdn[1], 1);

	//prep some variables
	int BCoffs = curFil->getBCoffset();
	bool checkBC2ndRd = curFil->checkBC2ndRd();
	bool dualBCs = curFil->doubleBarcodes();
	bool doBCsAtAll = curFil->doBarcodes();
	bool checkReversedRead = curFil->checkRevRd();

	int tagIdx(-2); int tagIdx2(-2);
	string presentBC(""); int c_err(0);
	bool isReversed(false);//was a reversion detected?

	if (MIDuse && tdn[2] != nullptr) {
		tagIdx = curFil->cutTag(tdn[2], presentBC, c_err, true); 
//		delete tdn[2]; 
		tdn[0]->setBCnumber(tagIdx, BCoffs);
	}

    if (checkBC2ndRd ) {
		if (!dualBCs) {
			bool revT = false;
			bool Pr1 = curFil->findPrimer(tdn[0], 0, false, 0);
			bool Pr2 = curFil->findPrimer(tdn[1], 0, false, 0);
			int chkRev1(-1), chkRev2(-1);
			tagIdx = curFil->findTag(tdn[0], presentBC, c_err, true, chkRev1);
			tagIdx2 = curFil->findTag(tdn[1], presentBC, c_err, true, chkRev2);
			if ( true &&checkReversedRead && (tagIdx2 < 0 && tagIdx < 0) ) {
				tdn[0]->reverse_transcribe(); tdn[1]->reverse_transcribe();
				Pr1 = curFil->findPrimer(tdn[0], 0, false, 0);
				Pr2 = curFil->findPrimer(tdn[1], 0, false, 0);
				tagIdx = curFil->findTag(tdn[0], presentBC, c_err, true, chkRev1);
				tagIdx2 = curFil->findTag(tdn[1], presentBC, c_err, true, chkRev2);
				revT = true;
			}
			if ((tagIdx2 >= 0 && tagIdx < 0 && !Pr1) || (Pr2 && !Pr1)) { //swap first & second read
				swap(tdn[0], tdn[1]);
				tdn[0]->constellationPairRev(true); tdn[1]->constellationPairRev(true);
				//revConstellation++;
			}
			/*else if (tagIdx2 < 0 && tagIdx < 0) {
				int x = 0;
			}*/
			if (revT) {
				tdn[0]->reverse_transcribe(); tdn[1]->reverse_transcribe();
			}
		}
		tagIdx2 = -2; tagIdx = -2;
		tdn[1]->setpairREV();		tdn[0]->setpairFWD();
	}

    
	//tdn[0]->reverse_transcribe();
	MD->analyzeDNA(tdn[0], -1, 0, tagIdx, curThread);
	//tdn[0]->matchSeqRev
	bool ch1(false); if (tdn[0] != NULL) { ch1 = tdn[0]->isGreenQual(); }
	bool ch2(false); bool ch2n(false);

	//this is all about barcodes..
	if (checkReversedRead  && tdn[0] != NULL && tagIdx < 0) {
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

	//test for reverse complemented reads (mohammad samples), when BC not found (NOT dual BC)
	//in that case, this is the first read
	if (false &&checkBC2ndRd && tagIdx < 0 && tdn[1] != NULL) {// && !tdn[0]->getBarcodeDetected() ) {
											   //dnaTemp2->reverse_transcribe();
		if (!MIDuse) { tagIdx = -2; }
//		curFil->sTotalMinus(0);
		MD->analyzeDNA(tdn[1], -1, 0, tagIdx, curThread);
		ch2n = tdn[1]->isGreenQual();
		if (!ch2n && checkReversedRead) {
			if (!MIDuse) { tagIdx = -2; }
//			curFil->sTotalMinus(0);
			tdn[1]->reverse_transcribe();
			MD->analyzeDNA(tdn[1], -1, 0, tagIdx, curThread);
			ch2n = tdn[1]->isGreenQual();
			isReversed = ch2n;
			if (!ch2n) { tdn[1]->reverse_transcribe(); }//reset to ori
		}
		if (ch2n) {//passed ch2 through BC filter, now really reverse
				   //1st, now 2nd pair_
			tdn[1]->setpairFWD();
			ch1 = ch2n;
			if (tdn[0] != NULL) {
				//tdn[0]->reverse_transcribe(); 
				tdn[0]->setpairREV();
				tdn[0]->reset();
				if (!dualBCs) { tagIdx2 = tdn[0]->getBarcodeNumber(); } // no 2nd BC, thus no BC search in 2nd read
				MD->analyzeDNA(tdn[0], -1, 1, tagIdx2, curThread);
				ch2 = tdn[0]->isGreenQual();
			}
			swap(tdn[0], tdn[1]);
			tdn[0]->reverse_transcribe(); tdn[1]->reverse_transcribe();
		}

	}

	//if ( ch1 ) {	cerr << cnt << " \n";	}
	//normal case for check 2nd read
	if (!ch2 && tdn[1] != NULL) { //ch1&&
								//dnaTemp2->setBCnumber(tdn[0]->getBarcodeNumber());
		if (doBCsAtAll && !dualBCs) { //only check in read1 for BC, if not dual BCing!!
			tagIdx2 = tdn[0]->getBarcodeNumber();  // no 2nd BC, thus no BC search in 2nd read
			if (tagIdx2 >= 0) {
				tagIdx2 -= BCoffs;
			}else if (tagIdx2 < -1) {//something wrong with BCoffs
				cerr << "tagidx2 wrongly truncated to " << tagIdx2 << endl;
			}
		}
		if (isReversed) { tdn[1]->reverse_transcribe(); }
		MD->analyzeDNA(tdn[1], -1, 1, tagIdx2, curThread);
		ch2 = tdn[1]->isGreenQual();
	}

	//set up BC in DNA header
	//remember that dual BCs are only valid after this step!
	if (dualBCs) {
		//tagIdx2 = -2; //reset just to be sure
		curFil->dblBCeval(tagIdx, tagIdx2, presentBC, tdn[0], tdn[1]);
		c_err = -1;

		//check a second time that barcode was correctly identified, just to be double sure...
		if (tagIdx != tagIdx2 || tdn[0]->getBarcodeNumber() != tdn[1]->getBarcodeNumber()) {
			cerr << "Unequal BC numbers:" << tagIdx << " : " << tagIdx2 << "; in object: " << tdn[0]->getBarcodeNumber() << " : " << tdn[1]->getBarcodeNumber() << endl;
			cerr << "In read:" << tdn[0]->getId() << endl;
			exit(835);
		}
	}
	else if (tagIdx >= 0) {
		if (MIDuse&&ch1) { curFil->BCintoHead(tagIdx, tdn[0], presentBC, c_err, true); }
		else { curFil->setBCdna(tagIdx, tdn[0]); }
		if (ch2) { curFil->BCintoHead(tagIdx, tdn[1], presentBC, c_err, true); }
	}

	
	if (tagIdx == -1 || tagIdx2 == -1) {
		tdn[0]->setBarcodeDetected(false);
		tdn[1]->setBarcodeDetected(false);
	}



	int idx1 = 1; int idx2 = 2;
	if (ch1 && !ch2) {
		idx1 = 3; idx2 = 4;
		if (tdn[1] != NULL) { tdn[1]->failed(); }
		//		delete dnaTemp2;
	}
	else if (ch2 && !ch1) {
		idx2 = 4; idx1 = 3;
		if (tdn[0] != NULL) { tdn[0]->failed(); }
		//		delete tdn[0];
	}
	else if (!ch1 && !ch2) { //nothing passes
		if (tdn[0] != NULL) { tdn[0]->failed(); }
		if (tdn[1] != NULL) { tdn[1]->failed(); }
		//		delete tdn[0]; delete dnaTemp2;
	}

	//pre-merge step
	if ( MD->mergeReads()) {MD->findSeedForMerge(tdn[0], tdn[1],curThread);	}

	//demultiplex write? do this first before DNA is deleted..
	//at this point the tagIDX *MUST* be correctly set + BCoffset (in the DNA object, tagIDX doesn;t matter)
	MD->write2Demulti(tdn[0], tdn[1], curFil->getBCoffset(), curThread);
	MD->dereplicateDNA(tdn[0], tdn[1]);
	MD->writeNonBCReads(tdn[0], tdn[1]);

    // Test
    /*if (OutStreamer->b_merge_pairs_derep_ && tdn[0]->merge_seed_pos_ > 0) {
        auto dna_merged = ReadMerger::merge(tdn[0], tdn[1]);
    }*/
	//testreadpair.unlock();
	//save for later .. and collect stats
	if (!MD->saveForWrite(tdn[0], idx1, curThread) || !MD->saveForWrite(tdn[1], idx2, curThread)) {
		return false;
	}
	return true;
}

struct job2 {
	bool inUse = false;
	future<bool> job;
};

bool read_paired(OptContainer& cmdArgs, shared_ptr<OutputStreamer> MD, 
	shared_ptr<InputStreamer> IS, bool MIDuse, int Nthreads) {
	
	DNAmap oldMIDs;
	bool fqHeadVer(true);
	cdbg( "Read paired routine\n");

	/*if (sync2pair && MIDuse) {
		cout << "Can not sync read pairs, while explicit MID sequences are being used! (not supported, sorry)\n";
		sync2pair = false;
	} */

	//bool syncedMID = false;
	
	//multithreading setup
	//int Nthrds = atoi(cmdArgs["-threads"].c_str()) -1 ;
	vector<job2> slots(Nthreads);
	int thrCnt = 0;
	bool doMC(false);
	if (Nthreads > 1) {
		doMC = true;
	}

	int DNAinMem(0);
	bool cont(true),cont2(true),cont3(true);
	bool keepPairedHD = IS->keepPairedHD();
	int revConstellation(0);
	int cnt(0); 
	bool switching(true); // important to keep track of this, to fix swapped read pairs

	vector<string>tmpLines2(4, "");
	vector<vector<string>> tmpLines(3, tmpLines2);
	vector<shared_ptr<DNA>> tdn(3,nullptr);

	while ( cont ) {
		//bool sync = false;
		//tests of different ways to read files..
		if (false) {
			tdn = IS->getDNAMC();//
		}
		else if (true) {
			for (uint i = 0; i < 3; i++) {
				if (!MIDuse && i == 2) {
					continue;
				}
				IS->getDNAlines(tmpLines[i], i);
			}
			
		
		} else {
			//tdn.resize(3, nullptr);
			tdn[0] = IS->getDNA(0);
			tdn[1] = IS->getDNA( 1);
			if (MIDuse) {
				tdn[2] = IS->getDNA(2);
				if (tdn[2] != nullptr) {
					tdn[2]->setMIDseq(true);
				}
			}
		}
		qual_score fastqVer = IS->fastQscore();

		if (fqHeadVer) { //some things just need to be done
			if (tdn[0] == nullptr) {
				tdn[0] = str2DNA(tmpLines[0], keepPairedHD, fastqVer, 0);
			}
			MD->checkFastqHeadVersion(tdn[0]); 
			fqHeadVer = false; 
			tdn[0] = nullptr;//delete again, this is a oneoff..
		}

		cnt++;
		//decide on MC or single core submission route
		if (true && doMC) {
//			if (tdn[0] == nullptr) { cont = false;  break; }
//			if (fqHeadVer) { MD->checkFastqHeadVersion(tdn[0]); fqHeadVer = false; }

			//work with threadpool instead 
			/*
			if (thrCnt >= Nthreads) { thrCnt = 0; }
			pool->enqueue([tdn, MIDuse, MD, thrCnt]
				{ read_paired_DNAready(tdn, MIDuse, MD, 0); }
			);
			thrCnt++;
			/**/
			
			bool notSubm(true);
			while (notSubm) {//go over possible submission slots
				if (thrCnt >= Nthreads) {thrCnt = 0;}
				if (slots[thrCnt].inUse == true && 
					slots[thrCnt].job.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
					slots[thrCnt].inUse = false;
					cont = slots[thrCnt].job.get();
				}
				if (slots[thrCnt].inUse == false) {
					//cdbg("submit readPairRdy ");
//					slots[thrCnt].job = async(std::launch::async, read_paired_DNAready,
//						tdn, MIDuse, MD, thrCnt);
					slots[thrCnt].job = async(std::launch::async, read_paired_STRready,
						tmpLines, MIDuse, MD, thrCnt, keepPairedHD, fastqVer);
					slots[thrCnt].inUse = true;
					notSubm = false;
				}
				thrCnt++;
			}
			if (!cont) { break; }
			/**/

		} else {
			if (0) {
				if (tdn[0] == nullptr) { cont = false;  break; }
				//if (fqHeadVer) { MD->checkFastqHeadVersion(tdn[0]); fqHeadVer = false; }
				cont = read_paired_DNAready(tdn, MIDuse, MD, 0);
				if (tdn[0]->isConstellationPairRev()) { revConstellation++; }
			}
			else {
				cont = read_paired_STRready(tmpLines, MIDuse, MD, 0, 
					keepPairedHD, fastqVer);
				//if (tdn[0]->isConstellationPairRev()) { revConstellation++; }
			}
		}
	}

	//get all slots
	for (int x = 0; x < slots.size(); x++){
		if (slots[x].inUse == true && slots[x].job.wait_for(std::chrono::milliseconds(1)) == std::future_status::ready) {
			slots[x].inUse = false;
			cont = slots[x].job.get();
		}
	}

	
	//close shop
	MD->revConstellationCnts(revConstellation);
	MD->closeOutStreams();
	return true;
}


bool readCmdArgs(int argc, char* argv[],OptContainer& cmdArgs){
	if (argc%2!=1){
		cerr<<"It seems command line arguments were not passed in pairs. Aborting.\n";
		exit(666);
	}
	for (int i=1; i<argc; i+=2){ //parsing of cmdline args
		string theNxtSt = string(argv[i+1]);
		if (theNxtSt[0] != '-'){
			cmdArgs[string(argv[i])] = theNxtSt;
		} else {
			cmdArgs[string(argv[i])] = "T";
		}
	}
	if (cmdArgs.find("-i_MID_fastq") == cmdArgs.end()) {
		cmdArgs["-i_MID_fastq"] = "";
	}
	//set to default (empty)
	if (cmdArgs.find("-OTU_fallback") == cmdArgs.end()){
		cmdArgs["-OTU_fallback"] = "";
	}
	if (cmdArgs.find("-otu_matrix") == cmdArgs.end()) {
		cmdArgs["-otu_matrix"] = "";
	}
	if (cmdArgs.find("-derepPerSR") == cmdArgs.end()) {
		cmdArgs["-derepPerSR"] = "0";
	}
	if (cmdArgs.find("-ucAdditionalCounts") == cmdArgs.end()) {//.ADD
		cmdArgs["-ucAdditionalCounts"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts1") == cmdArgs.end()) {//.REST
		cmdArgs["-ucAdditionalCounts1"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts_refclust") == cmdArgs.end()) {//.ADDREF
		cmdArgs["-ucAdditionalCounts_refclust"] = "";
	}
	if (cmdArgs.find("-optimalRead2Cluster_ref") == cmdArgs.end()) {//.ADDREF
		cmdArgs["-optimalRead2Cluster_ref"] = "";
	}
	//just for debuggin purposes: write out all seqs, where no BC can be detected..
	if (cmdArgs.find("-o_fastq_noBC") == cmdArgs.end()) {
		cmdArgs["-o_fastq_noBC"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts_refclust1") == cmdArgs.end()) {//.RESTREF
		cmdArgs["-ucAdditionalCounts_refclust1"] = "";
	}
	if (cmdArgs.find("-ucAdditionalCounts_refclust1") == cmdArgs.end()) {//.RESTREF
		cmdArgs["-ucAdditionalCounts_refclust1"] = "";
	}
	if (cmdArgs.find("-XfirstReads") == cmdArgs.end()) {
		cmdArgs["-XfirstReads"] = "";
	}
			//filter sequence file for a specific subset of sequences
	//these arguments can only occur together
	if (cmdArgs.find("-specificReads") == cmdArgs.end()) {
		cmdArgs["-specificReads"] = "";
	} else if (cmdArgs.find("-excludeFile") == cmdArgs.end()) {
		cmdArgs["-excludeFile"] = "";
	}
	if (cmdArgs.find("-onlyPair") == cmdArgs.end()) {
		cmdArgs["-onlyPair"] = "";
	}
	if (cmdArgs.find("-pairedDemulti") == cmdArgs.end()) {
		cmdArgs["-pairedDemulti"] = "0"; // by default: do not only report proper pairs
	}
	
	if (cmdArgs.find("-uparseVer") == cmdArgs.end()) {
		cmdArgs["-uparseVer"] = "";
	}


	if (cmdArgs.find("-i_path")  == cmdArgs.end()){ //ok files are not given in mapping file
		//check if dna and qual_ are passed
		if (cmdArgs.find("-i") != cmdArgs.end()) {
			string fmt = detectSeqFmt(cmdArgs["-i"]);
			if (fmt == "empty") {
				//exit(0);
				cerr << "Only empty input files\n";
			}
			cmdArgs[fmt] = cmdArgs["-i"];
		}
		if (cmdArgs.find("-i_fastq")  == cmdArgs.end()){ // fasta + quality format
			if (cmdArgs.find("-i_fna")  == cmdArgs.end()){
				cerr<<"You did not supply a fasta file. \nPlease give the path to your fasta file as command line argument:\n  -i_fna <yourFastaFile>\n";
				exit(2);
			}
			if (cmdArgs.find("-i_qual")  == cmdArgs.end()){
				string newQ = cmdArgs["-i_fna"];
				int pos = (int)newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual_");
				fstream fin;
				fin.open(newQ.c_str(),ios::in);
				if( fin.is_open() )	{
					cerr<<"Using quality file: "<<newQ <<endl;
				} else if ((cmdArgs.find("-number")!=  cmdArgs.end() && cmdArgs["-number"] =="T")||
					(cmdArgs["-specificReads"] != "")) {
					cmdArgs["-i_qual"] = "";
				} else {
					cerr<<"You did not supply a quality file. \nPlease give the path to your quality file as command line argument:\n  -i_qual <PathToQualityFile>\n";
					newQ = "";
					//fin.close();	exit(2);
				}
				fin.close();
				cmdArgs["-i_qual"] = newQ;
			}
		}
		//auto create output file name
		if (cmdArgs.find("-o_fna")  == cmdArgs.end()){
			if (cmdArgs.find("-o_fastq")  == cmdArgs.end()){
				//cmdArgs["-o_fna"] = cmdArgs["-i_fna"]+string(".sdm");
				//cerr<<"Writing output fasta into "<<cmdArgs["-o_fna"]<<endl;
				cerr << "No output file will be written\n";
			}
		} else {
			if (cmdArgs.find("-o_fastq")  != cmdArgs.end()){
				cerr<<"\"-o_fna\" was over-writen by \"-o_fastq\"\n";
				cmdArgs["-o_fna"] = "";
			}
		}
	} else {
		if (cmdArgs.find("-o_fna")  == cmdArgs.end() && cmdArgs.find("-o_fastq")  == cmdArgs.end()){
			cerr<<"Please give an output file (\"-o_fna\" || \"-o_fastq\") if you use sdm \"-i_path\" option.\n  Aborting..\n";
			exit(2);
		}
	}

	/*	if (cmdArgs.find("-base_map")  == cmdArgs.end()){
	cerr<<"You did not supply a mapping file. \nPlease give the path to your mapping file as command line argument:\n  -base_map <PathToMappingFile>\n";
	exit(2);
	}  */
	if (cmdArgs.find("-o_qual")  == cmdArgs.end()){
		cmdArgs["-o_qual"] = "";
	} else {
		if (cmdArgs.find("-o_fastq")  != cmdArgs.end()){
			cerr<<"\"-o_qual\" was over-writen by \"-o_fastq\"\n";
			cmdArgs["-o_qual"] = "";
		}
	}
	if (cmdArgs.find("-options")  == cmdArgs.end()){
		cmdArgs["-options"] = string("sdm_options.txt");
	}
	if (cmdArgs.find("-threads")  == cmdArgs.end()){
		cmdArgs["-threads"] = "1";
	}
	if (cmdArgs.find("-log")  == cmdArgs.end()){
		string ofile1 = cmdArgs["-o_fna"];
		if (ofile1==""){ofile1 = cmdArgs["-o_fastq"];}
		vector<string> tvec = splitByComma(ofile1,false); 
		ofile1 = tvec[0];
		//remove file ending
		size_t pos = ofile1.find_last_of(".");
		if (pos != string::npos){ofile1 = ofile1.substr(0,pos);	}
		if (tvec.size()==2){
			ofile1+= "_" + getFileNoPath(tvec[1]);
			pos =  ofile1.find_last_of(".");
			if (pos != string::npos){ofile1 = ofile1.substr(0,pos);	}
		}
		cmdArgs["-log"] = ofile1 + string(".log");
	}
	string ofile1 = cmdArgs["-log"];
	//ofile1.find_last_of(".log");
	size_t logPos = ofile1.find_last_of(".");
	if (logPos != std::string::npos){
		ofile1 = ofile1.substr(0,logPos);
	}
	if (cmdArgs.find("-length_hist")  == cmdArgs.end()){
		cmdArgs["-length_hist"]  = ofile1 + string("_lenHist.txt");
	}
	if (cmdArgs.find("-qual_hist") == cmdArgs.end()) {
		cmdArgs["-qual_hist"] = ofile1 + string("_qualHist.txt");
	}
	if (cmdArgs.find("-merg_readpos") == cmdArgs.end()) {
		cmdArgs["-merg_readpos"] = ofile1 + string("_mergRpos.txt");
	}
	if (cmdArgs.find("-	qual_readpos") == cmdArgs.end()) {
		cmdArgs["-qual_readpos"] = ofile1 + string("_qualRpos.txt");
	}

		//-length_hist   -qual_hist

	if (cmdArgs.find("-sample_sep")  == cmdArgs.end()){
		cmdArgs["-sample_sep"] = DEFAULT_BarcodeNameSep;
	} else 	if (cmdArgs["-sample_sep"]==""){
		cerr<<"Invalid sample separator (empty).\nAborting..\n";exit(82);
	}


	if (cmdArgs.find("-o_qual_offset") == cmdArgs.end()) {
		cmdArgs["-o_qual_offset"] = DEFAULT_output_qual_offset;
	}
	if (cmdArgs.find("-pairedRD_HD_out") == cmdArgs.end()) {
		cmdArgs["-pairedRD_HD_out"] = DEFAULT_pairedRD_HD_out;
	}

	if (cmdArgs.find("-ignore_IO_errors") == cmdArgs.end()) {
		cmdArgs["-ignore_IO_errors"] = DEFAULT_ignore_IO_errors;
	} else if (cmdArgs["-ignore_IO_errors"] != "0" && cmdArgs["-ignore_IO_errors"] != "1") {
		cerr << "Argument \"ignore_IO_errors\" can only be \"1\" or \"0\". Instead it has value: " << cmdArgs["-ignore_IO_errors"] << endl;
		exit(323);
	}
	if (cmdArgs.find("-o_dereplicate") == cmdArgs.end()) {
		cmdArgs["-o_dereplicate"] = "";
	}
	if (cmdArgs.find("-derep_map") == cmdArgs.end()) {
		cmdArgs["-derep_map"] = "";
	}

	

	//if (cmdArgs.count("-i_fna")==0){}

	return true;
}


/*******************************************
*				read_fasta   			   *
*******************************************

void openOutFiles(string files, string fmt, string xtr){
	ofstream fnaOut;

	vector<string> tfnaout(0);
	if (files.find(",") != string::npos){
		tfnaout = splitByCommas(files);
	} else {
		tfnaout.push_back(files);
	}
	bool multiple = tfnaout.size() > 1;
	string xtr2 = "";
	if (multiple){xtr2 = "paired ";}
	for (uint i =0; i< tfnaout.size(); i++){
		fnaOut.open ( tfnaout[i].c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open "<<xtr2<<xtr<<fmt<<" output file "<<i<<": "<<tfnaout[i]<<endl;exit(4);	}
		fnaOut.close();
		if (multiple){//also singletonfiles
			string tmp = tfnaout[0]+SingletonFileDescr;
			fnaOut.open(tmp.c_str(),ios_base::out);
			if (!fnaOut){	cerr<<"Could not open Singleton "<<xtr<<fmt<<" output file "<<i<<": "<<tmp<<endl;exit(4);	}
			fnaOut.close();
		}
	}

}

void prepareOutFiles(OptContainer& cmdArgs){
	ofstream fnaOut;
	
	//additional output files (secondary filtering)
	if (cmdArgs.find("-o_fastq2")  != cmdArgs.end() && cmdArgs["-o_fastq2"] != ""){
		openOutFiles(cmdArgs["-o_fastq2"],"fastq","add ");
	}
	if (cmdArgs.find("-o_fna2")  != cmdArgs.end() && cmdArgs["-o_fna2"] != ""){
		openOutFiles(cmdArgs["-o_fna2"],"fna","add ");
	}

	//fastq
	if (cmdArgs["-o_fna"]=="" && cmdArgs["-o_fastq"] != ""){
		openOutFiles(cmdArgs["-o_fastq"],"fastq","");
		return;
	}

	//fasta output
	vector<string> tfnaout = splitByComma(cmdArgs["-o_fna"],false);
	for (unsigned int i=0; i<tfnaout.size();i++){
		fnaOut.open(tfnaout[i].c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open Fasta "<< i<<" output file "<<tfnaout[0]<<endl;exit(4);	}
		fnaOut.close();
	}
	if (tfnaout.size()==2){//PE - singleton file
		string tmp = tfnaout[0]+SingletonFileDescr;
		fnaOut.open(tmp.c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open Singleton Fasta output file "<<tmp<<endl;exit(4);	}
		fnaOut.close();
		tmp = tfnaout[1]+SingletonFileDescr;
		fnaOut.open(tmp.c_str(),ios_base::out);
		if (!fnaOut){	cerr<<"Could not open Singleton Fasta output file "<<tmp<<endl;exit(4);	}
		fnaOut.close();
	}
	if (cmdArgs["-o_qual"] != ""){
		vector<string> tqout = splitByComma(cmdArgs["-o_qual"],false);
		for (unsigned int i=0; i<tqout.size();i++){
			fnaOut.open (tqout[0].c_str() ,ios_base::out);
			if (!fnaOut){			cerr<<"Could not open Quality "<<i<<" output file "<<tqout[0]<<endl;			exit(4);		}
			fnaOut.close();
		}
		if (tqout.size()==2){//PE - singleton file
			string tmp = tqout[0]+SingletonFileDescr;
			fnaOut.open(tmp.c_str(),ios_base::out);
			if (!fnaOut){	cerr<<"Could not open Singleton Quality output file "<<tmp<<endl;exit(4);	}
			fnaOut.close();
			tmp = tqout[1]+SingletonFileDescr;
			fnaOut.open(tmp.c_str(),ios_base::out);
			if (!fnaOut){	cerr<<"Could not open Singleton Quality output file "<<tmp<<endl;exit(4);	}
			fnaOut.close();
		}
	}
}
*/

//manages read in of several input files and associated primers / tags to each file
void separateByFile(Filters* mainFilter, OptContainer& cmdArgs){
#ifdef DEBUG
	cerr << "separateByFile"<<endl;
#endif

//	mainFilter->ini_filestruct(cmdArgs);
	
	vector<string> FastaF = mainFilter->getFastaFiles();
	vector<string> QualF = mainFilter->getQualFiles();
	vector<string> FastqF = mainFilter->getFastqFiles();
	vector<string> MIDfq = mainFilter->getMIDfqFiles();
	vector<string> fastXtar;
	
	// Indicates if FASTQ files were submitted
	bool isFastq = true;
	//setup once at start
	vector<ReadMerger*> merger;

	
	//prepareOutFiles(cmdArgs);
	string path="";

	//set up some log structures
	string deLog("");//dereplication main log
	string logF = cmdArgs["-log"], logFA = cmdArgs["-log"].substr(0, cmdArgs["-log"].length() - 3) + "add.log";

	
	// Set folder path
	if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() > 2){
		path=cmdArgs["-i_path"] + string("/");
	}

    // Set up b_derep_as_fasta_ or fastq way and save file vector in tar in case it is zipped
	if (FastaF.size() > 0) { // If b_derep_as_fasta_ vector contains elements
		fastXtar = FastaF;
        isFastq = false; // Set boolean Fastq to false
	} else { // If no FastaF present assume there are Fastq files
		fastXtar = FastqF;
		if (FastqF.size()==0){ // no Fasta and no Fastq files -> abort
			cerr<<"No FastQ or Fasta file given.\n  Aborting..\n";
			exit(12);
		}
	}

	// Unique Fas initialized with first element of tar (can be b_derep_as_fasta_ and fastq)
	// Contains all unique b_derep_as_fasta_ or fastq files from the mapping file
	//vector<string> uniqueFastxFiles(1, fastXtar[0]);
	unordered_map<string, int> uniqueFastxFiles;

	// idx content: [ [0] ]
	// idx contains one row (vector) for each unique string in tar
	// This vector then contains the indices at which this string occurs in tar
	vector < vector<int> > idx(0);
	//idx.push_back(vector<int> (1,0));


	// We dont know if it is a tar yet, but we call it tar
	//this routine is important for managing the blocks of files to be read together
	for (unsigned int i=0; i<fastXtar.size(); i++){
		bool suc = false;
		auto XX = uniqueFastxFiles.find(fastXtar[i]);
		if (XX == uniqueFastxFiles.end()) {//no entry for this fastq yet
			uniqueFastxFiles[fastXtar[i]] = (int) uniqueFastxFiles.size();
			idx.push_back(vector<int>(1, i));
		} else {//exists already..
			idx[XX->second].push_back(i);
		}

		/*for (unsigned int j=0; j < uniqueFastxFiles.size(); j++){
			if (fastXtar[i] == uniqueFastxFiles[j]){ //the same
				idx[j].push_back(i);
				suc=true;
				break;
			}
		}
		if (!suc){
			uniqueFastxFiles.push_back(fastXtar[i]);
			idx.push_back(vector<int> (1,i));
		}*/
	}
	

	//unique Fas files set up.. check for their existence
	shared_ptr<InputStreamer> testFiles = 
		make_shared<InputStreamer>(!isFastq,mainFilter->getuserReqFastqVer(),"1","1",1);
	// For each unique Fa file
//	for (unsigned int i = 0; i < uniqueFastxFiles.size(); i++) {
	for (auto uFX: uniqueFastxFiles) {
			int tarID = idx[uFX.second][0]; string tmp;
		string x = testFiles->setupInput(path, tarID, uFX.first, FastqF, FastaF,
			QualF, MIDfq, mainFilter->isPaired(), cmdArgs["-onlyPair"], tmp, true);
	}
	mainFilter->SRessentials((int)uniqueFastxFiles.size());
//	delete testFiles;

    // mainFile for processing
	string mainFile = "";
	string outFile = cmdArgs["-o_fna"];
	
	//special sdm functions ini
	UClinks *ucl = nullptr; // Seed extension ?
	shared_ptr<ReadSubset> RDSset (nullptr);  //read subset filtered out?
	shared_ptr<Dereplicate> dereplicator (nullptr); //dereplication of b_derep_as_fasta_ input?
	ucl = mainFilter->ini_SeedsReadsDerep(ucl, RDSset, dereplicator); // actual prep
	//needs to attach to existing file sometimes
	std::ios_base::openmode writeStatus = ios_base::out;
	bool shortStats = false;
	string shrtLog = "";

	// main loop that goes over different files
	int maxReads = mainFilter->getXreads(); // should program stop after having written a certain amount of reads?
	int totalReadsRead(0);
	uint accumBPwrite(0), accumBPwriteMerg(0);
	string lastSRblock ("");  //set up SequencingRun blocks to track

    //---------------------------------
    // Multithreading setup
    //------------------------------------
    bool multithreading = true;
    //ThreadPool *pool = nullptr;
    int threads = 1;
    if (multithreading) {
        if (cmdArgs.find("-threads") != cmdArgs.end()) {
            threads = stoi(cmdArgs["-threads"]);
        }
		/*if ( threads > 1) {
			pool = new ThreadPool(threads);
		}*/
		cerr << "Run with " << threads << " cores.";
    }

	//set up a read merger for each thread..
	merger.resize(threads, nullptr);
	for (int x = 0; x < threads; x++) {
		merger[x] = new ReadMerger(true);
	}


    //---------------------------------
	//uniqueFastxFiles: unique input fastx's, can be demultiplexed fastas, or un-demultiplexed (several) fastx
	//loop that goes over several input fxs
	//------------------------------------
//	for (unsigned int i = 0; i < uniqueFastxFiles.size(); i++ ) {
	for (auto uFX : uniqueFastxFiles) {
		cdbg("Unique file " + uFX.first + "\n");
		uint i = uFX.second;
		if (maxReads > 0 && maxReads - totalReadsRead <= 0) { break; }
		if (idx[i].size() == 0) {	cerr << "fastXtar vector for " << uFX.first << " is empty" << endl;	exit(10);		}

		//create subset of BCs for the currently processed fastq's (only relevant BCs)
		cdbg("new filter in round " + itos(i) + "\n");
//		Filters* filter = mainFilter->filterPerBCgroup(idx[i]);
		Filters* filter = mainFilter->filterPerBCgroup(idx[i]);
		cdbg("Setting up threads\n");
		filter->setThreads(threads);
		int tarID = idx[i][0];//just points to one file in uFX group..


		//initialize object to handle all input file combinations
		//main input of fastx, handles input IO
		cdbg("Ini InputStream");
		shared_ptr<InputStreamer> IS = make_shared<InputStreamer>(
			!isFastq, mainFilter->getuserReqFastqVer(),cmdArgs["-ignore_IO_errors"],
			cmdArgs["-pairedRD_HD_out"], threads);

		// there is an entry in tar for each barcode for this file. If idx[i].size() is 1 there is only one barcode
		if (idx[i].size() == 1 && uniqueFastxFiles.size() > 1) {
			IS->atFileYofX(i + 1, (unsigned int)uniqueFastxFiles.size(), 1);
		}

		// check first if derep block will change
		if (totalReadsRead>0 && mainFilter->doDereplicate() && dereplicator->DerepPerSR() 
			&& lastSRblock != mainFilter->SequencingRun[tarID]) {
				string deLogLocal = dereplicator->writeDereplDNA(mainFilter, lastSRblock);
				if (dereplicator->DerepPerSR()) {
					deLog += "Dereplication of SamplingRun " + lastSRblock + ":\n";
				}
				deLog += deLogLocal;
				if (cmdArgs["-log"] != "nolog") {
					dereplicator->writeLog(logF.substr(0, logF.length() - 3) + "dere", deLogLocal);
				}
				//in this case also needs to recheck merger prob
				if (merger[0] != NULL) {
					ReadMerger * cMerg = new ReadMerger();
					for (size_t x = 0; x < merger.size(); x++) {
						cMerg->addRMstats(merger[x]);
					}
					cMerg->printMergeHisto(subfile(cmdArgs["-merg_readpos"], lastSRblock));
					cMerg->printQualHisto(subfile(cmdArgs["-qual_readpos"], lastSRblock));
					cMerg->restStats();
					delete cMerg;
				}


			//continue?

				dereplicator->reset();//and reset..
				accumBPwrite = 0; accumBPwriteMerg = 0;
				//OutStreamer->resetDemultiBPperSR();
		}
		lastSRblock = mainFilter->SequencingRun[tarID];

		string mainFileShort;
		mainFile = IS->setupInput(path, tarID, uFX.first, FastqF, FastaF, QualF, MIDfq, filter->isPaired(),
			cmdArgs["-onlyPair"], mainFileShort, false);
		if (!IS->qualityPresent()) {
			filter->deactivateQualFilter();
			cerr << "\n*********\nWarning:: Quality file is not present.\nRecommended to abort demultiplexing.\n*********\n\n";
		}

		filter->BarcodePreStats();filter->checkDoubleBarcode();filter->checkDoubleSampleIDHead();


		if (mainFilter->doOptimalClusterSeq()) {
			ucl->findSeq2UCinstruction(IS, isFastq, mainFilter);
			continue;// after this is only qual_ filter, not required from this point on
		}
		cdbg("Setting up output\n");



		//OutputStreamer OutStreamer = OutputStreamer(&filter, cmdArgs, writeStatus, RDSset);

		//OutputStreamer also contains subfilters for MC processing and logging of reads
		shared_ptr<OutputStreamer> OutStreamer = make_shared<OutputStreamer>(filter, cmdArgs, 
			writeStatus, RDSset, "", threads);
		OutStreamer->attachDereplicator(dereplicator);
		OutStreamer->attachReadMerger(merger);
		OutStreamer->setBPwrittenInSR(accumBPwrite);
		OutStreamer->setBPwrittenInSRmerg(accumBPwriteMerg);	
		filter->setMultiDNA(OutStreamer);
		if (maxReads > 0) {
			OutStreamer->setReadLimit(maxReads - totalReadsRead);
		}
		writeStatus = ofstream::app;
		//prepare for BC checking (rev/fwd)
		if (filter->doDemultiplex()) {
			OutStreamer->setBCfixed(false, true);
			if (OutStreamer->isPEseq() == 2) {
				OutStreamer->setBCfixed(false, false);
			}
		}

		//only pull out a subset of sequences
		if (mainFilter->doSubselReads()) {
			if (cmdArgs.find("-mocatFix") != cmdArgs.end()) {
				cerr << "MOCAT fix appplies\n";
				RDSset->findMatches(IS, OutStreamer, true);
			}	else {
				RDSset->findMatches(IS, OutStreamer, false);
			}
			//delete OutStreamer;
			continue;
		}

		cdbg("Processing reads");

		// If multhreading is switched on and thread count is > 1  call multihtreaded functions.
		// if multithreading = false ignore thread option
		multithreading &= threads > 1;
		//**********************
		//heavy reading routine
		//**********************
		
		if (OutStreamer->isPEseq() == 2) {
			read_paired(cmdArgs, OutStreamer, IS, IS->hasMIDseqs(), threads);
		}else {
			read_single(cmdArgs, OutStreamer, IS);
			/*if (threads == 1) {
			}	else if (threads >= 2) {
				readSingleTp(cmdArgs, OutStreamer, IS, pool);
			}*/
		}
		//debug
		//std::cout << "asdlkjsdlkjsad> " << filter->statistics_[0].main_read_stats_[0]->total << std::endl;
		
		//here all subfilter can be merged (to get stats right)
		OutStreamer->mergeSubFilters();


		cdbg("All reads processed in uFX");
		outFile = OutStreamer->leadOutFile();
		OutStreamer->closeOutStreams(true);
		OutStreamer->closeOutFilesDemulti();
		accumBPwrite = OutStreamer->getBPwrittenInSR();
		accumBPwriteMerg = OutStreamer->getBPwrittenInSRmerg();
		cdbg("OutStreamer deleted for current uFX");

		//stats
		// at this point subfilters have already been merged in clouseOutStreams
		filter->prepStats();
		if (IS->getCurFileN() == 0) {
			filter->printStats(cerr, mainFile, outFile, true);
		}
		else {
			cerr << filter->shortStats(""); shortStats = true;
		}

		totalReadsRead += filter->totalAccepts();

		shrtLog += filter->shortStats(mainFileShort);
		//		delete IS;
				//write log file
		if (uniqueFastxFiles.size() > 1) {//only print sub log if neccessary
			ofstream log;
			string logF = cmdArgs["-log"] + string("0") + itos(i);
			log.open(logF.c_str(), ios_base::out);
			filter->printStats(log, mainFile, outFile, true);
			log.close();
		}

        mainFilter->addStats(filter, idx[i]);
		delete filter;
        // Print out merged read count (move to stort stats)
		if (OutStreamer->total_read_preMerge_) {
			std::cerr << "merged reads: " << OutStreamer->merged_counter_ << "/" << OutStreamer->total_read_preMerge_ << " (" << (double)OutStreamer->merged_counter_ / OutStreamer->total_read_preMerge_ << ")" << std::endl;
		}


    }//uFX
	cdbg("End of UFx loop");


    //-------------------------------------------------
    //reading of input fq's is done
    //postprocessing, write reads, dereplicates etc
    // ------------------------------------------------

//write log files
	cdbg("Prep final logging");
    if (uniqueFastxFiles.size() > 1){ mainFile = "several";}
    ofstream log; 
    //different logfile for SEED extension
    if (mainFilter->doOptimalClusterSeq()){
        //finish up dereplication file (creating pseudo seeds with counts)
        ucl->finishMAPfile();
        if (cmdArgs["-ucAdditionalCounts"] != ""){
            ucl->set2UC();
            ucl->finishUCfile(mainFilter, cmdArgs["-ucAdditionalCounts"], true);//with smplHead (.mid)
            ucl->finishUCfile(mainFilter, cmdArgs["-ucAdditionalCounts1"], false);//without smplHead (.rest)
        }
        if (cmdArgs["-ucAdditionalCounts_refclust"] != ""){
            //reference based clustering has some high qual_ seqs (no replacement with reads..)
            //takeOver even found high qual_ hits with these default seeds..
            shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true,
				mainFilter->getuserReqFastqVer(), cmdArgs["-ignore_IO_errors"], 
				cmdArgs["-pairedRD_HD_out"],1);
            //this reads in the SLV fna's & creates matrix entries for these
            FALL->setupFna(cmdArgs["-OTU_fallback_refclust"]);
            ucl->setRefMode();
            ucl->addDefSeeds(FALL, mainFilter);
            ucl->set2UC();
            //mapping from ref OTU clustering
            ucl->finishUCfile(mainFilter, cmdArgs["-optimalRead2Cluster_ref"], false);
            //mid / rest mappings
            ucl->finishUCfile(mainFilter, cmdArgs["-ucAdditionalCounts_refclust"], true);//with smplHead (.rest)
            ucl->finishUCfile(mainFilter, cmdArgs["-ucAdditionalCounts_refclust1"], false);//without smplHead (.rest)

        }
        if (cmdArgs["-log"] != "nolog") {
            log.open(logF.c_str(), ios_base::out);
            ucl->printStats(cerr);
            ucl->printStats(log);
            log.close();
        }

        //everything done on DNA? Then write & delete
        if (cmdArgs["-otu_matrix"] != "") {
            ucl->writeOTUmatrix(cmdArgs["-otu_matrix"]);
        }
        //polished OTU seeds need to be written after OTU matrix (renaming scheme)
        shared_ptr<OutputStreamer> MDx = make_shared<OutputStreamer>(mainFilter, cmdArgs, 
			ios::app, RDSset);
		vector<ReadMerger*> DerepM = vector<ReadMerger*>(1,NULL);
		DerepM[0] = new ReadMerger(false);
		MDx->attachReadMerger(DerepM);
        mainFilter->setMultiDNA(MDx);
        ucl->writeNewSeeds(MDx, mainFilter, false);
        //new fastas also need to be written..
        MDx.reset(new OutputStreamer(mainFilter, cmdArgs, 
			ios::app, RDSset, ".ref", 1));//force fna output
        mainFilter->setMultiDNA(MDx );
        ucl->writeNewSeeds(MDx, mainFilter, true, true);
        //delete MDx;
		delete DerepM[0];
        return;
    } else if (mainFilter->doDereplicate()) {
		//this is either the last time dereplicate is written (SRblocks),
		//or the main point at which to write the dereplicator for good
#ifdef DEBUG
        cerr << "write Dereplicated DNA" << endl;
#endif
		
        string deLogLocal = dereplicator->writeDereplDNA(mainFilter,lastSRblock);
		if (dereplicator->DerepPerSR()) {
			deLog += "Dereplication of SamplingRun " + lastSRblock + ":\n";
		}
		deLog += deLogLocal;
		if (cmdArgs["-log"] != "nolog") {
			dereplicator->writeLog(logF.substr(0, logF.length() - 3) + "dere" , deLogLocal);
		}
		dereplicator->finishMap();

		//last time merger stats to write
		if (merger[0] != nullptr) {
			ReadMerger* cMerg = new ReadMerger();
			for (size_t x = 0; x < merger.size(); x++) {
				cMerg->addRMstats(merger[x]);
			}
			cMerg->printMergeHisto(subfile(cmdArgs["-merg_readpos"], lastSRblock));
			cMerg->printQualHisto(subfile(cmdArgs["-qual_readpos"], lastSRblock));
			delete cMerg;
		}


#ifdef DEBUG
        cerr << "done write dereplicator" << endl;
#endif
    }
#ifdef DEBUG
    cerr << "Logging almost finished" << endl;
#endif

    if (cmdArgs["-log"] == "nolog") {
        return;
    }
    if (shortStats) {
        mainFilter->printStats(std::cerr, mainFile, outFile, true);
    }
#ifdef DEBUG
    cerr << "other logs start" << endl;
#endif
    //per sample success rate
    string logPS = logF.substr(0, logF.length() - 3) + "acceptsPerSample.log";
    log.open(logPS.c_str(), ios_base::out);
    mainFilter->SmplSpecStats(log);
    log.close();
    log.open(logF.c_str(), ios_base::out);
    mainFilter->printStats(log, mainFile, outFile, true);
    log.close();

    string logFs = logF.substr(0, logF.length() - 3) + "acceptsPerFile.log";
    log.open(logFs.c_str(), ios_base::out);
    log << shrtLog;
    log.close();

    string logFGC = logF.substr(0, logF.length() - 3) + "GC.txt";
    log.open(logFGC.c_str(), ios_base::out);
    mainFilter->printGC(log, mainFilter->isPaired());
    log.close();
#ifdef DEBUG
    cerr << "other logs end" << endl;
#endif


//for additional files
    if (mainFilter->secondaryOutput()){
        log.open (logFA.c_str() ,ios_base::out);
        mainFilter->printStats(log, mainFile, outFile, false);
        log.close();
    }


    //length histogram
    logF = cmdArgs["-length_hist"];
    log.open (logF.c_str() ,ios_base::out);
    mainFilter->printHisto(log, 0);
    log.close();
    //quality histogram
    logF = cmdArgs["-qual_hist"];
    log.open (logF.c_str() ,ios_base::out);
    mainFilter->printHisto(log, 1);
    log.close();

#ifdef DEBUG
    cerr << "separateByFile finished" << endl;
#endif

    // multithreading
    //delete pool;
	//ReadMerger no longer needed
	for (size_t x = 0; x < merger.size(); x++) {
		delete merger[x];
	}

}

void rewriteNumbers(OptContainer& cmdArgs){
    //no renumbering asked for
    if (!(cmdArgs.find("-number")  != cmdArgs.end() && cmdArgs["-number"]=="T")){
        return;
    }
    string prefix="";
    if (cmdArgs.find("-prefix")  != cmdArgs.end()){
        prefix = cmdArgs["-prefix"];
    }
    //read fasta & write with new headers
    int cnt=0;

    string line;
    //ofstream qualOut,fnaOut;
    ifstream fna;
    ofstream ofna;

    //        rerwite input fasta file
    string tname="",tseq="";
    fna.open(cmdArgs["-i_fna"].c_str(),ios::in);
    ofna.open(cmdArgs["-o_fna"].c_str(),ios::out);
    while (getline(fna,line,'\n')){

        if (line[0]=='$'){ //$ marks comment
            continue;
        }
        if(line[0] == '>'){ //fasta description
            if (cnt!=0){
                tname = ">"+prefix+itos(cnt)+"\n";
                ofna << tname << tseq;
            }
            cnt++;tseq="";
            continue;
        }
        tseq += line+"\n";
    }
    tname = ">"+prefix+itos(cnt)+"\n";
    ofna << tname << tseq;
    ofna.close(); fna.close();

    exit(0);
}

void Announce_sdm(){
    cerr << endl << "This is sdm (simple demultiplexer) " << sdm_version << " " << sdm_status << ".\n" << endl ;
}
void help_head(){
    cout <<"------------------------------\nThis is sdm version "<<sdm_version <<" "<< sdm_status <<" help print\n------------------------------\n"<<endl;
}
void general_help(){
    help_head();
    cout<<"sdm (simple demultiplexer) is a fast, memory efficient program to demultiplex fasta and fastq files or simply do quality filterings on these.\n";
#ifdef _gzipread
    cout<<"Compiled with gzip support\n";
#else
    cout << "No gzip support compiled\n";
#endif
#ifdef _THREADED
    cout<<"Multithreading not supported"
#else
    cout << "The compiled version does not support multithreading\n";
#endif
    cout<<"Select further help topics by typing:\nsdm -help_options : print help on configuring options files\nsdm -help_commands : help on command arguments for sdm\nsdm -help_map : base_map files and the keywords for barcodes etc.\n------------------------------\n";
    cout<<"Author: falk.hildebrand@gmail.com"<<endl;

}
void printCmdsHelp(){
    help_head();
    string def_sep = DEFAULT_BarcodeNameSep;
    cout << "Usage:\n./sdm\n  -i_path <path to several fastq / fasta files>\n------OR------\n -i <input sequence file, will autodetect fna/fastq>\n------OR------\n -i_fastq <fastQ file>\n------OR------\n -i_fna <your fasta input file> (required)\n -i_qual <corresponding quality file> (required, unless quality file is \"xx1.qual_\" and fasta is \"xx1.yy\")\n\n -base_map <mapping file in Qiime format> (optional)\n -o_fna <file to write output fasta> (optional)\n -o_qual <file to write corresponding quality values> (optional)\n -o_fastq <fastQ output file (overrides -o_qual & -o_fna)\n";
    cout << " -options <sdm option file>(optional)\n -log <file to save demultiplex log in>(optional). Set to \"nolog\" to deactivate alltogether.\n \n-sample_sep \"X\" string X is used to delimit samples and id_ (optional, default:\"" << def_sep << "\")\n -paired 1/2/3 (input is paired end sequenced(2), assumes two input files delimited by \',\'. 1=singleton (default); 3=paired end (R1,R3) + one file with MID (R2))\n";
    cout << " -o_demultiplex [path] write input into single, demultiplexed files\n";
    cout << " -onlyPair [1/2] consider only read pair_ 1 or 2. Useful when streamlining inputs (LotuS) or considering double barcoding.\n -i_MID_fastq fastq file with only MID sequences; if paired reads are supplied with -i_fna/-i_fastq and the MID identifier via -i_MID_fastq, paired has to be set to 2. If e.g. merged reads are supplied + mids, paired has to be set to 1.\n";
    cout << " -pairedDemulti [0/1] Only write complete pairs in demultiplexing paired input.\n";
    cout << " -oneLineFastaFormat [0/1] write Fasta and Quality file sequence string in one line, opposed to default 80 characters per line.\n -o_dereplicate <output fasta file> of dereplicated DNA reads (with size in header)\n -dere_size_fmt [0/1] either (0) usearch format \"size=X;\" or (1) \"_X\"\n -min_derep_copies only print seq if at least X copies present. Can be complex terms like \"10:1,3:3\" -> meaning at least 10x in 1 sample or 3x in 3 different samples.\n";
    cout << " -SyncReadPairs [T/F] sdm can check, if read pairs occur in the same (correct) order in the input files, and correct this in case not (T).\n";
    cout << " -maxReadsPerOutput number of filtered reads in output files. If more reads, a new file is created. Only works with -o_fna\n -mergedPairs <1/0> 1: paired sequences were merged externally, important for assumption that read quality is detoriating.\n -OTU_fallback <file>: Fallback fasta sequences for OTU's, only used in SEED extension mode\n";
    cout << " -i_qual_offset [0-64] fastq offset for quality values. Set this to \'0\' or \'auto\' if you are unsure which fastq version is being used (default: read from sdm option file)\n -o_qual_offset [0-64] set quality offset for fastq outfile. Default: 33\n";
    cout << " -ignore_IO_errors [0/1]: 1=Errors in fastq reads are ignored, with sdm trying to sync reads pairs after corrupted single reads (default: 0)\n";
    //-binomialFilterBothPairs [1/0]
    //-count_chimeras [T/F]
    // ucAdditionalCounts_refclust -OTU_fallback_refclust -optimalRead2Cluster_ref
    cout<<"\nMinimal Example:\n./sdm -i test.fna -base_map mapping.txt (assuming quality file is \"test.qual_\")\n";
    // further options (undocumented) :
    //-length_hist   -qual_hist
    //-suppressOutput[0/1]
    //
}
void printOptionHelp(){
    help_head();
    cout<<"The option file, specified via the \"-options\" argument, provides more specific control over filtering, barcode handling, and sequencing technologies, among others. A reference option file is printed below.\n\n";
    string helpOptionFile="";
    /*helpOptionFile += "minSeqLength - minimal accepted Sequence Length\nmaxSeqLength - maximal Length of Sequence\nminAvgQuality - minimal average Quality\nmaxAmbiguousNT - max number of Ambigous nt's in sequence\nQualWindowThreshhold - Q threshold where seq is rejected\nQualWindowWidth - average quality in this windows is used for QualWindowThreshhold\n";
    helpOptionFile += string("TrimWindowThreshhold - Q value below which sequence is 3' trimmed\nTrimWindowWidth - window size used for TrimWindowThreshhold\nmaxBarcodeErrs - max accepted barcode errors\nmaxPrimerErrs - max accepted Primer errors\nkeepBarcodeSeq - leave Barcode sequence_ on read? (0/1)\n");
    helpOptionFile += "keepPrimerSeq - keep Primer attached to seq? (0/1)\nmaxHomonucleotide - sequences with a homonucleotide run longer will be rejected\nmaxAccumulatedError - if P is surpassed, sequence is trimmed at that point\nTechnicalAdapter - sequence of the technical adapter (will be removed, if found 5')\nPEheaderPairFmt - ?\nTrimStartNTs - trim X nucleotides from the start of the sequence ";
    helpOptionFile += "fastqVersion - 1 = ";*/

    helpOptionFile += "#--- Example ---\n#copy into new file\n#sequence length refers to sequence length AFTER removal of Primers, Barcodes and trimming. this ensures that downstream analyis tools will have appropiate sequence information\nminSeqLength	250\nmaxSeqLength	1000\nminAvgQuality	25\n\n";
    helpOptionFile += "#Ambiguous bases in Sequence - uclust only supports 0 ambiguous nucleotides\nmaxAmbiguousNT	0\n\n#Homonucleotide Runs.. this should normally be filtered by sequencer software\nmaxHomonucleotide	8\n\n";
    helpOptionFile += "#Filter whole sequence if one window of quality scores is below average\nQualWindowWidth	50\nQualWindowThreshhold	25\n\n#Trim the end of a sequence if a window falls below quality threshhold. Useful for removing low qulaity trailing ends of sequence\n\nTrimWindowWidth	20\nTrimWindowThreshhold	25\n\n#Max number of accumulated P for a mismatch. After this length, the rest of the sequence will be deleted. Complimentary to TrimWindowThreshhold. (-1) deactivates this option.\nmaxAccumulatedError	1\n\n";
    helpOptionFile += "#Barcode Errors - currently this can only be 0; \nmaxBarcodeErrs	0\nmaxPrimerErrs	0\n\n#keep Barcode / Primer Sequence in the output fasta file - in a normal 16S analysis this should be deactivated (0) for Barcode and de-activated (0) for primer\nkeepBarcodeSeq	0\nkeepPrimerSeq	0\n\n";
    helpOptionFile += "#set fastqVersion to 1 if you use Sanger, Illumina 1.8+ or NCBI SRA files. Set fastqVersion to 2, if you use Illumina 1.3+ - 1.7+ or Solexa fastq files.\n\nfastqVersion	1\n\n#if one or more files have a technical adapter still included (e.g. TCAG 454) this can be removed by setting this option\n\nTechnicalAdapter	TCAG\n\n#delete X NTs (e.g. if the first 5 bases are known to have strange biases)\n\nTrimStartNTs	0\n";
    helpOptionFile += "#truncate total Sequence length to X (length after Barcode, Adapter and Primer removals)\nTruncateSequenceLength	200\n";
    helpOptionFile += "#correct PE header format (1/2) this is to accomodate the illumina miSeq paired end annotations 2=\"@XXX 1:0:4\" instead of 1=\"@XXX/1\". Note that the format will be automatically detected\nPEheaderPairFmt	1\n\n#sets if sequences without match to reverse primer will be accepted (T=reject ; F=accept all); default=F\nRejectSeqWithoutRevPrim	F\n";
    helpOptionFile += "#sets if sequences without a forward (LinkerPrimerSequence) primer will be accepted (T=reject ; F=accept all); default=T\nRejectSeqWithoutFwdPrim	T\n\n";
    helpOptionFile += "#checks if the whole amplicon was reverse-transcribed sequenced (Default = F)\nCheckForReversedSeqs	F\n\n";
    helpOptionFile += "#this option should be \"T\" if your amplicons are possibly shorter than a read in a paired end sequencing run (e.g. amplicon of 300 in 250x2 miSeq is \"T\")\nAmpliconShortPE	T\n\n";
    //CheckForMixedPairs CheckForReversedSeqs
    cout<<helpOptionFile<<endl;
}
void printMapHelp(){
    help_head();
    cout<<"The mapping file, specified via the \"-base_map\" argument, contains all information neccessary to demultiplex a sequencer fasta or fastq output file. If left out, the sequences are only quality checked.\n";
    cout<<"The Mapping file can contain comments, specified by \"#\" at the start of the comment line. Similarly, the header (required) has to start with \"#\"\n. Processed header fields are:\n";
    cout << "SampleID - sample_id_ Identifier, has to be unique for each Barcode\n";
    cout << "CombineSamples - combines sample_id_ Identifiers into one sample\n";
    cout << "BarcodeSequence - The Barcode (MID) tag assigned to each sample. Can contain IUPAC redundant nucleotides\n";
    cout << "Barcode2ndPair - in case of dual indexed reads, use this column to specify the BC on the 2nd read pair_\n";
    cout<<"ForwardPrimer (previously LinkerPrimerSequence) - Sequence used for 16S amplification, usually (unless paired end mode) is after the Barcode\n";
    cout<< "ReversePrimer - Reverse Primer Sequence (IUPAC code).\n";
    cout<< "fastqFile - if used in -i_path mode, gives relative location of fastq file, such that [-i_path][fastqFile] gives the absolute path to fastq file.\n";
    cout<<"fnaFile - see fastqFile. However, fasta formated file instead of fastq format.\n";
    cout<< "qualFile - see fastqFile. However, quality file corresponding to fasta file instead of fastq format.\n";
    cout<< "MIDfqFile - fastq file containing ONLY the MID sequence, corresponding fastq input file(s)\n";
    cout<< "SampleIDinHead - id_ in header of fasta/fastq file, that identifies sample_id_ - replaces Barcode (MID) scanning.\n";
    cout << "derepMin - this column can contain a single number X, only accepting reads in specific sample that have more than X dereplicated copies. Useful if a subset of samples was e.g. sampled 100x as deep as other samples.\nCombined sample_id_ WILL NOT be taken into account.\n";
    cout<<"\n";

}
void printVersion(){
    cout << "sdm " << sdm_version << " " << sdm_status << endl;
}


//////////////////////

/*
//manages read in of several input files and associated primers / tags to each file
void separateByFile(Filters* mainFil,OptContainer& cmdArgs){
//void separateByFile(Filters* mainFil,OptContainer& cmdArgs){
#ifdef DEBUG
    cerr << "separateByFile"<<endl;
#endif

//	mainFil->ini_filestruct(cmdArgs);
    
    vector<string> FastaF = mainFil->getFastaFiles();
    vector<string> QualF = mainFil->getQualFiles();
    vector<string> FastqF = mainFil->getFastqFiles();
    vector<string> MIDfq = mainFil->getMIDfqFiles();
    vector<string> tar;
    vector < vector<int> > idx(0);
    bool bFASTQ = true;
    //prepareOutFiles(cmdArgs);
    string path="";
    if (cmdArgs.find("-i_path")  != cmdArgs.end() && cmdArgs["-i_path"].length() > 2){
        path=cmdArgs["-i_path"] + string("/");
    }
    
    
    if (FastaF.size()>0){ //fasta way
        tar = FastaF;
        bFASTQ = false;
    } else { // fastq way
        tar = FastqF;
        if (FastqF.size()==0){
            cerr<<"No FastQ or Fasta file given.\n  Aborting..\n";
            exit(12);
        }
    }
    
    vector<string> uniqueFas(1,tar[0]);
    idx.push_back(vector<int> (1,0));
    
    for (unsigned int i=1; i<tar.size(); i++){
        bool suc = false;
        for (unsigned int j=0; j<uniqueFas.size(); j++){
            if (tar[i] == uniqueFas[j]){ //the same
                idx[j].push_back(i);
                suc=true;
                break;
            }
        }
        if (!suc){
            uniqueFas.push_back(tar[i]);
            idx.push_back(vector<int> (1,i));
        }
    }
    
    //unique Fas files set up.. check for their existence
    shared_ptr<InputStreamer> testFiles =
            make_shared<InputStreamer>(!bFASTQ, mainFil->getuserReqFastqVer(), "1","1");
    for (unsigned int i = 0; i < uniqueFas.size(); i++) {
        int tarID = idx[i][0]; string tmp;
        string x = testFiles->setupInput(path, i, tarID, uniqueFas, FastqF, FastaF, QualF, MIDfq, mainFil->isPaired(), cmdArgs["-onlyPair"], tmp, true);
    }
//	delete testFiles;
    
    string mainFile = "", outFile = cmdArgs["-o_fna"];
    
    //special sdm functions ini
    UClinks *ucl = NULL; // Seed extension ?
    shared_ptr<ReadSubset> RDSset;  //read subset filtered out?
    shared_ptr<Dereplicate> Dere ; //dereplication of fasta input?
    mainFil->ini_SeedsReadsDerep(ucl, RDSset, Dere); // actual prep
    
    
    //needs to attach to existing file sometimes
    std::ios_base::openmode writeStatus = ios_base::out;
    bool shortStats = false;
    string shrtLog = "";
    // main loop that goes over different files
    int maxRds = mainFil->getXreads();
    int totReadsRead(0);
    //---------------------------------
    //uniqueFas: unique input fastx's, can be demultiplexed fastas, or un-demultiplexed (several) fastx
    //loop that goes over several input fq's
    //------------------------------------
    for (uint i=0; i<uniqueFas.size();i++ ){
#ifdef DEBUG
        cerr << "new filter in round "<<i << endl;
#endif
        if (maxRds>0 && maxRds - totReadsRead <= 0) { break; }
        Filters* fil = make_shared<Filters>(mainFil, idx[i][0]);
        unsigned int tarSi = (unsigned int) idx[i].size();
        fil->allResize(tarSi);
        int tarID=-1;
        bool BC2mode = mainFil->doubleBarcodes();
        //int readsRead(0);
        
        
        for (unsigned int j=0; j<tarSi;j++){ //fill in filter
            tarID = idx[i][j];
            if (mainFil->PrimerIdx[tarID]>-1) {
                fil->addPrimerL(mainFil->PrimerL[mainFil->PrimerIdx[tarID]], j);
            }
            if (mainFil->doReversePrimers() && mainFil->PrimerIdxRev[tarID]>-1) {
                fil->addPrimerR(mainFil->PrimerR[mainFil->PrimerIdxRev[tarID]], j);
            }
            fil->Barcode[j] = mainFil->Barcode[tarID];
            if ( BC2mode ) {
                fil->Barcode2[j] = mainFil->Barcode2[tarID];
            }
            fil->SampleID[j] = mainFil->SampleID[tarID];
            fil->SampleID_Combi[j] = mainFil->SampleID_Combi[tarID];
            fil->HeadSmplID[j] = mainFil->HeadSmplID[tarID];
            if (fil->Demulti2Fls()) {
                fil->demultiSinglFiles[j] = mainFil->demultiSinglFiles[tarID];
                fil->demultiSinglFilesF[j] = mainFil->demultiSinglFilesF[tarID];
                
                //closing of ofstreams is only handled on the main object
                //mainFil->demultiSinglFiles[tarID] = vector<ofstream*> (2,NULL);
            }
        }
        //sanity check no double barcodes..
        fil->checkDoubleBarcode();
        
        if (tarID==-1){cerr<<"tar == -1. abort.\n";exit(10);}
        
        //initialize object to handle all input file combinations
        //main input of fastx, handles input IO
        shared_ptr<InputStreamer> IS = make_shared<InputStreamer>(!bFASTQ, mainFil->getuserReqFastqVer(),
                                                                  cmdArgs["-ignore_IO_errors"], cmdArgs["-pairedRD_HD_out"]);
        if (tarSi < 2 && uniqueFas.size() > 1) {
            IS->atFileYofX(i + 1, (uint)uniqueFas.size(), tarSi);
        }
        string mainFileShort = "";
        mainFile = IS->setupInput(path, i, tarID, uniqueFas, FastqF, FastaF, QualF, MIDfq, fil->isPaired(), cmdArgs["-onlyPair"], mainFileShort, false);
        if (!IS->qualityPresent()) {
            fil->deactivateQualFilter();
            cerr << "\n*********\nWarning:: Quality file is not present.\nRecommended to abort demultiplexing.\n*********\n\n";
        }
        fil->BarcodePreStats();
        fil->checkDoubleBarcode();
        fil->checkDoubleSampleIDHead();
        
        
        if (mainFil->doOptimalClusterSeq()){
            ucl->findSeq2UCinstruction(IS,bFASTQ, mainFil);
            continue;// after this is only qual_ filter, not required from this point on
        }


#ifdef DEBUG
        cerr << "Setting up output" << endl;
#endif
        //OutputStreamer OutStreamer = OutputStreamer(&fil, cmdArgs, writeStatus, RDSset);
        shared_ptr<OutputStreamer> OutStreamer = make_shared<OutputStreamer>(fil, cmdArgs, writeStatus, RDSset);
        fil->setMultiDNA(OutStreamer);
        if (maxRds > 0) { OutStreamer->setReadLimit(maxRds - totReadsRead); }
        writeStatus = ofstream::app;
        //prepare for BC checking (rev/fwd)
        if (fil->doDemultiplex()){
            OutStreamer->setBCfixed(false, true);
            if (OutStreamer->isPEseq() == 2) { OutStreamer->setBCfixed(false, false); }
        }
        if (cmdArgs.find("-oneLineFastaFormat") != cmdArgs.end() && cmdArgs["-oneLineFastaFormat"] == "1") {
            OutStreamer->setOneLinerFastaFmt(true);
        }
        //cout << Dere->Nms_size() << " DEBCs\n";
        OutStreamer->attachDereplicator(Dere);
        //only pull out a subset of sequences
        if (mainFil->doSubselReads()) {
            if (cmdArgs.find("-mocatFix") != cmdArgs.end()) {
                cerr << "MOCAT fix appplies\n";
                RDSset->findMatches(IS, OutStreamer, true);
            }else{
                RDSset->findMatches(IS, OutStreamer, false);
            }
            //delete OutStreamer;
            continue;
        }


#ifdef DEBUG
        cerr << "Processing reads" << endl;
#endif
        
        
        int threads = 1;
        if (cmdArgs.find("-threads") != cmdArgs.end()) {
            threads = stoi(cmdArgs["-threads"]);
        }
        
        
        //**********************
        //heavy reading routine
        //**********************
        if (OutStreamer->isPEseq() == 2){
            //read_paired(cmdArgs,OutStreamer,IS);
            while ( !read_paired(cmdArgs, OutStreamer, IS, IS->hasMIDseqs()) ) {
                //reset output files to previous state
                OutStreamer->resetOutFilesAndFilter();
            }
        } else {
            
            if (threads == 1) {
                cout << "single_threaded" << endl;
                read_single(cmdArgs,OutStreamer,IS);
                cout << "processed: " << counter << endl;
            } else if (threads == 2) {
                cout << "multi_threaded" << endl;
                read_single_threaded2(cmdArgs, OutStreamer, IS);
            } else if (threads > 2) {
                cout << "multi_threaded_mixed" << endl;
                read_single_mixed(cmdArgs, OutStreamer, IS);
            }
            //read_single_threaded(cmdArgs, OutStreamer, IS, 6);
        }
        wait.printResults();
        put_bm.printResults();
        get_dna.printResults();


#ifdef DEBUG
        cerr << "All read processed" << endl;
#endif
        outFile = OutStreamer->leadOutFile();
//		delete OutStreamer;
#ifdef DEBUG
        cerr << "OutStreamer deleted" << endl;
#endif
        
        //stats
        fil->prepStats();
        if (IS->getCurFileN() == 0) {
            fil->printStats(cerr, mainFile, outFile, true);
        } else {
            cerr<<fil->shortStats(""); shortStats = true;
        }
        
        totReadsRead += fil->totalAccepts();
        
        shrtLog += fil->shortStats( mainFileShort);
//		delete IS;
        //write log file
        if (uniqueFas.size() > 1){//only print sub log if neccessary
            ofstream log;
            string logF = cmdArgs["-log"] + string("0") + itos(i);
            log.open (logF.c_str() ,ios_base::out);
            fil->printStats(log,mainFile,outFile,true);
            log.close();
        }
        mainFil->addStats(fil,idx[i]);
#ifdef DEBUG
        cerr << "Delete tmp filter" << endl;
#endif
    }
    
    
    //-------------------------------------------------
    //reading of input fq's is done
    //postprocessing, write reads, dereplicates etc
    // ------------------------------------------------


#ifdef DEBUG
    cerr << "Prep final logging" << endl;
#endif
//write log files
    if (uniqueFas.size() > 1){
        mainFile = "several";
    }
    
    ofstream log; string deLog("");
    string logF = cmdArgs["-log"], logFA = cmdArgs["-log"].substr(0, cmdArgs["-log"].length()-3) + "add.log";
    
    //different logfile for SEED extension
    if (mainFil->doOptimalClusterSeq()){
        //finish up dereplication file (creating pseudo seeds with counts)
        ucl->finishMAPfile();
        if (cmdArgs["-ucAdditionalCounts"] != ""){
            ucl->set2UC();
            ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts"], true);//with smplHead (.mid)
            ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts1"], false);//without smplHead (.rest)
        }
        if (cmdArgs["-ucAdditionalCounts_refclust"] != ""){
            //reference based clustering has some high qual_ seqs (no replacement with reads..)
            //takeOver even found high qual_ hits with these default seeds..
            shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true,
                                                                        mainFil->getuserReqFastqVer(), cmdArgs["-ignore_IO_errors"], cmdArgs["-pairedRD_HD_out"]);
            //this reads in the SLV fna's & creates matrix entries for these
            FALL->setupFna(cmdArgs["-OTU_fallback_refclust"]);
            ucl->setRefMode();
            ucl->addDefSeeds(FALL, mainFil);
            ucl->set2UC();
            //mapping from ref OTU clustering
            ucl->finishUCfile(mainFil, cmdArgs["-optimalRead2Cluster_ref"], false);
            //mid / rest mappings
            ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts_refclust"], true);//with smplHead (.rest)
            ucl->finishUCfile(mainFil, cmdArgs["-ucAdditionalCounts_refclust1"], false);//without smplHead (.rest)
            
        }
        if (cmdArgs["-log"] != "nolog") {
            log.open(logF.c_str(), ios_base::out);
            ucl->printStats(cerr);
            ucl->printStats(log);
            log.close();
        }
        
        //everything done on DNA? Then write & delete
        if (cmdArgs["-otu_matrix"] != "") {
            ucl->writeOTUmatrix(cmdArgs["-otu_matrix"]);
        }
        //needs to be written after OTU matrix (renaming scheme)
        shared_ptr<OutputStreamer> OutStreamer = make_shared<OutputStreamer>(mainFil, cmdArgs, ios::out, RDSset);
        mainFil->setMultiDNA(OutStreamer);
        ucl->writeNewSeeds(OutStreamer, mainFil,false);
        //delete OutStreamer;
        //new fastas also need to be written..
        OutStreamer.reset(new OutputStreamer(mainFil, cmdArgs, ios::out, RDSset, ".ref", 1));//force fna output
        mainFil->setMultiDNA( OutStreamer );
        ucl->writeNewSeeds(OutStreamer, mainFil,true,true);
        //delete OutStreamer;
        
        
        return;
    } else if (mainFil->doDereplicate()) {
#ifdef DEBUG
        cerr << "write Dereplicated DNA" << endl;
#endif
        deLog = Dere->writeDereplDNA(mainFil);
#ifdef DEBUG
        cerr << "done write Dere" << endl;
#endif
    }
#ifdef DEBUG
    cerr << "Logging almost finished" << endl;
#endif
    
    if (cmdArgs["-log"] == "nolog") {
//		delete Dere;
        return;
    }
#ifdef DEBUG
    cerr << "DereLog start" << endl;
#endif
    if (mainFil->doDereplicate()) {
        string dereLog = logF.substr(0,logF.length()-3) + "dere";
        Dere->writeLog(dereLog, deLog);
//		delete Dere;
    }

#ifdef DEBUG
    cerr << "DereLog end" << endl;
#endif
    if (shortStats) {
        mainFil->printStats(std::cerr, mainFile, outFile, true);
    }
#ifdef DEBUG
    cerr << "other logs start" << endl;
#endif
    //per sample success rate
    string logPS = logF.substr(0, logF.length() - 3) + "acceptsPerSample.log";
    log.open(logPS.c_str(), ios_base::out);
    mainFil->SmplSpecStats(log);
    log.close();
    log.open(logF.c_str(), ios_base::out);
    mainFil->printStats(log,mainFile,outFile,true);
    log.close();
    
    string logFs = logF.substr(0, logF.length() - 3) + "acceptsPerFile.log";
    log.open(logFs.c_str(), ios_base::out);
    log << shrtLog;
    log.close();
    
    string logFGC = logF.substr(0, logF.length() - 3) + "GC.txt";
    log.open(logFGC.c_str(), ios_base::out);
    mainFil->printGC(log, mainFil->isPaired());
    log.close();
#ifdef DEBUG
    cerr << "other logs end" << endl;
#endif


//for additional files
    if (mainFil->secondaryOutput()){
        log.open (logFA.c_str() ,ios_base::out);
        mainFil->printStats(log,mainFile,outFile,false);
        log.close();
    }
    
    
    //length histogram
    logF = cmdArgs["-length_hist"];
    log.open (logF.c_str() ,ios_base::out);
    mainFil->printHisto(log,0);
    log.close();
    //quality histogram
    logF = cmdArgs["-qual_hist"];
    log.open (logF.c_str() ,ios_base::out);
    mainFil->printHisto(log,1);
    log.close();
    mainFil->closeOutFilesDemulti();
#ifdef DEBUG
    cerr << "separateByFile finished" << endl;
#endif

}*/

//bool readCmdArgs(int argc, char* argv[],base_map<char*, char*, lstr>& cmdArgs);
