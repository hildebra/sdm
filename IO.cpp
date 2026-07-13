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
#include "ThreadPool.h"
#include "IO.h"
#include "OutputStreamer.h"
#include "Filters.h"
#include "Benchmark.h"
#include <future>
#include <utility>


void multi_read_paired_STRget(shared_ptr<InputStreamer> IS, OutputStreamer* MD, int curThread, 
	int tmpBlockSize, bool MIDuse, bool keepPairHd, qual_score FastqVer) {


    std::unique_ptr<multi_tmp_lines> tmpLines = std::make_unique<multi_tmp_lines>(tmpBlockSize);
	bool cont(true);
	while (cont) {
        {
			Instr::ScopedTimer st(Instr::GETDNALINES);
			cont = IS->getDNAlines(tmpLines.get(), tmpBlockSize, MIDuse, true);
		}
		multi_read_paired_STRready(tmpLines.get(), MIDuse, MD, curThread, keepPairHd, FastqVer);
	}

	// unique_ptr will clean up tmpLines

}


bool multi_read_paired_STRready(multi_tmp_lines* tmpLines,
	bool MIDuse, OutputStreamer* MD, int curThread,
	bool keepPairHd, qual_score FastqVer) {

	//for loop over tmplines, each entry is a read pair
	Instr::ScopedTimer st_multi(Instr::MULTI_READ_READY);
	// measure async lag if submit_ts set
	if (tmpLines->submit_ts != 0) {
		long long lag = Instr::now_us() - tmpLines->submit_ts;
		Instr::add_us(Instr::ASYNC_LAG, lag);
		tmpLines->submit_ts = 0;
	}
	bool isOK(true);
	vector<shared_ptr<DNA>> ret(3, nullptr);

	for (size_t k = 0; k < tmpLines->size(); k++) {
		ret[0].reset();
		ret[1].reset();
		ret[2].reset();
		//cdbg("rpS");
            {
				Instr::ScopedTimer st(Instr::STR2DNA);
				ret[0] = str2DNA(std::move(tmpLines->tmp[k][0]), keepPairHd, FastqVer, 0);
			}
            {
				Instr::ScopedTimer st(Instr::STR2DNA);
				ret[1] = str2DNA(std::move(tmpLines->tmp[k][1]), keepPairHd, FastqVer, 1);
			}
		if (MIDuse) {
			ret[2].reset();
                {
					Instr::ScopedTimer st(Instr::STR2DNA);
					ret[2] = str2DNA(std::move(tmpLines->tmp[k][2]), keepPairHd, FastqVer, 2);
				}
		}
		if (ret[0] == nullptr) {
			continue;
			//let's just make sure that there are not trailing objects..
		}

		//GoldenAxe requires some extra rounds of checking DNA
		if (MD->getFilters(curThread)->isGoldAxe()) {
			vector<shared_ptr<DNA>> multDNA = MD->getFilters(curThread)->GoldenAxe(ret);
			int iniBCsGA(-2);
			if (multDNA.size()) {
				iniBCsGA = multDNA[0]->getBCnumber();
			}
			vector<shared_ptr<DNA>> ret2(3, nullptr);
			for (size_t xi = 0; xi < multDNA.size(); xi++) {
				ret2[1].reset();
				ret2[2].reset();
				ret2[0] = multDNA[xi];
				isOK=process_DNA(ret2, MIDuse, MD, curThread);
				if (ret2[0]->getBCnumber() != iniBCsGA) {//DEBUG
					int BCnewWrong = ret2[0]->getBCnumber();
					cerr << "Problem in GA read \n" << tmpLines->tmp[k][0][0]<<"\n"<< tmpLines->tmp[k][0][1] << "\n:Different BCs detected\n";
				}
			}
		} else {
			isOK = process_DNA(ret, MIDuse, MD, curThread);
		}
		//if (!isOK) {break;}
	}
	//delete tmpLines;
	//no longer needed..
	return isOK;
}

//is called from a while loop, that reads the DNA pairs
bool process_DNA(vector< shared_ptr<DNA>>& tdn,
	bool MIDuse, OutputStreamer* MD, int curThread) {

	if (tdn[0] == nullptr) {
	    return false;
	} //|| tdn->length()==0
    Instr::ScopedTimer st_proc(Instr::PROCESS_DNA);
	//DNA objects as they should be??
	if (MIDuse && tdn[2] == nullptr) {
		//cerr << "Missing MID read pair.\n";		exit(4);
		throw std::runtime_error("Missing MID read pair.");

	}
	if (tdn[1] == nullptr && tdn[0] != nullptr && MD->isPEseq() == 2) {
		//cerr << "Second provided file has not the same number of entries as first file.\n";exit(5);
		throw std::runtime_error("Second provided file has not the same number of entries as first file.");
	}
	//MD->checkFastqHeadVersion(tdn[0]);
	//testreadpair.lock();

	//cdbg("process_DNA"+ tdn[0]->getId()+"\n");

	Filters* curFil = MD->getFilters(curThread);
	//register read at all with stat counter:
	 curFil->sTotalPlus(0);
	if (tdn[1] != nullptr || curFil->isPaired() == 2) {
		curFil->sTotalPlus(1);
	}
	//collect some info on general run parameters
	curFil->preFilterSeqStat(tdn[0], 0);
	curFil->preFilterSeqStat(tdn[1], 1);



	//prep some variables
	int BCoffs = curFil->getBCoffset();
	bool checkSwitchedRdPairs = curFil->checkSwitchedRdPairs();
	bool dualBCs = curFil->doubleBarcodes();
	bool doBCsAtAll = curFil->doBarcodes();
	int pairedRd = curFil->isPaired();
	//bool checkReversedRead = curFil->checkRevRd();

	int tagIdx(-2); int tagIdx2(-2);
	string presentBC(""); int c_err(0);
	//bool isReversed(false);//was a reversion detected?

	if (MIDuse && tdn[2] != nullptr) {
		tagIdx = curFil->findTag(tdn[2], presentBC, c_err, true,0,true,false); 
		tdn[0]->setBCnumber(tagIdx, BCoffs);
	}
	//routine checks, and reverses/swaps DNA objects
	//bool reversedDNA = curFil->swapReverseDNApairs(tdn);


	//now split up into multiple PB reads, if GoldenAxe (or Kinnect)
	


	//manage reversed / swapped reads by detecting where the fwd primer is
	bool wasReversed;
	if (pairedRd == 2) {
		wasReversed= curFil->swapReverseDNApairs(tdn);
	} else if (tdn[1] == nullptr) {
		wasReversed = curFil->isReversedAmplicon(tdn[0]);
	}


    /*if (checkBC2ndRd) {
		if (!dualBCs) {
			//use this routine to set the BC inside the DNA object
			bool revT = false;
			bool Pr1 = curFil->findPrimer(tdn[0], 0, false, 0);
			bool Pr2 = curFil->findPrimer(tdn[1], 0, false, 0);
			//int chkRev1(-1), chkRev2(-1);
			tagIdx = curFil->findTag(tdn[0], presentBC, c_err, true, -1);
			tagIdx2 = curFil->findTag(tdn[1], presentBC, c_err, true, -1);
			if ((tagIdx2 >= 0 && tagIdx < 0 && !Pr1) || (Pr2 && !Pr1)) { //swap first & second read
				swap(tdn[0], tdn[1]);
				tdn[0]->constellationPairRev(true); tdn[1]->constellationPairRev(true);
				//revConstellation++;
			}
			//else if (tagIdx2 < 0 && tagIdx < 0) {
			//	int x = 0;
			//}
			if (revT) {
				tdn[0]->reverse_compliment(); tdn[1]->reverse_compliment();
			}
		}
		tagIdx2 = -2; tagIdx = -2;
		tdn[1]->setpairREV();		tdn[0]->setpairFWD();
	}
	*/

    
	//tdn[0]->reverse_compliment();
	//actually important routine
	MD->analyzeDNA(tdn[0], -1, 0, tagIdx, curThread);
	//tdn[0]->matchSeqRev
	bool ch1(false); if (tdn[0] != NULL) { ch1 = tdn[0]->isGreenQual(); }
	bool ch2(false); bool ch2n(false);



	//if ( ch1 ) {	cerr << cnt << " \n";	}
	//normal case for check 2nd read
	if (!ch2 && tdn[1] != nullptr) { //ch1&&
								//dnaTemp2->setBCnumber(tdn[0]->getBCnumber());
		if (doBCsAtAll && !dualBCs) { //only check in read1 for BC, if not dual BCing!!
			tagIdx2 = tdn[0]->getBCnumber();  // no 2nd BC, thus no BC search in 2nd read
			if (tagIdx2 >= 0) {
				tagIdx2 -= BCoffs;
			}else if (tagIdx2 < -1) {//something wrong with BCoffs
				cerr << "tagidx2 wrongly truncated to " << tagIdx2 << endl;
			}
		}
		//if (isReversed) { tdn[1]->reverse_compliment(); }
		MD->analyzeDNA(tdn[1], -1, 1, tagIdx2, curThread);
		ch2 = tdn[1]->isGreenQual();
	}

	//set up BC in DNA header
	//remember that dual BCs are only valid after this step!
	if (dualBCs && pairedRd == 2) {
		//in case of single reads (PB) with double BCs
		//if (dualBCs && tagIdx2 < 0 && pairedRd == 1) {
		//	tagIdx2 = curFil->findTag2(tdn[0], presentBC, c_err, false, -1);
		//}

		curFil->dblBCeval(tagIdx, tagIdx2, presentBC, tdn[0], tdn[1]);
		c_err = -1;

		//check a second time that barcode was correctly identified, just to be double sure...
		if ( tagIdx != tagIdx2 || tdn[0]->getBCnumber() != tdn[1]->getBCnumber()) {
			//cerr << "Unequal BC numbers:" << tagIdx << " : " << tagIdx2 << "; in object: " << tdn[0]->getBCnumber() << " : " << tdn[1]->getBCnumber() << endl;
			//cerr << "In read:" << tdn[0]->getId() << endl; exit(835);
			throw std::runtime_error("Unequal BC numbers:" + itos(tagIdx) + " : " + itos(tagIdx2) + "; in object: " + itos(tdn[0]->getBCnumber()) + " : " + itos(tdn[1]->getBCnumber()) + "\nIn read:" + tdn[0]->getId());			
		}
	}
	

	
	if (tagIdx == -1 || tagIdx2 == -1) {
		tdn[0]->setBarcodeDetected(false);
		if (tdn[1] != nullptr) {
			tdn[1]->setBarcodeDetected(false);
		}
	}	else  {
		curFil->setBCdna(tagIdx, tdn[0]);
		if (ch2) { curFil->setBCdna(tagIdx, tdn[1]); }
		if (MIDuse && ch1) {
			curFil->BCintoHead(tagIdx, tdn[0], presentBC, c_err, true);
			if (ch2) { curFil->BCintoHead(tagIdx, tdn[1], presentBC, c_err, true); }
		}
	}



	int idx1 = 1; int idx2 = 2;
	bool read2notNull = tdn[1] != nullptr;

	if (MD->isPEseq() == 2 ) {
		if (read2notNull && tdn[0] != nullptr &&
			(tdn[0]->isGreenYellowQual() && tdn[1]->isGreenYellowQual()) ) {
			//only this case allows for paired reads out
			idx1 = 1;  idx2 = 2;
		}
		else {
			idx1 = 3; idx2 = 4;
		}
		/*if (ch1 && !ch2) {
			idx1 = 3; idx2 = 4;
			//if (tdn[1] != NULL) { tdn[1]->failed(); }
			//		delete dnaTemp2;
		}
		else if (ch2 && !ch1) {
			idx2 = 4; idx1 = 3;
			//if (tdn[0] != NULL) { tdn[0]->failed(); }
			//		delete tdn[0];
		}
		/*else if (!ch1 && !ch2) { //nothing passes
			if (tdn[0] != NULL) { tdn[0]->failed(); }
			if (tdn[1] != NULL) { tdn[1]->failed(); }
			//		delete tdn[0]; delete dnaTemp2;
		}*/
	}
	



	//pre-merge step
	//if ( MD->mergeReads()) {MD->findSeedForMerge(tdn[0], tdn[1],curThread);	}
	//needed to eval qual of seed
	if (read2notNull) {//check that this is really correct read pair
		if (!tdn[0]->sameHeadPosFree(tdn[1])) {
			//cerr << "Error: invalid paired read pair headers:\n  read1 head: " << tdn[0]->getShortId() << "\n  read2 head: " << tdn[1]->getShortId()<<endl;exit(152);
			throw std::runtime_error("Error: invalid paired read pair headers:\n  read1 head: " + tdn[0]->getShortId() + "\n  read2 head: " + tdn[1]->getShortId());

		}
	}
	if (MD->isPEseq() == 2 && ( MD->mergeReads() || MD->doDerepMrgSrch()) && read2notNull) {
		MD->findSeedForMerge(tdn[0], tdn[1], curThread);
	}


	//demultiplex write? do this first before DNA is deleted..
	//at this point the tagIDX *MUST* be correctly set + BCoffset (in the DNA object, tagIDX doesn;t matter)
	//cdbg("w2D ");
	if (MD->Demulti2Fls()) {
		MD->write2Demulti(tdn[0], tdn[1], curFil->getBCoffset(), curThread);
	}
	if (MD->doWriteNonBCrds()) {
		MD->writeNonBCReads(tdn[0], tdn[1]);
	}
	if (MD->doDeriplicate()) {
		MD->dereplicateDNA(tdn[0], tdn[1]);
	}

    // Test
    /*if (OutStreamer->b_merge_pairs_derep_ && tdn[0]->merge_seed_pos_ > 0) {
        auto dna_merged = ReadMerger::merge(tdn[0], tdn[1]);
    }*/
	//testreadpair.unlock();
	//save for later .. and collect stats
	
	//writing of pairs needs to be mutex locked 
	int Cstr1(100), Cstr2(100);
	//cdbg("sFW\n");
    bool owr1;
	bool owr2;
	{
		Instr::ScopedTimer st_out(Instr::OUTPUT_WRITE);
		owr1 = MD->saveForWrite(tdn[0], idx1, curThread, Cstr1,false);
		owr2 = MD->saveForWrite(tdn[1], idx2, curThread, Cstr2,false);
		//mutex version to ensure read pairs together..
		MD->writeForWrite(tdn[0],idx1,Cstr1,tdn[1],idx2,Cstr2);
	}

	if (!owr1 || (!owr2 && read2notNull)) {
		return false;
	}
	return true;
}




struct job2 {
	bool inUse = false;
	future<bool> job;
	//thread job;
};

bool read_sequences(OutputStreamer* MD, 
	shared_ptr<InputStreamer> IS, int Nthreads) {
	OptContainer* cmdArgs = getGlobalOptContainer();
	if (cmdArgs == nullptr) {
		throw std::runtime_error("Global option container is not initialized.");
	}
	const OptContainer::RuntimeOptions& runtimeOpts = cmdArgs->runtime();
	if (Nthreads < 1) {
		Nthreads = 1;
	}
	
	DNAmap oldMIDs;
	bool MIDuse = IS->hasMIDseqs();
	bool fqHeadVer(true);
	cdbg( "Read paired routine\n");

	/*if (sync2pair && MIDuse) {
		cout << "Can not sync read pairs, while explicit MID sequences are being used! (not supported, sorry)\n";
		sync2pair = false;
	} */

	//bool syncedMID = false;
	
	//multithreading setup
	//int Nthrds = atoi((*cmdArgs)["-threads"].c_str()) -1 ;
	vector<job2> slots(Nthreads);
	int thrCnt = 0;
	bool doMC(false);
	if (Nthreads > 1) {
		doMC = true;
	}

	int DNAinMem(0);
	bool cont(true),cont2(true),cont3(true);
	bool keepPairedHD = IS->keepPairedHD();
	//int revConstellation(0);
	//int cnt(0); 
	//bool switching(true); // important to keep track of this, to fix swapped read pairs

	//vector<string>tmpLines2(4, "");
	//vector<vector<string>> tmpLines(3, tmpLines2);


 //ini blocks of empty strings for later reuse..
	int tmpBlockSize = runtimeOpts.iniBlockSize;

   qual_score fastqVer = 33;

	// Probe one block first and adapt batch size for short reads to reduce task scheduling overhead.
	// This keeps behavior unchanged for longer reads while increasing reads/task when records are tiny.
	int taskBatchMultiplier = 1;
	string empty(""); empty.reserve(151);
	{
		cdbg("task batch size probe");
		multi_tmp_lines probe(tmpBlockSize);
		bool probeCont = IS->getDNAlines(&probe, tmpBlockSize, MIDuse, true);
		if (probeCont && probe.size() > 0 && probe.tmp[0].size() > 0 && probe.tmp[0][0].size() > 1) {
			const size_t probeReadLen = probe.tmp[0][0][1].size();
			if (probeReadLen > 0) {
				if (probeReadLen <= 60) {
					taskBatchMultiplier = 32;
				}
				else if (probeReadLen <= 110) {
					taskBatchMultiplier = 16;
				}
				else if (probeReadLen <= 160) {
					taskBatchMultiplier = 8;
				}
				else if (probeReadLen <= 350) {
					taskBatchMultiplier = 4;
				}
			}

			empty.reserve(probeReadLen + 1);

			fastqVer = IS->fastQscore();
			if (fqHeadVer && probe.tmp[0][0].size() > 0) {
				shared_ptr<DNA> tmp = str2DNA(probe.tmp[0][0], keepPairedHD, fastqVer, 0);
				MD->checkFastqHeadVersion(tmp);
				fqHeadVer = false;
			}
			cdbg("FASTQversion: " + to_string(fastqVer) + ". ");
		}
		IS->allStreamReset();
		cdbg( " is " + to_string(taskBatchMultiplier) + ". ");

	}
	if (taskBatchMultiplier > 1 && tmpBlockSize <= (INT_MAX / taskBatchMultiplier)) {
		tmpBlockSize *= taskBatchMultiplier;
	}
	cdbg("read_seq check3 ");
	vector<multi_tmp_lines*> tmpStrHolders(Nthreads, nullptr);	
	vector<string>tmpLines2(4, empty);
	vector<vector<string>> tmpLines(3, tmpLines2);
	for (size_t ii = 0; ii < Nthreads; ii++) {
		multi_tmp_lines* tmpO = DBG_NEW multi_tmp_lines(); // will be deleted inside multi_read_paired_STRready function
		tmpO->tmp.resize(tmpBlockSize, tmpLines);
		tmpStrHolders[ii] = tmpO;
	}

    //object holding converted DNA.. not longer used
	//vector<shared_ptr<DNA>> tdn(3, nullptr);
	cdbg("read_seq check4 ");
	//doMC = false;//DEBUG
	while ( cont ) {
		//decide on MC or single core submission route
		if (doMC) { // multithreading
			// step 0: find first available slot (free or already completed)
			size_t slotIdx = static_cast<size_t>(thrCnt);
			bool hasFreeSlot = false;
			for (size_t probe = 0; probe < slots.size(); ++probe) {
				size_t idx = (slotIdx + probe) % slots.size();
				if (!slots[idx].inUse) {
					slotIdx = idx;
					hasFreeSlot = true;
					break;
				}
				if (slots[idx].job.valid() &&
					slots[idx].job.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
					bool doneOk = slots[idx].job.get();
					(void)doneOk;
					slots[idx].inUse = false;
					slotIdx = idx;
					hasFreeSlot = true;
					break;
				}
			}
			if (!hasFreeSlot) {
				if (slots[slotIdx].job.valid()) {
					bool doneOk = slots[slotIdx].job.get();
					(void)doneOk;
				}
				slots[slotIdx].inUse = false;
				hasFreeSlot = true;
			}
			thrCnt = static_cast<int>(slotIdx);

            //step 1: read block from input into tmp buffer
			{
				Instr::ScopedTimer st(Instr::GETDNALINES);
				cont = IS->getDNAlines(tmpStrHolders[thrCnt], tmpBlockSize, MIDuse);
			}

			// If we read any records into the tmp buffer, submit them for processing.
			// If cont==false but buffer contains data (partial final block), still process it.
			if (tmpStrHolders[thrCnt]->size() > 0) {
                if (slots[thrCnt].inUse == false) {
					// record submit timestamp and measure async creation cost
					long long submit_ts = Instr::now_us();
					tmpStrHolders[thrCnt]->submit_ts = submit_ts;
					{
						Instr::ScopedTimer st_async(Instr::ASYNC_SUBMIT);
						slots[thrCnt].job = ThreadPool::instance().submit(multi_read_paired_STRready,
							tmpStrHolders[thrCnt], MIDuse, MD, thrCnt, keepPairedHD, fastqVer);
						}
					slots[thrCnt].inUse = true;
				}
			}

			// If no data was read and cont is false, EOF reached without data -> stop.
			if (!cont && tmpStrHolders[thrCnt]->size() == 0) { break; }
			thrCnt++;
			if (thrCnt >= Nthreads) { thrCnt = 0; }
           if (!cont) { break; }
		} else {//single thread
           cont = IS->getDNAlines(tmpStrHolders[thrCnt], tmpBlockSize, MIDuse);
			// If we read any records into the tmp buffer, process them even if cont==false (partial final block).
			if (tmpStrHolders[thrCnt]->size() > 0) {
				//fastqVer = IS->fastQscore();
				bool cont2 = multi_read_paired_STRready(tmpStrHolders[thrCnt], MIDuse, 
					MD, 0, keepPairedHD, fastqVer);
			}
			if (!cont && tmpStrHolders[thrCnt]->size() == 0) { break; }
				//if (tdn[0]->isConstellationPairRev()) { revConstellation++; }
			//}
		}
	}
	cdbg("read_seq check5 ");

	//last jobs 


 // get all slots (single wrap-around pass starting at thrCnt)
	if (!slots.empty()) {
		for (size_t i = 0; i < slots.size(); ++i) {
			size_t x = (static_cast<size_t>(thrCnt) + i) % slots.size();
			if (slots[x].inUse == true) {
				if (slots[x].job.valid()) {
					cont = slots[x].job.get();
				}
				slots[x].inUse = false;
			}
		}
	}

	
	//close shop
	//MD->revConstellationCnts(revConstellation);
	MD->closeOutStreams();

	for (size_t ii = 0; ii < Nthreads; ii++) {
		delete tmpStrHolders[ii];
	}
	cdbg("read_seq check6 ");

	// Print instrumentation results
	Instr::printResults(std::cerr);

	return true;
}


bool readCmdArgs(int argc, char* argv[]){
	OptContainer* cmdArgs = getGlobalOptContainer();
	if (cmdArgs == nullptr) {
		throw std::runtime_error("Global option container is not initialized.");
	}
	if (argc%2!=1){
//		cerr<<"It seems command line arguments were not passed in pairs. Aborting.\n";exit(666);
		throw std::runtime_error("It seems command line arguments were not passed in pairs. Aborting.");
	}
	auto setDefault = [&](const string& key, const string& value) {
		cmdArgs->setDefault(key, value);
	};
	for (int i=1; i<argc; i+=2){ //parsing of cmdline args
		string theNxtSt = string(argv[i+1]);
		if (theNxtSt.empty() || theNxtSt[0] != '-'){
			cmdArgs->set(string(argv[i]), theNxtSt);
		} else {
			cmdArgs->set(string(argv[i]), "T");
		}
	}
	
	
	setDefault("-illuminaClip", "0");
	setDefault("-logLvsQ", "");
	setDefault("-GoldenAxe", "0");
	setDefault("-GoldenAxeMinAmpli", "0");
	setDefault("-GoldenAxeMaxAmpli", "0");
	setDefault("-derepSrchLen", "150");
	setDefault("-DNAseedSelectionClassic", "0");
	setDefault("-i_MID_fastq", "");
	//set to default (empty)
	setDefault("-OTU_fallback", "");
	setDefault("-otu_matrix", "");
	setDefault("-derepPerSR", "0");
	setDefault("-ucAdditionalCounts", "");
	setDefault("-ucAdditionalCounts1", "");
	setDefault("-ucAdditionalCounts_refclust", "");
	setDefault("-optimalRead2Cluster_ref", "");
	//just for debuggin purposes: write out all seqs, where no BC can be detected..
	setDefault("-o_fastq_noBC", "");
	setDefault("-ucAdditionalCounts_refclust1", "");
	
	setDefault("-XfirstReadsWritten", "");
	setDefault("-XfirstReadsRead", "");
	setDefault("-iniBlockSize", "120");
	//filter sequence file for a specific subset of sequences
	//these arguments can only occur together
	if (cmdArgs->find("-specificReads") == cmdArgs->end()) {
		cmdArgs->set("-specificReads", "");
	} else if (cmdArgs->find("-excludeFile") == cmdArgs->end()) {
		cmdArgs->set("-excludeFile", "");
	}
	setDefault("-onlyPair", "");
	setDefault("-onlyCanonicalNTs", "0");
	setDefault("-strictPositiveQualsOut", "0");
	setDefault("-pairedDemulti", "0"); // by default: do not only report proper pairs
	
	setDefault("-uparseVer", "");


	if (cmdArgs->find("-i_path")  == cmdArgs->end()){ //ok files are not given in mapping file
		//check if dna and qual_ are passed
		if (cmdArgs->find("-i") != cmdArgs->end()) {
			string fmt = detectSeqFmt((*cmdArgs)["-i"]);
			if (fmt == "empty") {
				//exit(0);
				cerr << "Only empty input files\n";
			}
			(*cmdArgs)[fmt] = (*cmdArgs)["-i"];
		}
		if (cmdArgs->find("-i_fastq")  == cmdArgs->end()){ // fasta + quality format
			if (cmdArgs->find("-i_fna")  == cmdArgs->end()){
				//cerr<<"Fasta not supplied. \nPlease give the path to your fasta file as command line argument:\n  -i_fna <yourFastaFile>\n";exit(2);
				throw std::runtime_error("Fasta not supplied. \nPlease give the path to your fasta file as command line argument:\n  -i_fna <yourFastaFile>");
			}
			if (cmdArgs->find("-i_qual")  == cmdArgs->end()){
				string newQ = (*cmdArgs)["-i_fna"];
				int pos = (int)newQ.find_last_of(".");
				newQ = newQ.substr(0,pos);
				newQ += string(".qual_");
				fstream fin;
				fin.open(newQ.c_str(),ios::in);
				if( fin.is_open() )	{
					cerr<<"Using quality file: "<<newQ <<endl;
				} else if ((cmdArgs->find("-number")!=  cmdArgs->end() && (*cmdArgs)["-number"] =="T")||
					((*cmdArgs)["-specificReads"] != "")) {
					(*cmdArgs)["-i_qual"] = "";
				} else {
					cerr<<"You did not supply a quality file. \nPlease give the path to your quality file as command line argument:\n  -i_qual <PathToQualityFile>\n";
					newQ = "";
					//fin.close();	exit(2);
				}
				fin.close();
				(*cmdArgs)["-i_qual"] = newQ;
			}
		}
		//auto create output file name
		if (cmdArgs->find("-o_fna")  == cmdArgs->end()){
			if (cmdArgs->find("-o_fastq")  == cmdArgs->end()){
				//(*cmdArgs)["-o_fna"] = (*cmdArgs)["-i_fna"]+string(".sdm");
				//cerr<<"Writing output fasta into "<<(*cmdArgs)["-o_fna"]<<endl;
				cerr << "No output file will be written\n";
			}
		} else {
			if (cmdArgs->find("-o_fastq")  != cmdArgs->end()){
				cerr<<"\"-o_fna\" was over-writen by \"-o_fastq\"\n";
				(*cmdArgs)["-o_fna"] = "";
			}
		}
	} else {
		if (cmdArgs->find("-o_fna")  == cmdArgs->end() && cmdArgs->find("-o_fastq")  == cmdArgs->end()){
			//cerr<<"Please give an output file (\"-o_fna\" || \"-o_fastq\") if you use sdm \"-i_path\" option.\n  Aborting..\n";exit(2);
			throw std::runtime_error("Please give an output file (\"-o_fna\" || \"-o_fastq\") if you use sdm \"-i_path\" option.\n  Aborting..");
		}
	}

	/*	if (cmdArgs->find("-map")  == cmdArgs->end()){
	cerr<<"You did not supply a mapping file. \nPlease give the path to your mapping file as command line argument:\n  -map <PathToMappingFile>\n";
	exit(2);
	}  */
	if (cmdArgs->find("-o_qual")  == cmdArgs->end()){
		cmdArgs->set("-o_qual", "");
	} else {
		if (cmdArgs->find("-o_fastq")  != cmdArgs->end()){
			cerr<<"\"-o_qual\" was over-writen by \"-o_fastq\"\n";
			(*cmdArgs)["-o_qual"] = "";
		}
	}
	setDefault("-options", string("sdm_options.txt"));
	setDefault("-threads", "1");
	setDefault("-threadIO", "1");
	if (cmdArgs->find("-log")  == cmdArgs->end()){
		string ofile1 = (*cmdArgs)["-o_fna"];
		if (ofile1==""){ofile1 = (*cmdArgs)["-o_fastq"];}
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
		(*cmdArgs)["-log"] = ofile1 + string(".log");
	}
	string ofile1 = (*cmdArgs)["-log"];
	//ofile1.find_last_of(".log");
	size_t logPos = ofile1.find_last_of(".");
	if (logPos != std::string::npos){
		ofile1 = ofile1.substr(0,logPos);
	}
	if (cmdArgs->find("-length_hist")  == cmdArgs->end()){
		cmdArgs->set("-length_hist", ofile1 + string("_lenHist.txt"));
	}
	if (cmdArgs->find("-qual_hist") == cmdArgs->end()) {
		cmdArgs->set("-qual_hist", ofile1 + string("_qualHist.txt"));
	}
	if (cmdArgs->find("-merg_readpos") == cmdArgs->end()) {
		cmdArgs->set("-merg_readpos", ofile1 + string("_mergRpos.txt"));
	}
	if (cmdArgs->find("-	qual_readpos") == cmdArgs->end()) {
		cmdArgs->set("-qual_readpos", ofile1 + string("_qualRpos.txt"));
	}

		//-length_hist   -qual_hist

	if (cmdArgs->find("-sample_sep")  == cmdArgs->end()){
		(*cmdArgs)["-sample_sep"] = DEFAULT_BarcodeNameSep;
	} else 	if ((*cmdArgs)["-sample_sep"]==""){
		//cerr<<"Invalid sample separator (empty).\nAborting..\n";exit(82);
		throw std::runtime_error("Invalid sample separator (empty).\nAborting..");
	}


	setDefault("-o_qual_offset", DEFAULT_output_qual_offset);
	setDefault("-pairedRD_HD_out", DEFAULT_pairedRD_HD_out);
	setDefault("-5PR1cut", DEFAULT_5PR1cut);
	setDefault("-5PR2cut", DEFAULT_5PR2cut);

	

	if (cmdArgs->find("-ignore_IO_errors") == cmdArgs->end()) {
		cmdArgs->set("-ignore_IO_errors", DEFAULT_ignore_IO_errors);
	} else if ((*cmdArgs)["-ignore_IO_errors"] != "0" && (*cmdArgs)["-ignore_IO_errors"] != "1") {
		//cerr << "Argument \"ignore_IO_errors\" can only be \"1\" or \"0\". Instead it has value: " << (*cmdArgs)["-ignore_IO_errors"] << endl;exit(323);
		throw std::runtime_error("Argument \"ignore_IO_errors\" can only be \"1\" or \"0\". Instead it has value: " + (*cmdArgs)["-ignore_IO_errors"]);
	}
	setDefault("-o_dereplicate", "");
	setDefault("-derep_map", "");

	

	//if (cmdArgs->count("-i_fna")==0){}

	setDNAseedSelectionClassic((*cmdArgs)["-DNAseedSelectionClassic"] == "1");

	return true;
}




//manages read in of several input files and associated primers / tags to each file
void separateByFile(Filters* mainFilter, Benchmark* sdm_benchmark){
	OptContainer* cmdArgs = getGlobalOptContainer();
	if (cmdArgs == nullptr) {
		throw std::runtime_error("Global option container is not initialized.");
	}
#ifdef DEBUG
	cerr << "separateByFile"<<endl;
#endif

//	mainFilter->ini_filestruct(cmdArgs);
	//setup once at start
	//vector<ReadMerger*> merger(0);
	//contains info on input files, format, order etc
	filesStr files;

	//handles testing for input files, resorting them, input file format..
	mainFilter->FileEssentials(files, cmdArgs);

    // mainFile for processing
	string mainFile = "";
	string outFile = (*cmdArgs)["-o_fna"];
	
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
	int maxReadsWr = mainFilter->getXreadsWr(); // should program stop after having written a certain amount of reads?
	int maxReadsRd = mainFilter->getXreadsRd();
	int totalReadsRead(0), totalReadsWrite(0);
	auto printConvertedSummary = [&]() {
		if ((*cmdArgs)["-o_fna"] != "") {
			cerr << "Converted sequences written: " << totalReadsWrite << "\n";
		}
	};
	uint accumBPwrite(0), accumBPwriteMerg(0);
	vector<string> lastSRblock (1,"");  //set up SequencingRun blocks to track
	string lastProcessedSRblock = "";

	const OptContainer::RuntimeOptions& runtimeOpts = cmdArgs->runtime();

    //---------------------------------
    // Multithreading setup
    //------------------------------------
    bool multithreading = true;
    //ThreadPool *pool = nullptr;
    int threads = 1;
    if (multithreading) {
		threads = runtimeOpts.threads;
		if (threads < 1) {
			threads = 1;
		}
		ThreadPool::configure(static_cast<size_t>(threads));
		/*if ( threads > 1) {
			pool = DBG_NEW ThreadPool(threads);
		}*/
		cerr << "Run with " << threads << " cores.\n";
    }
	if (maxReadsWr > 0) {
		cerr << "Only printing " << maxReadsWr << " first reads from file(s).\n";
	}
	if (maxReadsRd > 0) {
		cerr << "Only reading " << maxReadsRd << " first reads from file(s).\n";
	}

	//set up a read merger for each thread..
	bool detailedMergeStats(false);


    //---------------------------------
	//files.uniqueFastxFiles: unique input fastx's, can be demultiplexed fastas, 
	// or un-demultiplexed (several) fastx /loop that goes over several input fxs
	// Should be brought in right order (SRblock)
	//------------------------------------
//	for (unsigned int i = 0; i < files.uniqueFastxFiles.size(); i++ ) {
	for (auto &uFX : files.uniqFxFls) {
		uint i = uFX.second;
		cdbg("Unique file " + uFX.first + "(" + itos(i) + ")\n");

		// Sample memory at the start of processing each file
		Instr::g_memoryTracker.sample();

		if (maxReadsWr > 0 && maxReadsWr - totalReadsWrite <= 0) { cerr << "Skipping file " << uFX.first << "due to firstXreadsWrite .."; break; }
		if (maxReadsRd > 0 && maxReadsRd - totalReadsRead <= 0) { cerr << "Skipping file " << uFX.first << "due to firstXreadsRead .."; break; }
		if (files.idx[i].size() == 0) {	
			throw std::runtime_error("fastXtar vector for " + uFX.first + " is empty");
			//cerr << "fastXtar vector for " << uFX.first << " is empty" << endl;	exit(10);		
		}

		//create subset of BCs for the currently processed fastq's (only relevant BCs)
		cdbg("new filter in round " + itos(i) + "\n");
		Filters* filter = mainFilter->newFilterPerBCgroup(files.idx[i]);
		cdbg("Setting up threads\n");
		//filter->setThreads(threads);
		int tarID = files.idx[i][0];//just points to one file in uFX group..


		//initialize object to handle all input file combinations
		//main input of fastx, handles input IO
		cdbg("Ini InputStream");
		shared_ptr<InputStreamer> IS = make_shared<InputStreamer>(
			!files.isFastq, mainFilter->getuserReqFastqVer(), runtimeOpts.ignoreIOErrors,
			runtimeOpts.pairedRDHDOut, threads);
		IS->setTIO(runtimeOpts.threadIO);
		IS->setGlobalRdsRead(totalReadsRead); IS->setMaxRdsRead(maxReadsRd);


		// there is an entry in tar for each barcode for this file. If files.idx[i].size() is 1 there is only one barcode
		if (files.idx[i].size() == 1 && files.uniqueFastxFiles.size() > 1) {
			IS->atFileYofX(i + 1, (unsigned int)files.uniqueFastxFiles.size(), 1);
		}

		const string currentSRblock = mainFilter->SequencingRun[tarID];

		string mainFileShort;
		
		mainFile = IS->setupInput(files.path, tarID, uFX.first, files, filter->setPaired(),
			runtimeOpts.onlyPair, mainFileShort, false);
		if (mainFile.empty() || !IS->checkInFileStatus()) {
			cerr << "Warning: input file set is empty or unreadable; skipping " << uFX.first << "\n";
			delete filter;
			continue;
		}

		// check if derep block will change; only after input was confirmed readable
		if (totalReadsRead>0 && mainFilter->doDereplicate() && dereplicator->DerepPerSR()
			&& lastSRblock.back() != currentSRblock) {
			//sanity check
			for (uint ii = 0; ii < (lastSRblock.size()-1); ii++) {
				if (lastSRblock[ii] == lastSRblock.back()) {
					cerr << "\n\nWrong block order!!\n" << lastSRblock[ii] << ":" << lastSRblock.back() << endl;
				}
			}

			string deLogLocal = dereplicator->writeDereplDNA(mainFilter, lastSRblock.back());
			if (dereplicator->DerepPerSR()) {
				files.deLog += "Dereplication of SequencingRun " + lastSRblock.back() + ":\n";
			}
			files.deLog += deLogLocal;
			dereplicator->printMergeStats(subfile((*cmdArgs)["-merg_readpos"], lastSRblock.back()),
				subfile((*cmdArgs)["-qual_readpos"], lastSRblock.back()));

			dereplicator->reset();//and reset..
			accumBPwrite = 0; accumBPwriteMerg = 0;
			//OutStreamer->resetDemultiBPperSR();
		}
		if (currentSRblock != lastSRblock.back()) {
			lastSRblock.push_back(currentSRblock);
		}
		lastProcessedSRblock = currentSRblock;
		if (!IS->qualityPresent()) {
			filter->deactivateQualFilter();
			cerr << "\n*********\nWarning:: Quality file is not present.\nRecommended to abort demultiplexing.\n*********\n\n";
		}

		filter->BarcodePreStats();filter->checkDoubleBarcode();filter->checkDoubleSampleIDHead();


		if (mainFilter->doOptimalClusterSeq()) {
			ucl->findSeq2UCinstruction(IS, files.isFastq, mainFilter);
			delete filter;
			continue;// after this is only qual_ filter, not required from this point on
		}

		cdbg("Setting up output\n");



		//OutputStreamer OutStreamer = OutputStreamer(&filter, cmdArgs, writeStatus, RDSset);

		//OutputStreamer also contains subfilters for MC processing and logging of reads
		OutputStreamer* OutStreamer = DBG_NEW OutputStreamer(filter, cmdArgs, 
			writeStatus, RDSset, threads,"");
		OutStreamer->attachDereplicator(dereplicator);OutStreamer->activateReadMerger(threads);
		OutStreamer->setBPwrittenInSR(accumBPwrite);OutStreamer->setBPwrittenInSRmerg(accumBPwriteMerg);	
		OutStreamer->attachBenchmark(sdm_benchmark);filter->setMultiDNA(OutStreamer);
		if (maxReadsWr > 0) {
			OutStreamer->setReadLimit(maxReadsWr - totalReadsWrite);
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
			if (cmdArgs->find("-mocatFix") != cmdArgs->end()) {
				cerr << "MOCAT fix appplies\n";
				RDSset->findMatches(IS, OutStreamer, true);
			}	else {
				RDSset->findMatches(IS, OutStreamer, false);
			}
			//delete OutStreamer;
			delete filter;
			continue;
		}

		cdbg("Processing reads");

		// If multhreading is switched on and thread count is > 1  call multihtreaded functions.
		// if multithreading = false ignore thread option
		multithreading &= threads > 1;
		
		
		



		
		
		//**********************
		//heavy reading, demultiplexing, dereplicating routine
		//**********************
		read_sequences(OutStreamer, IS, threads);










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
		filter->addMergeStats(OutStreamer);
		filter->prepStats();
		if (IS->getCurFileN() == 0) {
			filter->printStats(cerr, mainFile, outFile, true);
		}
		else {
			cerr << filter->shortStats(""); shortStats = true;
		}

		//totalReadsRead += filter->getLocalAcceptRds();

		shrtLog += filter->shortStats(mainFileShort);
		totalReadsRead += filter->getLocalRdsRead();
		//		delete IS;
				//write log file
		if (files.uniqueFastxFiles.size() > 1) {//only print sub log if neccessary
			ofstream log;
			string logF2 = (*cmdArgs)["-log"] + string("0") + itos(i);
			log.open(logF2.c_str(), ios_base::out);
			filter->printStats(log, mainFile, outFile, true);
			log.close();
		}

        mainFilter->addStats(filter, files.idx[i]);
        // Print out merged read count (move to stort stats)

		//OutStreamer.reset();
		totalReadsWrite += OutStreamer->getReadsWritten();
		
		delete OutStreamer;
		delete filter;
    }//uFX
	cdbg("End of UFx loop");


    //-------------------------------------------------
    //reading of input fq's is done
    //postprocessing, write reads, dereplicates etc
    // ------------------------------------------------

//write log files
	cdbg("Prep final logging (derep etc)");
    if (files.uniqueFastxFiles.size() > 1){ mainFile = "several";}
    ofstream log; 
    //different logfile for SEED extension
    if (mainFilter->doOptimalClusterSeq()){
        //finish up dereplication file (creating pseudo seeds with counts)
        ucl->finishMAPfile();
        if ((*cmdArgs)["-ucAdditionalCounts"] != ""){
            ucl->set2UC();
            ucl->finishUCfile(mainFilter, (*cmdArgs)["-ucAdditionalCounts"], true);//with smplHead (.mid)
            ucl->finishUCfile(mainFilter, (*cmdArgs)["-ucAdditionalCounts1"], false);//without smplHead (.rest)
        }


        if ((*cmdArgs)["-ucAdditionalCounts_refclust"] != ""){
            //reference based clustering has some high qual_ seqs (no replacement with reads..)
            //takeOver even found high qual_ hits with these default seeds..
            shared_ptr<InputStreamer> FALL = make_shared<InputStreamer>(true,
				mainFilter->getuserReqFastqVer(), (*cmdArgs)["-ignore_IO_errors"], 
				(*cmdArgs)["-pairedRD_HD_out"],1);
            //this reads in the SLV fna's & creates matrix entries for these
            FALL->setupFna((*cmdArgs)["-OTU_fallback_refclust"]);
            ucl->setRefMode();
            ucl->addDefSeeds(FALL, mainFilter);
            ucl->set2UC();
            //mapping from ref OTU clustering
            ucl->finishUCfile(mainFilter, (*cmdArgs)["-optimalRead2Cluster_ref"], false);
            //mid / rest mappings
            ucl->finishUCfile(mainFilter, (*cmdArgs)["-ucAdditionalCounts_refclust"], true);//with smplHead (.rest)
            ucl->finishUCfile(mainFilter, (*cmdArgs)["-ucAdditionalCounts_refclust1"], false);//without smplHead (.rest)

        }
        if ((*cmdArgs)["-log"] != "nolog") {
            log.open(files.logF.c_str(), ios_base::out);
            ucl->printStats(cerr);
            ucl->printStats(log);
            log.close();
        }

        //everything done on DNA? Then write & delete
        if ((*cmdArgs)["-otu_matrix"] != "") {
            ucl->writeOTUmatrix((*cmdArgs)["-otu_matrix"]);
        }
        //polished OTU seeds need to be written after OTU matrix (renaming scheme)
        OutputStreamer* MDx = DBG_NEW OutputStreamer(mainFilter, cmdArgs, 
			ios::out, RDSset,0);
		//ReadMerger* DerepM = DBG_NEW ReadMerger(false);
		MDx->activateReadMerger(1);
        mainFilter->setMultiDNA(MDx);
        ucl->writeNewSeeds(MDx, mainFilter, false);
        
		//not used any longer
		/*
		//new fastas also need to be written..
        MDx.reset(DBG_NEW OutputStreamer(mainFilter, cmdArgs, 
			ios::app, RDSset,1,".ref"));//force fna output
        mainFilter->setMultiDNA(MDx );
        ucl->writeNewSeeds(MDx, mainFilter, true, true);
		*/
        delete MDx;
		//delete DerepM;
		printConvertedSummary();
        return;
    } else if (mainFilter->doDereplicate()) {
		//this is either the last time dereplicate is written (SRblocks),
		//or the main point at which to write the dereplicator for good
#ifdef DEBUG
        cerr << "write Dereplicated DNA" << endl;
#endif
		
		const string finalSRblock = (dereplicator->DerepPerSR() && !lastProcessedSRblock.empty()) ? lastProcessedSRblock : lastSRblock.back();
		string deLogLocal = dereplicator->writeDereplDNA(mainFilter, finalSRblock);
		if (dereplicator->DerepPerSR()) {
			files.deLog += "Dereplication of SamplingRun " + finalSRblock + ":\n";
		}
		files.deLog += deLogLocal;
		if ((*cmdArgs)["-log"] != "nolog") {
			dereplicator->writeLog(files.logF.substr(0, files.logF.length() - 3) + "dere" , files.deLog);
		}
		dereplicator->finishMap();

		//last time merger stats to write

		dereplicator->printMergeStats(subfile((*cmdArgs)["-merg_readpos"], finalSRblock),
			subfile((*cmdArgs)["-qual_readpos"], finalSRblock));
		/*if (merger[0] != nullptr) {
			ReadMerger* cMerg = DBG_NEW ReadMerger();
			for (size_t x = 0; x < merger.size(); x++) {
				cMerg->addRMstats(merger[x]);
			}
			cMerg->printMergeHisto(subfile((*cmdArgs)["-merg_readpos"], lastSRblock.back()));
			cMerg->printQualHisto(subfile((*cmdArgs)["-qual_readpos"], lastSRblock.back()));
			delete cMerg;
		}*/


#ifdef DEBUG
        cerr << "done write dereplicator" << endl;
#endif
    }
#ifdef DEBUG
    cerr << "Logging almost finished" << endl;
#endif

    if ((*cmdArgs)["-log"] == "nolog") {
		printConvertedSummary();
        return;
    }
    if (shortStats) {
        mainFilter->printStats(std::cerr, mainFile, outFile, true);
    }

#ifdef DEBUG
    cerr << "other logs start" << endl;
#endif
    //per sample success rate
    string logPS = files.logF.substr(0, files.logF.length() - 3) + "acceptsPerSample.log";
    log.open(logPS.c_str(), ios_base::out);
    mainFilter->SmplSpecStats(log);
    log.close();
    log.open(files.logF.c_str(), ios_base::out);
    mainFilter->printStats(log, mainFile, outFile, true);
    log.close();
	string primerVariantReport = files.logF.substr(0, files.logF.length() - 3) + "primerVars.txt";
#ifdef _gzipread
	primerVariantReport += ".gz";
#endif
	mainFilter->writePrimerVariantReport(primerVariantReport);

    string logFs = files.logF.substr(0, files.logF.length() - 3) + "acceptsPerFile.log";
    log.open(logFs.c_str(), ios_base::out);
    log << shrtLog;
    log.close();

    string logFGC = files.logF.substr(0, files.logF.length() - 3) + "GC.txt";
    log.open(logFGC.c_str(), ios_base::out);
    mainFilter->printGC(log, mainFilter->isPaired());
    log.close();


	if (mainFilter->isGoldAxe()) {
		string logFGA = files.logF.substr(0, files.logF.length() - 3) + "GoldAxe.txt";
		log.open(logFGA.c_str(), ios_base::out);
		mainFilter->GAstatistics->printBCtabs(log);
		log.close();
		logFGA = files.logF.substr(0, files.logF.length() - 3) + "GoldAxeAmpliN.txt";
		log.open(logFGA.c_str(), ios_base::out);
		mainFilter->GAstatistics->printBCAmpliNdistribution(log);
		log.close();
		logFGA = files.logF.substr(0, files.logF.length() - 3) + "GoldAxeAmpliL.txt";
		log.open(logFGA.c_str(), ios_base::out);
		mainFilter->GAstatistics->printBCAmpliLdistribution(log);
		log.close();
	}

#ifdef DEBUG
    cerr << "other logs end" << endl;
#endif


//for additional files
    if (mainFilter->secondaryOutput()){
        log.open (files.logFA.c_str() ,ios_base::out);
        mainFilter->printStats(log, mainFile, outFile, false);
        log.close();
    }


    //length histogram
    string logF2 = (*cmdArgs)["-length_hist"];
    log.open (logF2.c_str() ,ios_base::out);
    mainFilter->printHisto(log, 0);
    log.close();
    //quality histogram
    logF2 = (*cmdArgs)["-qual_hist"];
    log.open (logF2.c_str() ,ios_base::out);
    mainFilter->printHisto(log, 1);
    log.close();

	//length vs qual log file (large!!)
	logF2 = (*cmdArgs)["-logLvsQ"];
	if (logF2 != ""){
		if (logF2 == "1") {
			logF2 = files.logF.substr(0, files.logF.length() - 3) + "Len.Qual.txt";
		}
	log.open(logF2.c_str(), ios_base::out);
	mainFilter->printLenVsQual(log);
	log.close();
	}

#ifdef DEBUG
    cerr << "separateByFile finished" << endl;
#endif

    // multithreading
    //delete pool;
	//ReadMerger no longer needed
	/*for (size_t x = 0; x < merger.size(); x++) {
		if (merger[x] != nullptr) {
			delete merger[x];
		}
	}*/
	printConvertedSummary();

}

void rewriteNumbers(){
	OptContainer& cmdArgs = globalOptContainer();
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

        if (line.empty() || line[0]=='$'){ //$ marks comment
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
    cerr << endl << "This is sdm (simple demultiplexer) " << 
		sdm_version << " " << sdm_status << ".\n" << endl ;
}
void help_head(){
    cout <<"------------------------------\nThis is sdm version "<<sdm_version <<" "<< sdm_status <<" help print\n------------------------------\n"<<endl;
}
void general_help(){
    help_head();
    cout<<"sdm (simple demultiplexer) is a fast, memory efficient program to demultiplex fasta and fastq files or simply do quality filterings on these.\n";
#ifdef _gzipread
    cout<<"Compiled with zlib gzip support\n";
#else
    cout << "No gzip support compiled\n";
#endif
#ifdef _isa1gzip
	cout << "Compiled with ISA-1 gzip support\n";

#endif
#ifdef _THREADED
    cout<<"Multithreading not supported"
#else
    cout << "The compiled version does not support multithreading\n";
#endif
	cout << "Select further help topics by typing:\nsdm -help_options : print help on configuring options files\n";
	cout << "sdm -help_flags : help on command line (flags) arguments for sdm\n";
	cout << "sdm -help_map : base_map files and the keywords for barcodes etc.";
	cout << "sdm -v : version information";
	cout << "\n------------------------------\n"; 
	cout << "For more information, please visit the sdm github page: https://github.com/hildebra/sdm\n";
	cout << "If you find bugs, please report them there.\n";
	cout << "If you have questions, please contact me at: ";
	cout<<" falk.hildebrand@gmail.com"<<endl;

}
void printCmdsHelp(){//actually these are flags
    help_head();
    string def_sep = DEFAULT_BarcodeNameSep;
	cout << "Usage: ./sdm <option> <value> ...\n\n";
	cout << "Input (choose -i_path, -i, -i_fastq, or -i_fna):\n"
		 << " -i_path <path>                 Directory or path containing input files\n"
		 << " -i <file>                      Input sequence file; format is autodetected\n"
		 << " -i_fastq <file>                FASTQ input\n"
		 << " -i_fna <file>                  FASTA input\n"
		 << " -i_qual <file>                 Quality input for FASTA; defaults to <fasta>.qual_\n"
		 << " -i_MID_fastq <file>            FASTQ file containing MID sequences\n"
		 << " -i_qual_offset <0|auto|n>      Input FASTQ quality offset\n"
		 << " -paired <1|2|3>                1: singleton; 2: paired; 3: R1/R3 plus MID R2\n"
		 << " -onlyPair <1|2>                Process only the selected read of a pair\n"
		 << " -mergedPairs <0|1>             Input pairs were merged externally\n"
		 << " -SyncReadPairs <T|F>           Synchronize paired input reads\n"
		 << " -ignore_IO_errors <0|1>        Continue after corrupt FASTQ reads (default: 0)\n\n";
	cout << "Mapping and output:\n"
		 << " -map <file>                    Tab-delimited sample mapping file\n"
		 << " -options <file>                sdm option file (default: sdm_options.txt)\n"
		 << " -o_fna <file> [-o_qual <file>] FASTA and optional quality output\n"
		 << " -o_fastq <file>                FASTQ output; overrides -o_fna and -o_qual\n"
		 << " -o_fna2/-o_qual2/-o_fastq2     Output files for a second read\n"
		 << " -o_demultiplex <path>          Write individual demultiplexed files\n"
		 << " -o_demultiplex_gz <0|1>        Gzip demultiplexed output\n"
		 << " -o_fastq_noBC <file>           Write reads for which no barcode was detected\n"
		 << " -o_qual_offset <n>             Output FASTQ quality offset (default: 33)\n"
		 << " -oneLineFastaFormat <0|1>      Write FASTA/quality sequences on one line\n"
		 << " -sample_sep <text>             Sample/ID delimiter (default: \"" << def_sep << "\")\n"
		 << " -log <file|nolog>              Demultiplex log destination\n"
		 << " -maxReadsPerOutput <n>         Split FASTA output after n filtered reads\n"
		 << " -pairedDemulti <0|1>           Write only complete pairs when demultiplexing\n"
		 << " -pairedRD_HD_out <0|1>         Include paired-read identifiers in output headers\n\n";
	cout << "Filtering, selection, and diagnostics:\n"
		 << " -5PR1cut/-5PR2cut <n>          Remove the first n nucleotides from read 1/read 2\n"
		 << " -onlyCanonicalNTs <0|1>        Retain only canonical nucleotides\n"
		 << " -strictPositiveQualsOut <0|1>  Require positive output quality values\n"
		 << " -specificReads <file>          Include named reads; combine with -excludeFile <file>\n"
		 << " -number <T> [-prefix <text>]   Renumber FASTA headers and exit\n"
		 << " -XfirstReadsRead <n>           Read only the first n reads\n"
		 << " -XfirstReadsWritten <n>        Write only the first n reads\n"
		 << " -length_hist/-qual_hist <file> Write length/quality histograms\n"
		 << " -merg_readpos/-qual_readpos <file> Write merge/quality position reports\n"
		 << " -logLvsQ <file|1>              Write length-versus-quality data\n"
		 << " -suppressOutput <0|1>          Suppress sequence output\n"
		 << " -iniBlockSize <n>              Initial read block size (default: 120)\n"
		 << " -threads <n>                   Worker thread count (default: 1)\n"
		 << " -threadIO <0|1>                Enable threaded I/O (default: 1)\n\n";
	cout << "Dereplication, merging, and OTU support:\n"
		 << " -o_dereplicate <file>          Write dereplicated FASTA reads\n"
		 << " -dere_size_fmt <0|1>           0: size=X; format; 1: _X format\n"
		 << " -min_derep_copies <rule>       Minimum copies, e.g. 10:1,3:3\n"
		 << " -derep_map <file>              Dereplication mapping output\n"
		 << " -derep_reserve <n>             Reserve capacity for dereplication\n"
		 << " -derepPerSR <0|1>              Dereplicate per sample/read\n"
		 << " -derep_format <text>           Dereplication output format\n"
		 << " -derepSrchLen <n>              Dereplication search length (default: 150)\n"
		 << " -merge_pairs_filter/-merge_pairs_seed/-merge_pairs_demulti <0|1> Enable pair merging\n"
		 << " -merge_pairs_derep <0|1>       Merge pairs during dereplication\n"
		 << " -noSearchWithMerge <value>     Disable search when pair merging is enabled\n"
		 << " -OTU_fallback/-OTU_fallback_refclust <file> Fallback OTU FASTA sequences\n"
		 << " -otu_matrix <file>             OTU matrix output\n"
		 << " -ucAdditionalCounts/-ucAdditionalCounts1/-ucAdditionalCounts_refclust/-ucAdditionalCounts_refclust1 <file> Additional UC counts\n"
		 << " -optimalRead2Cluster/-optimalRead2Cluster_ref <file> Read-to-cluster assignments\n\n";
	cout << "Technology and advanced options:\n"
		 << " -GoldenAxe <0|1>               Process PacBio concatenated reads\n"
		 << " -GoldenAxeMinAmpli <n>         Minimum accepted GoldenAxe amplicons (default: 0)\n"
		 << " -GoldenAxeMaxAmpli <n>         Maximum accepted GoldenAxe amplicons\n"
		 << " -illuminaClip <0|1>            Enable Illumina adapter clipping\n"
		 << " -DNAseedSelectionClassic <0|1> Use classic DNA seed selection\n"
		 << " -uparseVer <text>              UPARSE compatibility version\n"
		 << " -normRdsToFiveNTs <0|1>        Normalize reads to five nucleotides\n"
		 << " -mocatFix <0|1>                Enable MOCAT compatibility handling\n"
		 << " -minQueryCov <n>               Minimum query coverage\n"
		 << " -DemultiBPperSR <n>            Demultiplex base pairs per sample/read\n"
		 << " -binomialFilterBothPairs <0|1> Apply binomial filtering to both pairs\n"
		 << " -count_chimeras <T|F>          Count chimeric reads\n";
	cout<<"\nMinimal Example:\n./sdm -i test.fna -map mapping.txt (assuming quality file is \"test.qual_\")\n";
}
void printOptionHelp(){
    help_head();
	cout << "The option file specified by -options controls filtering, barcode and primer handling, paired-read checks, and adapter sequences. "
			"Each non-comment line must contain an option name and value separated by a tab. Lines beginning with # are comments.\n\n";
	string helpOptionFile = "";
	helpOptionFile += "# A leading * enables an alternate setting for the modified output stream.\n";
	helpOptionFile += "# Supported for: minSeqLength, minAvgQuality, maxAmbiguousNT,\n";
	helpOptionFile += "# QualWindowThreshhold, TrimWindowThreshhold, maxPrimerErrs,\n";
	helpOptionFile += "# keepPrimerSeq, maxAccumulatedError, RejectSeqWithoutRevPrim, and\n";
	helpOptionFile += "# RejectSeqWithoutFwdPrim.\n\n";
	helpOptionFile += "# Length and whole-read quality filtering. Length is evaluated after\n";
	helpOptionFile += "# barcode, primer, adapter removal, and quality trimming.\n";
	helpOptionFile += "minSeqLength\t250\n";
	helpOptionFile += "maxSeqLength\t1000\n";
	helpOptionFile += "minAvgQuality\t25\n";
	helpOptionFile += "maxAmbiguousNT\t0\n";
	helpOptionFile += "maxHomonucleotide\t12\n";
	helpOptionFile += "trimHomonucleotide\t12\n\n";
	helpOptionFile += "# Reject a read when a quality window falls below this average. A threshold\n";
	helpOptionFile += "# of 0 disables this check.\n";
	helpOptionFile += "QualWindowWidth\t50\n";
	helpOptionFile += "QualWindowThreshhold\t0\n\n";
	helpOptionFile += "# Trim a low-quality 3' end using this window and threshold.\n";
	helpOptionFile += "TrimWindowWidth\t15\n";
	helpOptionFile += "TrimWindowThreshhold\t20\n\n";
	helpOptionFile += "# Maximum accumulated expected error; use -1 to disable this trimming.\n";
	helpOptionFile += "maxAccumulatedError\t-1\n\n";
	helpOptionFile += "# Binomial error model. Alpha must be -1 (disabled) or between 0 and 1.\n";
	helpOptionFile += "BinErrorModelMaxExpError\t2.5\n";
	helpOptionFile += "BinErrorModelAlpha\t-1\n\n";
	helpOptionFile += "# Barcode and primer handling. Set keep*Seq to 1 to retain the matched\n";
	helpOptionFile += "# sequence in output, or 0 to remove it.\n";
	helpOptionFile += "maxBarcodeErrs\t0\n";
	helpOptionFile += "maxPrimerErrs\t1\n";
	helpOptionFile += "# Optional primer-variant warnings. Reporting is disabled while\n";
	helpOptionFile += "# PrimerVariantMinCount is 0. Variants are substitutions only.\n";
	helpOptionFile += "PrimerVariantMaxMismatches\t3\n";
	helpOptionFile += "PrimerVariantMinCount\t0\n";
	helpOptionFile += "PrimerVariantMinPercent\t1\n";
	helpOptionFile += "# Warn in the main log when one variant reaches this percentage.\n";
	helpOptionFile += "PrimerVariantWarningPercent\t50\n";
	helpOptionFile += "keepBarcodeSeq\t0\n";
	helpOptionFile += "keepPrimerSeq\t0\n";
	helpOptionFile += "ExtensivePrimerChecks\tF\n";
	helpOptionFile += "RejectSeqWithoutFwdPrim\tT\n";
	helpOptionFile += "RejectSeqWithoutRevPrim\tF\n\n";
	helpOptionFile += "# Read formatting and trimming.\n";
	helpOptionFile += "# fastqVersion accepts auto or a numeric FASTQ version understood by sdm.\n";
	helpOptionFile += "fastqVersion\tauto\n";
	helpOptionFile += "PEheaderPairFmt\t0\n";
	helpOptionFile += "TrimStartNTs\t0\n";
	helpOptionFile += "TruncateSequenceLength\t-1\n\n";
	helpOptionFile += "# Paired-read processing. Use T to enable and F to disable.\n";
	helpOptionFile += "AmpliconShortPE\tF\n";
	helpOptionFile += "CheckForMixedPairs\tF\n";
	helpOptionFile += "CheckForReversedSeqs\tF\n";
	helpOptionFile += "SyncReadPairs\tT\n\n";
	helpOptionFile += "# Dereplication. If omitted, derepSrchLen defaults to minSeqLength.\n";
	helpOptionFile += "derepSrchLen\t250\n\n";
	helpOptionFile += "# Set TechnicalAdapter to an adapter sequence to remove it. Do not include\n";
	helpOptionFile += "# this line when no technical adapter should be removed.\n";
	helpOptionFile += "# TechnicalAdapter\tTCAG\n\n";
	helpOptionFile += "# Adapter sequences used when -illuminaClip 1 is supplied. Do not include\n";
	helpOptionFile += "# unused adapter lines.\n";
	helpOptionFile += "# illuminaFwd\t<forward paired-end adapter>\n";
	helpOptionFile += "# illuminaRev\t<reverse paired-end adapter>\n";
	helpOptionFile += "# illuminaSngUni\t<single-end universal adapter>\n";
	helpOptionFile += "# illuminaSngIdx\t<single-end indexed adapter>\n";
	cout << helpOptionFile << endl;
}
void printMapHelp(){
    help_head();
	cout << "A mapping file, supplied with -map <file>, assigns input reads or files to samples. Without -map, sdm only performs quality filtering.\n\n";
	cout << "Format rules:\n"
		 << " - The file is tab-delimited. Every data row must contain the same number of columns.\n"
		 << " - The first row is the header; it may begin with #. Lines beginning with # after the header are comments.\n"
		 << " - Header names are case-sensitive and may appear in any order. Unrecognized columns are ignored.\n"
		 << " - SampleID values must be unique. Barcode values, or dual-barcode pairs, must be unique.\n\n";
	cout << "Sample identification modes:\n"
		 << " - BarcodeSequence: identify samples from a read barcode (IUPAC bases are supported).\n"
		 << " - BarcodeSequence plus Barcode2ndPair: identify samples using dual indexes.\n"
		 << " - SampleIDinHead: identify samples from a matching string in the FASTA/FASTQ header.\n"
		 << " - In -i_path mode, unique fastqFile or fnaFile entries can identify one sample per file.\n"
		 << " Use one identification mode consistently across the mapping file.\n\n";
	cout << "Recognized columns:\n"
		 << " Sample identity\n"
		 << "  SampleID             Unique output sample identifier.\n"
		 << "  CombineSamples       Combine this sample into the specified sample identifier.\n"
		 << "  SampleIDinHead       Header string that identifies the sample instead of barcode scanning.\n"
		 << "  SequencingRun        Sequencing-run identifier for the sample.\n"
		 << " Barcode assignment\n"
		 << "  BarcodeSequence      Primary barcode/MID sequence; IUPAC bases are supported.\n"
		 << "  Barcode2ndPair       Barcode sequence on the second read for dual-indexed input.\n"
		 << " Primers and spacers\n"
		 << "  ForwardPrimer        Forward amplification primer sequence.\n"
		 << "  LinkerPrimerSequence Legacy alias for ForwardPrimer.\n"
		 << "  ReversePrimer        Reverse amplification primer sequence.\n"
		 << "  HetSpacerFwd         Forward heterogeneity spacer sequence.\n"
		 << "  HetSpacerRev         Reverse heterogeneity spacer sequence. Both HetSpacer columns are required to enable explicit spacer handling.\n"
		 << " Input files (used with -i_path)\n"
		 << "  fastqFile            Relative FASTQ path.\n"
		 << "  fnaFile              Relative FASTA path.\n"
		 << "  qualFile             Relative FASTA-quality path; if absent, sdm derives <fasta-base>.qual_.\n"
		 << "  MIDfqFile            Relative FASTQ path containing MID sequences.\n"
		 << " Dereplication\n"
		 << "  derepMin             Empty or a non-negative minimum copy count for this sample.\n\n";
	cout << "When using -i_path, the mapping file must contain either fastqFile or fnaFile. File paths are resolved relative to -i_path.\n\n";
	cout << "Examples (tabs are shown as spaces for readability; use real tab characters in the file):\n"
		 << "# Single barcode\n"
		 << "#SampleID    BarcodeSequence    ForwardPrimer\n"
		 << "sample_1     ACGTACGT           CCTACGGG\n"
		 << "sample_2     TGCATGCA           CCTACGGG\n\n"
		 << "# Dual-indexed paired reads\n"
		 << "#SampleID    BarcodeSequence    Barcode2ndPair\n"
		 << "sample_1     ACGTACGT           GATCAAGT\n\n"
		 << "# One sample per FASTQ file with -i_path reads\n"
		 << "#SampleID    fastqFile\n"
		 << "sample_1     sample_1.fastq.gz\n";

}
void printVersion(){
    cout << "sdm " << sdm_version << " " << sdm_status << endl;
	cout << "compiled on: " << __DATE__ << " at " << __TIME__ << endl;
	cout << "compiled with C++ v: " << __cplusplus << endl;

}

