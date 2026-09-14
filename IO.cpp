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

	auto& read1 = tdn[0];
	auto& read2 = tdn[1];
	auto& midRead = tdn[2];

	if (read1 == nullptr) {
	    return false;
	} //|| tdn->length()==0
    Instr::ScopedTimer st_proc(Instr::PROCESS_DNA);
	//DNA objects as they should be??
	if (MIDuse && midRead == nullptr) {
		//cerr << "Missing MID read pair.\n";		exit(4);
		throw std::runtime_error("Missing MID read pair.");

	}
	const bool read2notNull = read2 != nullptr;
	if (!read2notNull && MD->isPEseq() == 2) {
		//cerr << "Second provided file has not the same number of entries as first file.\n";exit(5);
		throw std::runtime_error("Second provided file has not the same number of entries as first file.");
	}
	//MD->checkFastqHeadVersion(tdn[0]);
	//testreadpair.lock();

	//cdbg("process_DNA"+ tdn[0]->getId()+"\n");

	Filters* curFil = MD->getFilters(curThread);
	const int mdPeSeq = MD->isPEseq();
	const int pairedRd = curFil->isPaired();
	const bool doBCsAtAll = curFil->doBarcodes();
	const bool dualBCs = curFil->doubleBarcodes();
	const int BCoffs = curFil->getBCoffset();
	//register read at all with stat counter:
	 curFil->sTotalPlus(0);
	if (read2notNull || pairedRd == 2) {
		curFil->sTotalPlus(1);
	}
	//collect some info on general run parameters
	curFil->preFilterSeqStat(read1, 0);
	curFil->preFilterSeqStat(read2, 1);



	//prep some variables
	bool checkSwitchedRdPairs = curFil->checkSwitchedRdPairs();
	(void)checkSwitchedRdPairs;
	//bool checkReversedRead = curFil->checkRevRd();

	int tagIdx(-2); int tagIdx2(-2);
	string presentBC(""); int c_err(0);
	//bool isReversed(false);//was a reversion detected?

	if (MIDuse && midRead != nullptr) {
		tagIdx = curFil->findTag(midRead, presentBC, c_err, true,0,true,false); 
		read1->setBCnumber(tagIdx, BCoffs);
	}
	//routine checks, and reverses/swaps DNA objects
	//bool reversedDNA = curFil->swapReverseDNApairs(tdn);


	//now split up into multiple PB reads, if GoldenAxe (or Kinnect)
	


	//manage reversed / swapped reads by detecting where the fwd primer is
	bool wasReversed = false;
	if (pairedRd == 2) {
		wasReversed= curFil->swapReverseDNApairs(tdn);
	} else if (!read2notNull) {
		wasReversed = curFil->isReversedAmplicon(read1);
	}
	(void)wasReversed;


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
	MD->analyzeDNA(read1, -1, 0, tagIdx, curThread);
	//tdn[0]->matchSeqRev
	const bool ch1 = read1->isGreenQual();
	bool ch2(false);



	//if ( ch1 ) {	cerr << cnt << " \n";	}
	//normal case for check 2nd read
	if (read2notNull) { //ch1&&
								//dnaTemp2->setBCnumber(tdn[0]->getBCnumber());
		if (doBCsAtAll && !dualBCs) { //only check in read1 for BC, if not dual BCing!!
			tagIdx2 = read1->getBCnumber();  // no 2nd BC, thus no BC search in 2nd read
			if (tagIdx2 >= 0) {
				tagIdx2 -= BCoffs;
			}else if (tagIdx2 < -1) {//something wrong with BCoffs
				cerr << "tagidx2 wrongly truncated to " << tagIdx2 << endl;
			}
		}
		//if (isReversed) { tdn[1]->reverse_compliment(); }
        if (!dualBCs && !read2->sampleBarcode().empty()) {
            tagIdx2 = curFil->detectCutBC(read2, true);
        }
		MD->analyzeDNA(read2, -1, 1, tagIdx2, curThread);
		ch2 = read2->isGreenQual();
	}

	//set up BC in DNA header
	//remember that dual BCs are only valid after this step!
	if (dualBCs && pairedRd == 2) {
		//in case of single reads (PB) with double BCs
		//if (dualBCs && tagIdx2 < 0 && pairedRd == 1) {
		//	tagIdx2 = curFil->findTag2(tdn[0], presentBC, c_err, false, -1);
		//}

		curFil->dblBCeval(tagIdx, tagIdx2, presentBC, read1, read2);
		c_err = -1;

		//check a second time that barcode was correctly identified, just to be double sure...
		if ( tagIdx != tagIdx2 || read1->getBCnumber() != read2->getBCnumber()) {
			//cerr << "Unequal BC numbers:" << tagIdx << " : " << tagIdx2 << "; in object: " << tdn[0]->getBCnumber() << " : " << tdn[1]->getBCnumber() << endl;
			//cerr << "In read:" << tdn[0]->getId() << endl; exit(835);
			throw std::runtime_error("Unequal BC numbers:" + itos(tagIdx) + " : " + itos(tagIdx2) + "; in object: " + itos(read1->getBCnumber()) + " : " + itos(read2->getBCnumber()) + "\nIn read:" + read1->getId());			
		}
	}
	

	
	if (tagIdx == -1 || tagIdx2 == -1) {
		read1->setBarcodeDetected(false);
		if (read2notNull) {
			read2->setBarcodeDetected(false);
		}
	}	else  {
		curFil->setBCdna(tagIdx, read1);
		if (ch2) { curFil->setBCdna(tagIdx, read2); }
		if (MIDuse && ch1) {
			curFil->BCintoHead(tagIdx, read1, presentBC, c_err, true);
			if (ch2) { curFil->BCintoHead(tagIdx, read2, presentBC, c_err, true); }
		}
	}



	int idx1 = 1; int idx2 = 2;

	if (mdPeSeq == 2 ) {
		if (read2notNull &&
			(read1->isGreenYellowQual() && read2->isGreenYellowQual()) ) {
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
		MD->dereplicateDNA(tdn[0], tdn[1], curThread);
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
		owr1 = MD->saveForWrite(read1, idx1, curThread, Cstr1,false);
		owr2 = MD->saveForWrite(read2, idx2, curThread, Cstr2,false);
		//mutex version to ensure read pairs together..
		MD->writeForWrite(read1,idx1,Cstr1,read2,idx2,Cstr2);
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

   qual_score fastqVer = IS->fastQscore();

	// Probe one block first and adapt batch size for short reads to reduce task scheduling overhead.
	// This keeps behavior unchanged for longer reads while increasing reads/task when records are tiny.
	int taskBatchMultiplier = 1;
	string empty(""); empty.reserve(151);
	{
		cdbg("task batch size probe");
		multi_tmp_lines probe(tmpBlockSize);
		bool probeCont = IS->getDNAlines(&probe, tmpBlockSize, MIDuse, true);
        fastqVer = IS->fastQscore() == 0 ? 33 : IS->fastQscore();
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
	vector<unique_ptr<multi_tmp_lines>> tmpStrHolders(Nthreads);	
	vector<string>tmpLines2(4, empty);
	vector<vector<string>> tmpLines(3, tmpLines2);
	for (size_t ii = 0; ii < Nthreads; ii++) {
		tmpStrHolders[ii] = make_unique<multi_tmp_lines>();
		tmpStrHolders[ii]->tmp.resize(tmpBlockSize, tmpLines);
	}

    //object holding converted DNA.. not longer used
	//vector<shared_ptr<DNA>> tdn(3, nullptr);
	cdbg("read_seq check4 ");
	//doMC = false;//DEBUG
	std::exception_ptr failure;
	try {
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
				cont = IS->getDNAlines(tmpStrHolders[thrCnt].get(), tmpBlockSize, MIDuse);
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
							tmpStrHolders[thrCnt].get(), MIDuse, MD, thrCnt, keepPairedHD, fastqVer);
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
           cont = IS->getDNAlines(tmpStrHolders[thrCnt].get(), tmpBlockSize, MIDuse);
			// If we read any records into the tmp buffer, process them even if cont==false (partial final block).
			if (tmpStrHolders[thrCnt]->size() > 0) {
				//fastqVer = IS->fastQscore();
				bool cont2 = multi_read_paired_STRready(tmpStrHolders[thrCnt].get(), MIDuse, 
					MD, 0, keepPairedHD, fastqVer);
			}
			if (!cont && tmpStrHolders[thrCnt]->size() == 0) { break; }
				//if (tdn[0]->isConstellationPairRev()) { revConstellation++; }
			//}
		}
	}
	}
	catch (...) {
		failure = std::current_exception();
	}

	// A packaged-task future does not join on destruction. Keep the input buffers
	// and OutputStreamer alive until every submitted worker has finished, even
	// when input parsing, submission, or another worker has already failed.
	for (size_t i = 0; i < slots.size(); ++i) {
		const size_t x = (static_cast<size_t>(thrCnt) + i) % slots.size();
		if (slots[x].job.valid()) {
			try {
				(void)slots[x].job.get();
			}
			catch (...) {
				if (!failure) { failure = std::current_exception(); }
			}
		}
		slots[x].inUse = false;
	}
	if (failure) { std::rethrow_exception(failure); }

	//close shop
	//MD->revConstellationCnts(revConstellation);
	MD->closeOutStreams();

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
		// Arguments are key/value pairs; retain negative numeric values such as -1.
		cmdArgs->set(string(argv[i]), theNxtSt);
	}
	
	
	setDefault("-illuminaClip", "0");
	setDefault("-logLvsQ", "");
	setDefault("-GoldenAxe", "0");
	setDefault("-GoldenAxeMinAmpli", "0");
	setDefault("-GoldenAxeMaxAmpli", "0");
	setDefault("-DNAseedSelectionClassic", "0");
	setDefault("-i_MID_fastq", "");
    setDefault("-cramRef", "");
    setDefault("-htsBarcode", "inline");
    if ((*cmdArgs)["-htsBarcode"] != "inline" && (*cmdArgs)["-htsBarcode"] != "BC")
        throw std::runtime_error("-htsBarcode must be inline or BC");
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
        for (const string key : {"-i_sam", "-i_bam", "-i_cram"}) {
            if (cmdArgs->contains(key)) {
                if (cmdArgs->contains("-i_fastq") || cmdArgs->contains("-i_fna"))
                    throw std::runtime_error("Specify only one input format");
                cmdArgs->set("-i_fastq", cmdArgs->get(key));
            }
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
	// Position reports are optional: collect only when an output is requested.
	setDefault("-merg_readpos", "");
	setDefault("-qual_readpos", "");

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
	unique_ptr<UClinks> ucl; // Own seed state through every exit path.
	shared_ptr<ReadSubset> RDSset (nullptr);  //read subset filtered out?
	shared_ptr<Dereplicate> dereplicator (nullptr); //dereplication of b_derep_as_fasta_ input?
	ucl.reset(mainFilter->ini_SeedsReadsDerep(nullptr, RDSset, dereplicator)); // actual prep
	//needs to attach to existing file sometimes
	std::ios_base::openmode writeStatus = ios_base::out;
	bool shortStats = false;
	string shrtLog = "";

	// main loop that goes over different files
	uint64_t maxReadsWr = mainFilter->getXreadsWr(); // should program stop after having written a certain amount of reads?
	uint64_t maxReadsRd = mainFilter->getXreadsRd();
	uint64_t totalReadsRead(0);
	uint64_t totalReadsWrite(0);
	auto printConvertedSummary = [&]() {
		if ((*cmdArgs)["-o_fna"] != "") {
			cerr << "Converted sequences written: " << totalReadsWrite << "\n";
		}
	};
	uint64_t accumBPwrite(0), accumBPwriteMerg(0);
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
	const bool detailedMergeStats = !(*cmdArgs)["-merg_readpos"].empty() || !(*cmdArgs)["-qual_readpos"].empty();
	ReadMerger positionStats;
	auto printPositionStats = [&] {
		if (!(*cmdArgs)["-merg_readpos"].empty()) { positionStats.printMergeHisto((*cmdArgs)["-merg_readpos"]); }
		if (!(*cmdArgs)["-qual_readpos"].empty()) { positionStats.printQualHisto((*cmdArgs)["-qual_readpos"]); }
	};


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

		if (maxReadsWr > 0 && totalReadsWrite >= maxReadsWr) { cerr << "Skipping file " << uFX.first << "due to firstXreadsWrite .."; break; }
		if (maxReadsRd > 0 && totalReadsRead >= maxReadsRd) { cerr << "Skipping file " << uFX.first << "due to firstXreadsRead .."; break; }
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
		OutStreamer->attachDereplicator(dereplicator);OutStreamer->activateReadMerger(threads, detailedMergeStats && !mainFilter->doDereplicate());
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
        mainFilter->copyOutputState(*filter);
        // Print out merged read count (move to stort stats)

		//OutStreamer.reset();
		totalReadsWrite += OutStreamer->getReadsWritten();
		
		if (detailedMergeStats && !mainFilter->doDereplicate()) { OutStreamer->addMergeStatsTo(positionStats); }
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
		MDx->activateReadMerger(1, detailedMergeStats);
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
        if (detailedMergeStats) {
            ucl->addMergeStatsTo(positionStats);
            MDx->addMergeStatsTo(positionStats);
            printPositionStats();
        }
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



#ifdef DEBUG
        cerr << "done write dereplicator" << endl;
#endif
    }
#ifdef DEBUG
    cerr << "Logging almost finished" << endl;
#endif

    if (detailedMergeStats && !mainFilter->doDereplicate()) { printPositionStats(); }

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
namespace {
void printBuildFeatures() {
#ifdef _gzipread
    cout << "Gzip support: enabled (zlib)\n";
#else
    cout << "Gzip support: disabled\n";
#endif
    // The std::thread worker pool is part of every current build.
    const OptContainer::RuntimeOptions defaults;
    cout << "ONT amplicon end matching: enabled (-ontMode 1)\n";
    cout << "Multithreading support: enabled (std::thread; -threads <n>, default: "
         << defaults.threads << ")\n";
    const string htsVersion = HtsInput::version();
    cout << "HTSlib: " << htsVersion;
    if (htsVersion == "disabled") {
        cout << " (SAM/BAM/CRAM unavailable; rebuild with make WITH_HTS=1)\n";
    } else {
        cout << " (SAM/BAM/CRAM input enabled)\n";
    }
}
}

void help_head(){
    cout << "sdm " << sdm_version << " " << sdm_status << " - help\n";
    printBuildFeatures();
    cout << endl;
}
void general_help(){
    help_head();
    cout << "sdm (simple demultiplexer) filters DNA reads, assigns samples, removes\n"
            "barcodes and primers, merges paired reads, and dereplicates sequences.\n"
            "FASTA and FASTQ input are supported; alignment input requires HTSlib.\n\n"
            "Usage: ./sdm <option> <value> ...\n"
            "Example: ./sdm -i reads.fastq -o_fastq filtered.fastq -options sdm_options.txt -threads 4\n\n"
            "Help topics (no value required):\n"
            "  sdm -help_options  Filtering options, defaults, and example settings\n"
            "  sdm -help_flags    Command-line inputs, outputs, and processing controls\n"
            "  sdm -help_map      Sample mapping columns, barcodes, and alignment files\n"
            "  sdm -v / -version  Version and compiled features\n\n"
            "Options files and mapping files are tab-delimited. If the options file\n"
            "cannot be opened, sdm continues with filtering disabled.\n\n"
            "Documentation and bug reports: https://github.com/hildebra/sdm\n"
            "Contact: falk.hildebrand@gmail.com\n";
}
void printCmdsHelp(){//actually these are flags
    help_head();
    string def_sep = DEFAULT_BarcodeNameSep;
	cout << "Usage: ./sdm <option> <value> ...\nValues are required except for help/version commands. Quote paths with spaces.\n\n";
	cout << "Input (choose one primary input mode):\n"
		 << " -i_path <path>                 Input directory; filenames come from -map\n"
		 << " -i <file>                      Autodetect FASTA/FASTQ/SAM/BAM/CRAM input\n"
		 << " -i_fastq <file>                FASTQ input\n"
         << " -i_sam/-i_bam/-i_cram <file>   Alignment input (requires HTSlib); -i detects these too\n"
         << " -cramRef <reference.fa>        Optional CRAM reference; default: none supplied\n"
         << " -htsBarcode <inline|BC>        Inline barcode (default), or exact BC:Z tag matching\n"
		 << " -i_fna <file>                  FASTA input\n"
		 << " -i_qual <file>                 FASTA quality file; e.g. reads.fa -> reads.qual_\n"
		 << " -i_MID_fastq <file>            FASTQ file containing MID sequences\n"
		 << " -i_qual_offset <0|auto|n>      FASTQ ASCII offset; 0/auto detects (not an HTS setting)\n"
		 << " -paired <1|2|3>                1: singleton; 2: paired; 3: FASTQ R1/MID/R3; default: infer\n"
		 << " -onlyPair <1|2>                Select a mate; HTS input still needs complete sorted pairs\n"
		 << " -mergedPairs <0|1>             Input pairs were merged externally\n"
		 << " -ignore_IO_errors <0|1>        Try FASTQ error recovery; HTS errors remain fatal (default: 0)\n\n";
	cout << "Paired FASTQ uses comma-separated R1,R2 filenames; -paired 3 uses R1,MID,R3.\n"
         << "Each HTS input set is one local SAM/BAM/CRAM file. Singleton records may be in\n"
         << "any order; paired primary records must be adjacent and query-name sorted.\n"
         << "Secondary/supplementary alignments are skipped. Reverse-strand reads are restored\n"
         << "to sequencing orientation. Missing sequence/quality and corrupt records are fatal.\n"
         << "HTS qualities use numeric Phred scores, bypassing FASTQ offset detection.\n"
         << "-cramRef applies to command-line and mapping inputs; use an indexed matching FASTA.\n"
         << "BC mode matches BC:Z tags exactly against BarcodeSequence/Barcode2ndPair and\n"
         << "does not cut bases. Separate MID FASTQ streams cannot accompany HTS input.\n\n";
    cout << "ONT amplicons (opt-in; mapping required):\n"
         << " -ontMode <0|1>                 Single-end FASTQ/FASTQ.gz end matching (default: 0)\n"
         << " -barcodeSearchWindow <bp>     End-relative start offsets 0..bp-1 (1..10000; ONT default: 200)\n"
         << " -ontBarcodeEnds <either|both> Unique consistent end evidence (default: either)\n"
         << "ONT reuses maxBarcodeErrs as edit distance (0..3), including indels.\n"
         << "Primers use maxPrimerErrs substitutions; barcode-only reads retain input strand.\n"
         << "Barcode2ndPair, paired/MID, HTS and GoldenAxe modes are incompatible.\n\n";
	cout << "Mapping and output:\n"
		 << " -map <file>                    Tab-delimited sample mapping file\n"
		 << " -options <file>                sdm option file (default: sdm_options.txt)\n"
		 << " -o_fna <file> [-o_qual <file>] FASTA and optional quality output\n"
		 << " -o_fastq <file>                FASTQ output; overrides -o_fna and -o_qual\n"
		 << " -o_fna2/-o_qual2/-o_fastq2 <file> Additional filtering tier, configured by * options\n"
		 << " -o_demultiplex <path>          Per-sample FASTQ output; directory must already exist\n"
		 << " -o_demultiplex_gz <0|1>        Gzip demultiplexed output\n"
		 << " -o_fastq_noBC <file>           No-barcode reads; paired: one interleaved file or R1,R2\n"
		 << " -o_qual_offset <n>             Output FASTQ quality offset (default: 33)\n"
		 << " -oneLineFastaFormat <0|1>      Write FASTA/quality sequences on one line\n"
		 << " -sample_sep <text>             Sample/ID delimiter (default: \"" << def_sep << "\")\n"
		 << " -log <file|nolog>              Demultiplex log destination\n"
		 << " -maxReadsPerOutput <n>         Split ordinary FASTA/FASTQ output after n record groups\n"
		 << " -pairedDemulti <0|1>           Write only complete pairs when demultiplexing\n"
		 << " -pairedRD_HD_out <0|1>         Include paired-read identifiers in output headers\n\n";
	cout << "Paired output uses comma-separated R1,R2 paths in -o_fastq or -o_fna/-o_qual.\n"
         << "The suffix 2 selects an additional filtering tier, not mate 2. Mates passing\n"
         << "different tiers are written together to the less stringent configured tier.\n"
         << "Record-group limits count one single read, pair or merged read as one group;\n"
         << "0 disables a limit. They do not randomly subsample input.\n\n";
	cout << "Filtering, selection, and diagnostics:\n"
		 << " -5PR1cut/-5PR2cut <n>          Remove the first n nucleotides from read 1/read 2\n"
		 << " -onlyCanonicalNTs <0|1>        Derep output: replace noncanonical bases with A, quality 0\n"
		 << " -strictPositiveQualsOut <0|1>  Derep FASTQ output: clamp negative qualities to 0\n"
		 << " -specificReads <file>          Extract named reads; -excludeFile receives unselected reads\n"
		 << " -number <T> [-prefix <text>]   Renumber FASTA headers and exit\n"
		 << " -XfirstReadsRead <n>           Read at most n input reads (mates counted individually)\n"
		 << " -XfirstReadsWritten <n>        Write at most n record groups (a pair counts as one)\n"
		 << " -length_hist/-qual_hist <file> Write length/quality histograms\n"
		 << " -merg_readpos/-qual_readpos <file> Enable detailed merge/quality position reports\n"
		 << " -logLvsQ <file|1>              Write length-versus-quality data\n"
		 << " -suppressOutput <0|1|2|3>      Suppress ordinary output: 0=none, 1=primary, 2=additional, 3=both\n"
		 << " -iniBlockSize <n>              Initial read block size (default: 120)\n"
		 << " -threads <n>                   Worker thread count (default: 1)\n"
		 << " -threadIO <0|1>                Enable threaded I/O (default: 1)\n\n";
	cout << "-threads controls SDM processing workers. -threadIO controls the FASTA/FASTQ\n"
         << "I/O backend, independently of worker count. HTS input decoding is currently\n"
         << "sequential; decoded reads still use the worker pool for processing.\n\n";
	cout << "Dereplication, merging, and OTU support:\n"
		 << " -o_dereplicate <file>          Dereplicated sequences (default format: fa)\n"
		 << " -dere_size_fmt <0|1>           0: size=X; format; 1: _X format\n"
		 << " -min_derep_copies <rule>       Minimum copies, e.g. 10:1,3:3\n"
		 << " -derep_map <file>              Dereplication mapping output\n"
		 << " -derep_reserve <n>             Reserve capacity for dereplication\n"
		 << " -derepPerSR <0|1>              Dereplicate per mapping-file SequencingRun (default: 0)\n"
		 << " -derep_format <fa|fq>          Dereplicated FASTA or FASTQ output (default: fa)\n"
		 << " -derepSrchLen <n>              Explicit dereplication search length (-1: full logical sequence)\n"
		 << " -derepPrefix <auto|0|1>        Auto-select, disable, or force 100% prefix matching\n"
         << " -derepIdentity <percent>      Streaming substitution clusters below 100 (default: 100)\n"
         << " -derepSeedHomopolymer <2..20> Exclude search seeds containing this run length (default: 4)\n"
         << " -derepReassign <0|1>          Final nearest-representative reassignment; retains variants/samples (default: 0)\n"
         << " -derepStoreDiffs <0|1>        Write nucleotide .diff: SDMDIFF2 single, SDMDIFF3 pairs (default: 0)\n"
         << " -derepStoreQuals <0|1>        Retain selected qualities per exact sequence or pair (default: 0)\n"
         << "                              Single: .subclusters.fq; paired: .subclusters.1.fq/.subclusters.2.fq\n"
         << "                              Complete unmerged reads with qualities; both mate qualities selected together.\n"
         << "                              Pair size counts fragments once; both mates must meet derepIdentity.\n"
         << "                              Implies in-memory variants, not .diff export; works with derepReassign.\n"
         << "                              Disable merge_pairs_derep, merge_pairs_filter, and merge_pairs_demulti.\n"
         << " -derepSubclusterFasta <0|1>   Export exact single/pair FASTA; requires derepStoreDiffs=1 (default: 0)\n"
		 << " -merge_pairs_filter/-merge_pairs_seed/-merge_pairs_demulti <0|1> Enable pair merging\n"
		 << " -merge_pairs_derep <0|1>       Merge pairs during dereplication\n"
		 << " -noSearchWithMerge <value>     Use R1 rather than merged derep search; presence enables, even value 0\n"
		 << " -OTU_fallback/-OTU_fallback_refclust <file> Fallback OTU FASTA sequences\n"
		 << " -otu_matrix <file>             OTU matrix output\n"
		 << " -ucAdditionalCounts/-ucAdditionalCounts1/-ucAdditionalCounts_refclust/-ucAdditionalCounts_refclust1 <file> Additional UC counts\n"
		 << " -optimalRead2Cluster/-optimalRead2Cluster_ref <file> Read-to-cluster assignments\n\n";
	cout << "A copy rule 10:1,3:3 means >=10 copies in one sample OR >=3 copies in each\n"
         << "of three samples. derepPrefix=auto uses exact keys for guaranteed fixed lengths\n"
         << "or an explicit positive derepSrchLen; otherwise it uses 100% prefix matching.\n"
         << "Forced prefix mode overrides search length; TruncateSequenceLength still caps\n"
         << "clustering keys. See -help_options for length and output interactions.\n\n";
	cout << "Technology and advanced options:\n"
		 << " -GoldenAxe <0|1>               Process PacBio concatenated reads\n"
		 << " -GoldenAxeMinAmpli <n>         Minimum accepted amplicons (default: 0; -1 disables)\n"
		 << " -GoldenAxeMaxAmpli <n>         Maximum accepted amplicons (0/-1: no upper bound)\n"
		 << " -illuminaClip <0|1>            Enable Illumina adapter clipping\n"
		 << " -DNAseedSelectionClassic <0|1> Use classic DNA seed selection\n"
		 << " -uparseVer <text>              UPARSE compatibility version\n"
		 << " -normRdsToFiveNTs <value>      Legacy placeholder: not implemented; omit this option\n"
		 << " -mocatFix <0|1>                Enable MOCAT compatibility handling\n"
		 << " -minQueryCov <fraction>        Seed/OTU query coverage (default: 0.8)\n"
		 << " -DemultiBPperSR <n>            Demultiplex bases per sequencing run (both mates; 0=unlimited)\n"
		 << " -binomialFilterBothPairs <0|1> Apply binomial filtering to both pairs\n"
		 << " -id <percent>                 Seed/OTU identity threshold (default: 97)\n"
		 << " -count_chimeras <T|F>          Count chimeric reads\n";
	cout << "\nExamples:\n"
         << "./sdm -i reads.fastq -o_fastq filtered.fastq -options sdm_options.txt -threads 4\n"
         << "./sdm -i reads.by-name.bam -paired 2 -o_fastq R1.fq,R2.fq -options sdm_options.txt\n"
         << "./sdm -i_path reads -map samples.tsv -o_fastq filtered.fq -cramRef reference.fa\n"
         << "\nHelp/version (no value): -help_flags, -help_options, -help_map, -v, -version\n"
         << "Documentation: docs/options.md, docs/mapping.md, docs/hts.md, docs/output.md\n";
}
void printOptionHelp(){
    help_head();
    cout << R"SDM(The file selected by -options configures filtering, trimming, primers and barcodes.
Use case-sensitive names without a leading dash, separated from values by a tab.
Lines beginning with # are comments. Keep the spellings Threshhold and Homonucleotide.
Sample names, barcode sequences and primers belong in -map; filenames, threads and
HTS controls (-cramRef, -htsBarcode) are command-line options.

Without -options, sdm looks for sdm_options.txt in the current directory. If the
file cannot be opened, processing continues with filtering disabled. The settings
below show defaults when an options file is opened; adjust them for your data.
Copy only the settings and # comments into your options file.

# Length and quality. Minimum length is checked after removal and trimming.
# minSeqLength >1 is an absolute length; 0..1 is relative to observed read length.
# A negative maxAmbiguousNT disables that check. Homopolymer settings use 0 to disable.
minSeqLength	250
maxSeqLength	1000
minAvgQuality	25
maxAmbiguousNT	0
maxHomonucleotide	12
trimHomonucleotide	12

# Reject low-quality windows; a threshold of 0 disables the check.
QualWindowWidth	50
QualWindowThreshhold	0
# Trim the low-quality 3' end using this window and threshold.
TrimWindowWidth	15
TrimWindowThreshhold	20
# Accumulated expected-error trimming; -1 disables it.
maxAccumulatedError	-1
# Binomial alpha: -1 disables the model; otherwise use 0..1.
# -binomialFilterBothPairs 1 applies the model to both mates.
BinErrorModelMaxExpError	2.5
BinErrorModelAlpha	-1

# Barcode and primer mismatches; 0/1 removes/retains matched sequence.
# BC:Z tag mode uses exact matching and never cuts barcode bases from reads.
maxBarcodeErrs	0
maxPrimerErrs	1
keepBarcodeSeq	0
keepPrimerSeq	0
ExtensivePrimerChecks	F
RejectSeqWithoutFwdPrim	F
RejectSeqWithoutRevPrim	F

# Primer substitution-variant reports: positive MinCount enables reporting.
# MaxMismatches accepts 1..3; percentages accept 0..100.
# These report variants without changing the primer matching threshold.
PrimerVariantMaxMismatches	3
PrimerVariantMinCount	0
PrimerVariantMinPercent	0
PrimerVariantWarningPercent	50

# FASTQ version selectors: auto, 1 (offset 33), 2 (offset 64), 3 (legacy offset 59).
# Use -i_qual_offset for a raw offset. If both are given, this option-file setting
# is applied after -i_qual_offset. HTS input always bypasses FASTQ offset selection.
fastqVersion	auto
# Paired headers: 0 preserves formatting; 1 converts Illumina headers to /1 or /2.
# Conversion to mode 2 is not implemented; use 0 or 1.
PEheaderPairFmt	0
TrimStartNTs	0
# -1 disables the cap on ordinary reads and unmerged dereplication output.
# A smaller positive cap also lowers an absolute ordinary-read minimum.
# Merged output and HQ mates retain full length; seed extension ignores this cap.
TruncateSequenceLength	-1

# Use T/F for paired-read checks. IDs are always checked for synchronization.
# AmpliconShortPE checks reverse primers on R1; MixedPairs checks swapped mates;
# ReversedSeqs checks reads in the reverse orientation.
AmpliconShortPE	F
CheckForMixedPairs	F
CheckForReversedSeqs	F

# Dereplication: command-line values override the corresponding file settings.
# auto: exact keys for guaranteed fixed lengths or an explicit positive search
# length; otherwise 100% identity prefix matching at the same sequence start.
# 0: legacy length-matched keys. 1: force prefix matching, overriding search length.
# TruncateSequenceLength caps both R1 and merged clustering keys in prefix modes,
# without shortening stored merged output or HQ mates. No gaps/mismatches are allowed.
derepPrefix	auto
# Below 100: streaming same-start substitution clustering over the shorter key.
# Fixed initial anchors; one high-quality representative and all assigned counts.
# This mode takes precedence over derepPrefix and still respects search/truncation caps.
derepIdentity	100
# Search seeds containing this many consecutive identical bases are excluded.
derepSeedHomopolymer	4
# Optional nucleotide-only .diff: SDMDIFF2 single sequences; SDMDIFF3 linked pairs.
# Paired storage is selected by paired input together with derepStoreQuals=1.
# Neither format stores qualities, per-read names, or sample IDs.
# Setting this to 1 also selects the streaming mode at identity 100.
derepStoreDiffs	0
# Optional selected qualities per exact sequence or (R1,R2) variant; default off.
# Single-end output: .subclusters.fq; paired: .subclusters.1.fq and .subclusters.2.fq.
# Both mate quality vectors come from one fragment; ;size= counts fragments once.
# Requires complete, unmerged reads with qualities. Disable merge_pairs_derep,
# merge_pairs_filter, and merge_pairs_demulti. No additional paired retention flag.
# Retains nucleotide variants and selects streaming mode at identity 100;
# does not enable binary .diff output. Can be combined with derepReassign=1.
derepStoreQuals	0
# Optional exact subcluster FASTA with ;size= counts; requires derepStoreDiffs=1.
# Single: .subclusters.fna; paired quality mode: .subclusters.1.fna and .subclusters.2.fna.
derepSubclusterFasta	0
# Optional final reassignment to strictly closer qualifying representatives; default off.
# Retains variants/sample counts in memory; does not enable sidecar output.
# Paired quality mode moves whole variants, both quality vectors, and fragment counts.
# Both mates must qualify; combined mismatch rate ranks targets. Ties keep the original.
derepReassign	0
# Optional derepSrchLen: positive bases or -1 for the full logical sequence.
# An explicit positive value selects length-matched keys in auto mode.
# Legacy default: effective absolute minSeqLength, or 150 when unavailable.
# derepSrchLen	250

# Optional technical adapter; omit when unused.
# TechnicalAdapter	TCAG
# Illumina adapters also require -illuminaClip 1. Supply only relevant sequences.
# illuminaFwd	<forward paired-end adapter>
# illuminaRev	<reverse paired-end adapter>
# illuminaSngUni	<single-end universal adapter>
# illuminaSngIdx	<single-end indexed adapter>

Second filtering tier:
A leading * sets an alternate threshold. Put the primary setting before its
alternate and specify the alternate thresholds you need explicitly. For example:
minAvgQuality	25
*minAvgQuality	20

Supported alternates: minSeqLength, minAvgQuality, maxAmbiguousNT,
QualWindowThreshhold, TrimWindowThreshhold, maxPrimerErrs, keepPrimerSeq,
maxAccumulatedError, RejectSeqWithoutRevPrim, RejectSeqWithoutFwdPrim.
Use -o_fastq2, or -o_fna2/-o_qual2, for additional-tier output. If alternate settings
are present and these paths are omitted, names are derived from the primary output.
For paired output, supply R1,R2 paths within each tier; mates with different passing
tiers stay together in the less stringent configured output.

Examples:
./sdm -i reads.fq -o_fastq filtered.fq -options sdm_options.txt
./sdm -i reads.fq -o_fastq primary.fq -o_fastq2 additional.fq -options sdm_options.txt

More detail: docs/options.md and docs/derep-prefix-evaluation.md.
)SDM";
}
void printMapHelp(){
    help_head();
	cout << "A mapping file, supplied with -map <file>, assigns input reads or files to samples. Without -map, sdm only performs quality filtering.\n\n";
	cout << "Format rules:\n"
		 << " - The file is tab-delimited. Every data row must contain the same number of columns.\n"
		 << " - Blank lines and comments may precede the header, which may begin with #. Later # lines are comments.\n"
		 << " - Header names are case-sensitive and may appear in any order. Unrecognized columns are ignored.\n"
		 << " - SampleID values must be unique. Barcode values, or dual-barcode pairs, must be unique.\n\n";
	cout << "Sample identification modes:\n"
		 << " - BarcodeSequence: identify samples from a read barcode (IUPAC bases are supported).\n"
		 << " - BarcodeSequence plus Barcode2ndPair: identify samples using dual indexes.\n"
		 << " - SampleIDinHead: identify samples from a matching string in the FASTA/FASTQ header.\n"
		 << " - In -i_path mode, unique fastqFile, alignmentFile or fnaFile entries can identify one sample per file.\n"
		 << " - -htsBarcode BC: match BarcodeSequence (and Barcode2ndPair) to BC:Z tags exactly.\n"
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
		 << "  fastqFile            Relative FASTQ path; comma-separated R1,R2 for paired reads.\n"
		 << "  alignmentFile        One SAM/BAM/CRAM path; leave fastqFile empty for this sample.\n"
		 << "  fnaFile              Relative FASTA path.\n"
		 << "  qualFile             Relative FASTA-quality path; if absent, sdm derives <fasta-base>.qual_.\n"
		 << "  MIDfqFile            Relative FASTQ path containing MID sequences.\n"
		 << " Dereplication\n"
		 << "  derepMin             Empty or a non-negative minimum copy count for this sample.\n\n";
	cout << "When using -i_path, the mapping file must contain fastqFile, alignmentFile or fnaFile. File paths are resolved relative to -i_path.\n\n";
	cout << "With alignmentFile present, each row must fill exactly one of fastqFile or alignmentFile. Both columns may appear in a mixed FASTQ/alignment map.\n";
	cout << "Alignment rows cannot use MIDfqFile or qualFile. Use a separate run for FASTA inputs.\n";
	cout << "Paired alignment files must be query-name sorted; singleton files may have any order. Use -cramRef for an external CRAM reference.\n\n";
	cout << "Preserve empty cells, including a final tab when the last column is blank.\n"
         << "The FASTQ column is fastqFile (not fqFile); alignment-only maps may omit it.\n"
         << "Mixed FASTQ/alignment runs need the same singleton/paired layout throughout.\n"
         << "Use indexed reference.fa with -cramRef for externally referenced CRAM.\n"
         << "BC:Z:ACGT matches BarcodeSequence=ACGT; BC:Z:ACGT-TGCA matches that barcode\n"
         << "plus Barcode2ndPair=TGCA. Matching is exact and does not trim the biological read.\n"
         << "A tag on one mate applies to both; missing, invalid or conflicting tags are errors.\n"
         << "Repeat an alignment path across barcode rows to demultiplex that file once.\n"
         << "See docs/mapping.md and docs/hts.md for complete examples.\n\n";
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
		 << "sample_1     sample_1.fastq.gz\n\n"
		 << "# One sample per SAM/BAM/CRAM file with -i_path reads\n"
		 << "#SampleID    alignmentFile\n"
		 << "sample_1     sample_1.bam\n\n"
         << "# Shared BAM with -htsBarcode BC\n"
         << "#SampleID    alignmentFile    BarcodeSequence\n"
         << "sample_1     shared.bam       ACGT\n"
         << "sample_2     shared.bam       TTAA\n";

}
void printVersion(){
    cout << "sdm " << sdm_version << " " << sdm_status << endl;
	cout << "compiled on: " << __DATE__ << " at " << __TIME__ << endl;
	cout << "compiled with C++ v: " << __cplusplus << endl;
    printBuildFeatures();

}

