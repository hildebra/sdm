#include "Filters.h"
#include "OutputStreamer.h"



//*******************************************
//*        FILTERS OBJECT
//*******************************************


Filters::Filters(OptContainer* cmdArgs1) :
	PrimerL(0), PrimerR(0), PrimerL_RC(0), PrimerR_RC(0), PrimerIdx(0),
	Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
	HeadSmplID(0),
	hetPrimer(2, vector<string>(0)),
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



	bool alt_bRequireRevPrimSet = false;

	string optF("");
	if (cmdArgs->find("-options") != cmdArgs->end()) {
		optF = (*cmdArgs)["-options"];
	}

	iniSpacer = (*cmdArgs)["-sample_sep"];
	//***************************************
	//default options
	int maxAmb(0), PrimerErrs(1), TagErrs(0);
	float minQual(25);
	float minL(250.f);
	int maxL(1000);
	int QualWinWidth = 50;
	float QualWinThr = 0;
	int EndWinWidth = 15;
	float EndWinThr = 20;
	int maxHomoNT(12); int trimHomoNT(12);
	bool keepTag(false), keepPrimer(false);
	bool addModConf = false;

	//set up some basic objects
	if ((*cmdArgs).find("-paired") != cmdArgs->end()) {
		pairedSeq = atoi((*cmdArgs)["-paired"].c_str()); //fakeEssentials();
		if (pairedSeq < 1 || pairedSeq>3) { cerr << "Argument \"-paired\" supplied with unknown parameter. Aborting.\n"; exit(28); }
		if ((*cmdArgs)["-onlyPair"] == "1" || (*cmdArgs)["-onlyPair"] == "2") {
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
		demultiBPperSR = (uint)d;
	}
	//important for fastq format
	if (cmdArgs->find("-i_qual_offset") != cmdArgs->end()) {
		if ((*cmdArgs)["-i_qual_offset"] == "auto") {
			userReqFastqVer = 0;
		}
		else {
			userReqFastqVer = atoi((*cmdArgs)["-i_qual_offset"].c_str());
		}
	}

	cut5PR1 = atoi((*cmdArgs)["-5PR1cut"].c_str());
	cut5PR2 = atoi((*cmdArgs)["-5PR2cut"].c_str());
	//cerr<<(*cmdArgs)["-o_qual_offset"]<<endl;
	userReqFastqOutVer = atoi((*cmdArgs)["-o_qual_offset"].c_str());
	//statistic tracker
		//do new SEED sequence selection?
	if (cmdArgs->find("-optimalRead2Cluster") != cmdArgs->end()) {
		b_optiClusterSeq = true;
	}
	//do selection of specific reads?
	if ((*cmdArgs)["-specificReads"] != "") {
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
	opt.open(optF.c_str(), ios::in);
	if (!opt || optF == "") {
		cerr << "NO filtering will be done on your reads (just rewriting / log files created)." << endl;
		b_doFilter = false;
		return;
	}
	string line;
	while (getline(opt, line, '\n')) {

		if (line.length() <= 1 || line.substr(0, 1) == "#") {
			continue;
		}

		bool addMod = false;
		if (line.substr(0, 1) == "*") {
			addMod = true;
			line = line.substr(1);
		}
		string segs;
		string segs2;
		stringstream ss;
		ss << line;
		getline(ss, segs, '\t');
		getline(ss, segs2, '\t');

		if (strcmp(segs.c_str(), "minSeqLength") == 0) {
			if (addMod) {
				float tmp = (float)atof(segs2.c_str());
				if (tmp > 1.f) {
					alt_min_l = (int)tmp;
				}
				else {
					alt_min_l = -1;
					alt_min_l_p = tmp;
				}
				if (alt_min_l != minL) { addModConf = true; }
			}
			else {
				minL = (float)atof(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "maxSeqLength") == 0) {
			maxL = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "minAvgQuality") == 0) {
			if (addMod) {
				alt_min_q = (float)atof(segs2.c_str());
				if (alt_min_q != minQual) { addModConf = true; }
			}
			else {
				minQual = (float)atof(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "maxAmbiguousNT") == 0) {
			if (addMod) {
				alt_MaxAmb = atoi(segs2.c_str());
				if (MaxAmb != alt_MaxAmb) { addModConf = true; }
			}
			else {
				maxAmb = atoi(segs2.c_str());

			}
		}
		else if (strcmp(segs.c_str(), "QualWindowThreshhold") == 0) {
			if (addMod) {
				alt_FQWthr = (float)atof(segs2.c_str());
				if (alt_FQWthr != QualWinThr) { addModConf = true; }
			}
			else {
				QualWinThr = (float)atof(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "QualWindowWidth") == 0) {
			QualWinWidth = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "BinErrorModelMaxExpError") == 0) {
			BinFilErr = (float)atof(segs2.c_str());
			if (BinFilErr < 0) {
				cerr << "BinErrorModelMaxExpError was set to <0. Set to 0 instead.\n";
				BinFilErr = 0;
			}
		}
		else if (strcmp(segs.c_str(), "BinErrorModelAlpha") == 0) {
			BinFilP = (float)atof(segs2.c_str());
			if (BinFilP != -1.f && (BinFilP < 0.f || BinFilP>1.f)) {
				cerr << "BinErrorModelAlpha has to be between 0 and 1 (or -1 to deactivate).\nAborting..\n";
				exit(542);
			}

		}
		else if (strcmp(segs.c_str(), "TrimWindowWidth") == 0) {
			EndWinWidth = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "TrimWindowThreshhold") == 0) {
			if (addMod) {
				alt_EWthr = (float)atof(segs2.c_str());
				if (alt_EWthr != EndWinThr) { addModConf = true; }
			}
			else {
				EndWinThr = (float)atof(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "maxBarcodeErrs") == 0) {
			TagErrs = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "maxPrimerErrs") == 0) {
			if (addMod) {
				alt_PrimerErrs = atoi(segs2.c_str());
				if (alt_PrimerErrs != PrimerErrs) { addModConf = true; }
			}
			else {
				PrimerErrs = atoi(segs2.c_str());
			}
		}
		else if (strcmp(segs.c_str(), "keepBarcodeSeq") == 0) {
			atoi(segs2.c_str()) == 0 ? keepTag = false : keepTag = true;
		}
		else if (strcmp(segs.c_str(), "keepPrimerSeq") == 0) {
			if (addMod) {
				atoi(segs2.c_str()) == 0 ? alt_BcutPrimer = false : alt_BcutPrimer = true;
				if (alt_BcutPrimer != keepPrimer) { addModConf = true; }
			}
			else {
				atoi(segs2.c_str()) == 0 ? keepPrimer = false : keepPrimer = true;
			}
		}
		else if (strcmp(segs.c_str(), "maxHomonucleotide") == 0) {
			maxHomoNT = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "trimHomonucleotide") == 0) {
			trimHomoNT = atoi(segs2.c_str());
		}
		else if (strcmp(segs.c_str(), "maxAccumulatedError") == 0) {
			if (addMod) {
				alt_maxAccumQP = double(atof(segs2.c_str()));
				if (alt_maxAccumQP != maxAccumQP) { addModConf = true; }
			}
			else {
				maxAccumQP = double(atof(segs2.c_str()));
			}
		}
		else if (strcmp(segs.c_str(), "TechnicalAdapter") == 0) {
			tAdapter = segs2.c_str();
			transform(tAdapter.begin(), tAdapter.end(), tAdapter.begin(), ::toupper);
			tAdapterLength = (int)tAdapter.length();
			removeAdapter = true;
		}
		else if (segs == "PEheaderPairFmt") {
			PEheaderVerWr = atoi(segs2.c_str());
		}
		else if (segs == "TrimStartNTs") {
			TrimStartNTs = atoi(segs2.c_str());
		}
		else if (segs == "fastqVersion") {
			if (segs2 == "auto") {
				userReqFastqVer = 0;
			}
			else {
				userReqFastqVer = FastqVerMod(atoi(segs2.c_str()));
			}
		}
		else if (segs == "ExtensivePrimerChecks") {
			if (segs2 == "T") {
				BextensivePrimerChecks = true;
			}
		}
		else if (segs == "RejectSeqWithoutRevPrim") {
			if (addMod) {
				alt_bRequireRevPrimSet = true;
				if (segs2 == "T") {
					alt_bRequireRevPrim = true;
				}
				else { alt_bRequireRevPrim = false; }
				if (alt_bRequireRevPrim != bRequireRevPrim) { addModConf = true; }
			}
			else {
				if (segs2 == "T") {
					bRequireRevPrim = true;
				}
				else { bRequireRevPrim = false; }
			}
		}
		else if (segs == "RejectSeqWithoutFwdPrim") {
			if (addMod) {
				alt_bRequireFwdPrim = true;
				if (segs2 == "T") {
					alt_bRequireFwdPrim = true;
				}
				else { alt_bRequireFwdPrim = false; }
				if (alt_bRequireFwdPrim != bRequireFwdPrim) { addModConf = true; }
			}
			else {
				if (segs2 == "T") {
					bRequireFwdPrim = true;
				}
				else { bRequireFwdPrim = false; }
			}
		}
		else if (segs == "TruncateSequenceLength") {
			TruncSeq = atoi(segs2.c_str());
			if (TruncSeq != -1 && TruncSeq < (int)minL) { minL = (float)TruncSeq; }
		}
		else if (segs == "AmpliconShortPE") {
			if (segs2 == "T") {
				bShortAmplicons = true;
			}
			else { bShortAmplicons = false; }
		}
		else if (segs == "CheckForMixedPairs") {
			if (segs2 == "T") {
				b2ndRDBcPrimCk = true;
			}
			else { b2ndRDBcPrimCk = false; }
		}
		else if (segs == "CheckForReversedSeqs") {
			if (segs2 == "T") {
				bRevRdCk = true;
			}
			else { bRevRdCk = false; }
		}
		else if (segs == "SyncReadPairs") {
			if (segs2 == "T") {
				bChkRdPrs = true;
			}
			else { bChkRdPrs = false; }
		}
		else if (segs == "illuminaFwd") {
			illuPEfwd = segs2;
		}
		else if (segs == "illuminaRev") {
			illuPErev = segs2;
		}
		else if (segs == "illuminaSngUni") {
			illuSEuni = segs2;
		}
		else if (segs == "illuminaSngIdx") {
			illuSEidx = segs2;
		}
	}

	//report some non-std options
	if (bShortAmplicons) {
		cerr << "Checking for reverse primers on 1st read.\n";
	}
	if (b2ndRDBcPrimCk) {
		cerr << "Checking for switched pairs.\n";
	}

	opt.close();
	//set in filter object
	this->setSeqLength(minL, maxL);
	this->setPrimerErrs(PrimerErrs);
	this->setTagErrs(TagErrs);
	this->removePrimer(!keepPrimer);
	this->removeTag(!keepTag);
	this->setMaxAmb(maxAmb);
	this->setAvgMinQual(minQual);
	this->setFloatingQWin(QualWinWidth, QualWinThr);
	this->setFloatingEWin(EndWinWidth, EndWinThr);
	this->setMaxHomo(maxHomoNT);
	this->setTrimHomo(trimHomoNT);

	//alternative options (mid qual filtering)
	if (addModConf) {
		if (!alt_bRequireRevPrimSet) { alt_bRequireRevPrim = bRequireRevPrim; }
		if (cmdArgs->find("-o_fna") != cmdArgs->end() && (*cmdArgs)["-o_fna"].length() > 1) {
			if (cmdArgs->find("-o_fna2") == cmdArgs->end()) {
				(*cmdArgs)["-o_fna2"] = additionalFileName((*cmdArgs)["-o_fna"]);
				//(*cmdArgs)["-o_fna2"] = (*cmdArgs)["-o_fna"].substr(0,(*cmdArgs)["-o_fna"].length()-4)+".add.fna";
			}
		}
		else if (cmdArgs->find("-o_fastq") != cmdArgs->end() && (*cmdArgs)["-o_fastq"].length() > 1) {
			if (cmdArgs->find("-o_fastq2") == cmdArgs->end()) {
				(*cmdArgs)["-o_fastq2"] = additionalFileName((*cmdArgs)["-o_fastq"]);
			}
		}
		bAdditionalOutput = true;
	}
}


void Filters::addDNAtoCStats(const shared_ptr<DNA>& d, int Pair) {
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
	}
	else {
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

    if (d->isGreenQual() || d->isYellowQual()) {
		countBCdetected(d->getBCnumber(), easyPair, d->isYellowQual());
		//and register as success
	}
	else {
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



Filters::Filters(Filters* of, int BCnumber, bool takeAll, size_t threads)
	:
	PrimerL(0, ""), PrimerR(0, ""),
	PrimerL_RC(0, ""), PrimerR_RC(0, ""),
	PrimerIdx(of->PrimerIdx),
	Barcode(0), revBarcode(0), Barcode2(0), revBarcode2(0),
	HeadSmplID(0), hetPrimer(2, vector<string>(0)),
	collectStatistics(2, nullptr), statAddition(2, nullptr),
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

	SequencingRun(0), cmdArgs(of->cmdArgs), passed_interval_reads(0)
{
	cdbg("New Filter object from copy\n");
	ReadsWritten = of->writtenReads();
	OFileIncre = of->getFileIncrementor();
	BCdFWDREV[0].reset(); BCdFWDREV[1].reset();
	//collectStatistics.resize(2); statAddition.resize(2);
	//csMTX[0] = new mutex(); csMTX[1] = new mutex();

	collectStatistics[0] = make_shared<collectstats>();
	collectStatistics[1] = make_shared<collectstats>();
	collectStatistics[0]->setbLvsQlogsPreFilt(of->collectStatistics[0]->getbLvsQlogsPreFilt());
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
	//cdbg("Done\n");
}


Filters* Filters::newFilterPerBCgroup(const vector<int> idxi) {

	if (idxi.size() < 1) {
		return nullptr;
	}
	cdbg("newFilterPerBCgroup::start : " + itos(idxi[0]) + "\n");

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
		if (!d->getFwdPrimDetect() && bRequireFwdPrim) {
			d->setYellowQual(true);
		}
		//conditions for failing read on not finding rev primer
		if (d2 != nullptr && bRequireRevPrim && !d->getRevPrimDetect()) {//failed to find reverse primer
			d->setYellowQual(true);
		}
	}

}

//service function to ini OTU Seed extension 
UClinks* Filters::ini_SeedsReadsDerep(UClinks* ucl, shared_ptr<ReadSubset>& RDSset,
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
				this->getuserReqFastqVer(), (*cmdArgs)["-ignore_IO_errors"], (*cmdArgs)["-pairedRD_HD_out"], 1);
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
bool Filters::setcmdArgsFiles() {

	if (FastqF.size() == 0 && QualF.size() == 0 && FastaF.size() > 0) {
		//fasta entry but no qual entries

		string path = "";
		if (cmdArgs->find("-i_path") != cmdArgs->end() && (*cmdArgs)["-i_path"].length() >= 1) {
			path = (*cmdArgs)["-i_path"] + string("/");
		}

		QualF.resize(FastaF.size());
		for (unsigned int i = 0; i < FastaF.size(); i++) {
			string newQ = FastaF[i];
			int pos = (int)newQ.find_last_of(".");
			newQ = newQ.substr(0, pos);
			newQ += string(".qual_");
			fstream fin;
			string fullQ = path + newQ;
			fin.open(fullQ.c_str(), ios::in);
			if (fin.is_open()) {
				cerr << "Using quality file: " << fullQ << endl;
			}
			else if (cmdArgs->find("-number") != cmdArgs->end() && (*cmdArgs)["-number"] == "T") {
				;
			}
			else {
				cerr << "You did not supply a quality file for" << path + FastaF[i] << ". \nPlease give the path to your quality file as command line argument:\n  -i_qual <PathToQualityFile>\n";
				newQ = "";
				//fin.close();return false;
			}
			fin.close();
			QualF[i] = newQ;
		}
	}

	int fileSiz = (int)Barcode.size();
	//instead of max:
	if (!bDoMultiplexing) { fileSiz = 1; }

	if (FastaF.size() == 0 && FastqF.size() == 0) {
		//set up fasta/fastq vector specific to corresponding BC (that should be in this file)
		if (cmdArgs->find("-i_fastq") == cmdArgs->end()) {
			FastaF.resize(fileSiz);
			QualF.resize(fileSiz);
			for (unsigned int i = 0; i < FastaF.size(); i++) {
				FastaF[i] = (*cmdArgs)["-i_fna"];
				QualF[i] = (*cmdArgs)["-i_qual"];
			}
		}
		else {//fastq input
			vector<string> fqTmp(1, (*cmdArgs)["-i_fastq"]);
			if ((*cmdArgs)["-i_fastq"].find(";") != string::npos) {//";" denotes several files
				if (fileSiz == 1) {//no BC, 
					fqTmp = splitByCommas((*cmdArgs)["-i_fastq"], ';');
					this->allResize((uint)fqTmp.size());
					fileSiz = (int)fqTmp.size();
					cerr << "Detected " << fileSiz << " input files (pairs)." << endl;
					FastqF = fqTmp;
				}
				else {
					cerr << "Fastq string contains symbol \";\". Not allowed in input string"; exit(32);
				}
			}
			else {
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
void Filters::reverseTS_all_BC() {
	//	for (int i=0; i<Barcode.size();i++){
	//		reverseTS(Barcode[i]);
	//	}
	Barcode = revBarcode;
	revBarcode.resize(0);
	barcodes1_.clear();
	for (uint i = 0; i < Barcode.size(); i++) {
		barcodes1_[Barcode[i]] = i;
		barcodeLengths1_[i] = (int)Barcode[i].length();
	}

}
void Filters::reverseTS_all_BC2() {
	//	for (int i=0; i<Barcode.size();i++){
	//		reverseTS(Barcode[i]);
	//	}
	Barcode2 = revBarcode2;
	revBarcode2.resize(0);
	barcodes2_.clear();
	for (uint i = 0; i < Barcode2.size(); i++) {
		barcodes2_[Barcode2[i]] = i;
		barcodeLengths2_[i] = (int)Barcode2[i].length();
	}
}

bool Filters::isReversedAmplicon(shared_ptr<DNA> tdn) {
	if ((!checkRevRd())) {
		return false;
	}


	//BC should be already cut at this point..
	//int tagIdx = 0;//just try

	//method 1: just check if primer is found reversed, most basic and seems to work fine..
	//simple check if fwd rev primer is in reverse position: then reverse transcribe

	bool fwdRd1Primer = checkIfPrimerHits(tdn, 0, 0);
	bool fwdRd1PrimerRev = checkIfRevPrimerHits(tdn, 0, 0);
	if (fwdRd1Primer) {
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
		idx = detectCutBC(tdn[0], true);
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
		if (limitF2 > limitR) {
			limitF = limitF2;
		}
		else {
			limitF = dn->matchSeq(PrimerL[PrimerIdx[idx]], PrimerErrs, SearchL + limitR, limitR);
		}
		if (limitF < 0) { break; }
		limitF2 = dn->matchSeq(PrimerL[PrimerIdx[idx]], PrimerErrs, SearchL + limitF + 10, limitF + 10);
		//records reverse-searched primer positions
		limitR = dn->matchSeq(PrimerR_RC[PrimerIdx[idx]], PrimerErrs, SearchL + limitF, limitF + 1);

		if (limitF2 > 0 && limitF2 < limitR) {//something went wrong.. couldn't detect correct reverse primer?
			limitF2 = -1;
			isProblem.push_back(true);
		}
		else {
			isProblem.push_back(false);
		}

		if (limitR < 0) { break; }
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
            dn->getDNAsubseq((int)posF[i], (int)(posR[i] + PrimerR_RC[PrimerIdx[idx]].length()),
                dn->getId() + "_" + itos((int)i))
		);
		BCintoHead(idx, retDNA.back(), "", -1, false, true);

		if (i > 0) {
			int disGAs = posF[i] - posR[i - 1] - (int)PrimerR[PrimerIdx[idx]].length();
			if (disGAs > 12) {
				//cerr << ("disGA:"+itos(disGA));
				//std::cout << "disGA" << " ";
				missedGAs++;
			}
		}
		else if (posF[i] > 100) {
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
		(GoldAxeMaxAmpli != -1 && retDNA.size() > GoldAxeMaxAmpli)) {
		tdn[0]->failed();
		retDNA.resize(0);
	}
	//this->addGAstats(dn, retDNA);
	GAstatistics->addBaseGAStats(dn, retDNA, missedGAs);
	//GAstatistics->addMissedGAs(missedGAs);

	return retDNA;

}

bool Filters::swapReverseDNApairs(vector< shared_ptr<DNA>>& tdn) {
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
	if (checkIfRevPrimerHits(tdn[0], 0, 0)) {
		tdn[0]->reverse_compliment();
		collectStatistics[0]->reversedRds++;
		if (tdn[1] != nullptr) {
			tdn[1]->reverse_compliment();
			//take stats on this
			collectStatistics[1]->reversedRds++;
			return true;
		}
	}
	else if (checkSwitchedRdPairs() && isPaired() == 2 && tdn[1] != nullptr) {
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
		}
		else if (checkIfPrimerHits(tdn[1], 0, 0)) {//switched pairs only
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

void Filters::preFilterSeqStat(const shared_ptr<DNA>& d, int pair) {
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
			min_l = (int)((float)maxReadLength * min_l_p);
		}
	}
}
void Filters::setSeqLength(float minL, int maxL) {
	if (minL > 1.f) {
		min_l = (int)minL;
		min_l_p = -1.f;
	}
	else {
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
	shared_ptr<DNAunique> ref, float ever_best,
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

	return whoIsBetter(d1, d1->getPair(), d1->getMerge(),
		ref, ref->getPair(), ref->getMerge(), ever_best, true);

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
bool Filters::checkYellowAndGreen(const shared_ptr<DNA>& d, int pairPre,
	int& tagIdx, bool doSeeding) {
	int pair = max(0, pairPre);//corrects for -1 (undefined pair_) to set to 0

	// We trim for poly-G / homonucleotide tail here, before filtering for length, and before trimming adapter sequences etc. 
	if (trimHomonucleotide != 0) {
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
	if (tagIdx == -2) {
		if (bDoBarcode2 && pair == 1) {
			tagIdx = detectCutBC(d, false); //barcode 2nd part
		}
		else if (bDoBarcode && pair == 0) {
			tagIdx = detectCutBC(d, true); //barcode
		}
		else {
			tagIdx = 0;
		}
	}
	if ((bDoBarcode || bDoBarcode2) && tagIdx < 0) {
		d->QualCtrl.TagFail = true;
		return false;
	}

	if ((BextensivePrimerChecks || BcutPrimer || Bcheck4illuAdapts) &&
		(tagIdx < 0 || (size_t)tagIdx >= PrimerIdx.size() || (size_t)tagIdx >= PrimerIdxRev.size())) {
		d->QualCtrl.TagFail = true;
		d->failed();
		return false;
	}

	if (BextensivePrimerChecks) {
		if (checkIfRevPrimerHits(d, PrimerIdx[tagIdx], 0, bShortAmplicons)) {
			d->reverse_compliment(false);
		}
	}
	//bShortAmplicons checks for reverse primer on 1st read
	if (BcutPrimer || Bcheck4illuAdapts) {
       if (pair == 1) {//pair_ == 1, check for fwd primer in pair_ 2 (rev-compl)
			bool revCheck = false;
			cutPrimerRev(d, PrimerIdxRev[tagIdx], revCheck, false);
			if (bShortAmplicons) {//also check other end of primer..
				cutPrimer(d, PrimerIdx[tagIdx], revCheck, pair);
			}
		}
		else {//pair 0 or -1
			bool fwdRC = false;
			cutPrimer(d, PrimerIdx[tagIdx], fwdRC, pair, BextensivePrimerChecks);
			if (bShortAmplicons) {//also check other end of primer.. also use for PacBio amplicons
               bool revRC = true;
				cutPrimerRev(d, PrimerIdxRev[tagIdx], revRC, BextensivePrimerChecks);
				//case for 1) long read 2) look for both primers 3) rev required
			}
		}
		//conditions for failing read on not finding fwd primer
		if (pair != 1 && !d->getFwdPrimDetect() && bRequireFwdPrim) {
			if (alt_bRequireFwdPrim) {
				d->QualCtrl.PrimerFwdFail = true;
				d->failed(); return false;
			}
			d->setYellowQual(true);
		}
		//conditions for failing read on not finding rev primer
		if (bRequireRevPrim && !d->getRevPrimDetect()
			&& (pair != 0 || isPaired() == 1)  //only pair2 or 1-read-PacBio will fail
			) {//failed to find reverse primer
			if (alt_bRequireRevPrim) {
				d->QualCtrl.PrimerRevFail = true;
				d->failed(); return false;
			}
			else {
				d->setYellowQual(true);
			}
		}

	}


	if (doSeeding) {
		//cut off low qual, hard limits
		d->qualWinPos(EWwidth, EWthr);
		return true;
	}

	//if seq needs to be cut, than here
	if (TruncSeq > 0) {
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
				// Mark as mid-quality (yellow) when trimming due to accumulated error
				d->setYellowQual(true);  d->QualCtrl.minLqualTrim = true;//sMinQTrim(pair_);
				cP = d->qualAccumulate(alt_maxAccumQP);
				if (check_lengthXtra(d, 0, cP)) {//check if passes alt
					d->failed(); return false;
				}
			}
			else {
				d->qualAccumTrim(maxAccumQP);//) {  AccErrTrim = true; }
			}
		}


        const bool usePrimaryQualWin = (min_q > 0 || FQWthr > 0);
		const bool useAltQualWin = (alt_min_q > 0 || alt_FQWthr > 0);
		int rea(2), rea2(2);
		float avgQ(min_q);
		float avgQalt(alt_min_q);
		if (usePrimaryQualWin) {
			avgQ = d->qualWinfloat(FQWwidth, FQWthr, rea);
		}
		if ((avgQ < min_q || rea == 1) && useAltQualWin) {
			avgQalt = d->qualWinfloat(FQWwidth, alt_FQWthr, rea2);
		}
		if (avgQ < min_q) {
			d->QualCtrl.AvgQual = true; //sAvgQual(pair_);
			if (avgQalt < alt_min_q) {
				d->QualCtrl.AvgQual = true; //statAddition.AvgQual++;
				d->failed(); return false;
			}
			else {
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
		if (b_BinFilBothPairs || pair == 0) {
			float ExpErr = d->binomialFilter((int)BinFilErr, BinFilP);
			if (ExpErr > BinFilErr) {
				d->QualCtrl.BinomialErr = true;
				d->failed(); return false;
			}
		}
	}
	// Fused single-pass: ambiguity count + homonucleotide run detection
	if (MaxAmb >= 0 || maxHomonucleotide != 0) {
		int ambNTs = 0;
		bool homoOK = d->scanSequenceChecks(maxHomonucleotide, ambNTs);
		if (!homoOK) {
			d->QualCtrl.HomoNT = true;
			d->failed(); return false;
		}
		if (MaxAmb >= 0 && ambNTs > MaxAmb) {
			d->QualCtrl.MaxAmb = true;
			if (alt_MaxAmb != -1 && ambNTs >= alt_MaxAmb) {
				d->QualCtrl.MaxAmb = true;
				d->failed(); return false;
			}
			else {
				d->setYellowQual(true);
			}
		}
	}

	//adapter removed, quality filtering done. If no map is provided, that is all that is needed
	if (!bDoMultiplexing) {
		if (TrimStartNTs > 0) {
           const size_t readLen = d->length();
			if (readLen <= (unsigned int)TrimStartNTs) {
				d->QualCtrl.minL = true;
				d->failed();
				return false;
			}
         if (readLen - TrimStartNTs > max_l) {//length check
				d->QualCtrl.maxL = true; //sMaxLength(pair_);
				d->failed();
				return false;
			}
			//remove start NTs
			d->cutSeq(0, TrimStartNTs);

		}
		d->setPassed(true);
		return true;
	}

	if (!d->isYellowQual()) {
		d->setPassed(true);
	}

	return true;
}
void Filters::noMapMode() {
	string noMapTxt = "sdm run in No Map Mode.";
	if (cmdArgs->find("-paired") != cmdArgs->end() && ((*cmdArgs)["-paired"] == "2" || (*cmdArgs)["-paired"] == "2")) {
		pairedSeq = 2; //fakeEssentials();
		noMapTxt += " Using paired end sequencing files.";
	}

	BcutPrimer = false; bDoBarcode = false; bDoBarcode2 = false; bDoBarcode2Rd1 = false;
	removeAdapter = false; bDoMultiplexing = false;
	bDoHeadSmplID = false;
	minBCLength1_ = 0; minBCLength2_ = 0; maxBCLength1_ = 0; maxBCLength2_ = 0; minPrimerLength_ = 0;
	cerr << noMapTxt << endl;

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
void Filters::fakeEssentials(bool all) {
	//create fake entries
	Barcode.push_back("NA");
	barcodeLengths1_.push_back(0);
	barcodeLengths2_.push_back(0);
	SequencingRun.push_back("");
	SequencingRun2id[""] = vector<int>(1, 0);
	SampleID.push_back("NA"); SampleID_Combi.push_back("NA");
	HeadSmplID.push_back(""); bDoHeadSmplID = false;
	collectStatistics[0]->BarcodeDetected.push_back(-1);
	collectStatistics[1]->BarcodeDetected.push_back(-1);
	collectStatistics[0]->BarcodeDetectedFail.push_back(-1);
	collectStatistics[1]->BarcodeDetectedFail.push_back(-1);
	if (all) {
		PrimerIdx.push_back(0); PrimerL.push_back(""); PrimerL_RC.push_back("");
	}

}
void Filters::allResize(unsigned int x) {
	//cerr<<"resize "<<x<<endl;
	PrimerIdx.resize(x, 0);
	PrimerIdxRev.resize(x, 0);
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

bool Filters::remove_adapter(const shared_ptr<DNA>& d) { //technical adapter
	//allows for 0 errors, no shifts
	if (d->getTA_cut()) {
		return true;
	}
	const string& se = d->getSequence();
 if (se.size() < tAdapterLength) {
		return false;
	}
	for (unsigned int i = 0; i < tAdapterLength; i++) {
		if (se[i] != tAdapter[i]) {
			return false;
		}
	}
	d->cutSeq(0, tAdapterLength);
	d->setTA_cut(true);
	return true;
}
//only identifies based on dual BCding
void Filters::dblBCeval(int& tagIdx, int& tagIdx2, string& presentBC,
	shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2) {
	//bool BCfail = false;// , BCfail2 = false;

	if (tagIdx < 0 || tagIdx2 < 0 || !tdn->getBarcodeDetected() ||
		(tdn2 != nullptr && !tdn2->getBarcodeDetected())) {
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != nullptr) {
			tdn->setPassed(false); /*BCfail = true; */
			tdn->setBCnumber(tagIdx, BCoffset); tdn->setYellowQual(false);
		}
		if (tdn2 != nullptr) { tdn2->setPassed(false); tdn2->setYellowQual(false); tdn2->setBCnumber(tagIdx2, BCoffset); }

		collectStatistics[0]->dblTagFail++;
		return;
	}
	string BC1 = Barcode[tagIdx];
	string BC2 = Barcode2[tagIdx2];
	bool hit(false);
	if (tagIdx == tagIdx2) {
		hit = true;
	}
	else {
		//this routine finds two matching barcodes (as several combinations are possible)
		for (uint i = 0; i < Barcode.size(); i++) {
			if (Barcode[i] == BC1 && Barcode2[i] == BC2) {
				tagIdx = i; tagIdx2 = i; hit = true; break;
			}
		}
	}

	if (!hit) {
		//no BC, useless
		tagIdx = -1; tagIdx2 = -1;
		if (tdn != nullptr) { tdn->setPassed(false); tdn->setYellowQual(false); tdn->setBCnumber(tagIdx, BCoffset); }
		if (tdn2 != nullptr) { tdn2->setPassed(false); tdn2->setYellowQual(false); tdn2->setBCnumber(tagIdx2, BCoffset); }
		return;
	}
	presentBC = BC1 + "|" + BC2;
	//add new BC info to DNA
	//also reset BC in DNA
	if (tdn != NULL) {
		BCintoHead(tagIdx, tdn, presentBC, -1, false, true);
		//already done in BCintoHead
	}
	if (tdn2 != NULL) {
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
int Filters::findTag(const shared_ptr<DNA>& d, string& presentBC, int& c_err,
	bool isPair1, int revChecks, bool cutBC, bool endCheck) {

	//cout << "FIND TAG DO HEAD " << bDoHeadSmplID << endl;
	if (bDoHeadSmplID) {
		for (unsigned int i = 0; i < HeadSmplID.size(); i++) {
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

	scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1, false, endCheck);

	if (!BCdFWDREV[!isPair1].b_BCdirFix) {
		if (start == -1) {//check reverse transcription
			//d->reverseTranscribe();
			scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1, true, endCheck);
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
		scanBC(d, start, stop, idx, c_err, scanRegion, presentBC, isPair1, true, endCheck);
		if (idx >= 0) {
			revChecks = 0;
		}
	}
	if (start == -1) {
		idx = -1;
	}
	if (idx != -1) {
		if (cutBC && BcutTag && !d->isMIDseq()) {
			//remove tag from DNA
			if (endCheck) {
				d->cutSeq(start, -1);
			}
			else {
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
int Filters::detectCutBC(const shared_ptr<DNA>& d, bool isPair1) {
	//seq too short for BC
	if (d->length() < minBCLength1_) {
		return -1;
	}
	//already detected barcode
	if (d->getBarcodeCut()) {// && !scndBC) {
		return d->getBCnumber() - BCoffset;
	}
	if ((isPair1 && !bDoBarcode) || (!isPair1 && !bDoBarcode2)) {
		d->setBCnumber(0, BCoffset);
		return BCoffset; //not failed, just not requested
	}

	//ok, really start looking for BC in seq
	int idx(-1);
	if (bDoHeadSmplID) {
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
	}
	else if (bOneFileSample) {
		// if theres only one sample per file then there is no need for demultiplexing
		d->setBCnumber(0, this->currentBCnumber());

		// needs to be locked when multithreading
		BCintoHead(0, d, "FileName", isPair1, false);

		return 0;
	}




	//int start(-1);int stop(-1);
	string presentBC; int c_err(0);
	bool useBC1 = isPair1; bool useBC2 = !useBC1;
	idx = findTag(d, presentBC, c_err, useBC1, 0, true, false);

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
		if (idx != idxX) {
			cerr << "(2) Unequal BC numbers:" << idx << " : " << idxX << "; in object: " << d->getBCnumber() << endl;
			cerr << "In read:" << d->getId() << endl;
			exit(835);
		}

		//cout << "scanBCrev: " << start << "," << stop << "," << idx << "," << c_err << "," << presentBC << endl;
	}


	d->setBCnumber(idx, BCoffset);

	return idx;
}

void Filters::BCintoHead(int idx, const shared_ptr<DNA>& d, std::string_view presentBC,
	const int c_err, bool pair1, bool atEnd) {
	vector<string>* locBC = (!pair1) ? &Barcode2 : &Barcode;
	const string& oldId = d->getId();
	const size_t idStop = atEnd ? oldId.size() : oldId.find_first_of(' ', 0);
	const size_t idLen = idStop == string::npos ? oldId.size() : idStop;

	const string& samplePrefix = bDoCombiSamples ? SampleID_Combi[idx] : SampleID[idx];
	string nID;
	nID.reserve(samplePrefix.size() + iniSpacer.size() + idLen + (*locBC)[idx].size() + presentBC.size() + 40);
	nID.append(samplePrefix);
	nID.append(iniSpacer);
	nID.append(oldId.data(), idLen);
	nID += " orig_bc=";
	nID += (*locBC)[idx];
	if (!presentBC.empty() && (c_err > 0 || atEnd)) {//atEnd: dbl barcode
		nID += " new_bc=";
		nID.append(presentBC.data(), presentBC.size());
		nID += " bc_diffs=";
		nID += std::to_string(c_err);
	}
	d->setNewID(nID);
	d->setBCnumber(idx, BCoffset);
}

void Filters::SampleIntoHead(const int idx, shared_ptr<DNA> d, const size_t pos) {
    std::string_view on = d->getId();
	size_t pos2 = on.find_first_of(' ', pos);
	// on2 = on.substr(0, pos) + on.substr(pos2 + 1)
	size_t part2_start = (pos2 == std::string_view::npos) ? on.size() : pos2 + 1;
	size_t part2_len = (pos2 == std::string_view::npos) ? 0 : (on.size() - part2_start);
	// on3 = on.substr(pos, pos2)  (note: original used pos2 as length param)
	std::string_view on3 = (pos2 == std::string_view::npos) ? on.substr(pos) : on.substr(pos, pos2);

	string nID = bDoCombiSamples ? SampleID_Combi[idx] : SampleID[idx];
	nID.reserve(nID.size() + iniSpacer.size() + pos + part2_len + 1 + 12 + on3.size());
	nID.append(iniSpacer);
	// append on.substr(0, pos)
	nID.append(on.data(), pos);
	// append on.substr(pos2 + 1)
	if (part2_len > 0) {
		nID.append(on.data() + part2_start, part2_len);
	}
	nID += ' ';
	nID += "orig_hdPart=";
	nID.append(on3.data(), on3.size());
	d->setNewID(nID);
	d->setBCnumber(idx, BCoffset);
}
void Filters::countBCdetected(int BC, int Pair, bool MidQ) {

	if (!bDoMultiplexing) { return; }
    const int bcIdx = BC - BCoffset;
	if (bcIdx < 0) {
		return;
	}
	if (Pair < 0) { Pair = 0; }
    if (Pair >= (int)collectStatistics.size()) { return; }
	if (!MidQ) {
      if (bcIdx < (int)collectStatistics[Pair]->BarcodeDetected.size()) {
			collectStatistics[Pair]->BarcodeDetected[bcIdx]++;
		}
	}
	else {
       if (Pair < (int)statAddition.size() && statAddition[Pair] && bcIdx < (int)statAddition[Pair]->BarcodeDetected.size()) {
			statAddition[Pair]->BarcodeDetected[bcIdx]++;
		}
	}
}
bool Filters::eval_reversingBC(bool fwd) {
	if (!fwd && !bDoBarcode2) {
		return true;
	}
	/*BCdecide lbcd(BCdFWD);
	if ( !fwd ) {
		lbcd = BCdREV;
	}*/
	if (BCdFWDREV[!fwd].b_BCdirFix) { return true; }
	BCdFWDREV[!fwd].b_BCdirFix = true; lMD->setBCfixed(true, fwd);
	if (BCdFWDREV[!fwd].BCrevhit > BCdFWDREV[!fwd].BChit * 8) {//use reversed_ BC ..
		BCdFWDREV[!fwd].reversedBCs = true;
		if (fwd) {
			reverseTS_all_BC();
		}
		else {
			reverseTS_all_BC2();
		}
		if (BCdFWDREV[!fwd].BChit > 0) {
			restartSet = true;
			return false;
		}
	}
	else if (BCdFWDREV[!fwd].BCrevhit > 0) {
		restartSet = true;
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
           revBarcodes1_[revBarcode[i]] = (int)i;
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
	int scanRegion, string& presentBC, bool fwdStrand, bool revBC, bool endScan) {
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
	vector<int>* localBarcodeLengths;


	// was locked with an omp pragma
	unsigned int maxBCLength;
	unsigned int minBCLength;
	if (fwdStrand) {
		if (revBC) {
			localBarcodes = &revBarcodes1_;
		}
		else {
			localBarcodes = &barcodes1_;
		}
		localBarcodeLengths = &barcodeLengths1_;
		maxBCLength = maxBCLength1_;
		minBCLength = minBCLength1_;
	}
	else {
		if (revBC) {
			localBarcodes = &revBarcodes2_;
		}
		else { localBarcodes = &barcodes2_; }
		localBarcodeLengths = &barcodeLengths2_;
		maxBCLength = maxBCLength2_;
		minBCLength = minBCLength2_;
	}
	if (d->length() < maxBCLength) {
		return;
	}

	int seqLen = d->mem_length();

	//    cout << "minBCLength: " << minBCLength << endl;
	//    cout << "maxBCLength: " << maxBCLength << endl;
	//    cout << "barcodeErrs: " << barcodeErrors_ << endl;
	//

		//bool found = false;


    // Use a thread-local scratch string to avoid repeated allocations for temporary substring
	thread_local std::string scratch;
	for (int start2 = 0; start2 < scanRegion; start2++) {
		start = endScan ? (seqLen - start2 - maxBCLength) : start2;
		if (start < 0) { break; }
		// get a non-owning view of the substring
		std::string_view test = d->getSubSeq(start, (int)maxBCLength);
		// assign into scratch (reuses capacity across iterations/threads) for map lookup
		scratch.assign(test.data(), test.size());
		auto barcodeIterator = localBarcodes->find(scratch);

		if (barcodeIterator != localBarcodes->end()) {// set index if found
			idx = (*barcodeIterator).second;
			stop = start + (int)(*localBarcodeLengths)[idx];
			presentBC.assign(test.data(), test.size()); // copy into output string
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
					stop = start + (int)(*localBarcodeLengths)[idx];
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
				stop = start + (int)(*localBarcodeLengths)[idx];
				presentBC = d->getSubSeq(start, stop);
			}
		}


	}

	if (start == -1) {
		idx = -1;
	}


	return;
}

bool Filters::checkIfRevPrimerHits(shared_ptr<DNA> d, int primerID, int pair, bool twoPrimers) {
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0) { return true; }
	if (d->getFwdPrimDetect()) {
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
	if (pair != 1 || twoPrimers) {
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit);
	}

	if (start != -1) {
		return true;
	}

	return false;

}
bool Filters::passedReads(int n) {
	if (collectStatistics[0]->total - passed_interval_reads > (uint)n) {
		passed_interval_reads = collectStatistics[0]->total;
		return true;
	}
	return false;
}

bool Filters::checkIfPrimerHits(shared_ptr<DNA> d, int primerID, int pair) {
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0) { return true; }
	if (d->getFwdPrimDetect()) {
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
bool Filters::cutPrimer(shared_ptr<DNA> d, int primerID, bool& RC, int pair, bool extensivePrimerCheck) {
	//only adapted to singular BC
	if (PrimerL.size() == 0 || PrimerL[0].length() == 0) { return true; }
	if (d->getFwdPrimDetect()) {
		return true;
	}
	int start(-1), stop(-1);
	int limit(50), startSearch(0), limit2(50);
	if (Bcheck4illuAdapts) {
		limit = 40;
	}
	else	if (!d->getBarcodeCut() && maxBCLength1_ > 0) {
		limit = maxBCLength1_ + 4;
	}
	else { limit = 22; }//in this case nothing is known about 5' end
	if (extensivePrimerCheck) {
		int QS = (int)d->length(); limit = min(QS - 1, max(((int)((float)d->length() * 0.75)), 250));
		limit2 = (d->length() - 100);;
	}

	if (!BcutTag) {
		//Tag was not cut out of sequence_, take this into account
		startSearch = minBCLength1_ - 2;
		limit += (maxBCLength1_ - minBCLength1_) + 4;
	}
	if (!RC) {
		start = d->matchSeq(PrimerL[primerID], PrimerErrs, limit, startSearch);
		stop = start + (int)PrimerL[primerID].length();
	}
	else {
		//if (1 && start == -1){
		int QS = d->length(); int limit = max(int(QS * 0.75), QS - 250); stop = QS;
		start = d->matchSeqRev(PrimerL_RC[primerID], PrimerErrs, limit, startSearch);
		//		start = d->matchSeq(PrimerL_RC[primerID], PrimerErrs, tolerance, startSearch);
		if (0 && start != -1) {
			int y = 0;//debug
		}
	}
	if (start < 0) {//failed to match primer
		//this should only be set if it led to failing to read
		//d->QualCtrl.PrimerFwdFail = true;
		//sPrimerFail(pair_);// max(0, (int)d->getReadMatePos()));
		if (alt_PrimerErrs != 0 && PrimerErrs < alt_PrimerErrs) {
			if (!RC) {
				start = d->matchSeq(PrimerL[primerID], alt_PrimerErrs, limit, startSearch);
				stop = start + (int)PrimerL[primerID].length();
			}
			else {
				int QS = d->length();  stop = QS;
				start = d->matchSeqRev(PrimerL_RC[primerID], alt_PrimerErrs, limit2, startSearch);
			}
		}
		if (start < 0) {
			//statAddition.PrimerFail++;
			return false;
		}
		else if (pair != 1) { //2nd read shouldnt be affected by fwd primer (but still checked in short read mode)
			d->setYellowQual(true);
		}
	}
	//remove everything before/after primer cut
	if (!BcutPrimer) {
		if (!RC) { d->cutSeq(0, start); }
		else { d->cutSeq(stop, -1); }
		return true;
	}
	else {
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
bool Filters::cutPrimerRev(shared_ptr<DNA> d, int primerID, bool& RC, bool extraCheck) {
	//const string& se = d->getSequence();
	if (d->getRevPrimDetect() || PrimerR.size() == 0) {
		return true;
	}

	int start(-1), stop(d->length());
	int QS = d->length();
	int limit = max(QS >> 1, QS - 250);
	int limit2 = 50;
	if (true) {
		limit = min(QS - 1, max(((int)((float)d->length() * 0.75)), 250));
		limit2 = (d->length() - 100);
	}

	//int limit = QS>>1;

	if (RC) {
		start = d->matchSeqRev(PrimerR_RC[primerID], PrimerErrs, limit, extraCheck);
	}
	else {
		start = d->matchSeq(PrimerR[primerID], PrimerErrs, limit2, 0, extraCheck);
		stop = start + (int)PrimerR[primerID].length();
		if (start >= 0) {
			RC = false;
		}
	}

	if (start < 0) {//failed to match primer
		//d->QualCtrl.PrimerRevFail = true;
		return false;
	}

	if (!BcutPrimer) { //found it, but no cut
		if (!RC) { d->cutSeq(0, start); }
		else { d->cutSeq(stop, -1); }
		return true;
	}
	else {
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
bool Filters::readMap() {//core routine to read map info
	if (cmdArgs->find("-map") == cmdArgs->end()) {
		this->noMapMode();
		return true;
	}
	string MapF = (*cmdArgs)["-map"];
	if (cmdArgs->find("-optimalRead2Cluster") != cmdArgs->end()) {
		bDoSeedExtension = true;
	}



	string path = ""; bool pathMode = false;
	if (cmdArgs->find("-i_path") != cmdArgs->end() && (*cmdArgs)["-i_path"].length() >= 1) {
		path = (*cmdArgs)["-i_path"] + string("/");
		pathMode = true;//check later if mapping file contains fasta/fastq
	}

	minBCLength1_ = 100000; minBCLength2_ = 1000000; maxBCLength1_ = 0; maxBCLength2_ = 0; minPrimerLength_ = 100000;
	string line;
	ifstream in(MapF.c_str());
	if (!in) {
		cerr << "Could not find " << MapF << " mapping file. Exiting.\n"; exit(2);
	}
	int ini_ColPerRow(0), cnt(0), skips(0);

	//check MAP format
	//while(getline(in,line,'\n')) {
	while (!safeGetline(in, line).eof()) {
		if (line.substr(0, 1) == "#" || line.length() == 0) { skips++; continue; }
		string segments;
		int ColsPerRow = 0; // Initialize counter.
		stringstream ss;
		ss << line;
		while (getline(ss, segments, '\t')) {
			ColsPerRow++;
		}
		if (segments == "") { ColsPerRow++; }

		if (cnt == 0) {
			ini_ColPerRow = ColsPerRow;
		}
		else {
			if (ColsPerRow != ini_ColPerRow) {
				cerr << "Number of columns on line " << cnt + skips << " is " << ColsPerRow << ". Expected " << ini_ColPerRow << " columns.\n";
				return false;
			}
		}
		cnt++;
	}
	if (ini_ColPerRow == 0) {
		cerr << "Mapping File exists, but appears to be badly formated (0 columns detected). Exiting\n"; exit(2);
	}
	if (cnt == 0) {
		cerr << "Mapping File exists, but appears to be badly formated (0 lines detected). Exiting\n"; exit(2);
	}


	PrimerIdx.resize(cnt, -1); PrimerIdxRev.resize(cnt, -1);
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
	cnt = 0;
	bool hasQualityColumn = false;
	vector<string> terms(16); terms[0] = "SampleID";
	terms[1] = "BarcodeSequence"; terms[2] = "LinkerPrimerSequence";
	terms[3] = "ReversePrimer"; terms[4] = "fastqFile";
	terms[5] = "fnaFile"; terms[6] = "qualFile";
	terms[7] = "SampleIDinHead"; terms[8] = "MIDfqFile";
	terms[9] = "CombineSamples"; terms[10] = "ForwardPrimer";
	terms[11] = "Barcode2ndPair";
	terms[12] = "HetSpacerFwd";	terms[13] = "HetSpacerRev";
	terms[14] = "derepMin"; terms[15] = "SequencingRun";
	bool hetOneSide = false;

	vector<int> termIdx(terms.size(), -1);

	while (!safeGetline(in, line).eof()) {
		//	while(getline(in,line,'\n')) {
		if (cnt != 0 && line.substr(0, 1) == "#") { continue; }
		if (line.length() < 10) { continue; }
		if (cnt == 0) {
			line = line.substr(1);
		}

		string segments;
		stringstream ss;
		ss << line;
		int tbcnt = 0;

		//(*cmdArgs)["-i_MID_fastq"]
		while (getline(ss, segments, '\t')) {
			trim(segments);
			if (cnt == 0) { //search for header
				//Primer, BarcodeSequence, LinkerPrimerSequence
				//PrLCol(-1), PrRCol(-1), BCCol(-1), SIDCol(-1);
				for (unsigned int i = 0; i < terms.size(); i++) {
					if (segments == terms[i]) {
						if (i == 6) { hasQualityColumn = true; }
						if (i == 12) { if (hetOneSide) { doHetPrimerExplicit = true; } else { hetOneSide = true; } }
						if (i == 13) { if (hetOneSide) { doHetPrimerExplicit = true; } else { hetOneSide = true; } }
						if (termIdx[i] != -1) {
							cerr << "Header contains ambiguous entries: " << segments << "\n";
							exit(9);
						}
						termIdx[i] = tbcnt;
					}
				}
			}
			else { //read data into entries
				for (uint k = 0; k < terms.size(); k++) {
					if (termIdx[k] == tbcnt) {
						extractMap(k, cnt - 1, tbcnt, segments, hasQualityColumn);
					}
				}
			}
			tbcnt++;
		}
		cnt++;
	}
	//check some prerequisites 
	if (HeadSmplID.size() == 0 && Barcode.size() == 0) {
		cerr << "Could not find in mapping file either (1) valid Barcodes or (2) valid SampleID's in Sequence header. Exiting\n";
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
	if (termIdx[4] != -1 && termIdx[5] != -1 && termIdx[6] != -1) {
		checkDoubleBarcode();
	}
	if (termIdx[9] != -1) { bDoCombiSamples = true; }
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
	}
	else {
		this->removeSRs();
	}

	this->BarcodePreStats();

	return true;
}
void Filters::decideHeadBC() {
	bDoMultiplexing = true;
	if (HeadSmplID[0].length() > 0 && Barcode[0].length() == 0) {
		bDoBarcode = true; bDoHeadSmplID = true; minPrimerLength_ = 0; return;
	}
	else if (HeadSmplID[0].length() == 0 && Barcode[0].length() > 0 && Barcode2[0].length() > 0) {
		bDoBarcode = true; bDoHeadSmplID = false; bDoBarcode2 = true;  return;
	}
	else if (HeadSmplID[0].length() == 0 && Barcode[0].length() > 0) {
		bDoBarcode = true; bDoHeadSmplID = false; bDoBarcode2 = false;  return;
	}
	else if (HeadSmplID[0].length() == 0 && Barcode[0].length() == 0) {
		//simply check if each filename is different
		if (FastqF.size() > 0) {
			for (uint i = 0; i < FastqF.size(); i++) {
				for (uint j = i + 1; j < FastqF.size(); j++) {
					if (FastqF[i] == FastqF[j]) {
						cerr << "File names " << i << " and " << j << " are equal - no identification by filename possible.\n   Aborting..\n"; exit(55);
					}
				}
			}
			bOneFileSample = true; bDoBarcode = true; bDoHeadSmplID = false;
			return;
		}
		if (FastaF.size() > 0) {
			for (uint i = 0; i < FastaF.size(); i++) {
				for (uint j = i + 1; j < FastaF.size(); j++) {
					if (FastaF[i] == FastaF[j]) {
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


void Filters::checkDoubleBarcode() {
	if (!bDoBarcode || bDoHeadSmplID || Barcode.size() == 0) { return; }
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	if (bDoBarcode2) {
		if (Barcode.size() != Barcode2.size()) {
			cerr << "Unequal Barcode vector sizes in dual barcoding controls. Exiting.."; exit(45);
		}
		for (unsigned int i = 0; i < Barcode.size(); i++) {
			for (unsigned int j = i + 1; j < Barcode.size(); j++) {
				if (strcmp(Barcode[i].c_str(), Barcode[j].c_str()) == 0 && strcmp(Barcode2[i].c_str(), Barcode2[j].c_str()) == 0) {
					empty[0] = i; empty[1] = j; doubles.push_back(empty);
				}
			}
		}
		if (doubles.size() > 0) {
			for (uint x = 0; x < doubles.size(); x++) {
				int i = doubles[x][0]; int j = doubles[x][1];
				cerr << "Duplicate dual Barcode detected: Barcode1 " << i + 1 << " (" << Barcode[i] << ") and " << j + 1 << " (" << Barcode[j] << ")  as well as Barcode1 " << i + 1 << " (" << Barcode2[i] << ") and " << j + 1 << " (" << Barcode2[j] << ") are equal.\n";
			}
			exit(8);
		}
	}
	else {
		for (unsigned int i = 0; i < Barcode.size(); i++) {
			for (unsigned int j = i + 1; j < Barcode.size(); j++) {
				if (strcmp(Barcode[i].c_str(), Barcode[j].c_str()) == 0) {
					empty[0] = i; empty[1] = j; doubles.push_back(empty);
				}
			}
		}
		if (doubles.size() > 0) {
			for (uint x = 0; x < doubles.size(); x++) {
				int i = doubles[x][0]; int j = doubles[x][1];
				cerr << "Duplicate Barcode detected: Barcode " << i + 1 << " (" << Barcode[i] << ") and " << j + 1 << " (" << Barcode[j] << ") are equal.\n";
			}
			exit(8);
		}

	}
}
void Filters::checkDoubleSampleIDHead() {
	if (!bDoHeadSmplID) { return; }
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	for (unsigned int i = 0; i < HeadSmplID.size(); i++) {
		for (unsigned int j = i + 1; j < HeadSmplID.size(); j++) {
			if (strcmp(HeadSmplID[i].c_str(), HeadSmplID[j].c_str()) == 0) {
				empty[0] = i; empty[1] = j; doubles.push_back(empty);
			}
		}
	}
	if (doubles.size() > 0) {
		for (uint x = 0; x < doubles.size(); x++) {
			int i = doubles[x][0]; int j = doubles[x][1];
			cerr << "Duplicate Header2split detected: pattern " << i + 1 << " and " << j + 1 << " are equal.\n";
		}
		exit(8);
	}

}

void Filters::checDoubleSampleID() {
	vector<vector<int>> doubles(0);	vector<int> empty(2, 0);
	for (unsigned int i = 0; i < SampleID.size(); i++) {
		for (unsigned int j = i + 1; j < SampleID.size(); j++) {
			if (strcmp(SampleID[i].c_str(), SampleID[j].c_str()) == 0) {
				empty[0] = i; empty[1] = j; doubles.push_back(empty);
			}
		}
	}
	if (doubles.size() > 0) {
		for (uint x = 0; x < doubles.size(); x++) {
			int i = doubles[x][0]; int j = doubles[x][1];
			cerr << "Duplicate SampleID detected: SampleID " << i + 1 << " and " << j + 1 << " are equal.\n";
		}
		exit(8);
	}
	if (!bDoCombiSamples) {
		return;
	}
	bDoCombiSamples = false;
	if (SampleID_Combi.size() <= 1) {
		return;
	}

	string prevCSID = SampleID_Combi[0];
	for (unsigned int i = 1; i < SampleID_Combi.size(); i++) {
		if (SampleID_Combi[i] != prevCSID) {
			bDoCombiSamples = true;
		}
		if (SampleID_Combi[i] == "") {
			SampleID_Combi[i] = SampleID[i];
		}
	}
}

void Filters::removeSRs(void) {
	SequencingRun2id.clear();
	SequencingRun2id[""] = vector(0, 0);
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
			}
			else {
				FQ = FastqF[J];
			}
			auto XX = files.find(FQ);
			if (XX == files.end()) {
				files[FQ] = curBlk;
			}
			else if (XX->second != curBlk) {
				cerr << "Sample " << SampleID[J] << "(" << FQ << ") has SequencingRun " << SequencingRun[J] << " already found for other files. This is not allowed.\n";
				exit(323);
			}
		}
	}
}

void Filters::BarcodePreStats() {
	minBCLength1_ = 100000; maxBCLength1_ = 0;
	for (unsigned int i = 0; i < Barcode.size(); i++) {
		if (Barcode[i].length() < minBCLength1_) minBCLength1_ = (unsigned int)Barcode[i].length();
		if (Barcode[i].length() > maxBCLength1_) maxBCLength1_ = (unsigned int)Barcode[i].length();
		//initialize Barcodes
		barcodes1_[Barcode[i]] = i;
		barcodeLengths1_[i] = (int)Barcode[i].length();
	}
	if (minBCLength1_ == maxBCLength1_) {
		bBarcodeSameSize = true;
	}
	minBCLength2_ = 100000; maxBCLength2_ = 0;
	for (unsigned int i = 0; i < Barcode2.size(); i++) {
		//create index
		barcodes2_[Barcode2[i]] = i;
		if (Barcode2[i].length() < minBCLength2_) minBCLength2_ = (unsigned int)Barcode2[i].length();
		if (Barcode2[i].length() > maxBCLength2_) maxBCLength2_ = (unsigned int)Barcode2[i].length();
		barcodeLengths2_[i] = (int)Barcode2[i].length();
	}
	//fix empty last column specifically for derepMinNum
	if (derepMinNum.size() > 0) {
		derepMinNum.resize(Barcode.size(), -1);
	}

}
void Filters::resetStats() {
	for (size_t i = 0; i < 2; i++) {
		collectStatistics[i]->reset();
        if (statAddition[i]) {
			statAddition[i]->reset();
		}
	}

}

void Filters::failedStats2(const shared_ptr<DNA>& d, int pair) {
	int pa = max(pair, 0);
  if (pa >= (int)collectStatistics.size()) {
		return;
	}
	if (bDoMultiplexing) {
		int idx = d->getBCnumber() - BCoffset;
		if (bOneFileSample) {
            if (!collectStatistics[pa]->BarcodeDetectedFail.empty()) {
				collectStatistics[pa]->BarcodeDetectedFail[0]++;
			}
		}
		else if (idx >= 0) {


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

void Filters::debugVerifyStats(const char* context) const {
#ifdef _DEBUG
	cerr << "[Filters::debugVerifyStats] " << context << ": Verifying statistics consistency...\n";
	for (size_t i = 0; i < 2; i++) {
		if (!collectStatistics[i]) {
			cerr << "[Filters::debugVerifyStats] " << context << " pair " << i << ": collectStatistics is null\n";
			continue;
		}
     const collectstats& cs = *collectStatistics[i];
		const std::vector<std::string> warnings = get_stats_invariant_warnings(cs, "pair " + itos((int)i));
		for (const auto& w : warnings) {
			cerr << "[Filters::debugVerifyStats] " << context << " " << w << "\n";
		}

		if (bAdditionalOutput) {
			if (i >= statAddition.size() || !statAddition[i]) {
				cerr << "[Filters::debugVerifyStats] " << context << " pair " << i << ": statAddition is null\n";
				continue;
			}
          const collectstats& sa = *statAddition[i];
			const std::vector<std::string> addWarnings = get_stats_invariant_warnings(sa, "add pair " + itos((int)i));
			for (const auto& w : addWarnings) {
				cerr << "[Filters::debugVerifyStats] " << context << " " << w << "\n";
			}
		}
	}
#else
	(void)context;
#endif
}

void Filters::addMergeStats(OutputStreamer* out) {
	mergeStats->BPwritten = (uint)out->getBPwrittenInSR();
	mergeStats->BPmergeWritte = (uint)out->getBPwrittenInSRmerg();
	mergeStats->total_read_preMerge_ = (int)out->total_read_preMerge_;
	mergeStats->merged_counter_ = (int)out->merged_counter_;



}
void Filters::prepStats() {
    debugVerifyStats("prepStats:before");
	for (size_t i = 0; i < 2; i++) {
     float remSeqs = float(collectStatistics[i]->total - collectStatistics[i]->totalRejected);
		collectStatistics[i]->PostFilt.calcSummaryStats(remSeqs, min_l, min_q);
		if (bAdditionalOutput) {
           remSeqs = float(statAddition[i]->total - statAddition[i]->totalRejected);
			statAddition[i]->PostFilt.calcSummaryStats(remSeqs, min_l, min_q);//yellow

		}
		collectStatistics[i]->PreFilt.calcSummaryStats(1, min_l, min_q);
	}
	debugVerifyStats("prepStats:after");

}


void Filters::addPrimerL(string segments, int cnt) {
	int used = -1;
	for (unsigned int i = 0; i < PrimerL.size(); i++) {
		if (segments == PrimerL[i]) { used = i; }
	}
	if (used == -1) {
		PrimerL.push_back(segments);
		PrimerL_RC.push_back(reverseTS2(segments));
		PrimerIdx[cnt] = (int)PrimerL.size() - 1;
		if (segments.length() < minPrimerLength_) minPrimerLength_ = (unsigned int)segments.length();
	}
	else {
		PrimerIdx[cnt] = used;
	}
}
void Filters::addPrimerR(string segments, int cnt) {
	bPrimerR = true;
	int used = -1;
	for (unsigned int i = 0; i < PrimerR.size(); i++) {
		if (segments == PrimerR[i]) { used = i; }
	}
	if (used == -1) {
		PrimerR.push_back(segments);
		PrimerR_RC.push_back(reverseTS2(segments));
		PrimerIdxRev[cnt] = (int)PrimerR.size() - 1;
	}
	else {
		PrimerIdxRev[cnt] = used;
	}
}


void Filters::extractMap(int k, int cnt, int tbcnt, string& segments,
	bool hasQualityColumn) {
	//terms[4] = "fastqFile";terms[5] = "fnaFile"; terms[6] = "qualFile";
	switch (k)
	{
	case 4: // fastq file
		trim(segments);
		FastqF.push_back(segments);
		break;
	case 5: // fna file
		trim(segments);
		FastaF.push_back(segments);
		if (!hasQualityColumn) {//code to create artificial quality file name
			string newQ = segments;
			size_t pos = newQ.find_last_of(".");
			newQ = newQ.substr(0, pos);
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
		transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
		this->addPrimerL(segments, cnt);
		break;
	case 3: // right primer
		trim(segments);
		transform(segments.begin(), segments.end(), segments.begin(), ::toupper);
		this->addPrimerR(segments, cnt);
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
		}
		else if (!is_digits(segments)) {
			cerr << "Wrong map entry \"" << segments << "\". For header derepMin only number can be used.\n"; exit(313);
		}
		else {
			int nint = atoi(segments.c_str());
			derepMinNum.push_back(nint);
		}
		break;
	case 12://het primer fw
		if (!doHetPrimerExplicit) { break; }
		hetPrimer[0][cnt] = segments;
		break;
	case 13://het primer rv
		if (!doHetPrimerExplicit) { break; }
		hetPrimer[1][cnt] = segments;
		break;
	case 15: //SequencingRun
		trim(segments);
		SequencingRun[cnt] = segments;
		auto XX = SequencingRun2id.find(segments);
		if (XX != SequencingRun2id.end()) {
			XX->second.push_back(cnt);
		}
		else {
			SequencingRun2id[segments] = vector<int>(1, cnt);
		}
		break;
	}

	if (k == 6 || k == 5) {//a qual_ pushback was "" (was empty); takeOver
		if (QualF.size() == FastaF.size() && QualF.back() == "") {
			string newQ = FastaF.back();
			size_t pos = newQ.find_last_of(".");
			newQ = newQ.substr(0, pos);
			newQ += string(".qual_");
			QualF.back() = newQ;
		}
	}

}

void Filters::printLenVsQual(ostream& give) {
	collectStatistics[0]->PreFilt.printLvsQ(give);
}

void Filters::printHisto(ostream& give, int which, int set) {
	bool p2stat = pairedSeq > 1;

	if (set == 0) {
		vector<uint> colStats(collectStatistics[0]->PostFilt.get_rstat_Vmed(which));
		vector<size_t> ra(collectStatistics[0]->PostFilt.getVrange(which));

		if (which == 1) {
			give << "qual_\tFilterObs" << endl;
		}
		else {
			give << "Length\tFilterObs" << endl;
		}
		for (size_t i = ra[0]; i < ra[1]; i++) {
			give << i << "\t" << colStats[i] << endl;
		}
	}
	else if (set == 1) {
		vector<size_t> ra(2, 0), tra;
		vector<uint> stat;
		vector<bool> skips(6, false);
		vector<vector<uint>> matHist;
		if (which == 1) {
			give << "#qual_\t";
		}
		else { give << "#Length\t"; }
		if (p2stat && b_doFilter) { give << "FilteredP1\tFilteredP2\t"; }
		else if (b_doFilter) { give << "Filtered\t"; skips[1] = true; }
		else { skips[1] = true; skips[0] = true; }

		if (bAdditionalOutput && b_doFilter) {
			if (p2stat) { give << "AddFilterP1\tAddFilterP2\t"; }
			else { give << "AddFilter\t"; skips[3] = true; }
		}
		else {
			skips[2] = true; skips[3] = true;
		}
		if (p2stat) { give << "RawReadsP1\tRawReadsP2"; }
		else { give << "RawReads"; skips[5] = true; }
		give << endl;
		ra = collectStatistics[0]->PostFilt.getVrange(which);
		if (p2stat) { tra = collectStatistics[1]->PostFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
		if (!skips[2]) {
			tra = statAddition[0]->PostFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]);
           if (p2stat) { tra = statAddition[1]->PostFilt.getVrange(which); ra[0] = min(tra[0], ra[0]); ra[1] = max(tra[1], ra[1]); }
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
				if (i >= ra[1]) { break; }
				matHist[kk][i] = stat[i];
			}
		}
		for (size_t i = ra[0]; i < ra[1]; i++) {
			give << i;
			for (size_t kk = 0; kk < matHist.size(); kk++) {
				if (skips[kk]) { continue; }
				give << "\t" << matHist[kk][i];
			}
			give << endl;
		}

	}
}
void Filters::FileEssentials(filesStr& files, OptContainer* cmdArgs) {//unordered_map<string, int>& UFF) {

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
		pair<string, string> tmp(SR, uFX.first);
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

vector<int> Filters::combiSmplConvergeVec(const vector<string>& inNames) {
	vector<int> retV(inNames.size(), -1);
	unordered_map<string, int> smpl2combi;
	unordered_map<string, int>::iterator s2cIT;
	int cntGrps(-1);
	for (size_t i = 0; i < SampleID_Combi.size(); i++) {
		s2cIT = combiMapCollectGrp.find(SampleID_Combi[i]);
		if (s2cIT == combiMapCollectGrp.end()) {
			cntGrps++;
			combiMapCollectGrp[SampleID_Combi[i]] = cntGrps;
		}
		smpl2combi[SampleID[i]] = combiMapCollectGrp[SampleID_Combi[i]];
	}
	for (size_t i = 0; i < inNames.size(); i++) {
		s2cIT = smpl2combi.find(inNames[i]);
		if (s2cIT == smpl2combi.end()) {
			cerr << "Can't find SampleID " << inNames[i] << " in reference sample_id_ Names\n"; exit(113);
		}
		retV[i] = s2cIT->second;
	}
	return retV;
}


string Filters::shortStats(const string& file) {
	shared_ptr<collectstats> cst = collectStatistics[0];
	string ret("");
	if (file != "") {
		ret += file + "\n";
	}
	if (pairedSeq > 1) {
		ret += "Pair 1: ";
	}
	//char buffer[50];
   const float p1Total = static_cast<float>(cst->total);
	float tmp = p1Total > 0.f ? (100.f * float(cst->total - cst->totalRejected) / p1Total) : 0.f;
	float p1Trimmed = p1Total > 0.f ? (100.f * float(cst->Trimmed) / p1Total) : 0.f;
	ostringstream os;
  os << tmp << "% of " << cst->total << " reads accepted (" << p1Trimmed << "% end-trimmed)\n";
	//	sprintf(buffer, "%.3f%% of %d", tmp, cst->total); ret += buffer;
	//	sprintf(buffer," reads accepted (%.3f%% end-trimmed)\n", (100.f* float(cst->total - cst->Trimmed) / (float)cst->total)); ret += buffer;

	if (pairedSeq > 1) {
		shared_ptr<collectstats> cst = collectStatistics[1];
      const float p2Total = static_cast<float>(cst->total);
		float p2Accepted = p2Total > 0.f ? (100.f * float(cst->total - cst->totalRejected) / p2Total) : 0.f;
		float p2Trimmed = p2Total > 0.f ? (100.f * float(cst->Trimmed) / p2Total) : 0.f;
		os << "Pair 2: " << p2Accepted << "% of " << cst->total;
		os << " reads accepted (" << p2Trimmed << "% end-trimmed)\n";
		//		sprintf(buffer,"Pair 2: %.3f%% of %d", (100.f*float(cst->total - cst->totalRejected) / (float)cst->total), cst->total); ret += buffer;
		//		sprintf(buffer," reads accepted (%.3f%% end - trimmed)\n", (100.f* float(cst->total - cst->Trimmed) / (float)cst->total)); ret += buffer;
	}
	ret = os.str();
	return ret;
}
void Filters::printGC(ostream& os, int Npair) {
	os << "Subset\t\tOccurence\t\tAvg.Quality\n";
	os << "\tA\tT\tG\tC\tA\tT\tG\tC\n";
	os << "R1 pre-filter";
	collectStatistics[0]->PreFilt.printGCstats(os);
	if (Npair > 1) {
		os << "R2 pre-filter";
		collectStatistics[1]->PreFilt.printGCstats(os);
	}
	if (!b_doFilter) { return; }
	os << "R1 filtered";
	collectStatistics[0]->PostFilt.printGCstats(os);
	if (Npair > 1) {
		os << "R2 filtered";
		collectStatistics[1]->PostFilt.printGCstats(os);
	}
}

void Filters::printStats(ostream& give, string file, string outf, bool greenQualStats) {
    debugVerifyStats("printStats");
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
	if (file.length() > 0) {
		give << "Input File:  " << file << endl;
	}
	if (outf == "-") {
		give << "Output File: stdout\n";
	}
	else if (outf.length() > 0) {
		give << "Output File: " << outf << endl;
	}
	if (!b_doFilter) {
		give << "No valid Filter file provided; no filtering done on files\n";
		return;
	}
	/*if (!greenQualStats) {
		give << "Statistics of reads that passed the mid qual filter\n";
	} else {
		give << "Statistics of high quality reads\n";
	}*/
	float remSeqs = float(cst->total - cst->totalRejected);

	give << endl;

	string ReadTag = "Reads";
	if (isGoldAxe()) { ReadTag = "(Sub)Reads"; }

	if (!greenQualStats) {
		give << ReadTag << " not High qual_: " << intwithcommas((int)cst->totalRejected);
	}
	else {
		give << "Reads processed: " << intwithcommas((int)cst->total);
		if (p2stat) {
			give << "; " << intwithcommas((int)cst2->total) << " (pair 1;pair 2)";
		}
	}
	give << endl;
	if (cst->reversedRds > 0 || cst->swappedRds > 0) {
		if (cst->reversedRds > 0) {
			give << cst->reversedRds;
			if (p2stat) { give << ";" << cst2->reversedRds; }
			give << " reads reverse-translated";
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
	if (!greenQualStats) {
		give << "Rejected:" << intwithcommas((int)(cst->totalRejected)) << endl;
		give << "Accepted (High qual): " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst->Trimmed) << " end-trimmed)\n";
		give << "Accepted (Mid qual): " << intwithcommas(numMid) << endl;
	}
	else {
		if (p2stat) {
			int numMid2 = (int)cst2->totalMid;
			int numAccept2 = int(cst2->totalSuccess - cst2->totalMid);
			give << "Rejected: " << intwithcommas((int)cst->totalRejected) << "; " << intwithcommas((int)cst2->totalRejected) << endl;
			give << "Accepted (High qual): " << intwithcommas((int)numAccept) << "; " << intwithcommas((int)numAccept2) << " (" << intwithcommas((int)cst->Trimmed) << "; " << intwithcommas((int)cst2->Trimmed) << " end-trimmed)\n";
			give << "Accepted (Mid qual): " << intwithcommas(numMid) << "; " << intwithcommas(numMid2) << endl;
		}
		else {
			give << "Rejected: " << intwithcommas((int)cst->totalRejected) << endl;
			give << "Accepted (High qual): " << intwithcommas((int)numAccept) << " (" << intwithcommas((int)cst->Trimmed) << " end-trimmed)\n";
			give << "Accepted (Mid qual): " << intwithcommas(numMid) << endl;
		}
	}


	if (false && bPrimerR) { //confusing collectStatistics
		give << ", with rev. primer: " << intwithcommas((int)cst->RevPrimFound); if (p2stat) { give << "; " << intwithcommas((int)cst2->RevPrimFound); }
	}

	if (pairedSeq > 1) {
		give << "Singletons among these: " << intwithcommas((int)cst->singleton) << "; " << intwithcommas((int)cst2->singleton) << endl;
	}
	give << "Bad Reads recovered with dereplication: " << intwithcommas((int)cst->DerepAddBadSeq) << endl;

	if (bShortAmplicons) {
		give << "Short amplicon mode.\n";
	}
	if (checkSwitchedRdPairs()) {
		//give << "Looked for switched read pairs (" << intwithcommas(revConstellationN) << " detected)" << endl;
	}
	if (greenQualStats) {
		cst->PostFilt.printStats2(give, remSeqs, 0);
		collectStatistics[1]->PostFilt.printStats2(give, remSeqs, 1);
	}
	else {
		statAddition[0]->PostFilt.printStats2(give, remSeqs, 0);
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
	if (maxAccumQP > 0.0) {
		give << "  > (" << dval << ") acc. errors, trimmed seqs : " << spaceX(8 - digitsFlt(dval)) << intwithcommas((int)cst->AccErrTrimmed);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->AccErrTrimmed); } give << endl;
	}

	give << "  > (" << trimHomonucleotide << ") homo-nt trimmed  : " << spaceX(17 - digitsInt(maxHomonucleotide)) << intwithcommas((int)cst->HomoNTtrimmed);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->HomoNTtrimmed); } give << endl;

	give << "Rejected due to:\n";
	float val = (float)min_l;
	if (val == -1.f) { val = min_l_p; }
	if (!greenQualStats) { val = (float)alt_min_l; }

	give << "  < min Sequence length (" << val << ")  : " << spaceX(18 - digitsFlt(val)) << intwithcommas((int)cst->minL);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->minL); } give << endl;
	if (cst->minLqualTrim > 0) {//this is failed because seq was too short after trimming
		give << "       -after Quality trimming : " << spaceX(10) << intwithcommas((int)cst->minLqualTrim);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->minLqualTrim); } give << endl;
	}
	float valf = min_q;	if (!greenQualStats) { valf = alt_min_q; }
	give << "  < avg Quality (" << valf << ")  : " << spaceX(21 - digitsInt((int)min_q)) << intwithcommas((int)cst->AvgQual);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->AvgQual); } give << endl;
	give << "  < window (" << FQWwidth << " nt) avg. Quality (" << FQWthr << ")  : " << spaceX(5 - digitsInt(FQWwidth)) << intwithcommas((int)cst->QualWin);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->QualWin); } give << endl;
	give << "  > max Sequence length (" << max_l << ")  : " << spaceX(18 - digitsInt(max_l)) << intwithcommas((int)cst->maxL);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->maxL); } give << endl;
	give << "  > (" << maxHomonucleotide << ") homo-nt run  : " << spaceX(21 - digitsInt(maxHomonucleotide)) << intwithcommas((int)cst->HomoNT);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->HomoNT); } give << endl;
	int val2 = MaxAmb;	if (!greenQualStats) { val2 = alt_MaxAmb; }
	give << "  > (" << val2 << ") amb. Bases  : " << spaceX(22 - digitsInt(val2)) << intwithcommas((int)cst->MaxAmb);
	if (p2stat) { give << "; " << intwithcommas((int)cst2->MaxAmb); } give << endl;
	if (BinFilP >= 0.f) {
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
	if ((bDoMultiplexing || cst->PrimerFail > 0) || ((!greenQualStats && alt_bRequireFwdPrim) || bRequireFwdPrim)) {
		give << "  -With fwd Primer remaining (<= " << PrimerErrs << " mismatches";
		if ((!greenQualStats && alt_bRequireFwdPrim) || bRequireFwdPrim) {
			give << ", required) : ";
			give << spaceX(1 - digitsInt(PrimerErrs));
		}
		else {
			give << ") : " << spaceX(11 - digitsInt(PrimerErrs));
		}
		give << intwithcommas((int)cst->PrimerFail);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->PrimerFail) << endl; }
		else { give << endl; }
	}
	if (bPrimerR) {
		give << "  -With rev Primer remaining (<= " << PrimerErrs << " mismatches";
		if ((!greenQualStats && alt_bRequireRevPrim) || bRequireRevPrim) {
			give << ", required) : ";
			give << spaceX(1 - digitsInt(PrimerErrs));
		}
		else {
			give << ") : " << spaceX(11 - digitsInt(PrimerErrs));
		}
		give << intwithcommas((int)cst->PrimerRevFail);
		if (p2stat) { give << "; " << intwithcommas((int)cst2->PrimerRevFail) << endl; }
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

		}
		else if (bDoHeadSmplID) {
			give << endl << "SampleID\t";
			if (bDoCombiSamples) { give << "\tSampleGroup"; }
			give << "\tSampleID\tInstances\n";
			for (unsigned int i = 0; i < Barcode.size(); i++) {
				give << SampleID[i] << "\t";
				if (bDoCombiSamples) { give << SampleID_Combi[i] << "\t"; }
				give << HeadSmplID[i] << "\t" << intwithcommas((int)cst->BarcodeDetected[i]) << endl;
			}
		}
	}

}

//statistics for each single sample 
void Filters::SmplSpecStats(ostream& give) {
	shared_ptr<collectstats> cst = collectStatistics[0];
	shared_ptr<collectstats> cst2 = collectStatistics[1];
	bool p2stat = pairedSeq > 1;
	give << std::setprecision(3);
	give << "SampleID";
	if (bDoCombiSamples) {
		give << "\tSampleGroup";
	}
	if (bDoHeadSmplID) {
		give << "\tSampleID";

	}
	else {
		if (bDoBarcode2) { give << "\tBarcode\tBarcode2"; }
		else { give << "\tBarcode"; }
	}
	if (p2stat) {
		give << "\tRead1Accepted\tRead1Filtered\tRead1PassedFrac\tRead2Accepted\tRead2Filtered\tRead2PassedFrac\n";
	}
	else {
		give << "\tReadsAccepted\tReadsFailed\tPassed%\n";
	}
	for (unsigned int i = 0; i < Barcode.size(); i++) {
		give << SampleID[i] << "\t";
		if (bDoCombiSamples) { give << SampleID_Combi[i] << "\t"; }
		if (p2stat) {
			give << Barcode[i] << "\t";
			if (bDoBarcode2) { give << Barcode2[i] << "\t"; }
			give << cst->BarcodeDetected[i] << "\t" << cst->BarcodeDetectedFail[i] << "\t";
			float totSum = (float(cst->BarcodeDetected[i]) + float(cst->BarcodeDetectedFail[i]));
			if (totSum > 0) {
				give << float(cst->BarcodeDetected[i]) / totSum << "\t";
			}
			else { give << "NA\t"; }

			give << cst2->BarcodeDetected[i] << "\t" << cst2->BarcodeDetectedFail[i] << "\t";
			totSum = (float(cst2->BarcodeDetected[i]) + float(cst2->BarcodeDetectedFail[i]));
			if (totSum > 0) {
				give << float(cst2->BarcodeDetected[i]) / totSum << "";
			}
			else { give << "NA"; }
			give << endl;
		}
		else {
			give << Barcode[i] << "\t";
			if (bDoBarcode2) { give << Barcode2[i] << "\t"; }
			give << (int)cst->BarcodeDetected[i]
				<< "\t" << cst->BarcodeDetectedFail[i] << "\t";

			float totSum = (float(cst->BarcodeDetected[i]) + float(cst->BarcodeDetectedFail[i]));
			if (totSum > 0) {
				give << float(cst->BarcodeDetected[i]) / totSum << "";
			}
			else { give << "NA"; }
			give << endl;
		}
	}



}



void Filters::addStats(Filters* fil, vector<int>& idx) {
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
