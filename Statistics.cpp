//
// Created by fritsche on 10/11/2020.
//

#include "Statistics.h"



void collectstats::addStats(shared_ptr<collectstats> cs, vector<int>& idx) {
    if (BarcodeDetected.size() > (uint)10000) {
        cerr << "Unrealistic number of barcodes (>10000) in addStats\n"; exit(79);
    }
    int BCS = (int)BarcodeDetected.size();
    for (unsigned int i = 0; i < idx.size(); i++) {
        if (idx[i] >= BCS) { return; }
        //assert(idx[i] < BCS);
        BarcodeDetected[idx[i]] += cs->BarcodeDetected[i];
        BarcodeDetectedFail[idx[i]] += cs->BarcodeDetectedFail[i];
    }
    maxL += cs->maxL;	PrimerFail += cs->PrimerFail;
    AvgQual += cs->AvgQual; HomoNT += cs->HomoNT; HomoNTtrimmed += cs->HomoNTtrimmed; 
    totalMid += cs->totalMid;
    PrimerRevFail += cs->PrimerRevFail;
    minL += cs->minL; minLqualTrim += cs->minLqualTrim; TagFail += cs->TagFail;
    MaxAmb += cs->MaxAmb; QualWin += cs->QualWin;
    Trimmed += cs->Trimmed;
    AccErrTrimmed += cs->AccErrTrimmed;
    total.fetch_add(cs->total.load(), std::memory_order_relaxed);
    QWinTrimmed += cs->QWinTrimmed;
    totalRejected += cs->totalRejected;
    fail_correct_BC += cs->fail_correct_BC; suc_correct_BC += cs->suc_correct_BC;
    failedDNAread += cs->failedDNAread; adapterRem += cs->adapterRem;
    RevPrimFound += cs->RevPrimFound;
    singleton += cs->singleton;
    BinomialErr += cs->BinomialErr;
    dblTagFail += cs->dblTagFail;
    DerepAddBadSeq += cs->DerepAddBadSeq;
    total2 += cs->total2; totalSuccess += cs->totalSuccess;
    swappedRds += cs->swappedRds; reversedRds += cs->reversedRds;
    PostFilt.addRepStats(cs->PostFilt);
    PreFilt.addRepStats(cs->PreFilt);

}

void collectstats::reset() {
    singleton = 0;
    size_t BCsiz = BarcodeDetected.size();
    for (unsigned int i = 0; i < BCsiz; i++) {
        BarcodeDetected[i] = 0;
        BarcodeDetectedFail[i] = 0;
    }
    maxL = 0; PrimerFail = 0; AvgQual = 0; HomoNT = 0; HomoNTtrimmed = 0;
    PrimerRevFail = 0;
    minL = 0; TagFail = 0; MaxAmb = 0; QualWin = 0;
    Trimmed = 0; total = 0; totalRejected = 0;
    fail_correct_BC = 0; suc_correct_BC = 0; failedDNAread = 0;
    adapterRem = 0; RevPrimFound = 0; DerepAddBadSeq = 0;
    total2 = 0; totalSuccess = 0;

    //reportStats
    PostFilt.reset(); PreFilt.reset();
}



void ReportStats::addDNAStats(shared_ptr<DNA> d) {
    float avq = d->getAvgQual();
    float ace = (float)d->getAccumError();
    uint len = d->length();
    int median(0);
    if (bLvsQlogs) {
        median = d->getMedianQual();
    }
    stats_mutex.lock();
    //pretty fast
    addMeanStats(len, (int)avq, ace);
    //NT specific quality scores

//    d->ntSpecQualScores(QperNT, NTcounts); // Thread safe

    // Test this function
    //addNtSpecQualScores(d);

    //more memory intensive
    if (bLvsQlogs) {
        addLvsQlogs(len, avq, median);
    }
    stats_mutex.unlock();
    if (bMedianCalcs) {
        //quali
        addMedian2Histo(uint(avq + 0.5f), rstat_VQmed); // Thread safe (?)
        addMedian2Histo(len, rstat_VSmed); // Thread safe
    }
}


/*void ReportStats::mergeStats(ReportStats& data) {
    rstat_NTs += data.total_nts;
    rstat_totReads += data.total_reads;
    rstat_qualSum += data.qual_sum;
    rstat_accumError +=  data.accum_error;

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

    //listOfLengths.merge(data.listOfLengths);

}
*/

/*

 void ReportStats::addDNAStatsMT(shared_ptr<DNA> d, ReportStats *data){
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
*/



void ReportStats::reset() {
    rstat_totReads = 0; rstat_NTs = 0; rstat_qualSum = 0;
    rstat_Qmed = 0; rstat_Smed = 0;
    RSQS = 0.f; USQS = 0.f; rstat_accumError = 0.f;
    QperNT.resize(1000, 0); NTcounts.resize(1000, 0);
    std::fill(QperNT.begin(), QperNT.end(), 0);
    std::fill(NTcounts.begin(), NTcounts.end(), 0);
    rstat_VQmed.resize(0); rstat_VSmed.resize(0);
    listOfLengths.resize(0); listOfQuals.resize(0); listOfQualMeds.resize(0);
}
unsigned int ReportStats::lowest(const vector<uint>& in) {
    for (int i = 0; i < (int)in.size(); i++) {
        if (in[i] > 0) { return i; }
    }
    return(0);
}
unsigned int ReportStats::highest(const vector<uint>& in) {
    if (in.size() == 0) { return 0; }
    for (int i = (int)in.size() - 1; i >= 0; i--) {
        if (in[i] > 0) { return i; }
    }
    return 0;
}

void ReportStats::printLvsQ(ostream& give) {
    std::list<int>::iterator it1(listOfLengths.begin());
    std::list<float>::iterator it2(listOfQuals.begin());
    std::list<float>::iterator it3(listOfQualMeds.begin());
    give << "Length\tAvgQual\tMedianQual\n";

    for (; it1 != listOfLengths.end() && it2 != listOfQuals.end() && it3 != listOfQualMeds.end(); ++it1, ++it2, ++it3) {
        give << *it1 << "\t" << *it2 << "\t" << *it3 << "\n";
    }

}






void ReportStats::printGCstats(ostream& give) {
    //NT_POS['A'] = 0; NT_POS['T'] = 1; NT_POS['G'] = 2; NT_POS['C'] = 3;	NT_POS['N'] = 4;
    vector<string> NTs(6, "X"); NTs[0] = "A"; NTs[1] = "T";
    NTs[2] = "G"; NTs[3] = "C"; NTs[4] = "N";
    for (uint i = 0; i < 4; i++) {
        give << "\t" << NTcounts[i];
    }
    //give << endl;
    for (uint i = 0; i < 4; i++) {
        give << "\t" << float(QperNT[i]) / float(NTcounts[i]);
    }
    give << endl;
}
void ReportStats::printStats2(ostream& give, float remSeqs, int pair) {
    if (pair == 1) {
        return;//deactivate for now
    }
    if (pair == 0) {
        if (bMedianCalcs) {
            unsigned int minS = lowest(rstat_VSmed);
            unsigned int maxS = highest(rstat_VSmed);
            unsigned int minQ = lowest(rstat_VQmed);
            unsigned int maxQ = highest(rstat_VQmed);
            give << "Min/Avg/Max stats Pair 1";// -RSQS : "<<RSQS;
            if (remSeqs == 0) {
                give << "\n     - sequence Length : " << "0/0/0"
                    << "\n     - Quality :   " << "0/0/0";
            }
            else {
                give << "\n     - sequence Length : " << minS << "/" << float(rstat_NTs) / (float)rstat_totReads << "/" << maxS
                    << "\n     - Quality :   " << minQ << "/" << float(rstat_qualSum) / (float)rstat_totReads << "/" << maxQ;
            }
        }
        else {
            give << "Average Stats - RSQS : " << RSQS;
            if (remSeqs == 0) {
                give << "\n     - sequence Length : " << "0/0/0"
                    << "\n     - Quality :   " << "0/0/0";
            }
            else {
                give << "\n     - sequence Length : " << float(rstat_NTs) / (float)rstat_totReads
                    << "\n     - Quality :   " << float(rstat_qualSum) / (float)rstat_totReads;
            }
        }
    }
    else {
        give << "Pair 2 stats";
    }

    if (bMedianCalcs) {
        //give << "Median Stats Pair 1";// -USQS : " << USQS;
        give << "\n     - Median sequence Length : " << rstat_Smed << ", Quality : " << rstat_Qmed;// << "\n";
    }
    give << "\n     - Accum. Error " << (rstat_accumError / (float)rstat_totReads) << "\n";
}

//add the stats from a different ReportStats object
void ReportStats::addRepStats(ReportStats& RepoStat) {
    //report stats:
    rstat_NTs += RepoStat.rstat_NTs; rstat_totReads += RepoStat.rstat_totReads;
    rstat_qualSum += RepoStat.rstat_qualSum;
    rstat_accumError += RepoStat.rstat_accumError;
    for (uint i = 0; i < 6; i++) {
        QperNT[i] += RepoStat.QperNT[i];
        NTcounts[i] += RepoStat.NTcounts[i];
    }
    if (bMedianCalcs) {
        //vectors for median calcs
        if (rstat_VQmed.size() < RepoStat.rstat_VQmed.size()) {
            rstat_VQmed.resize(RepoStat.rstat_VQmed.size(), 0);
            assert(rstat_VQmed.size() < 10000);
        }
        for (unsigned int i = 0; i < RepoStat.rstat_VQmed.size(); i++) {
            rstat_VQmed[i] += RepoStat.rstat_VQmed[i];
        }
        if (rstat_VSmed.size() < RepoStat.rstat_VSmed.size()) {
            rstat_VSmed.resize(RepoStat.rstat_VSmed.size(), 0);
            //times change.. wrong assert here
            //assert(rstat_VSmed.size() < 10000);
        }
        for (unsigned int i = 0; i < RepoStat.rstat_VSmed.size(); i++) {
            rstat_VSmed[i] += RepoStat.rstat_VSmed[i];
        }
    }
    if (bLvsQlogs) {
        //l1.splice(l1.end(), l2);

        listOfLengths.splice(listOfLengths.end(), RepoStat.listOfLengths);
        listOfQuals.splice(listOfQuals.end(), RepoStat.listOfQuals);
        listOfQualMeds.splice(listOfQualMeds.end(), RepoStat.listOfQualMeds);
    }
}
//calculate median value from data stored as histogram-vector
// for median use perc = 0.5f
float ReportStats::calc_median(vector<uint>& in, float perc) {
    unsigned int sum = 0;
    for (unsigned int i = 0; i < in.size(); i++) {
        sum += in[i];
    }
    unsigned int threshold = (unsigned int)(((float)sum) * perc);
    sum = 0;
    for (unsigned int i = 0; i < in.size(); i++) {
        sum += in[i];
        if (sum >= threshold) {
            return (float)i;
        }
    }

    return 0.f;
}
void ReportStats::add_median2histo(vector<unsigned int>& in, vector<unsigned int>& histo)
{
    unsigned int max = *max_element(in.begin(), in.end());
    if (max > histo.size()) {
        if (max > 10000) { cerr << "max bigger 10000.\n"; exit(77); }
        histo.resize(max, 0);
    }
    for (unsigned int i = 0; i < histo.size(); i++) {
        histo[in[i]]++;
    }
}
void ReportStats::calcSummaryStats(float remSeqs, unsigned int min_l, float min_q) {
    if (remSeqs == 0) { return; }
    if (bMedianCalcs) {
        rstat_Smed = (int)calc_median(rstat_VSmed, 0.5f);
        rstat_Qmed = (int)calc_median(rstat_VQmed, 0.5f);
        USQS = 0.f;
    }
    RSQS = (((float(rstat_NTs) / remSeqs) / (float)min_l) +
        ((float(rstat_qualSum) / remSeqs) / min_q)) / 2.f;
}

vector<size_t> ReportStats::getVrange(int which) {
    if (which == 1) {
        return medVrange(rstat_VQmed);
    }
    else {
        return medVrange(rstat_VSmed);
    }
}
vector<size_t> ReportStats::medVrange(const vector<uint> x) {
    vector<size_t> ret(2, 0); ret[1] = 0; bool empty = true;
    for (size_t i = 0; i < x.size(); i++) {
        if (x[i] > 0) {
            if (empty) { ret[0] = i; empty = false; }
            ret[1] = i;
        }
    }
    ret[1]++;
    return ret;
}

void ReportStats::addMedian2Histo(const uint in, vector<unsigned int>& histo)
{
    if (in >= histo.size()) {
        stats_mutex.lock();
        histo.resize(in + 3, 0);
        stats_mutex.unlock();
        assert(in < 1000000);
    }
    histo[in]++;
}



void GAstats::addBaseGAStats(shared_ptr<DNA> dn, vector<shared_ptr<DNA>> GA, int missedGAs) {
    totalRds++; totalGAs += GA.size();
    totalRdLen += dn->length();
    missGAs += missedGAs;
    double totL(0.f); int incP(0);
    for (size_t i = 0; i < GA.size(); i++) {
        totL += (double)GA[i]->length();
        if (GA[i]->isYellowQual()) { incP++; }
        //if (GA[i]->length() > 2000) {cout << "TOO LONG:\n" << GA[i]->getSequence()<<endl;}
    }
    inccorrectPrimers+=incP;
    totalGALen += totL;
    int idx = dn->getBCnumber();
    if ((idx + 1) > GAperBC.size()) {
        GAperBC.resize(idx + 1, 0); CNTperBC.resize(idx + 1, 0);
        GALENperBC.resize(idx + 1, 0); rdLENperBC.resize(idx + 1, 0);
        GAfailsPerBC.resize(idx + 1, 0);
        inCorPrimPerBC.resize(idx+1, 0);
        missGAsPerBC.resize(idx+1, 0);
        GALDISperBC.resize(idx + 1, vector<int>(0));
        GALAMPLperBC.resize(idx + 1, vector<double>(0));
    }
    if (idx >= 0) {
        if (dn->isFailed()) {
            GAfailsPerBC[idx]++;
        } else {
            int GAL = GA.size();
            GAperBC[idx] += GAL;
            CNTperBC[idx]++;
            GALENperBC[idx] += totL;
            inCorPrimPerBC[idx] += incP;
            rdLENperBC[idx] += dn->length();
            if (GALDISperBC[idx].size() < (GAL+1)) {
                GALDISperBC[idx].resize(GAL + 1);
                GALAMPLperBC[idx].resize(GAL + 1);
            }
            GALDISperBC[idx][GAL] ++;
            GALAMPLperBC[idx][GAL] += totL;
            missGAsPerBC[idx] += missedGAs;
        }
    }
    else {

    }
}
void GAstats::setBCs(vector<string> SI, vector<string> B1, vector<string> B2) {
    assert(SI.size() >= GAperBC.size());
    SmplIDBC = SI; BC1 = B1; BC2 = B2;
    size_t idx = SI.size();
    GAperBC.resize(idx, 0); CNTperBC.resize(idx, 0);
    GALENperBC.resize(idx, 0); rdLENperBC.resize(idx, 0);
    GAfailsPerBC.resize(idx, 0);
    inCorPrimPerBC.resize(idx, 0);
    missGAsPerBC.resize(idx, 0);
    GALDISperBC.resize(idx , vector<int>(0));
    GALAMPLperBC.resize(idx , vector<double>(0));

}

void GAstats::addStats(shared_ptr<GAstats> o) {
    totalRds += o->totalRds;
    totalGAs += o->totalGAs;
    inccorrectPrimers += o->inccorrectPrimers;
    totalRdLen += o->totalRdLen;
    totalGALen += o->totalGALen;
    missGAs += o->missGAs;
    size_t idx = o->GAperBC.size();
    if (idx > GAperBC.size()) {
        GAperBC.resize(idx, 0); CNTperBC.resize(idx, 0);
        GALENperBC.resize(idx, 0); rdLENperBC.resize(idx, 0);
        GAfailsPerBC.resize(idx, 0);
        SmplIDBC.resize(idx, ""); BC1.resize(idx, ""); BC2.resize(idx, "");
        inCorPrimPerBC.resize(idx, 0);
        missGAsPerBC.resize(idx, 0);
        GALDISperBC.resize(idx, vector<int>(0));
        GALAMPLperBC.resize(idx, vector<double>(0));
    }
    for (idx = 0; idx < GAperBC.size(); idx++) {
        if (o->GAperBC.size() <= idx) { break; }
        GAperBC[idx] += o->GAperBC[idx];
        CNTperBC[idx] += o->CNTperBC[idx];
        GALENperBC[idx] += o->GALENperBC[idx];
        rdLENperBC[idx] += o->rdLENperBC[idx];
        GAfailsPerBC[idx] += o->GAfailsPerBC[idx];
        if (GALDISperBC[idx].size() < o->GALDISperBC[idx].size() ) {
            size_t oSze = o->GALDISperBC[idx].size();
            GALDISperBC[idx].resize(oSze );
            GALAMPLperBC[idx].resize(oSze );
        }
        for (size_t x = 0; x < GALDISperBC[idx].size(); x++) {
            GALDISperBC[idx][x] += o->GALDISperBC[idx][x];
            GALAMPLperBC[idx][x] += o->GALAMPLperBC[idx][x];
        }

    }

}


void GAstats::printBCAmpliNdistribution(ostream& give) {
    assert(GAperBC.size() == GALAMPLperBC.size()); assert(GAperBC.size() == GALDISperBC.size());
    int maxL(0);
    for (unsigned int i = 0; i < GALDISperBC.size(); i++) {
        if (GALDISperBC[i].size() > maxL) { maxL = GALDISperBC[i].size(); }
    }
    give << "SampleID";
    for (int i = 0; i < maxL; i++) { give << "\t#Amplis_" << i; }
    give << endl;

    for (unsigned int idx = 0; idx < GALDISperBC.size(); idx++) {
        give << SmplIDBC[idx];
        for (int i = 0; i < maxL; i++) {
            if (i < GALDISperBC[idx].size()) {
                give << "\t" << GALDISperBC[idx][i];
            }
            else {
                give << "\t0";
            }

        }
        give << endl;
    }

}

void GAstats::printBCAmpliLdistribution(ostream& give) {
    assert(GAperBC.size() == GALAMPLperBC.size()); assert(GAperBC.size() == GALDISperBC.size());
    int maxL(0);
    for (unsigned int i = 0; i < GALAMPLperBC.size(); i++) {
        if (GALAMPLperBC[i].size() > maxL) { maxL = (int)GALAMPLperBC[i].size(); }
    }
    give << "SampleID";
    for (uint i = 0; i < maxL; i++) { give << "\t#Amplis_" << i; }
    give << endl;

    for (unsigned int idx = 0; idx < GALAMPLperBC.size(); idx++) {
        give << SmplIDBC[idx];
        for (int i = 0; i < maxL; i++) {
            if (i < GALAMPLperBC[idx].size() && GALAMPLperBC[idx][i]>0) {
                give << "\t" << (GALAMPLperBC[idx][i]/(double)GALDISperBC[idx][i]/(double)i);
            }
            else {
                give << "\t0";
            }
        }
        give << endl;
    }

}



void GAstats::printBCtabs(ostream& give) {
    give << "SampleID";
    give << "\t#GA_Success\tGA_Fails\t#Ampli/Read\tmean_Read_Length\tmean_Ampli_Length\tIncorBCs\tMissedAmplics\n";
    assert(GAperBC.size() == CNTperBC.size()); assert(GAperBC.size() == GALENperBC.size());
    assert(GAperBC.size() == rdLENperBC.size()); assert(GAperBC.size() == SmplIDBC.size());
    assert(GAperBC.size() == GAfailsPerBC.size());
    for (unsigned int i = 0; i < GAperBC.size(); i++) {
        give << SmplIDBC[i] << "\t" << CNTperBC[i] << "\t" << GAfailsPerBC[i] << "\t";
        if (CNTperBC[i] > 0) {
            give << float(GAperBC[i]) / float(CNTperBC[i]) << "\t";
            give << int(float(rdLENperBC[i]) / float(CNTperBC[i]) + 0.5f) << "\t";
        } else { give <<"0\t0\t"; }
        if (GAperBC[i] > 0) {
            give << int(float(GALENperBC[i]) / float(GAperBC[i]) + 0.5f) << "\t";
        }else { give << "0\t"; }

        give<< inCorPrimPerBC[i] << "\t"<< missGAsPerBC[i];
        give << "\n";

    }

}

void GAstats::printSummary(ostream& give) {
    int totalFails = 0; for (size_t i = 0; i < GAfailsPerBC.size(); i++) { totalFails += GAfailsPerBC[i]; }
    give << "\nGoldenAxe Stats\n";
    give << "  Total reads; fails; inc. Primers; amplicons/read: " << totalRds << "; " << totalFails << "; "<< inccorrectPrimers<< "; " << (float)totalGAs / (float)totalRds << endl;
    give << "  Length reads; amplicon; missed Amplicons     : " << totalRdLen / (double)totalRds << "; " << totalGALen / (double)totalGAs << "; " << missGAs<< endl;
}




void Statistics::merge(Statistics &statistics, vector<int> &idx) {
    if (BarcodeDetected.size() > (uint) 10000) {
        cerr << "Unrealistic number of barcodes (>10000) in addStats\n";
        exit(79);
    }
    int BCS = (int) BarcodeDetected.size();
    for (unsigned int i = 0; i < idx.size(); i++) {
        if (idx[i] >= BCS) return;
        BarcodeDetected[idx[i]] += statistics.BarcodeDetected[i];
        BarcodeDetectedFail[idx[i]] += statistics.BarcodeDetectedFail[i];
    }
    maxL += statistics.maxL;
    PrimerFail += statistics.PrimerFail;
    AvgQual += statistics.AvgQual;
    HomoNT += statistics.HomoNT;
    HomoNTtrimmed += statistics.HomoNTtrimmed;
    PrimerRevFail += statistics.PrimerRevFail;
    minL += statistics.minL;
    minLqualTrim += statistics.minLqualTrim;
    TagFail += statistics.TagFail;
    MaxAmb += statistics.MaxAmb;
    QualWin += statistics.QualWin;
    Trimmed += statistics.Trimmed;
    AccErrTrimmed += statistics.AccErrTrimmed;
    total += statistics.total;
    QWinTrimmed += statistics.QWinTrimmed;
    totalRejected += statistics.totalRejected;
    fail_correct_BC += statistics.fail_correct_BC;
    suc_correct_BC += statistics.suc_correct_BC;
    failedDNAread += statistics.failedDNAread;
    adapterRem += statistics.adapterRem;
    RevPrimFound += statistics.RevPrimFound;
    singleton += statistics.singleton;
    BinomialErr += statistics.BinomialErr;
    dblTagFail += statistics.dblTagFail;
    DerepAddBadSeq += statistics.DerepAddBadSeq;
    total2 += statistics.total2;
    totalSuccess += statistics.totalSuccess;
}

void Statistics::reset() {
    //Former collectstats object
    singleton = 0;
    size_t BCsiz = BarcodeDetected.size();
    for (size_t i = 0; i < BCsiz; i++) {
        BarcodeDetected[i] = 0;
        BarcodeDetectedFail[i] = 0;
    }
    maxL = 0;
    PrimerFail = 0;
    AvgQual = 0;
    HomoNT = 0;
    HomoNTtrimmed = 0;
    PrimerRevFail = 0;
    minL = 0;
    TagFail = 0;
    MaxAmb = 0;
    QualWin = 0;
    Trimmed = 0;
    total = 0;
    totalRejected = 0;
    fail_correct_BC = 0;
    suc_correct_BC = 0;
    failedDNAread = 0;
    adapterRem = 0;
    RevPrimFound = 0;
    DerepAddBadSeq = 0;
    total2 = 0;
    totalSuccess = 0;

    // Former ReportStats object
    rstat_totReads = 0;
    rstat_NTs = 0;
    rstat_qualSum = 0;
    rstat_Qmed = 0;
    rstat_Smed = 0;
    RSQS = 0.f;
    USQS = 0.f;
    rstat_accumError = 0.f;
    QperNT.resize(1000, 0);
    NTcounts.resize(1000, 0);
    std::fill(QperNT.begin(), QperNT.end(), 0);
    std::fill(NTcounts.begin(), NTcounts.end(), 0);
    rstat_VQmed.resize(0);
    rstat_VSmed.resize(0);
}

float Statistics::GCcontent() {
    return float(NTcounts[2] + NTcounts[3]) / float(NTcounts[0] + NTcounts[1] + NTcounts[2] + NTcounts[3]);
}

void Statistics::addDNAStats(shared_ptr<DNA> d) {
    rstat_NTs += d->length();
    rstat_totReads++;
    rstat_qualSum += (int)(d->getAvgQual()+0.5f);
    rstat_accumError += (float)d->getAccumError();

    d->ntSpecQualScores(QperNT, NTcounts);

    if (bMedianCalcs) {
        //quali
        uint avq = (uint) (d->getAvgQual() + 0.5f);
        addMedian2Histo(avq, rstat_VQmed);
        addMedian2Histo(d->length(), rstat_VSmed); // Thread safe
    }
}

void Statistics::addMedian2Histo(unsigned int in, vector<unsigned int> &histo) {
    if (in >= histo.size()) {
        histo.resize(in + 3, 0);
        assert(in < 1e6);
    }
    histo[in] += 1;
}

inline float Statistics::calc_median(vector<unsigned int> &in, float perc) {
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
//
//void Statistics::addDNAStats(shared_ptr<DNA> d) {
//    //here should be the only place to count Barcodes!
//    int Pair = 1; // remove
//    int easyPair = Pair < 3 ? Pair - 1 : Pair - 3;
////    total2++;
////    if (d->isPassed() || d->isMidQual()) {
////        this->DNAstatLQ(d, easyPair, d->isMidQual());
////        totalSuccess++;
////    } else {
////        totalRejected++;
////    }
//
//    //some general stats that always apply:
//    if (d->QualCtrl.PrimerFwdFail) {
//        PrimerFail++;
//    }
//    if (d->QualCtrl.PrimerRevFail) {
//        PrimerRevFail++;
//    }
//    if (d->QualCtrl.minLqualTrim) {
//        minLqualTrim++;
//    }
//    if (d->QualCtrl.TagFail) {
//        TagFail++;
//    }
//    if (d->QualCtrl.fail_correct_BC) {
//        fail_correct_BC++;
//    }
//    if (d->QualCtrl.suc_correct_BC) {
//        suc_correct_BC++;
//    }
//    if (d->QualCtrl.RevPrimFound) {
//        RevPrimFound++;
//    }
//    if (d->QualCtrl.QWinTrimmed || d->QualCtrl.AccErrTrimmed) {
//        Trimmed++;
//    }
//    if (d->getTA_cut()) {
//        adapterRem++;
//    }
//
////exit(0);
//
//    if (d->isPassed() || d->isMidQual()) {
//        countBCdetected(d->getBCnumber(), easyPair, false);
//        //and register as success
//    } else {
//        if (d->getBarcodeDetected()) {
//            //DNA is no longer useful
//            failedStats2(d, easyPair);
//        }
//        //delete d;
//        if (d->QualCtrl.AvgQual) {
//            AvgQual++;
//        }
//        if (d->QualCtrl.minL) {
//            minL++;
//        }
//        if (d->QualCtrl.maxL) {
//            maxL++;
//        }
//        if (d->QualCtrl.HomoNT) {
//            HomoNT++;
//        }
//        if (d->QualCtrl.MaxAmb) {
//            MaxAmb++;
//        }
//        if (d->QualCtrl.BinomialErr) {
//            BinomialErr++;
//        }
//        if (d->QualCtrl.QualWin) {
//            QualWin++;
//        }
//    }
////    if (d->isDereplicated()) {
////        if (d->getBarcodeDetected() && !d->isPassed() && !d->isMidQual()) {
////            this->statAddDerepBadSeq(d->getBCnumber());
////        }
////    }
//}
