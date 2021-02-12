//
// Created by fritsche on 10/11/2020.
//

#pragma once


#include "DNAconsts.h"
#include "InputStream.h"

class Statistics {
public:
    Statistics() : maxL(0), PrimerFail(0), AvgQual(0), HomoNT(0),
                   PrimerRevFail(0), minL(0), minLqualTrim(0),
                   TagFail(0), MaxAmb(0), QualWin(0),
                   Trimmed(0), AccErrTrimmed(0), QWinTrimmed(0),
                   total(0), totalRejected(0),
                   fail_correct_BC(0), suc_correct_BC(0),
                   failedDNAread(0), adapterRem(0), RevPrimFound(0),
                   total2(0), totalSuccess(0),
                   DerepAddBadSeq(0), BinomialErr(0), dblTagFail(0),
                   singleton(0), BarcodeDetected(0), BarcodeDetectedFail(0) {}

    unsigned int maxL, PrimerFail, AvgQual, HomoNT;
    unsigned int PrimerRevFail; //Number of sequences, where RevPrimer was detected (and removed)
    unsigned int minL, minLqualTrim, TagFail, MaxAmb, QualWin;
    unsigned int Trimmed, AccErrTrimmed, QWinTrimmed, total, totalRejected;
    unsigned int fail_correct_BC, suc_correct_BC, failedDNAread;
    unsigned int adapterRem, RevPrimFound;
    uint totalSuccess;
    uint DerepAddBadSeq;
    //binomial error model
    unsigned int BinomialErr;
    uint dblTagFail;
    //recovered singletons within pairs
    unsigned int singleton;
    vector<int> BarcodeDetected;
    vector<int> BarcodeDetectedFail;

    // ReportStats
    bool bMedianCalcs;

    const vector<unsigned int> &get_rstat_Vmed(int x) {
        if (x == 1) { return rstat_VQmed; } else { return rstat_VSmed; }
    }

    //median
    unsigned int rstat_totReads, rstat_NTs, rstat_qualSum, rstat_Qmed, rstat_Smed;
    //means, Relative sample_id_ Quality Score (RSQS), Unifying sample_id_ Quality Score (USQS)
    float RSQS, USQS;
    float rstat_accumError;
    vector<long> QperNT, NTcounts;
    vector<unsigned int> rstat_VQmed, rstat_VSmed;


    void addDNAStats(shared_ptr<DNA> d);

    void calcSummaryStats(float remSeqs, unsigned int min_l, float min_q);

    static float calc_median(vector<unsigned int> &in, float perc);

    static void addMedian2Histo(unsigned int in, vector<unsigned int> &histo);

    unsigned int lowest(const vector<uint> &in);

    unsigned int highest(const vector<uint> &in);

    void printStats2(ostream &give, float remSeqs, int pair);

    void printGCstats(ostream &give);

    vector<size_t> medVrange(const vector<uint>);

    vector<size_t> getVrange(int which);

    float GCcontent();
    void merge(Statistics &, vector<int> &idx);
    void reset();
    uint total2;
};
