//
// Created by fritsche on 10/11/2020.
//

#pragma once


#include "DNAconsts.h"
#include "InputStream.h"
#include <atomic>

class Statistics {
public:
    Statistics() : maxL(0), PrimerFail(0), AvgQual(0), HomoNT(0), HomoNTtrimmed(0), 
                   PrimerRevFail(0), minL(0), minLqualTrim(0),
                   TagFail(0), MaxAmb(0), QualWin(0),
                   Trimmed(0), AccErrTrimmed(0), QWinTrimmed(0),
                   total(0), totalRejected(0),
                   fail_correct_BC(0), suc_correct_BC(0),
                   failedDNAread(0), adapterRem(0), RevPrimFound(0),
                   total2(0), totalSuccess(0),
                   DerepAddBadSeq(0), BinomialErr(0), dblTagFail(0),
                   singleton(0), BarcodeDetected(0), BarcodeDetectedFail(0),
                   rstat_totReads(0), rstat_NTs(0), rstat_qualSum(0), rstat_Qmed(0), rstat_Smed(0),
                   RSQS(0.f), USQS(0.f), rstat_accumError(0.f),
                   rstat_VQmed(0), rstat_VSmed(0){}

    unsigned int maxL, PrimerFail, AvgQual, HomoNT, HomoNTtrimmed;
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

    const vector<unsigned int>& get_rstat_Vmed(int x);

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

class MEstats {
public:
    MEstats() :total_read_preMerge_(0), merged_counter_(0) {}
    ~MEstats() {}
    void addStats(shared_ptr<MEstats> o);

    void print(ostream& give);
    //variables
    int total_read_preMerge_, merged_counter_;
    uint BPwritten, BPmergeWritte;


};



//reported stats on sequence properties
class ReportStats {
public:
    ReportStats(bool MedianDo = true) :
        bMedianCalcs(MedianDo), bLvsQlogs(false), rstat_totReads(0), rstat_NTs(0), rstat_qualSum(0),
        rstat_Qmed(0), rstat_Smed(0), RSQS(0.f), USQS(0.f), rstat_accumError(0.f),
        QperNT(6, 0), NTcounts(6, 0),
        rstat_VQmed(0), rstat_VSmed(0),
        listOfLengths(0), listOfQuals(0), listOfQualMeds(0)
    {
    }

    ~ReportStats() {
        cdbg("~ReportStats: N DNAs: " + to_string(rstat_totReads) + " ");
    }
    void reset();
    void addDNAStats(shared_ptr<DNA> d);
    //void mergeStats(data_MT &data);
    void setbLvsQlogs(bool b);
    bool getbLvsQlogs();
    //void addDNAStatsMT(shared_ptr<DNA> d, data_MT *data);
    void calcSummaryStats(float remSeqs, unsigned int min_l, float min_q);
    float calc_median(vector<uint>& in, float perc);
    void add_median2histo(vector<unsigned int>& in, vector<unsigned int>& histo);
    void addMedian2Histo(unsigned int in, vector<unsigned int>& histo);
    void addMeanStats(unsigned int NT, unsigned int Qsum, float AccErr);
    void addLvsQlogs(uint L, float avg, int med);
    // TEST IF PRODUCES SAME RESULTS
    void addNtSpecQualScores(shared_ptr<DNA> dna);

    unsigned int lowest(const vector<uint>& in);
    unsigned int highest(const vector<uint>& in);
    void printStats2(ostream& give, float remSeqs, int pair);
    void printGCstats(ostream& give);
    void printLvsQ(ostream& give);
    void addRepStats(ReportStats&);
    bool bMedianCalcs;
    bool bLvsQlogs;
    const vector<unsigned int>& get_rstat_Vmed(int x);
    //const vector<unsigned int> &get_rstat_VSmed(){return rstat_VSmed;}
    vector<size_t> getVrange(int which);

protected:

    //median
    vector<size_t> medVrange(const vector<uint>);
    unsigned long rstat_totReads, rstat_NTs, rstat_qualSum, rstat_Qmed, rstat_Smed;
    //means, Relative sample_id_ Quality Score (RSQS), Unifying sample_id_ Quality Score (USQS)
    float RSQS, USQS;
    double rstat_accumError;
    vector<long> QperNT, NTcounts;
    float GCcontent();

    //bin based median calculation's
    vector<unsigned int> rstat_VQmed, rstat_VSmed;
    std::list<int> listOfLengths;
    std::list<float> listOfQuals;
    std::list<float> listOfQualMeds;

private:
    std::mutex stats_mutex;
};


class collectstats {
public:
    collectstats() : maxL(0), PrimerFail(0), AvgQual(0), HomoNT(0), HomoNTtrimmed(0), 
        PrimerRevFail(0), minL(0), minLqualTrim(0),
        TagFail(0), MaxAmb(0), QualWin(0),
        Trimmed(0), AccErrTrimmed(0), QWinTrimmed(0),
        total(0), totalMid(0), totalRejected(0),
        fail_correct_BC(0), suc_correct_BC(0),
        failedDNAread(0), adapterRem(0), RevPrimFound(0),
        total2(0), totalSuccess(0),
        DerepAddBadSeq(0), BinomialErr(0),
        dblTagFail(0),
        reversedRds(0), swappedRds(0),
        singleton(0), BarcodeDetected(0), BarcodeDetectedFail(0),
        PostFilt(ReportStats()), PreFilt(ReportStats())
    {
        cdbg("Ini collectstats");
        //PostFilt = DBG_NEW ReportStats();
        //PreFilt = DBG_NEW ReportStats();


    }
    //collectstats(const collectstats&) = default;
    //collectstats& operator=(const collectstats&) = default;


    ~collectstats() {
        //delete PostFilt; PostFilt = nullptr;
        //delete PreFilt; PreFilt = nullptr;
    }

    void addPostFilt(shared_ptr<DNA> d);
    void addPreFilt(shared_ptr<DNA> d);

    unsigned int maxL, PrimerFail, AvgQual, HomoNT, HomoNTtrimmed;
    unsigned int PrimerRevFail; //Number of sequences, where RevPrimer was detected (and removed)
    unsigned int minL, minLqualTrim, TagFail, MaxAmb, QualWin;
    unsigned int Trimmed, AccErrTrimmed, QWinTrimmed, totalMid;
    std::atomic<unsigned int> total;
    std::atomic<unsigned int> totalRejected;
    unsigned int fail_correct_BC, suc_correct_BC, failedDNAread;
    unsigned int adapterRem, RevPrimFound;
    std::atomic<unsigned int> total2;
    unsigned int totalSuccess;
    uint DerepAddBadSeq;
    //binomial error model
    unsigned int BinomialErr;
    uint dblTagFail;
    //swapping/reversing reads
    uint reversedRds;
    uint swappedRds;
    //recovered singletons within pairs
    unsigned int singleton;
    vector<int> BarcodeDetected;
    vector<int> BarcodeDetectedFail;
    void addStats(shared_ptr<collectstats>, vector<int>& idx);
    void reset();
    //void ini_repStat(bool midQ) {}
    void ini_repStat(void);
    void setbLvsQlogsPreFilt(bool b);
    bool getbLvsQlogsPreFilt();

    ReportStats PostFilt;//green
    ReportStats PreFilt; //before filtering
};


class GAstats {
public:
    GAstats() :totalRds(0), totalGAs(0), inccorrectPrimers(0), missGAs(0),
        totalRdLen(0.f), totalGALen(0.f),
        GAperBC(0), CNTperBC(0), GALENperBC(0), rdLENperBC(0), GAfailsPerBC(0), 
        inCorPrimPerBC(0), missGAsPerBC(0),
        SmplIDBC(0), BC1(0), BC2(0)
    {
    }
    ~GAstats() {}
    void reset();
    void setBCs(vector<string> SI, vector<string> B1, vector<string> B2);
    void addStats(shared_ptr<GAstats> o);
    void addMissedGAs(int X);

    void printSummary(ostream& give);//summary of GA stats
    void printBCtabs(ostream& give); //stats per BC
    void printBCAmpliNdistribution(ostream& give); //distribution of lengths of BCs
    void printBCAmpliLdistribution(ostream& give); //distribution of lengths of BCs
    void addBaseGAStats(shared_ptr<DNA> dn, vector<shared_ptr<DNA>> GA, int missedGAs);
    uint totalRds, totalGAs; //total reads, total amplicons
    uint inccorrectPrimers, missGAs;
    double totalRdLen, totalGALen; //length in bp
    vector<int> GAperBC, CNTperBC, GALENperBC, rdLENperBC, GAfailsPerBC;
    vector<int> inCorPrimPerBC, missGAsPerBC;
    vector<string> SmplIDBC, BC1, BC2;
    vector<vector<int>> GALDISperBC;
    vector<vector<double>>  GALAMPLperBC;

};
