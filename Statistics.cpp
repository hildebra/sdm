//
// Created by fritsche on 10/11/2020.
//

#include "Statistics.h"

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
