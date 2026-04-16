#include "InputStream.h"

#include <fstream>

using namespace std;

// Recover missing implementations when InputStream.cpp tail is unavailable

bool InputStreamer::setupFastq_2(string p1, string p2, string midp) {
    string file_type = "";
    size_t bufS = INPUT_BUFFER_SIZE;

    if (!p1.empty()) {
        file_type = "fastq file 1";
        fastq_istreams[0] = DBG_NEW ifbufstream(p1, (size_t)round(bufS * 0.8), doTIO);
        if (fastq_istreams[0]->eof()) {
            cerr << "\nWarning: Could not open or empty " << file_type << " " << p1 << " — skipping this input.\n";
            fastq_istreams[0] = nullptr;
            return false;
        }
    }
    if (!p2.empty()) {
        file_type = "fastq file 2";
        fastq_istreams[1] = DBG_NEW ifbufstream(p2, (size_t)round(bufS * 1.2), doTIO);
        if (fastq_istreams[1]->eof()) {
            cerr << "\nWarning: Could not open or empty " << file_type << " " << p2 << " — skipping this input.\n";
            fastq_istreams[1] = nullptr;
            return false;
        }
    }
    if (!midp.empty()) {
        this->openMIDseqs("", midp);
    }
    return true;
}

string InputStreamer::setupInput(string path, int t, const string& uniqueFastxFile,
    filesStr& files, int& paired, string onlyPair,
    string& mainFilename, bool simulate) {
    string mainFilepath("");
    vector<string> fastqFiles = files.FastqF;
    vector<string> fastaFiles = files.FastaF;
    vector<string> qualityFiles = files.QualF;
    vector<string> midFiles = files.MIDfq;

    if (isFasta) {
        if (fastaFiles[t] != uniqueFastxFile) {
            cerr << "Error in matching FASTA target filenames.\n";
            exit(11);
        }
        this->setupFastaQual(path, fastaFiles[t], qualityFiles[t], paired, onlyPair);
        mainFilepath = path + fastaFiles[t];
        mainFilename = fastaFiles[t];
    }
    else {
        if (fastqFiles[t] != uniqueFastxFile) {
            cerr << "Error in matching target filenames.\n";
            exit(11);
        }
        this->setupFastq(path, fastqFiles[t], paired, onlyPair, simulate, !simulate);
        mainFilepath = path + fastqFiles[t];
        mainFilename = fastqFiles[t];
    }
    if ((size_t)t < midFiles.size()) {
        this->openMIDseqs(path, midFiles[t]);
    }
    return mainFilepath;
}

void InputStreamer::setupFna(string in) {
    int paired = 1;
    if (detectSeqFmt(in) == "-i_fna") {
        setupFastaQual("", in, "", paired, "", false);
    }
    else {
        setupFastq("", in, paired, "", false, false);
    }
}

int InputStreamer::auto_fq_version() {
    if (fastQver != 0) {
        return fastQver;
    }
    fastQver = (qual_score)auto_fq_version(minQScore, maxQScore);
    return fastQver;
}

void InputStreamer::maxminQualWarns_fq() {
    // no-op fallback; warnings are optional
}

ofbufstream::ofbufstream(const string IF, int mif, bool isMC, size_t bufferS)
    : file(IF), keeper(bufferS), keeperW(bufferS), modeIO(mif), used(0), usedW(0),
      coutW(false), isGZ(false), doMC(isMC), primary(nullptr), bufS(bufferS), hasKickoff(false) {
    if (bufS == 0) bufS = 20000;
    if (file == "" || file == "-") {
        coutW = true;
        return;
    }
    activate();
}

ofbufstream::~ofbufstream() {
    finishWrites();
    deactivate();
}

void ofbufstream::finishWrites() {
    if (hasKickoff) {
        writeKickoff.wait();
        hasKickoff = false;
    }
    emptyStream();
}

bool ofbufstream::operator! (void) {
    if (coutW) return false;
    if (!primary) return true;
    return !(*primary);
}

void ofbufstream::operator<< (const string& X) {
    if (X.empty()) return;
    if (coutW) {
        cout << X;
        return;
    }
    if (bufS == 0) {
        if (primary) (*primary) << X;
        return;
    }

    if (used + X.size() > keeper.size()) {
        emptyStream();
    }

    if (X.size() > keeper.size()) {
        if (primary) primary->write(X.data(), (streamsize)X.size());
        return;
    }

    memcpy(keeper.data() + used, X.data(), X.size());
    used += X.size();
}

void ofbufstream::emptyStream() {
    if (coutW) return;
    if (!primary || used == 0) return;
    primary->write(keeper.data(), (streamsize)used);
    used = 0;
}

void ofbufstream::activate() {
    if (coutW) return;
    if (primary) return;
    ios_base::openmode m = ios::out;
    if (modeIO == ios::app) m |= ios::app;
    else m |= ios::trunc;
    primary = std::make_unique<std::ofstream>(file.c_str(), m | ios::binary);
}

void ofbufstream::deactivate() {
    if (coutW) return;
    if (primary) {
        primary->flush();
        primary.reset();
    }
}

void ofbufstream::writeStream(bool doKickoff) {
    (void)doKickoff;
    emptyStream();
}

dualOfBufStream::dualOfBufStream(void)
    : buf1S(20000), buf2S(20000), bufs(2, ""), FileNames(2, ""), ostr(2), opened(2, false), active(false) {
}

dualOfBufStream::~dualOfBufStream(void) {
    emptyStreams(true);
    deactivate();
}
