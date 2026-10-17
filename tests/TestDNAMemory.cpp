// Reproducible production-object benchmark. Compatible with the pre-change API.
#include "DNA.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
uint64_t digest = 1469598103934665603ULL;
void consume(std::string_view value) {
    for (unsigned char c : value) { digest = (digest ^ c) * 1099511628211ULL; }
}
std::shared_ptr<DNA> read(size_t i, size_t length, int quality = 30) {
    std::string sequence(length, 'A');
    uint64_t state = i + 1234567;
    for (char& c : sequence) { state ^= state << 13; state ^= state >> 7; state ^= state << 17; c = "ACGT"[state & 3]; }
    std::vector<std::string> lines{"@instrument:flowcell:lane:tile:coordinate:" + std::to_string(i) + " long metadata retained in input headers", std::move(sequence), "+", std::string(length, char(33 + quality))};
    return std::make_shared<DNA>(std::move(lines), 33);
}
long peakRSS() {
    std::ifstream status("/proc/self/status"); std::string line;
    while (std::getline(status, line)) { if (line.starts_with("VmHWM:")) { return std::stol(line.substr(6)); } }
    return 0;
}
}
int main(int argc, char** argv) {
    if (argc != 4) { std::cerr << "Usage: TestDNAMemory <dna|stream|equal|variable|history|history_raw> <objects> <length>\n"; return 2; }
    const std::string mode = argv[1]; const size_t count = std::stoull(argv[2]), length = std::stoull(argv[3]);
    if (length < 20 || count == 0) { return 2; }
    ini_DNAconstants();
    const auto start = std::chrono::steady_clock::now();
    std::vector<std::shared_ptr<DNA>> reads;
    std::vector<std::shared_ptr<DNAunique>> unique;
    if (mode == "dna") { reads.reserve(count); }
    else if (mode != "stream") { unique.reserve(count); }
    for (size_t i = 0; i < count; ++i) {
        auto dna = read(i, length);
        if (mode == "dna") { reads.push_back(std::move(dna)); continue; }
        if (mode == "stream") {
            dna->cutSeq(0, 7); dna->cutSeq(10, 13);
            auto sub = dna->getDNAsubseq(0, 10, "subsequence");
            consume(sub->getSequence());
            if (!(i % 8)) { dna->setNewID("edited_" + std::to_string(i)); }
            dna->prepareWrite(33); consume(dna->writeFastQ());
            continue;
        }
        auto group = std::make_shared<DNAunique>(dna, 0); group->saveMem();
        if (mode == "equal" || mode == "variable") {
            for (size_t observation = 0; observation < 4; ++observation) {
                dna->resetTruncation();
                if (mode == "variable" && observation % 2) { dna->cutSeqPseudo(length - 7); }
                dna->setAllQual(qual_score(10 + observation * 10));
                group->sumQualities(dna); group->incrementSampleCounter(0);
            }
        } else if (mode == "history" || mode == "history_raw") {
            auto merged = read(i, length * 2);
            group->attachMerge(std::make_shared<DNAunique>(merged, -1));
            if (mode == "history") { group->getMerge()->sumQualities(merged); }
            group->takeOverDNA(dna, nullptr, nullptr);
        } else { throw std::invalid_argument("Unknown benchmark mode"); }
        unique.push_back(std::move(group));
    }
    for (const auto& dna : reads) { dna->prepareWrite(33); consume(dna->writeFastQ()); }
    for (const auto& group : unique) {
        auto output = group->cloneForOutput();
        if (mode == "history" || mode == "history_raw") {
            auto merged = read(0, length * 2, 40);
            output->takeOverDNA(read(0, length), nullptr, merged);
            auto consensus = output->getMerge();
            consensus->releaseRawQualitiesIfPromoted(); consensus->prepareWrite(33);
            consume(consensus->writeFastQ());
        } else {
            output->releaseRawQualitiesIfPromoted(); output->prepareWrite(33);
            consume(output->writeFastQ());
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "mode,objects,length,seconds,peak_rss_kib,checksum,dna_bytes,unique_bytes\n"
        << mode << ',' << count << ',' << length << ',' << std::setprecision(9) << seconds << ',' << peakRSS() << ',' << digest << ',' << sizeof(DNA) << ',' << sizeof(DNAunique) << '\n';
}
