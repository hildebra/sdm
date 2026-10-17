#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#if defined(__linux__)
#include <unistd.h>
#endif

#include "DNA.h"
#include "DNAconsts.h"
#include "Filters.h"
#include "InputStream.h"
#include "IO.h"
#include "OutputStreamer.h"
#include "ReadMerger.h"
#include "Statistics.h"
#include "containers.h"

using namespace std;

static_assert(is_same<decltype(declval<Filters&>().getReadsWritten()), uint64_t>::value,
    "filter read counters must expose a 64-bit value");
static_assert(is_same<decltype(declval<OutputStreamer&>().getReadsWritten()), uint64_t>::value,
    "output read counters must expose a 64-bit value");
static_assert(is_same<decltype(declval<collectstats&>().total.load()), uint64_t>::value,
    "cumulative statistics counters must expose a 64-bit value");
static_assert(is_same<sample_count_t, uint64_t>::value,
    "dereplication sample counts must use 64-bit storage");
static_assert(is_same<DNAunique::qual_accum_t, uint64_t>::value,
    "dereplication quality sums must use 64-bit storage");

namespace {

int failures = 0;

void check(bool condition, const string& message) {
    if (!condition) {
        cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

template <class Func>
void checkThrows(Func&& func, const string& message) {
    try {
        func();
        check(false, message + " (no exception was thrown)");
    }
    catch (const runtime_error&) {
    }
}

void writeFile(const string& path, const string& contents) {
    ofstream out(path, ios::binary);
    if (!out) {
        throw runtime_error("Could not create regression-test file: " + path);
    }
    out << contents;
}

string readFile(const string& path) {
    ifstream in(path, ios::binary);
    if (!in) {
        throw runtime_error("Could not read regression-test file: " + path);
    }
    return string(istreambuf_iterator<char>(in), istreambuf_iterator<char>());
}

bool testFileExists(const string& path) {
    ifstream in(path, ios::binary);
    return static_cast<bool>(in);
}

size_t countOccurrences(const string& text, const string& needle) {
    size_t count = 0;
    size_t pos = 0;
    while ((pos = text.find(needle, pos)) != string::npos) {
        ++count;
        pos += needle.size();
    }
    return count;
}

vector<string> fastaHeaders(const string& path) {
    ifstream in(path, ios::binary);
    if (!in) {
        throw runtime_error("Could not read regression-test file: " + path);
    }
    vector<string> headers;
    string line;
    while (getline(in, line)) {
        if (!line.empty() && line.front() == '>') {
            headers.push_back(line.substr(1));
        }
    }
    return headers;
}

class TestFilters : public Filters {
public:
    explicit TestFilters(OptContainer* options) : Filters(options) {}
    using Filters::cutPrimerRev;
};

shared_ptr<DNA> makeDNA(const string& sequence, const string& id = "read") {
    auto dna = make_shared<DNA>(sequence, id);
    dna->setQual(vector<qual_score>(sequence.size(), 30));
    check(dna->seal(), "test DNA should seal successfully");
    return dna;
}

OptContainer makeSeedOutputOptions(const string& r1Path, const string& r2Path);

void testFastaQualityParsing() {
    vector<string> valid = {">read", "ACGT", "10 20 30 40"};
    DNA dna(valid);
    check(dna.getQual().size() == 4, "FASTA quality count should match sequence length");
    check(dna.getQual()[0] == 10 && dna.getQual()[3] == 40,
        "FASTA quality values should be parsed exactly");

    vector<string> tooFew = {">short", "ACGT", "10 20 30"};
    checkThrows([&] { DNA invalid(tooFew); },
        "FASTA input with too few quality scores should be rejected");

    vector<string> invalidToken = {">bad", "ACGT", "10 20 nope 40"};
    checkThrows([&] { DNA invalid(invalidToken); },
        "FASTA input with a non-numeric quality should be rejected");

    vector<string> negative = {">negative", "AC", "-5 20"};
    checkThrows([&] { DNA invalid(negative); },
        "FASTA input with a negative Phred quality should be rejected");

    vector<string> tooHigh = {">too_high", "AC", "20 110"};
    checkThrows([&] { DNA invalid(tooHigh); },
        "FASTA input above the supported Phred range should be rejected");
}

void testMedianQualityDoesNotMutateRead() {
	DNA empty;
	check(empty.getMedianQual() == 0,
		"median quality of a read without qualities should be zero");

    DNA dna("ACGTA", "median_read");
    const vector<qual_score> qualities = {5, 40, 10, 30, 20};
    dna.setQual(qualities);
    check(dna.seal(), "median test DNA should seal successfully");

    check(dna.getMedianQual() == 20, "median quality should be calculated correctly");
    check(dna.getQual() == qualities,
        "getMedianQual must not reorder the read's per-base qualities");
    check(dna.getMedianQual() == 20 && dna.getQual() == qualities,
        "repeated median-quality queries must be stable and non-mutating");

    check(dna.cutSeqPseudo(3), "median test should accept a logical 3' trim");
    check(dna.getMedianQual() == 10,
        "median quality must exclude qualities beyond the logical sequence length");
    check(dna.getQual() == qualities,
        "logical-length median calculation must leave backing qualities untouched");
}

void testQualityCachesFollowMutations() {
    DNA dna("ACGT", "quality_cache_read");
    dna.setQual({40, 40, 0, 0});
    check(dna.seal(), "quality-cache test DNA should seal successfully");

    const double fullError = dna.getAccumError();
    check(std::fabs(dna.getAvgQual() - 20.f) < 0.001f,
        "average quality should include the complete logical read");

    check(dna.cutSeqPseudo(2), "quality-cache test should accept a logical trim");
    check(std::fabs(dna.getAvgQual() - 40.f) < 0.001f,
        "logical trimming must invalidate the cached average quality");
    check(dna.getAccumError() < fullError,
        "logical trimming must invalidate the cached accumulated error");

    dna.resetTruncation();
    dna.setAllQual(10);
    check(std::fabs(dna.getAvgQual() - 10.f) < 0.001f,
        "quality mutation must invalidate the cached average quality");
	const double expectedMutatedError = 4.0 * SAqualP[10];
    check(std::fabs(dna.getAccumError() - expectedMutatedError) < 1e-12,
        "quality mutation must invalidate the cached accumulated error");
}

void testQualityMathBoundsAndFractionalWindow() {
    DNA bounded("ACG", "bounded_quality_read");
    bounded.setQual({static_cast<qual_score>(-5), static_cast<qual_score>(110),
        static_cast<qual_score>(127)});
    check(bounded.seal(), "bounded-quality test DNA should seal successfully");
    const double expected = SAqualP[0] + (2.0 * SAqualP[maxSAqualP]);
    check(std::fabs(bounded.getAccumError() - expected) < 1e-12,
        "accumulated-error lookup must clamp every quality to its table domain");
    check(bounded.qualAccumulate(100.0) == 3,
        "quality accumulation must safely consume out-of-domain stored qualities");

    DNA windowed("AC", "fractional_window_read");
    windowed.setQual({0, 0});
    check(windowed.seal(), "quality-window test DNA should seal successfully");
    int reason = -1;
    check(windowed.qualWinfloat(2, 0.9f, reason) == 0.f && reason == 1,
        "fractional quality-window thresholds must not be truncated to zero");
}

void testLazyHeadersAndPhysicalTrimming() {
    auto dna = makeDNA("ACGTACGTACGT", string(90, 'h') + " description");
    const string original = dna->getOldId();
    dna->prepareWrite(33);
    check(dna->writeFastQ().starts_with("@" + original + "\n"), "Unedited lazy header is emitted by FASTQ writers");
    check(dna->writeSeq(true).starts_with(">" + original + "\n"), "Unedited lazy header is emitted by FASTA writers");
    dna->setNewID("");
    check(dna->writeFastQ().starts_with("@\n") && dna->getOldId() == original,
        "An explicitly empty new header is not the original header");
    dna->setHeader("new/1"); dna->normalizePositionFreeId();
    check(dna->writeFastQ().starts_with("@new\n"), "Normalization resets the lazy header");
    DNA copied(*dna); copied.setNewID("snapshot");
    check(dna->writeFastQ().starts_with("@new\n") && copied.writeFastQ().starts_with("@snapshot\n"),
        "Copied lazy headers remain independently editable");
    dna->cutSeqPseudo(10); dna->cutSeq(3, 6);
    check(dna->getSequence() == "ACGGTACGT" && dna->length() == 7 && dna->getQual().size() == 9,
        "Physical trimming preserves the backing and logical lengths");
    dna->resetTruncation();
    check(dna->length() == 9, "Pseudo-trimmed tail survives physical compaction");
    DNAunique saved(dna, 0); saved.saveMem(); saved.saveMem(); saved.Count2Head(true); saved.prepareWrite(33);
    check(saved.writeFastQ().starts_with("@new;size=1;\n"), "Repeated storage compaction retains the output header");
}

void testUniformCoverageAndQualityHistory() {
    auto initial = makeDNA(string(31, 'A'), "coverage"); initial->setAllQual(20);
    DNAunique group(initial, 0);
    vector<uint64_t> sums(31, 20), coverage(31, 1);
    // Includes negative input scores, zero-length observations, changing lengths,
    // raw representative replacement, and independent output snapshots.
    for (size_t observation = 0; observation < 80; ++observation) {
        const size_t n = observation < 8 ? 31 : (observation * 17) % 53;
        auto read = makeDNA(string(n, 'A'), "replacement");
        vector<qual_score> qualities(n);
        for (size_t i = 0; i < n; ++i) { qualities[i] = qual_score(int((i + observation * 7) % 100) - 5); }
        read->setQual(qualities);
        group.sumQualities(read);
        sums.resize(max(sums.size(), n), 0); coverage.resize(sums.size(), 0);
        for (size_t i = 0; i < n; ++i) { sums[i] += max<int>(0, qualities[i]); ++coverage[i]; }
        if (n && observation % 7 == 0) { group.takeOverDNA(read, nullptr, nullptr); }
        check(group.quality_sum_per_base_ == sums, "Compressed coverage retains every quality sum");
        for (size_t i = 0; i < coverage.size(); ++i) {
            check(group.qualityCoverage(i) == coverage[i], "Scalar and nonuniform coverage match an independent per-base oracle");
        }
        if (observation % 11 == 0) {
            auto snapshot = group.cloneForOutput(); snapshot->releaseRawQualitiesIfPromoted();
            snapshot->prepareWrite(33);
            string expected;
            for (size_t i = 0; i < snapshot->length(); ++i) {
                expected += static_cast<char>(std::clamp(int(round(double(sums[i]) / coverage[i])) + 33, 33, 126));
            }
            check(snapshot->writeFastQ().ends_with("\n+\n" + expected + "\n"),
                "Output snapshot uses exact per-position coverage and rounding");
        }
    }
    auto singleton = makeDNA(string(133, 'A'), "singleton");
    vector<qual_score> allScores(133);
    for (int i = 0; i < 133; ++i) { allScores[i] = static_cast<qual_score>(i - 5); }
    singleton->setQual(allScores);
    for (int offset : {33, 64}) {
        singleton->prepareWrite(offset);
        const string rawOutput = singleton->writeFastQ();
        DNAunique promoted(singleton, 0);
        promoted.releaseRawQualitiesIfPromoted(); promoted.prepareWrite(offset);
        check(promoted.writeFastQ() == rawOutput,
            "Singleton evidence preserves raw writer output for negative and high qualities at both offsets");
    }
    DNAunique attached(initial, 0);
    auto heldMerge = make_shared<DNAunique>(initial, -1);
    heldMerge->sumQualities(initial);
    attached.attachMerge(heldMerge);
    const auto heldSums = heldMerge->quality_sum_per_base_;
    attached.takeOverDNA(initial, nullptr, nullptr);
    check(heldMerge->quality_sum_per_base_ == heldSums,
        "Detaching history does not clear a former merge still owned by its caller");
    DNAQualityEvidence zero;
    zero.initializeQualityEvidence(vector<qual_score>(3, 0), 3);
    for (unsigned bit = 0; bit < 63; ++bit) { const auto copy = zero; zero.addQualityEvidence(copy); }
    check(zero.qualityCoverage(0) == (uint64_t(1) << 63), "Uniform coverage retains 64-bit depth");
    const auto huge = zero;
    checkThrows([&] { zero.addQualityEvidence(huge); }, "Uniform coverage overflow is rejected");
    // Trigger vector coverage, then verify its overflow path independently.
    zero.addQualityEvidence(vector<qual_score>(2, 0), 2);
    const auto variableHuge = zero;
    checkThrows([&] { zero.addQualityEvidence(variableHuge); }, "Nonuniform coverage overflow is rejected");

    auto merged = makeDNA("ACGTAC", "historical"); merged->setAllQual(10);
    auto replacementMerge = makeDNA("ACGTAC", "current"); replacementMerge->setAllQual(40);
    auto root = makeDNA("ACGT", "root");
    DNAunique unpromoted(root, 0);
    unpromoted.attachMerge(make_shared<DNAunique>(merged, -1));
    unpromoted.takeOverDNA(root, nullptr, nullptr);
    auto snapshot = unpromoted.cloneForOutput();
    unpromoted.matchedDNA(root, nullptr, replacementMerge, -1, false, true);
    unpromoted.takeOverDNA(root, nullptr, replacementMerge);
    snapshot->takeOverDNA(root, nullptr, replacementMerge);
    auto output = unpromoted.getMerge(); output->releaseRawQualitiesIfPromoted(); output->prepareWrite(33);
    check(output->writeFastQ().ends_with("\n+\n::::::\n"), "Raw detached merge evidence contributes once when promoted");
    output = snapshot->getMerge(); output->releaseRawQualitiesIfPromoted(); output->prepareWrite(33);
    check(output->writeFastQ().ends_with("\n+\nIIIIII\n"),
        "A raw-history snapshot keeps the existing replacement behaviour and is independent of later promotion");
}

void testDerepRepresentativeReplacementPreservesQualitySums() {
    auto first = makeDNA("ACGT", "first_representative");
    first->setAllQual(10);
    auto replacement = makeDNA("ACGT", "better_representative");
    replacement->setAllQual(30);

    DNAunique unique(first, 0);
	unique.attachPair(make_shared<DNAunique>(makeDNA("TGCA", "stale_pair"), -1));
	unique.attachMerge(make_shared<DNAunique>(makeDNA("ACGTTGCA", "stale_merge"), -1));
    unique.incrementSampleCounter(0);
    unique.sumQualities(replacement);
    unique.takeOverDNA(replacement, nullptr, nullptr, true);

    const vector<DNAunique::qual_accum_t> sums = unique.transferPerBaseQualitySum();
    check(sums.size() == 4,
        "representative replacement must retain the derep quality accumulator length");
    check(sums.size() == 4 && sums[0] == 40 && sums[3] == 40,
        "representative replacement must retain all previously accumulated qualities");
	check(unique.getPair() == nullptr && unique.getMerge() == nullptr &&
		unique.getMergeLength() == -1,
		"representative replacement without pair/merge data must clear stale attachments");

	auto shortRead = makeDNA("ACGTGG", "short_quality_member");
	shortRead->cutSeqPseudo(4);
	shortRead->setAllQual(10);
	auto longRead = makeDNA("ACGTAC", "long_quality_member");
	longRead->setAllQual(40);
	DNAunique variableLength(shortRead, 0);
	variableLength.incrementSampleCounter(0);
	variableLength.sumQualities(longRead);
	variableLength.takeOverDNA(longRead, nullptr, nullptr, true);
	variableLength.releaseRawQualitiesIfPromoted();
	variableLength.prepareWrite(33);
	const string variableFastq = variableLength.writeFastQ(false);
	check(variableFastq.find("\n+\n::::II\n") != string::npos,
		"variable-length derep qualities must average each base over reads that cover it");
}

void testDerepSetBestAbsorbsQualityEvidence() {
	auto dominantRead = makeDNA("ACGT", "dominant_subgroup");
	dominantRead->setAllQual(10);
	auto dominant = make_shared<DNAunique>(dominantRead, 0);
	dominant->incrementSampleCounter(0);
	dominant->sumQualities(dominantRead);

	auto absorbedRead = makeDNA("ACGT", "absorbed_subgroup");
	absorbedRead->setAllQual(30);
	auto absorbed = make_shared<DNAunique>(absorbedRead, 0);

	DNAuniqSet alternatives;
	alternatives[0] = dominant;
	alternatives[1] = absorbed;
	auto selected = alternatives.best(true);
	check(selected && selected->totalSum() == 3,
		"setBest must retain occurrence counts from every absorbed subgroup");
	selected->releaseRawQualitiesIfPromoted();
	selected->prepareWrite(33);
	const string fastq = selected->writeFastQ(false);
	check(fastq.find("\n+\n2222\n") != string::npos,
		"setBest must average qualities over the same absorbed population as its counts");
}

void testDerepSnapshotsRetainAllSubgroups() {
    auto read = makeDNA("ACGT", "unmerged");
    read->setAllQual(10);
    DNAuniqSet groups;
    groups.addNewDNAuniq(read, nullptr, nullptr, -1, 0, false);
    groups[-1]->incrementSampleCounterBy(0, 9);
    auto merged = makeDNA("ACGTAC", "merged");
    merged->setAllQual(30);
    groups.addNewDNAuniq(read, nullptr, merged, 10, 0, false);
    auto first = groups.best(true);
    check(first->totalSum() == 11 && groups[-1]->totalSum() == 10,
        "an unmerged winner must include all subgroups without modifying the tracker");
    groups[10]->incrementSampleCounterBy(1, 20);
    auto second = groups.best(true);
    check(second->totalSum() == 31 && second->getMerge() && first->totalSum() == 11,
        "later snapshots must reevaluate existing subgroups and preserve older snapshots");
    groups.addNewDNAuniq(read, nullptr, merged, 20, 2, false);
    check(groups.best(true)->totalSum() == 32 && groups.best(true)->totalSum() == 32,
        "new subgroups must appear once in every subsequent snapshot");
    second->getMerge()->cutSeqPseudo(2);
    check(groups[10]->getMerge()->length() == 6,
        "merged output trimming must not mutate live subgroup attachments");
}

void testMergedDerepQualityEvidence() {
    auto first = makeDNA("ACGT", "first"); first->setAllQual(20); first->setPassed(true);
    auto firstMerge = makeDNA("ACGTAC", "first_merge"); firstMerge->setAllQual(20);
    DNAunique group(first, 0);
    group.attachMerge(make_shared<DNAunique>(firstMerge, 0));
    auto second = makeDNA("ACGT", "second"); second->setAllQual(40); second->setPassed(true);
    auto secondMerge = makeDNA("ACGTACGG", "second_merge"); secondMerge->setAllQual(40);
    group.matchedDNA(second, nullptr, secondMerge, 0, false);
    group.takeOverDNA(second, nullptr, secondMerge);
    auto consensus = group.getMerge();
    consensus->releaseRawQualitiesIfPromoted(); consensus->prepareWrite(33);
    check(consensus->writeFastQ(false).find("\n+\n??????II\n") != string::npos,
        "merged qualities must survive replacement and average only covered positions");
    group.takeOverDNA(second, nullptr, nullptr);
    check(!group.getMerge(), "a missing accepted merge must clear the selected attachment");
    group.takeOverDNA(second, nullptr, secondMerge);
    consensus = group.getMerge();
    consensus->releaseRawQualitiesIfPromoted(); consensus->prepareWrite(33);
    check(consensus->writeFastQ(false).find("\n+\n??????II\n") != string::npos,
        "accepted consensus evidence must survive a temporary unmerged representative");

    DNAunique unmergedOutput(first, 0);
    unmergedOutput.attachMerge(make_shared<DNAunique>(firstMerge, 0));
    unmergedOutput.matchedDNA(second, nullptr, secondMerge, 0, false, false);
    check(unmergedOutput.getMerge()->quality_sum_per_base_.empty(),
        "ordinary FASTQ output must not accumulate unused consensus-quality arrays");

    DNAuniqSet sets;
    sets.addNewDNAuniq(first, nullptr, firstMerge, 1, 0, false);
    sets.addNewDNAuniq(second, nullptr, secondMerge, 2, 1, false);
    auto selected = sets.best(true)->getMerge();
    selected->releaseRawQualitiesIfPromoted(); selected->prepareWrite(33);
    check(selected->writeFastQ(false).find("\n+\n??????II\n") != string::npos,
        "subgroup aggregation must preserve merged quality sums and per-base coverage");
}

void testDerepQualityEvidenceAcrossLengthChanges() {
    auto longRead = makeDNA("ACGTAC", "long_member"); longRead->setAllQual(10);
    auto shortRead = makeDNA("ACGT", "short_member"); shortRead->setAllQual(30);
    DNAunique group(longRead, 0);
    group.incrementSampleCounter(0); group.sumQualities(shortRead);
    group.takeOverDNA(shortRead, nullptr, nullptr);
    group.prepareDerepQualities(33);
    group.releaseRawQualitiesIfPromoted();
    group.incrementSampleCounter(0); group.sumQualities(longRead);
    group.takeOverDNA(longRead, nullptr, nullptr);
    group.prepareDerepQualities(33);
    const string path = "tests/reg_quality_evidence.fq";
    { ofstream output(path); group.writeDerepFastQ(output, false); }
    check(readFile(path).find("\n+\n2222++\n") != string::npos,
        "derep writer must retain tail evidence across shorter representatives and count actual coverage");
    std::remove(path.c_str());

    DNAunique moved(longRead, -1);
    auto source = make_shared<DNAunique>(longRead, 0);
    source->incrementSampleCounter(0); source->sumQualities(shortRead);
    moved.takeOver(source); moved.releaseRawQualitiesIfPromoted(); moved.prepareWrite(33);
    check(moved.writeFastQ(false).find("\n+\n5555++\n") != string::npos,
        "transferring a representative must transfer quality coverage with its sums");
}

void testEnabledMergeStatistics() {
    const string ref = "ACGTTGCAAGTCCTGATCGTAGGCTAACGTTCAGTACCGATGGCATTCGATCCGTAAGCTGACTACGGTTCGACTAGCTTACCGGATCGTAGCTAGGCTAACGATCGTTCAGGCATCGATGCTACGTAGCATGAC";
    ReadMerger merger(true), finder;
    const string mergePath = "tests/reg_enabled_merge_stats.tsv";
    const string qualPath = "tests/reg_enabled_quality_stats.tsv";
    auto mergePair = [&](size_t firstLength, size_t secondLength, ReadMerger& seedFinder) {
        auto first = makeDNA(ref.substr(0, firstLength), "statistics/1");
        auto second = makeDNA(seedFinder.reverseComplement(ref.substr(ref.size() - secondLength)), "statistics/2");
        return seedFinder.findSeedForMerge(first, second) && merger.merge(first, second) != nullptr;
    };
    check(mergePair(40, 120, finder), "asymmetric reverse-only statistics fixture must merge");
    merger.printMergeHisto(mergePath);
    merger.printQualHisto(qualPath);
    check(readFile(mergePath).find("FracRev") != string::npos && readFile(qualPath).find("AvgQ_r2") != string::npos,
        "enabled reverse-only statistics must produce both reports without invalid memory access");
    merger.restStats();

    constexpr int threads = 12, perThread = 75;
    std::atomic<int> failures{0};
    vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t] {
            ReadMerger localFinder;
            for (int i = 0; i < perThread; ++i) {
                if (!mergePair(t % 2 ? 120 : 40, t % 2 ? 40 : 120, localFinder)) { ++failures; }
            }
        });
    }
    for (auto& worker : workers) { worker.join(); }
    check(failures == 0, "all concurrent statistics merges must succeed");
    merger.printMergeHisto(mergePath);
    double observations = 0;
    istringstream rows(readFile(mergePath)); string row;
    getline(rows, row);
    while (getline(rows, row)) {
        vector<string> fields; istringstream columns(row); string field;
        while (getline(columns, field, '\t')) { fields.push_back(field); }
        for (size_t column : {size_t(2), size_t(3), size_t(5), size_t(6)}) {
            if (column < fields.size() && !fields[column].empty()) { observations += stod(fields[column]); }
        }
    }
    check(observations == threads * perThread * (160 - ref.size()),
        "shared-merger statistics must retain every overlap observation across threads");
    ReadMerger aggregate;
    aggregate.addRMstats(&merger); aggregate.addRMstats(&aggregate);
    aggregate.printMergeHisto(qualPath);
    check(readFile(mergePath) == readFile(qualPath), "merger aggregation must preserve both orientation histograms");
    merger.restStats(); merger.printMergeHisto(mergePath); merger.printQualHisto(qualPath);
    const string emptyMerge = readFile(mergePath), emptyQual = readFile(qualPath);
    check(count(emptyMerge.begin(), emptyMerge.end(), '\n') == 1 && count(emptyQual.begin(), emptyQual.end(), '\n') == 1,
        "reset statistics must emit empty report headers without stale observations");
    std::remove(mergePath.c_str()); std::remove(qualPath.c_str());
}

void testOutputFinalizationFailureSurfaces() {
#if defined(__linux__)
	ofbufstream failing("/dev/full", ios::out, false, 32);
	failing << string(4096, 'X');
	checkThrows([&] { failing.deactivate(); },
		"an output error delayed until flush/close must be reported to the caller");

#if defined(_gzipread)
	const string gzipFull = "/tmp/sdm_regression_full.gz";
	::unlink(gzipFull.c_str());
	check(::symlink("/dev/full", gzipFull.c_str()) == 0,
		"parallel-gzip failure fixture should create a /dev/full symlink");
	{
		ofbufstream failingGzip(gzipFull, ios::out, true, 64);
		failingGzip << string(4096, 'Y');
		checkThrows([&] { failingGzip.deactivate(); },
			"parallel gzip must report a delayed compressed flush/close failure");
	}
	::unlink(gzipFull.c_str());
#endif
#endif
}

void testRepresentativeScoringDirections() {
    auto uniform = makeDNA(string(100, 'A'), "uniform_Q30");
    uniform->setAllQual(30);
    auto uneven = makeDNA(string(100, 'A'), "alternating_Q20_Q40");
    vector<qual_score> qualities(100);
    for (size_t i = 0; i < qualities.size(); ++i) { qualities[i] = i % 2 ? 40 : 20; }
    uneven->setQual(qualities);
    for (bool seed : {false, true}) {
        float best = -1.f;
        check(!whoIsBetter(uneven, nullptr, nullptr, uniform, nullptr, nullptr, best, seed, nullptr, 2, 2),
            "equal mean quality must not conceal a worse expected-error score");
        best = -1.f;
        check(whoIsBetter(uniform, nullptr, nullptr, uneven, nullptr, nullptr, best, seed, nullptr, 2, 2),
            "lower expected error must win when average quality and length are equal");
    }

	auto lowError = makeDNA("ACGTACGT", "low_error_candidate");
	auto highError = makeDNA("ACGTACGT", "high_error_reference");
	lowError->setAllQual(40);
	highError->setAllQual(10);
	float everBest = -1.f;
	check(whoIsBetter(lowError, nullptr, nullptr, highError, nullptr, nullptr,
		everBest, true, nullptr, 2, 2),
		"seed scoring must prefer lower accumulated error when other factors agree");
	everBest = -1.f;
	check(!whoIsBetter(highError, nullptr, nullptr, lowError, nullptr, nullptr,
		everBest, true, nullptr, 2, 2),
		"seed scoring must not reward higher accumulated error");

	auto closeToMean = makeDNA(string(100, 'A'), "close_merge_length");
	auto farFromMean = makeDNA(string(100, 'A'), "far_merge_length");
	closeToMean->setMergeLength(100);
	farFromMean->setMergeLength(110);
	auto closeMerged = makeDNA(string(100, 'A'), "close_merged");
	auto farMerged = makeDNA(string(100, 'A'), "far_merged");
	DNAunique history(makeDNA(string(100, 'A'), "merge_history"), -1);
	for (int i = 0; i < 50; ++i) { history.recordWhoIsBetterCompare(100, 100); }
	everBest = -1.f;
	check(whoIsBetter(closeToMean, nullptr, closeMerged, farFromMean, nullptr,
		farMerged, everBest, false, &history, 2, 2),
		"representative scoring must prefer merge lengths closer to the historical mean");
	everBest = -1.f;
	check(!whoIsBetter(farFromMean, nullptr, farMerged, closeToMean, nullptr,
		closeMerged, everBest, false, &history, 2, 2),
		"representative scoring must not reward merge-length outliers");

	// Keep the seed-extension output merged once an eligible merge exists.
	auto mergedSeed = make_shared<DNAunique>(makeDNA(string(100, 'A'), "merged_seed"), -1);
	mergedSeed->attachMerge(make_shared<DNAunique>(
		makeDNA(string(150, 'A'), "merged_seed_product"), -1));
	mergedSeed->setTempFloat(99.5f);
	auto unmergedSeed = make_shared<DNAunique>(makeDNA(string(100, 'A'), "unmerged_seed"), -1);
	unmergedSeed->setTempFloat(99.5f);

	betterSeedStruct startsMerged(mergedSeed);
	check(!startsMerged.whoIsBetter2(unmergedSeed) &&
		startsMerged.tdn1->getMerge() != nullptr,
		"an eligible unmerged candidate must not replace a merged OTU seed");
	betterSeedStruct startsUnmerged(unmergedSeed);
	check(startsUnmerged.whoIsBetter2(mergedSeed) &&
		startsUnmerged.tdn1->getMerge() != nullptr,
		"an eligible merged candidate must replace an unmerged OTU seed");
}

void testSequenceBounds() {
    auto allN = makeDNA("NNNN");
    allN->stripLeadEndN();
    check(allN->length() == 0 && allN->getSequence().empty(),
        "stripping terminal Ns from an all-N read should produce an empty read safely");

    int position = 0;
    check(allN->matchSeq_tot("ACGT", 1, 0, position) == -1,
        "primer matching should reject a primer longer than the read");

    auto trimmed = makeDNA("ACGTACGT", "physical_after_logical_trim");
    check(trimmed->cutSeqPseudo(6), "test read should accept an initial logical trim");
    check(trimmed->cutSeq(4), "physical tail removal after a logical trim should succeed");
    check(trimmed->length() == 4 && trimmed->mem_length() == 4,
        "physical tail removal must not subtract an already removed logical tail twice");

	DNA pseudo("ACGTAAAA", "pseudo_reverse");
	pseudo.setQual({1, 2, 3, 4, 5, 6, 7, 8});
	check(pseudo.seal() && pseudo.cutSeqPseudo(4),
		"logical reverse-complement test read should initialize and trim");
	check(pseudo.matchSeqRev("AAAA", 0, 0) == -1,
		"reverse-primer matching must ignore a pseudo-trimmed backing suffix");
	vector<uint64_t> qualityByNt(6, 0), countsByNt(6, 0);
	pseudo.ntSpecQualScores(qualityByNt, countsByNt);
	check(countsByNt[0] == 1 && countsByNt[1] == 1 && countsByNt[2] == 1 &&
		countsByNt[3] == 1 && qualityByNt.size() == 6,
		"nucleotide statistics must use only the logical read and fixed category bins");
	pseudo.reverse_compliment(false);
	check(pseudo.getSequence() == "ACGT" && pseudo.length() == 4 &&
		pseudo.mem_length() == 4 && pseudo.getQual() == vector<qual_score>({4, 3, 2, 1}),
		"reverse complementation must materialize and reverse only the logical prefix");
}

void testReversePrimerMatching() {
    OptContainer options;
    options["-sample_sep"] = "___";
    options["-o_qual_offset"] = "33";
    TestFilters filters(&options);
    filters.PrimerR = {"ACGT"};
    filters.PrimerR_RC = {"ACGT"};
    filters.setPrimerErrs(1);
    filters.removePrimer(true);

    bool reverseComplement = true;
    auto noPrimer = makeDNA(string(300, 'G'), "no_primer");
    const string original = noPrimer->getSequence();
    check(!filters.cutPrimerRev(noPrimer, 0, reverseComplement, true),
        "reverse-primer search must not accept a partial one-base match at the read end");
    check(noPrimer->getSequence() == original,
        "failed reverse-primer matching must leave the sequence unchanged");

    reverseComplement = true;
    auto withPrimer = makeDNA(string(296, 'G') + "ACGT", "with_primer");
    check(filters.cutPrimerRev(withPrimer, 0, reverseComplement, true),
        "a complete reverse primer near the read end should be found");
    check(withPrimer->length() == 296,
        "a confirmed reverse primer should be removed at its full-match start");
}

void testPairInput(const string& stem, const string& r1Contents,
    const string& r2Contents, int batchSize, bool shouldThrow) {
    const string r1Path = stem + "_r1.fastq";
    const string r2Path = stem + "_r2.fastq";
    writeFile(r1Path, r1Contents);
    writeFile(r2Path, r2Contents);

    bool threw = false;
    try {
        InputStreamer input(false, 33, "0", "1", 1);
        int pairs = 2;
        check(input.setupFastq("", r1Path + "," + r2Path, pairs, "", false, false),
            "paired FASTQ setup should succeed");
        multi_tmp_lines batch(batchSize);
        input.getDNAlines(&batch, batchSize, false, true);
    }
    catch (const runtime_error&) {
        threw = true;
    }

    std::remove(r1Path.c_str());
    std::remove(r2Path.c_str());
    check(threw == shouldThrow, shouldThrow
        ? "mismatched paired input should be rejected during raw record ingestion"
        : "matching paired input should be accepted");
}

void testPairedReadSynchronization() {
    const string oneR1 = "@read42/1\nACGT\n+\nIIII\n";
    const string oneR2 = "@read42/2\nTGCA\n+\nIIII\n";
    testPairInput("tests/reg_pair_ok", oneR1, oneR2, 1, false);

    const string wrongR2 = "@different/2\nTGCA\n+\nIIII\n";
    testPairInput("tests/reg_pair_id", oneR1, wrongR2, 1, true);

    const string extraR2 = oneR2 + "@extra/2\nAAAA\n+\nIIII\n";
    testPairInput("tests/reg_pair_eof", oneR1, extraR2, 2, true);

    const string legacyR1 = "@legacy#1:0\nACGT\n+\nIIII\n";
    const string legacyR2 = "@legacy#2:0\nTGCA\n+\nIIII\n";
    testPairInput("tests/reg_pair_legacy", legacyR1, legacyR2, 1, false);

}

void testOnlyPairSelectionUsesSingleInputTopology() {
	const string fastqR1 = "tests/reg_only_pair_r1.fastq";
	const string fastqR2 = "tests/reg_only_pair_r2.fastq";
	writeFile(fastqR1, "@selected/1\nAAAA\n+\nIIII\n");
	writeFile(fastqR2, "@selected/2\nCCCC\n+\nJJJJ\n");
	{
		InputStreamer input(false, 33, "0", "1", 1);
		int pairs = 2;
		check(input.setupFastq("", fastqR1 + "," + fastqR2, pairs,
			"2", false, false), "-onlyPair 2 FASTQ setup should succeed");
		const auto reads = input.getDNApairs();
		check(pairs == 1 && input.pairNum() == 1 && reads[0] && !reads[1] &&
			reads[0]->getSequence() == "CCCC",
			"-onlyPair 2 must expose R2 as one single-end input without a stale mate");
	}

	const string fastaR1 = "tests/reg_only_pair_r1.fna";
	const string fastaR2 = "tests/reg_only_pair_r2.fna";
	writeFile(fastaR1, ">selected/1\nAAAA\n");
	writeFile(fastaR2, ">selected/2\nGGGG\n");
	{
		InputStreamer input(true, 33, "0", "1", 1);
		int pairs = 2;
		check(input.setupFastaQual("", fastaR1 + "," + fastaR2, "", pairs,
			"2", false), "-onlyPair 2 FASTA setup should succeed");
		const auto reads = input.getDNApairs();
		check(pairs == 1 && input.pairNum() == 1 && reads[0] && !reads[1] &&
			reads[0]->getSequence() == "GGGG",
			"-onlyPair 2 must open only the selected FASTA mate");
	}

	std::remove(fastqR1.c_str());
	std::remove(fastqR2.c_str());
	std::remove(fastaR1.c_str());
	std::remove(fastaR2.c_str());
}

void testDirectPairedDNAInput() {
    const string r1Path = "tests/reg_direct_pair_r1.fastq";
    const string r2Path = "tests/reg_direct_pair_r2.fastq";
    writeFile(r1Path, "@direct#1:0\nACGT\n+\nIJKL\n");
    writeFile(r2Path, "@direct#2:0\nTGCA\n+\nLKJI\n");

    InputStreamer input(false, 33, "0", "1", 1);
    int pairs = 2;
    check(input.setupFastq("", r1Path + "," + r2Path, pairs, "", false, false),
        "direct paired FASTQ setup should succeed");
    const vector<shared_ptr<DNA>> reads = input.getDNApairs();
    check(reads.size() == 3 && reads[0] && reads[1] && !reads[2],
        "getDNApairs should return one synchronized R1/R2 pair");
    check(reads[0] && reads[0]->getSequence() == "ACGT" &&
        reads[1] && reads[1]->getSequence() == "TGCA",
        "direct paired input should preserve both sequences");
    check(reads[0] && reads[0]->getQual() == vector<qual_score>({40, 41, 42, 43}) &&
        reads[1] && reads[1]->getQual() == vector<qual_score>({43, 42, 41, 40}),
        "direct paired input should preserve per-base quality order");
    const vector<shared_ptr<DNA>> eof = input.getDNApairs();
    check(!eof[0] && !eof[1] && !eof[2],
        "getDNApairs should return an empty synchronized set at EOF");

    std::remove(r1Path.c_str());
    std::remove(r2Path.c_str());
}

void testFastaQualitySynchronization() {
    const string mismatchFasta = "tests/reg_fasta_qual_suffix.fna";
    const string mismatchQual = "tests/reg_fasta_qual_suffix.qual";
    writeFile(mismatchFasta, ">same/1\nACGT\n");
    writeFile(mismatchQual, ">same/2\n10 20 30 40\n");
    {
        InputStreamer input(true, 33, "0", "1", 1);
        int pairs = 1;
        check(input.setupFastaQual("", mismatchFasta, mismatchQual, pairs, "", false),
            "FASTA/QUAL suffix-mismatch fixture should open");
        checkThrows([&] { (void)input.getDNApairs(); },
            "FASTA and QUAL headers with different mate suffixes must not match");
    }
    std::remove(mismatchFasta.c_str());
    std::remove(mismatchQual.c_str());

    const string extraFasta = "tests/reg_fasta_qual_extra.fna";
    const string extraQual = "tests/reg_fasta_qual_extra.qual";
    writeFile(extraFasta, ">same\nACGT\n");
    writeFile(extraQual,
        ">same\n10 20 30 40\n"
        ">orphan\n11 21 31 41\n");
    {
        InputStreamer input(true, 33, "0", "1", 1);
        int pairs = 1;
        check(input.setupFastaQual("", extraFasta, extraQual, pairs, "", false),
            "extra-QUAL fixture should open");
        const vector<shared_ptr<DNA>> first = input.getDNApairs();
        check(first[0] && first[0]->getId() == "same",
            "the synchronized FASTA/QUAL record before an orphan should remain readable");
        checkThrows([&] { (void)input.getDNApairs(); },
            "an extra trailing QUAL record must be rejected");
    }
    std::remove(extraFasta.c_str());
    std::remove(extraQual.c_str());
}

void testInputFormatMismatchIsRejected() {
    const string fastaPath = "tests/reg_format_as_fastq.fna";
    const string fastqPath = "tests/reg_format_as_fasta.fastq";
    writeFile(fastaPath, ">format\nACGT\n");
    writeFile(fastqPath, "@format\nACGT\n+\nIIII\n");

    {
        InputStreamer input(false, 33, "0", "1", 1);
        int pairs = 1;
        checkThrows([&] {
            (void)input.setupFastq("", fastaPath, pairs, "", false, false);
        }, "FASTA data supplied through the FASTQ path must be rejected explicitly");
    }
    {
        InputStreamer input(true, 33, "0", "1", 1);
        int pairs = 1;
        checkThrows([&] {
            (void)input.setupFastaQual("", fastqPath, "", pairs, "", false);
        }, "FASTQ data supplied through the FASTA path must be rejected explicitly");
    }
	{
		InputStreamer input(true, 33, "0", "1", 1);
		int pairs = 2;
		checkThrows([&] {
			(void)input.setupFastaQual("", fastaPath, "", pairs, "", true);
		}, "paired FASTA setup must reject a missing mate during preflight");
	}

    std::remove(fastaPath.c_str());
    std::remove(fastqPath.c_str());
}

void testFilterOptionValidationDoesNotTruncateOutput() {
	OptContainer invalidQuota = makeSeedOutputOptions("", "");
	invalidQuota["-DemultiBPperSR"] = "-1";
	checkThrows([&] { Filters filters(&invalidQuota); },
		"negative demultiplex base quotas must be rejected");

	const string derepPath = "tests/reg_existing_derep.fna";
	const string sentinel = "existing output must survive preflight\n";
	writeFile(derepPath, sentinel);
	OptContainer options = makeSeedOutputOptions("", "");
	options["-o_dereplicate"] = derepPath;
	Filters filters(&options);
	check(filters.setcmdArgsFiles(), "dereplication option setup should succeed");
	check(readFile(derepPath) == sentinel,
		"option setup must not truncate dereplication output before input preflight");
	std::remove(derepPath.c_str());

	OptContainer onlyPairOptions = makeSeedOutputOptions("", "");
	onlyPairOptions["-paired"] = "2";
	onlyPairOptions["-onlyPair"] = "2";
	Filters onlyPairFilter(&onlyPairOptions);
	check(onlyPairFilter.isPaired() == 1,
		"no-map mode must preserve the effective single-end topology of -onlyPair");

}

void testDirectTripleInputRoutingAndReuse() {
    const string r1Path = "tests/reg_direct_triple_r1.fastq";
    const string midPath = "tests/reg_direct_triple_mid.fastq";
    const string r2Path = "tests/reg_direct_triple_r2.fastq";
    const string singlePath = "tests/reg_direct_triple_single.fastq";
    writeFile(r1Path, "@triple\nACGT\n+\nIIII\n");
    writeFile(midPath, "@triple\nGGGG\n+\nJJJJ\n");
    writeFile(r2Path, "@triple\nTGCA\n+\nKKKK\n");
    writeFile(singlePath, "@single\nCCCC\n+\nLLLL\n");

    InputStreamer input(false, 33, "0", "1", 1);
    int pairs = 3;
    check(input.setupFastq("", r1Path + "," + midPath + "," + r2Path,
        pairs, "", false, false),
        "direct R1/MID/R2 FASTQ setup should succeed");
    const vector<shared_ptr<DNA>> triple = input.getDNApairs();
    check(triple.size() == 3 && triple[0] && triple[1] && triple[2],
        "paired=3 direct input should return R1, R2, and MID records");
    check(triple[0] && triple[0]->getSequence() == "ACGT" &&
        triple[1] && triple[1]->getSequence() == "TGCA" &&
        triple[2] && triple[2]->getSequence() == "GGGG",
        "paired=3 file order must route file 1 to R1, file 2 to MID, and file 3 to R2");
    check(triple[2] && triple[2]->isMIDseq(),
        "the directly read MID record should retain its MID marker");
	checkThrows([&] { input.openMIDseqs("", singlePath); },
		"an inline MID stream must not be overwritten by a separately configured MID file");

    input.allStreamReset();
    input.setMaxRdsRead(2);
    multi_tmp_lines limitedTriple(1);
    check(input.getDNAlines(&limitedTriple, 1, true, true) &&
        limitedTriple.size() == 1 &&
        limitedTriple.tmp[0][0][0] == "@triple" &&
        limitedTriple.tmp[0][1][0] == "@triple" &&
        limitedTriple.tmp[0][2][0] == "@triple",
        "a two-read limit must admit one R1/MID/R2 record because MID is metadata, not a counted mate");

    pairs = 1;
    check(input.setupFastq("", singlePath, pairs, "", false, false),
        "an InputStreamer previously used for paired=3 should be reusable for single input");
    check(input.pairNum() == 1 && !input.hasMIDseqs(),
        "single-input reuse must clear paired and MID topology state");
    const vector<shared_ptr<DNA>> single = input.getDNApairs();
    check(single.size() == 3 && single[0] && !single[1] && !single[2] &&
        single[0]->getSequence() == "CCCC",
        "single-input reuse must not read stale R2 or MID streams");

    std::remove(r1Path.c_str());
    std::remove(midPath.c_str());
    std::remove(r2Path.c_str());
    std::remove(singlePath.c_str());
}

void expectDirectFastqFailure(const string& stem, const string& contents,
    const string& message) {
    const string path = stem + ".fastq";
    writeFile(path, contents);
    InputStreamer input(false, 33, "0", "1", 1);
    int pairs = 1;
    check(input.setupFastq("", path, pairs, "", false, false),
        "malformed direct FASTQ fixture should open before record validation");
    checkThrows([&] { (void)input.getDNApairs(); }, message);
    std::remove(path.c_str());
}

void testMalformedDirectFastqInput() {
    expectDirectFastqFailure("tests/reg_direct_truncated",
        "@truncated\nACGT\n+\n",
        "a truncated direct FASTQ record must throw");
    expectDirectFastqFailure("tests/reg_direct_bad_plus",
        "@bad_plus\nACGT\n-\nIIII\n",
        "a direct FASTQ record with an invalid plus line must throw");
    expectDirectFastqFailure("tests/reg_direct_bad_length",
        "@bad_length\nACGT\n+\nIII\n",
        "a direct FASTQ record with mismatched sequence and quality lengths must throw");
}

void testMalformedBatchedFastqInput() {
    const string stem = "tests/reg_batched_malformed";
    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-merge_pairs_seed"] = "0";
    options["-suppressOutput"] = "3";
    Filters filters(&options);
    OutputStreamer output(&filters, &options, ios::out, nullptr, 1);

    multi_tmp_lines records(1);
    records.tmp[0][0] = {"@batched/1", "ACGT", "+", "IIII"};
    records.tmp[0][1] = {"@batched/2", "TGCA", "-", "IIII"};
    checkThrows([&] {
        (void)multi_read_paired_STRready(&records, false, &output, 0, true, 33);
    }, "a malformed required R2 in a converted input batch must throw");
}

void testThreadedPipelineDrainsAfterWorkerFailure() {
    const string r1Path = "tests/reg_worker_failure_r1.fastq";
    const string r2Path = "tests/reg_worker_failure_r2.fastq";
    // Cover errors consumed when reusing a slot and during the final drain,
    // multiple failed jobs, and a reader error while workers are outstanding.
    for (int threadCount : {1, 2, 4, 12}) {
        for (int scenario = 0; scenario < 4; ++scenario) {
            string r1Contents;
            string r2Contents;
            for (int index = 0; index < 65; ++index) {
                const string id = "worker_failure_" + to_string(index);
                const bool malformed = (scenario == 0 && index == 32) ||
                    (scenario == 1 && index == 64) ||
                    (scenario == 2 && (index == 32 || index == 64));
                const string mateId = (scenario == 3 && index == 64) ? "wrong_mate" : id;
                r1Contents += "@" + id + "/1\nACGT\n+\nIIII\n";
                r2Contents += "@" + mateId + "/2\nTGCA\n" +
                    string(malformed ? "-\n" : "+\n") + "IIII\n";
            }
            writeFile(r1Path, r1Contents);
            writeFile(r2Path, r2Contents);

            OptContainer options = makeSeedOutputOptions(
                "tests/reg_worker_failure.1.fq", "tests/reg_worker_failure.2.fq");
            options["-merge_pairs_seed"] = "0";
            options["-suppressOutput"] = "3";
            options["-iniBlockSize"] = "1";
            Filters filters(&options);
            auto input = make_shared<InputStreamer>(false, 33, "0", "1", threadCount);
            int pairs = 2;
            check(input->setupFastq("", r1Path + "," + r2Path, pairs, "", false, false),
                "threaded worker-failure input should open successfully");

            OptContainer* previousOptions = getGlobalOptContainer();
            setGlobalOptContainer(&options);
            bool threw = false;
            try {
                OutputStreamer output(&filters, &options, ios::out, nullptr, threadCount);
                (void)read_sequences(&output, input, threadCount);
            }
            catch (const runtime_error&) {
                threw = true;
            }
            setGlobalOptContainer(previousOptions);
            check(threw,
                "a pipeline failure must be rethrown after all active workers are drained (" +
                to_string(threadCount) + " workers, scenario " + to_string(scenario) + ")");
        }
    }
    std::remove(r1Path.c_str());
    std::remove(r2Path.c_str());
}

void testAllStreamResetResetsLocalReadLimit() {
    const string path = "tests/reg_stream_reset.fastq";
    writeFile(path,
        "@reset_0\nACGT\n+\nIIII\n"
        "@reset_1\nTGCA\n+\nJJJJ\n");

    InputStreamer input(false, 33, "0", "1", 1);
    int pairs = 1;
    check(input.setupFastq("", path, pairs, "", false, false),
        "read-limit reset FASTQ setup should succeed");
    input.setMaxRdsRead(1);

    multi_tmp_lines first(1);
    check(input.getDNAlines(&first, 1, false, true) && first.size() == 1 &&
        first.tmp[0][0][0] == "@reset_0",
        "the local maximum-read counter should initially permit one record");
    multi_tmp_lines limited(1);
    check(!input.getDNAlines(&limited, 1, false, true) && limited.size() == 0,
        "the local maximum-read counter should stop the stream at its limit");

    input.allStreamReset();
    multi_tmp_lines replay(1);
    check(input.getDNAlines(&replay, 1, false, true) && replay.size() == 1 &&
        replay.tmp[0][0][0] == "@reset_0",
        "allStreamReset must reset the local maximum-read counter as well as the stream cursor");

    std::remove(path.c_str());
}

void testConcurrentFastaHeaderOrder() {
    const string pairedR1Path = "tests/reg_parallel_fasta_r1.fna";
    const string pairedR2Path = "tests/reg_parallel_fasta_r2.fna";
    const string singlePath = "tests/reg_parallel_fasta_single.fna";
    const string singleQualPath = "tests/reg_parallel_fasta_single.qual";
    constexpr int recordCount = 64;

    string pairedR1;
    string pairedR2;
    string singleFasta;
    string singleQual;
    for (int i = 0; i < recordCount; ++i) {
        pairedR1 += ">pair_" + to_string(i) + "/1\nACGT\n";
        pairedR2 += ">pair_" + to_string(i) + "/2\nTGCA\n";
        singleFasta += ">single_" + to_string(i) + "\nGATT\n";
        singleQual += ">single_" + to_string(i) + "\n" +
            to_string(i % 20 + 1) + " 10 20 30\n";
    }
    writeFile(pairedR1Path, pairedR1);
    writeFile(pairedR2Path, pairedR2);
    writeFile(singlePath, singleFasta);
    writeFile(singleQualPath, singleQual);

    InputStreamer pairedInput(true, 33, "0", "1", 2);
    InputStreamer singleInput(true, 33, "0", "1", 2);
    int paired = 2;
    int single = 1;
    check(pairedInput.setupFastaQual("", pairedR1Path + "," + pairedR2Path,
        "", paired, "", false),
        "concurrent paired FASTA setup should succeed");
    check(singleInput.setupFastaQual("", singlePath, singleQualPath,
        single, "", false),
        "concurrent single FASTA/QUAL setup should succeed");

    vector<pair<string, string>> pairedHeaders;
    vector<string> singleHeaders;
    vector<vector<qual_score>> singleQualities;
    exception_ptr pairedError;
    exception_ptr singleError;
    thread pairedReader([&] {
        try {
            for (;;) {
                const vector<shared_ptr<DNA>> reads = pairedInput.getDNApairs();
                if (!reads[0]) { break; }
                pairedHeaders.emplace_back(reads[0]->getId(), reads[1]->getId());
            }
        }
        catch (...) {
            pairedError = current_exception();
        }
    });
    thread singleReader([&] {
        try {
            for (;;) {
                const vector<shared_ptr<DNA>> reads = singleInput.getDNApairs();
                if (!reads[0]) { break; }
                singleHeaders.push_back(reads[0]->getId());
                singleQualities.push_back(reads[0]->getQual());
            }
        }
        catch (...) {
            singleError = current_exception();
        }
    });
    pairedReader.join();
    singleReader.join();

    check(!pairedError && !singleError,
        "concurrent paired FASTA and single FASTA/QUAL reads should not throw");
    check(pairedHeaders.size() == recordCount && singleHeaders.size() == recordCount,
        "concurrent FASTA readers should retain every record");
    for (size_t i = 0; i < pairedHeaders.size(); ++i) {
        check(pairedHeaders[i].first == "pair_" + to_string(i) + "/1" &&
            pairedHeaders[i].second == "pair_" + to_string(i) + "/2",
            "paired FASTA headers must remain synchronized and in input order");
    }
    for (size_t i = 0; i < singleHeaders.size(); ++i) {
        check(singleHeaders[i] == "single_" + to_string(i),
            "single FASTA headers must remain in input order under concurrent reading");
        check(i < singleQualities.size() && singleQualities[i].size() == 4 &&
            static_cast<int>(singleQualities[i][0]) == static_cast<int>(i % 20 + 1) &&
            static_cast<int>(singleQualities[i][3]) == 30,
            "single FASTA qualities must remain paired with their matching header");
    }

    std::remove(pairedR1Path.c_str());
    std::remove(pairedR2Path.c_str());
    std::remove(singlePath.c_str());
    std::remove(singleQualPath.c_str());
}

void testMappingRows() {
    const string mapPath = "tests/reg_mapping.tsv";
    writeFile(mapPath,
		"# generated regression-test map\n"
        "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\n"
        "sample_a\tACGT\tTTAA\n"
        "sample_b\tTGCA\tCCGG\n");

    OptContainer options;
    options["-sample_sep"] = "___";
    options["-o_qual_offset"] = "33";
    options["-map"] = mapPath;
    options["-paired"] = "1";
    Filters filters(&options);
    check(filters.readMap(), "mapping file should parse successfully");
    check(filters.SampleID.size() == 2 && filters.Barcode.size() == 2,
        "mapping parser should create exactly one sample per data row");
    check(filters.SampleID[0] == "sample_a" && filters.SampleID[1] == "sample_b",
        "mapping parser should not treat the header as a sample");
    std::remove(mapPath.c_str());
}

void testHighDepthDereplicationCountsAndQualities() {
    const string mapPath = "tests/reg_large_derep.map";
    std::remove(mapPath.c_str());
    auto unique = make_shared<DNAunique>(makeDNA("ACGT", "large_count"), -1);
    const uint64_t largeCount =
        static_cast<uint64_t>(std::numeric_limits<int>::max()) + 17ULL;
    unique->incrementSampleCounterBy(0, largeCount);
    check(unique->totalSum() == largeCount,
        "dereplication totals must preserve counts above signed 32-bit range");
    const vector<uint64_t> sorted = unique->getDerepMapSort(1);
    check(sorted.size() == 1 && sorted[0] == largeCount,
        "dereplication threshold sorting must preserve large unsigned counts");

    vector<uint64_t> perSample(1, 0);
    vector<uint64_t> combined;
    vector<int> combinedIds;
    vector<int> touched;
    {
        ofstream map(mapPath, ios::out | ios::trunc);
        unique->writeMap(map, "large_count", perSample, combinedIds, combined, touched);
    }
    check(perSample[0] == largeCount &&
        readFile(mapPath).find("0:" + to_string(largeCount)) != string::npos,
        "dereplication maps must write exact counts above INT_MAX");

    DNAunique::qual_accum_t qualityTotal =
        static_cast<DNAunique::qual_accum_t>(std::numeric_limits<uint32_t>::max()) - 2;
    unique->addToQualSum(qualityTotal, 10);
    check(qualityTotal ==
        static_cast<DNAunique::qual_accum_t>(std::numeric_limits<uint32_t>::max()) + 8,
        "per-base quality accumulation must not saturate at 32 bits");

    auto capped = make_shared<DNAunique>(makeDNA("ACGT", "capped_count"), -1);
    capped->incrementSampleCounterBy(0, std::numeric_limits<uint64_t>::max());
    checkThrows([&] { capped->incrementSampleCounter(0); },
        "per-sample dereplication overflow must throw instead of wrapping");
    std::remove(mapPath.c_str());
}

void testDereplicationUsesLogicalTrimmedSequence() {
    const string outputPath = "tests/reg_trimmed_derep.fna";
    const string base = "tests/reg_trimmed_derep";
    const vector<string> paths = {
        outputPath, outputPath + ".rest", base + ".map",
        base + ".1.hq.fq", base + ".2.hq.fq", base + ".merg.fna"
    };
    for (const string& path : paths) { std::remove(path.c_str()); }

    OptContainer options = makeSeedOutputOptions("", "");
    options["-paired"] = "1";
    options["-o_dereplicate"] = outputPath;
    options["-derep_format"] = "fa";
    options["-derepSrchLen"] = "-1";
    Filters filters(&options);
    filters.SampleID = {"sample"};
    filters.collectStatistics[0]->BarcodeDetected.assign(1, 2);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);

    {
        Dereplicate dereplicator(&options, &filters);
        dereplicator.BCnamesAdding(&filters);
        auto first = makeDNA(string(40, 'A') + string(20, 'C'), "trimmed_first");
        auto second = makeDNA(string(40, 'A') + string(20, 'G'), "trimmed_second");
        first->cutSeqPseudo(40);
        second->cutSeqPseudo(40);
        first->setBCnumber(0, 0);
        second->setBCnumber(0, 0);
        first->setPassed(true);
        second->setPassed(true);
        check(dereplicator.addDNA(first, nullptr) && dereplicator.addDNA(second, nullptr),
            "trimmed dereplication fixtures should both pass insertion");
        (void)dereplicator.writeDereplDNA(&filters, "");
    }

    check(fastaHeaders(outputPath).size() == 1,
        "reads with the same logical prefix but different trimmed tails must dereplicate together");
    for (const string& path : paths) { std::remove(path.c_str()); }
}

void testStatistics() {
    ReportStats stats;
    vector<uint64_t> histogram(11, 0);
    histogram[10] = 1;
    check(stats.calc_median(histogram, 0.5f) == 10.0f,
        "median of a one-observation histogram should be that observation");

    ReportStats percentileStats;
    for (int value = 1; value <= 20; ++value) {
        auto read = make_shared<DNA>(string(static_cast<size_t>(value), 'A'), "percentile_read");
        read->setQual(vector<qual_score>(static_cast<size_t>(value), static_cast<qual_score>(value)));
        check(read->seal(), "percentile test DNA should seal successfully");
        percentileStats.addDNAStats(read);
    }
    ostringstream percentileReport;
    percentileStats.printStats2(percentileReport, 20.0f, 0);
    const string percentileText = percentileReport.str();
    check(percentileText.find("5%/50%/95% quantiles Pair 1") != string::npos,
        "read-characteristic report should identify its percentile triplet");
    check(percentileText.find("sequence Length : 1/10/19") != string::npos,
        "read-length report should contain the nearest-rank 5th, 50th, and 95th percentiles");
    check(percentileText.find("Quality :   1/10/19") != string::npos,
        "quality report should contain the nearest-rank 5th, 50th, and 95th percentiles");
    check(percentileText.find("Min/Avg/Max") == string::npos,
        "read-characteristic report should no longer label percentiles as extrema and a mean");

    mergeStats merge;
    merge.logDistri(5, 1, 8, true);
    merge.logDistri(5, 1, 8, false);
    merge.prepStats();
    check(merge.percFwd.size() > 5 && fabs(merge.percFwd[5] - 0.5f) < 0.0001f,
        "merge mismatch fraction should be mismatches divided by all observations");

    mergeStats reverseOnly;
    reverseOnly.logDistri(1, 5, 8, false);
    reverseOnly.prepStats();
    check(reverseOnly.percRev.size() > 5 &&
        fabs(reverseOnly.percRev[5] - 1.0f) < 0.0001f,
        "reverse-only merge statistics must be prepared without forward observations");

    auto destination = make_shared<GAstats>();
    auto source = make_shared<GAstats>();
    destination->setBCs({"sample"}, {"ACGT"}, {""});
    source->setBCs({"sample"}, {"ACGT"}, {""});
    destination->GALDISperBC[0] = {5, 6, 7};
    destination->GALAMPLperBC[0] = {0.5, 0.6, 0.7};
    source->GALDISperBC[0] = {2};
    source->GALAMPLperBC[0] = {0.2};
    destination->addStats(source);
    check(destination->GALDISperBC[0][0] == 7 && destination->GALDISperBC[0].size() == 3,
        "GoldenAxe aggregation should iterate over the source distribution safely");
}

void testMergeIdentityCutoff() {
    ReadMerger merger;
    const string reference = "ACGTACGTACGT";
    string oneMismatch = reference;
    oneMismatch[5] = oneMismatch[5] == 'A' ? 'C' : 'A';

    const double identity = merger.percentIdentity(
        reference.data(), oneMismatch.data(), static_cast<int>(reference.size()), 1);
    check(identity > 0.9,
        "one mismatch in a 12-base overlap should pass the greater-than-90-percent cutoff");

    string twoMismatches = oneMismatch;
    twoMismatches[8] = twoMismatches[8] == 'A' ? 'C' : 'A';
    check(merger.percentIdentity(reference.data(), twoMismatches.data(),
        static_cast<int>(reference.size()), 1) < 0.0,
        "a mismatch count above the inclusive limit should still be rejected early");
    check(merger.percentIdentity(string_view(reference), string_view(oneMismatch), 1) > 0.9,
        "both percent-identity overloads should use the same inclusive mismatch semantics");
}

qual_score legacyMergedQuality(qual_score q1, qual_score q2, bool same) {
    const double p1 = pow(10, -1 * static_cast<int>(q1) / 10.f);
    const double p2 = pow(10, -1 * static_cast<int>(q2) / 10.f);
    double pNew = 0.0;
    if (same) {
        pNew = (p1 * p2 / 3) / (1 - p1 - p2 + (4 * p1 * p2 / 3));
    }
    else {
        if (p1 > p2) {
            pNew = p2 * (1 - p1 / 3) / (p2 + p1 - 4 * p2 * p1 / 3);
        }
        else {
            pNew = p1 * (1 - p2 / 3) / (p1 + p2 - 4 * p1 * p2 / 3);
        }
    }
    qual_score qNew = static_cast<qual_score>(-10 * log10(pNew) + 0.5);
    if (qNew >= 40) { qNew = 40; }
    return qNew;
}

void testMergedQualityLookupMatchesLegacyCalculation() {
    ReadMerger merger;
    for (int q1 = 1; q1 <= 40; ++q1) {
        for (int q2 = 1; q2 <= 40; ++q2) {
            for (const bool same : {false, true}) {
                const qual_score left = static_cast<qual_score>(q1);
                const qual_score right = static_cast<qual_score>(q2);
                check(merger.getMergedQual(left, right, same) ==
                    legacyMergedQuality(left, right, same),
                    "merged-quality lookup must match the legacy probability calculation");
            }
        }
    }
}

void testMergedQualityBoundsForUnusualInputs() {
    ReadMerger merger;
    const vector<qual_score> qualities = {
        static_cast<qual_score>(-5),
        static_cast<qual_score>(0),
        static_cast<qual_score>(1),
        static_cast<qual_score>(40),
        static_cast<qual_score>(93),
        static_cast<qual_score>(109),
        static_cast<qual_score>(110),
        static_cast<qual_score>(120),
        static_cast<qual_score>(127)
    };
    for (const qual_score q1 : qualities) {
        for (const qual_score q2 : qualities) {
            for (const bool same : {false, true}) {
                const int merged = static_cast<int>(merger.getMergedQual(q1, q2, same));
                check(merged >= 0 && merged <= 40,
                    "merged qualities must be clamped to the supported 0..40 range");
            }
        }
    }
    check(merger.getMergedQual(static_cast<qual_score>(-5),
        static_cast<qual_score>(-5), true) == 0,
        "invalid negative merged-quality inputs must conservatively produce Q0");
}

void testFailedSeedSearchClearsMergeMetadata() {
    ReadMerger merger;
    auto r1 = makeDNA(string(40, 'A'), "seed_state/1");
    auto r2 = makeDNA(string(40, 'T'), "seed_state/2");
    check(merger.findSeedForMerge(r1, r2),
        "reverse-complementary reads should establish merge metadata");
    check(r1->merge_seed_pos_ >= 0 && r1->merge_offset_ >= 0 &&
        r2->merge_seed_pos_ >= 0 && r2->merge_offset_ >= 0 && r2->reversed_merge_,
        "a successful seed search should populate both reads' merge metadata");

    r2->setSequence(string(40, 'C'));
    check(!merger.findSeedForMerge(r1, r2),
        "unrelated reads should fail a second seed search");
    check(r1->merge_seed_pos_ == -1 && r1->merge_offset_ == -1 &&
        r2->merge_seed_pos_ == -1 && r2->merge_offset_ == -1 &&
        !r1->reversed_merge_ && !r2->reversed_merge_,
        "a failed seed search must clear all stale merge positions, offsets, and orientation");

    r1->merge_seed_pos_ = 4;
    r1->merge_offset_ = 2;
    r1->reversed_merge_ = true;
    check(!merger.findSeedForMerge(r1, nullptr) &&
        r1->merge_seed_pos_ == -1 && r1->merge_offset_ == -1 && !r1->reversed_merge_,
        "a missing mate must still clear merge metadata on the available read");
}

void testMergePreservesInputsAndLogicalTrim() {
    ReadMerger merger;
    auto r1 = makeDNA(string(40, 'A') + string(20, 'C'), "trimmed_merge/1");
    auto r2 = makeDNA(string(40, 'T') + string(20, 'G'), "trimmed_merge/2");
    check(r1->cutSeqPseudo(40) && r2->cutSeqPseudo(40),
        "merge trim fixture should apply logical 3-prime trimming");
    const string r1Sequence = r1->getSequence();
    const string r2Sequence = r2->getSequence();
    const vector<qual_score> r1Qualities = r1->getQual();
    const vector<qual_score> r2Qualities = r2->getQual();

    check(merger.findSeedForMerge(r1, r2),
        "logically trimmed reverse-complementary reads should find a merge seed");
    const shared_ptr<DNA> first = merger.merge(r1, r2);
    check(first && first->length() <= 40 &&
        first->getSequence().find_first_not_of('A') == string::npos,
        "merge construction must exclude logically trimmed backing-sequence tails");
    check(r1->getSequence() == r1Sequence && r2->getSequence() == r2Sequence &&
        r1->getQual() == r1Qualities && r2->getQual() == r2Qualities,
        "reverse-complement merging must not mutate either input read");

    check(merger.findSeedForMerge(r1, r2),
        "a repeated merge seed search should remain stable");
    const shared_ptr<DNA> second = merger.merge(r1, r2);
    check(first && second && first->getSequence() == second->getSequence() &&
        first->getQual() == second->getQual(),
        "repeated merges of the same inputs must produce identical output");
    check(r1->getSequence() == r1Sequence && r2->getSequence() == r2Sequence &&
        r1->getQual() == r1Qualities && r2->getQual() == r2Qualities,
        "repeated merges must leave source sequence and quality orientation unchanged");
}


void testRejectedMergeClearsCachedConsensus() {
    auto first = makeDNA(string(80, 'A'), "first_representative");
    auto representative = make_shared<DNAunique>(first, 0);
    representative->attachMerge(make_shared<DNAunique>(first, 0));
    auto replacement = makeDNA(string(80, 'C'), "replacement_without_merge");
    representative->takeOverDNA(replacement, nullptr, nullptr, true);
    check(!representative->getMerge(),
        "a replacement without an accepted consensus must discard the previous cached merge");
}

void testMergeOriginalOverlapAndBoundaries() {
    const string reference = "ACGTTGCAAGTCGATCGGATCTAGCTAACGGTTCAGTACCGATGCTAGGCTAACGTAGCTTGACCGTACAGGTTCGATACCTGA";
    ReadMerger merger;
    for (bool reversed : {false, true}) {
        auto r1 = makeDNA(reference.substr(0, 60), "overlap/1");
        string mate = reference.substr(24);
        if (reversed) { reverseTS(mate); }
        auto r2 = makeDNA(mate, "overlap/2");
        r1->setAllQual(20);
        r2->setAllQual(35);
        r1->cutSeqPseudo(40);
        r2->cutSeqPseudo(40);
        const auto q1 = r1->getQual(), q2 = r2->getQual();
        check(merger.findSeedForMerge(r1, r2), "original overlap must survive mate-output truncation");
        const auto merged = merger.merge(r1, r2);
        check(merged && merged->getSequence() == reference,
            "both forward and reverse overlaps must use consistent full-read coordinates");
        if (merged) {
            check(merged->getQual()[0] == 20 && merged->getQual().back() == 35 &&
                merged->getQual()[30] == merger.getMergedQual(20, 35, true),
                "overhang and overlap qualities must follow their sequence coordinates");
        }
        check(r1->getQual() == q1 && r2->getQual() == q2 &&
            r1->length() == 40 && r2->length() == 40 && r2->getSequence() == mate,
            "merging must preserve mate qualities, orientation and logical trims");
    }
    auto r1 = makeDNA(reference + "TTGACCTAGGTC", "readthrough/1");
    auto r2 = makeDNA(reverseTS2(reference) + "GTCAGCTTACGA", "readthrough/2");
    r1->setFwdPrimCut(); r1->setBarcodeCut(); r1->setBCnumber(3, 0);
    r2->setRevPrimCut(); r2->setBarcodeCut();
    check(merger.findSeedForMerge(r1, r2), "read-through fixture must overlap");
    auto merged = merger.merge(r1, r2);
    check(merged && merged->getSequence() == reference && merged->has2PrimersDetected() &&
        merged->getBarcodeCut() && merged->getBCnumber() == 3,
        "read-through must not restore removed primers/barcodes and cleanup metadata must survive");

    // R1 can read past the cleaned R2 5' boundary even without a left dovetail.
    r1 = makeDNA(reference, "asymmetric/1");
    r2 = makeDNA(reverseTS2(reference.substr(20, 50)), "asymmetric/2");
    check(merger.findSeedForMerge(r1, r2), "unequal-length reverse mates must overlap");
    merged = merger.merge(r1, r2);
    check(merged && merged->getSequence() == reference.substr(0, 70),
        "the cleaned reverse mate must bound the merged 3' end");
    r1->merge_offset_ = 999;
    check(!merger.merge(r1, r2), "out-of-range merge offsets must fail safely");
    check(!merger.merge(nullptr, r2), "missing merge inputs must fail safely");

    std::atomic<int> bad{0};
    vector<thread> workers;
    for (int worker = 0; worker < 4; ++worker) {
        workers.emplace_back([&] {
            ReadMerger seeds;
            auto a = make_shared<DNA>(reference.substr(0, 60), "parallel/1");
            auto b = make_shared<DNA>(reverseTS2(reference.substr(24)), "parallel/2");
            a->setQual(vector<qual_score>(60, 30)); b->setQual(vector<qual_score>(60, 30));
            if (!seeds.findSeedForMerge(a, b)) { ++bad; return; }
            for (int repeat = 0; repeat < 100; ++repeat) {
                const auto result = merger.merge(a, b);
                if (!result || result->getSequence() != reference) { ++bad; }
            }
        });
    }
    for (auto& worker : workers) { worker.join(); }
    check(bad == 0, "a shared dereplication merger must keep workers' reverse buffers separate");
}

void testMergedFilteringAndCachedLength() {
    const string stem = "tests/reg_merge_filter";
    const string optionFile = stem + ".options";
    writeFile(optionFile,
        "minSeqLength\t50\nmaxSeqLength\t200\nTruncateSequenceLength\t60\n"
        "minAvgQuality\t0\nQualWindowWidth\t10\nQualWindowThreshhold\t0\n"
        "TrimWindowWidth\t10\nTrimWindowThreshhold\t25\nmaxAccumulatedError\t-1\n"
        "BinErrorModelAlpha\t-1\nmaxAmbiguousNT\t0\nmaxHomonucleotide\t0\ntrimHomonucleotide\t0\n");
    auto options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-options"] = optionFile;
    Filters filters(&options);
    check(filters.readMap(), "merged quality test must initialize no-map filtering");
    auto good = makeDNA(string(90, 'A'), "merged_good");
    check(filters.checkMergedDNA(good) && good->length() == 90,
        "merged consensus must keep its full length despite the ordinary-read truncation cap");
    auto shortRead = makeDNA(string(40, 'A'), "merged_short");
    check(!filters.checkMergedDNA(shortRead), "merged consensus shorter than minSeqLength must fail");
    auto trimmed = makeDNA(string(90, 'A'), "merged_quality_trim");
    vector<qual_score> qualities(90, 0);
    fill(qualities.begin(), qualities.begin() + 20, 40);
    trimmed->setQual(qualities);
    check(!filters.checkMergedDNA(trimmed) && trimmed->length() < 50,
        "minSeqLength must also be checked after consensus-quality trimming");
    auto ambiguous = makeDNA(string(70, 'A') + "N" + string(19, 'A'), "merged_ambiguous");
    check(!filters.checkMergedDNA(ambiguous), "merged consensus must respect maxAmbiguousNT beyond the ordinary-read cap");
    check(filters.collectStatistics[0]->total == 0 && filters.collectStatistics[1]->total == 0,
        "consensus filtering must not count a pair as an additional input read");

    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        auto source = make_shared<DNAunique>(makeDNA(string(90, 'A'), "cached_source"), -1);
        auto cached = make_shared<DNAunique>(makeDNA(string(90, 'A'), "cached_merge"), -1);
        cached->cutSeqPseudo(20);
        source->attachMerge(cached);
        check(output.saveForWrite_merge(source, "cached_merge", 0, false),
            "an eligible cached merge must be written");
        output.closeOutStreams();
    }
    const string written = readFile(stem + ".merg.fq");
    check(written.find("\n" + string(90, 'A') + "\n+\n") != string::npos,
        "cached merged seed output must restore full length and ignore the clustering cap");
    std::remove(optionFile.c_str());
    // The output helper is declared below this test.
    for (const string& suffix : {".1.fq", ".2.fq", ".1.singl.fq", ".2.singl.fq", ".merg.fq"}) {
        std::remove((stem + suffix).c_str());
    }
}

OptContainer makeSeedOutputOptions(const string& r1Path, const string& r2Path) {
    OptContainer options;
    options["-sample_sep"] = "___";
    options["-paired"] = "2";
    options["-onlyPair"] = "";
    options["-i_qual_offset"] = "33";
    options["-o_qual_offset"] = "33";
    options["-5PR1cut"] = "0";
    options["-5PR2cut"] = "0";
    options["-logLvsQ"] = "";
    options["-GoldenAxe"] = "0";
    options["-o_fastq"] = r1Path + "," + r2Path;
    options["-o_fna"] = "";
    options["-o_qual"] = "";
    options["-o_fastq2"] = "";
    options["-o_fastq_noBC"] = "";
    options["-suppressOutput"] = "0";
    options["-threadIO"] = "0";
    options["-merge_pairs_seed"] = "1";
    options["-uparseVer"] = "";
    options["-count_chimeras"] = "F";
    return options;
}

OptContainer makeFastaOutputOptions(const string& sequencePaths,
    const string& qualityPath, int paired) {
    OptContainer options = makeSeedOutputOptions("", "");
    options["-paired"] = to_string(paired);
    options["-o_fastq"] = "";
    options["-o_fna"] = sequencePaths;
    options["-o_qual"] = qualityPath;
    options["-merge_pairs_seed"] = "0";
    return options;
}

void removeSeedOutputFiles(const string& stem) {
	for (int generation = 0; generation < 5; ++generation) {
		std::remove(applyFileIT(stem + ".1.fq", generation).c_str());
		std::remove(applyFileIT(stem + ".2.fq", generation).c_str());
		std::remove(applyFileIT(stem + ".1.fq", generation, SingletonFileDescr).c_str());
		std::remove(applyFileIT(stem + ".2.fq", generation, SingletonFileDescr).c_str());
		std::remove(applyFileIT(stem + ".merg.fq", generation).c_str());
	}
}

void testNoBarcodeOutputArityAndSingleEnd() {
    const string singlePath = "tests/reg_no_barcode_single.fq";
    const string pairedPath = "tests/reg_no_barcode_paired.fq";
    const string unusedPath = "tests/reg_no_barcode_unused.fq";
    std::remove(singlePath.c_str());
    std::remove(pairedPath.c_str());
    std::remove(unusedPath.c_str());

    {
        OptContainer options = makeSeedOutputOptions("", "");
        options["-paired"] = "1";
        options["-o_fastq"] = "";
        options["-o_fastq_noBC"] = singlePath;
        Filters filters(&options);
        {
            OutputStreamer output(&filters, &options, ios::out, nullptr, 1);
            output.writeNonBCReads(makeDNA("ACGT", "single_no_bc"), nullptr);
        }
        const string singleRecords = readFile(singlePath);
        check(countOccurrences(singleRecords, "@single_no_bc\n") == 1,
            "single-end no-barcode output must accept one path and flush on destruction");
		check(countOccurrences(singleRecords, "\n+\n????\n") == 1,
			"no-barcode output must apply the configured FASTQ quality offset");
    }

    {
        OptContainer options = makeSeedOutputOptions("", "");
        options["-o_fastq"] = "";
        options["-o_fastq_noBC"] = pairedPath;
        Filters filters(&options);
        {
            OutputStreamer output(&filters, &options, ios::out, nullptr, 1);
            output.writeNonBCReads(
                makeDNA("ACGT", "pair_no_bc/1"), makeDNA("TGCA", "pair_no_bc/2"));
        }
        const string pairedRecords = readFile(pairedPath);
        check(countOccurrences(pairedRecords, "@pair_no_bc/1\n") == 1 &&
            countOccurrences(pairedRecords, "@pair_no_bc/2\n") == 1,
            "paired no-barcode output with one path must retain both mates exactly once");
    }

	{
		OptContainer options = makeSeedOutputOptions("", "");
		options["-paired"] = "1";
		options["-o_fastq"] = "";
		options["-o_fastq_noBC"] = "tests";
		Filters filters(&options);
		checkThrows([&] {
			OutputStreamer output(&filters, &options, ios::out, nullptr, 1);
		}, "no-barcode output must reject a path that cannot be opened");
	}

    {
        OptContainer options = makeSeedOutputOptions("", "");
        options["-paired"] = "1";
        options["-o_fastq"] = "";
        options["-o_fastq_noBC"] = singlePath + "," + unusedPath;
        Filters filters(&options);
        checkThrows([&] {
            OutputStreamer output(&filters, &options, ios::out, nullptr, 1);
        }, "single-end no-barcode output must reject two output paths");
    }

    std::remove(singlePath.c_str());
    std::remove(pairedPath.c_str());
    std::remove(unusedPath.c_str());
}

void testFallbackSeedUsesOnlySingletonOutput() {
    const string stem = "tests/reg_seed_fallback";
    const string fallbackPath = stem + ".fna";
    const string assignmentsPath = stem + ".paf";
    removeSeedOutputFiles(stem);
    writeFile(fallbackPath, ">OTU_4208\nACGTACGT\n");
    writeFile(assignmentsPath, "");

    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-optimalRead2Cluster"] = assignmentsPath;
    Filters filters(&options);

    {
        UClinks links(&options);
        auto fallback = make_shared<InputStreamer>(true, 33, "0", "1", 1);
        fallback->setupFna(fallbackPath);
        links.setupDefSeeds(fallback, {});

        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        output.activateReadMerger(1);
        links.writeNewSeeds(&output, &filters, false);
    }

    const string merged = readFile(stem + ".merg.fq");
    const string singleton = readFile(stem + ".1.singl.fq");
    check(merged.empty(),
        "an unmergeable fallback seed must not be written to merged output");
    check(countOccurrences(singleton, "@OTU_4208.1\n") == 1,
        "an unmergeable fallback seed must be written exactly once as an R1 singleton");

    std::remove(fallbackPath.c_str());
    std::remove(assignmentsPath.c_str());
    removeSeedOutputFiles(stem);
}

void testMergedSeedKeepsFullLength() {
    const string stem = "tests/reg_seed_full_length";
    const string optionFile = stem + ".options";
    writeFile(optionFile,
        "minSeqLength\t50\nmaxSeqLength\t200\nTruncateSequenceLength\t150\n"
        "minAvgQuality\t45\nQualWindowWidth\t10\nQualWindowThreshhold\t45\n"
        "TrimWindowWidth\t10\nTrimWindowThreshhold\t45\nmaxAccumulatedError\t0.005\n"
        "BinErrorModelAlpha\t-1\nmaxAmbiguousNT\t0\nmaxHomonucleotide\t0\n");
    auto options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-options"] = optionFile;
    Filters filters(&options);
    check(filters.readMap(), "full merged seed test must initialize filtering");
    auto merged = makeDNA(string(180, 'A'), "full_seed");
    merged->setAllQual(40);
    merged->cutSeqPseudo(60);
    check(filters.checkMergedDNA(merged, true) && merged->length() == 180,
        "processed merged seeds must restore full length without clustering truncation or quality trimming");
    check(filters.checkMergedDNA(merged, true) && merged->length() == 180,
        "preparing a cached merged seed for output must preserve its full length");
    auto shortSeed = makeDNA(string(40, 'A'), "short_seed");
    check(!filters.checkMergedDNA(shortSeed, true),
        "restoring a merged seed must still respect minSeqLength");
    filters.setSeqLength(50, 140);
    check(!filters.checkMergedDNA(merged, true),
        "restoring a merged seed must still respect maxSeqLength");
    filters.setSeqLength(160, 200);
    check(filters.checkMergedDNA(merged, true) && merged->length() == 180,
        "a clustering cap below minSeqLength must not reject an eligible full-length seed");
    filters.setSeqLength(190, 200);
    check(!filters.checkMergedDNA(merged, true),
        "full-length merged seeds below the actual minimum must still be rejected");
    std::remove(optionFile.c_str());
}

void testSeedExtensionDoesNotRefilterProcessedPairs() {
	const string stem = "tests/reg_seed_extension_merge";
	const string fallbackPath = stem + ".fna";
	const string assignmentsPath = stem + ".paf";
	const string r1Path = stem + ".input.1.fq";
	const string r2Path = stem + ".input.2.fq";
	removeSeedOutputFiles(stem);
	writeFile(fallbackPath, ">OTU_4208\n" + string(80, 'A') + "\n");
	writeFile(assignmentsPath,
		"processed_pair/1\t80\t0\t80\t+\tOTU_4208\t80\t0\t80\t80\t80\t60\n");
	writeFile(r1Path, "@processed_pair/1\n" + string(80, 'A') +
		"\n+\n" + string(80, '?') + "\n");
	writeFile(r2Path, "@processed_pair/2\n" + string(80, 'T') +
		"\n+\n" + string(80, '?') + "\n");

	OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
	options["-optimalRead2Cluster"] = assignmentsPath;
	// These cuts have already been applied upstream to derep.*.hq.fq. If seed
	// extension applies them again, both reads become too short to merge.
	options["-5PR1cut"] = "65";
	options["-5PR2cut"] = "65";
	Filters filters(&options);

	{
		UClinks links(&options);
		auto fallback = make_shared<InputStreamer>(true, 33, "0", "1", 1);
		fallback->setupFna(fallbackPath);
		links.setupDefSeeds(fallback, {});
		links.activateMerger();

		auto input = make_shared<InputStreamer>(false, 33, "0", "1", 2);
		int paired = 2;
		check(input->setupFastq("", r1Path + "," + r2Path, paired, "", false, false),
			"processed paired seed input should initialize");
		links.findSeq2UCinstruction(input, true, &filters);

		OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
		output.activateReadMerger(1);
		links.writeNewSeeds(&output, &filters, false);
		ostringstream stats;
		links.printStats(stats);
		check(stats.str().find("1/1 seed candidates") != string::npos,
			"seed-extension diagnostics must report the actual overlap success rate");
	}

	check(countOccurrences(readFile(stem + ".merg.fq"), "@OTU_4208\n") == 1,
		"already-processed HQ pairs must merge without a second trim/orientation pass");
	check(readFile(stem + ".1.singl.fq").empty(),
		"a successfully merged seed must not also be routed to singleton output");

	std::remove(fallbackPath.c_str());
	std::remove(assignmentsPath.c_str());
	std::remove(r1Path.c_str());
	std::remove(r2Path.c_str());
	removeSeedOutputFiles(stem);
}

void testUnmergeablePairUsesOnlyPairedOutput() {
    const string stem = "tests/reg_seed_pair";
    removeSeedOutputFiles(stem);
    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    Filters filters(&options);

    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        output.activateReadMerger(1);
        auto r1 = make_shared<DNAunique>(makeDNA(string(80, 'A'), "r1"), -1);
        auto r2 = make_shared<DNAunique>(makeDNA(string(80, 'C'), "r2"), -1);
        r1->attachPair(r2);
        r1->setPassed(true);
        r2->setPassed(true);
        r1->setNewID("OTU_9.1");
        r2->setNewID("OTU_9.2");

        check(!output.saveForWrite_merge(r1, "OTU_9", 0),
            "an unmergeable read pair must not report merged output");
        r1->merge_seed_pos_ = 3;
        r1->merge_offset_ = 1;
        r2->merge_seed_pos_ = 4;
        r2->merge_offset_ = 0;
        r2->reversed_merge_ = true;
        check(!output.saveForWrite_merge(r1, "OTU_9", 99) &&
            r1->merge_seed_pos_ == -1 && r1->merge_offset_ == -1 &&
            r2->merge_seed_pos_ == -1 && r2->merge_offset_ == -1 &&
            !r2->reversed_merge_,
            "an invalid merge worker must fail safely and clear stale pair metadata");
        int r1Stream = 100;
        int r2Stream = 100;
        output.saveForWrite(r1, 1, -1, r1Stream, false);
        output.saveForWrite(r2, 2, -1, r2Stream, false);
        output.writeForWrite(r1, 1, r1Stream, r2, 2, r2Stream);
        output.closeOutStreams();
    }

    check(readFile(stem + ".merg.fq").empty(),
        "an unmergeable read pair must leave merged output empty");
    check(countOccurrences(readFile(stem + ".1.fq"), "@OTU_9.1\n") == 1,
        "an unmergeable read pair must retain exactly one R1 record");
    check(countOccurrences(readFile(stem + ".2.fq"), "@OTU_9.2\n") == 1,
        "an unmergeable read pair must retain exactly one R2 record");
    removeSeedOutputFiles(stem);
}

void testMergedSeedWriterReportsActualOutput() {
    const string stem = "tests/reg_seed_merged";
    removeSeedOutputFiles(stem);
    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    Filters filters(&options);

    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        auto representative = make_shared<DNAunique>(makeDNA("ACGT", "representative"), -1);
        auto merged = make_shared<DNAunique>(makeDNA("ACGTTGCA", "merged"), -1);
        representative->attachMerge(merged);
        check(output.saveForWrite_merge(representative, "OTU_7", 0),
            "a prepared merged seed should report that it was written");
        check(output.getReadsWritten() == 1,
            "a successful merged seed must reserve one global output slot");
        output.closeOutStreams();
    }

    const string mergedOutput = readFile(stem + ".merg.fq");
    check(countOccurrences(mergedOutput, "@OTU_7\n") == 1,
        "a successful merged seed should produce exactly one merged record");
    check(mergedOutput.find("@OTU_7.1\n") == string::npos,
        "merged seed headers must not retain the R1 suffix");
    removeSeedOutputFiles(stem);

    OptContainer suppressedOptions = makeSeedOutputOptions(
        stem + ".suppressed.1.fq", stem + ".suppressed.2.fq");
    suppressedOptions["-suppressOutput"] = "3";
    Filters suppressedFilters(&suppressedOptions);
    OutputStreamer suppressedOutput(&suppressedFilters, &suppressedOptions,
        ios::out, nullptr, 0);
    auto representative = make_shared<DNAunique>(makeDNA("ACGT", "representative"), -1);
    representative->attachMerge(
        make_shared<DNAunique>(makeDNA("ACGTTGCA", "merged"), -1));
    check(!suppressedOutput.saveForWrite_merge(representative, "OTU_8", 0),
        "a merged seed must not report success when no merged output stream exists");
}

void testMergedOutputRotationAndLimit() {
    const string stem = "tests/reg_seed_merged_limit";
    removeSeedOutputFiles(stem);
    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-maxReadsPerOutput"] = "1";
    Filters filters(&options);

    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        output.setReadLimit(2);
        for (int index = 0; index < 3; ++index) {
            auto representative =
                make_shared<DNAunique>(makeDNA("ACGT", "representative"), -1);
            representative->attachMerge(make_shared<DNAunique>(
                makeDNA("ACGTTGCA", "merged_" + to_string(index)), -1));
            const bool written = output.saveForWrite_merge(
                representative, "OTU_limit_" + to_string(index), 0);
            check(written == (index < 2),
                "merged output must honor its global read limit exactly");
        }
        check(output.getReadsWritten() == 2,
            "merged output accounting must stop at the configured global limit");
        output.closeOutStreams();
    }

    check(countOccurrences(readFile(applyFileIT(stem + ".merg.fq", 0)), "\n+\n") == 1 &&
        countOccurrences(readFile(applyFileIT(stem + ".merg.fq", 1)), "\n+\n") == 1,
        "merged records must rotate at maxReadsPerOutput without loss or duplication");
    check(!testFileExists(applyFileIT(stem + ".merg.fq", 2)),
        "hitting the merged global limit must not create an empty output generation");
    removeSeedOutputFiles(stem);
}

void testConcurrentPairedOutput() {
    const string stem = "tests/reg_parallel_output";
    removeSeedOutputFiles(stem);
    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-merge_pairs_seed"] = "0";
    Filters filters(&options);

    constexpr int threadCount = 4;
    constexpr int recordsPerThread = 40;
    vector<vector<pair<shared_ptr<DNA>, shared_ptr<DNA>>>> records(threadCount);
    for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
        records[threadIndex].reserve(recordsPerThread);
        for (int recordIndex = 0; recordIndex < recordsPerThread; ++recordIndex) {
            const string id = "parallel_" + to_string(threadIndex) + "_" +
                to_string(recordIndex);
            records[threadIndex].push_back({
                makeDNA("ACGT", id + "/1"), makeDNA("TGCA", id + "/2")});
        }
    }

    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, threadCount);
        vector<thread> writers;
        for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
            writers.emplace_back([&, threadIndex] {
                for (const auto& record : records[threadIndex]) {
                    output.writeForWrite(record.first, 1, 0, record.second, 2, 0);
                }
            });
        }
        for (thread& writer : writers) { writer.join(); }
        check(output.getReadsWritten() == threadCount * recordsPerThread,
            "concurrent paired output should reserve exactly one counter slot per pair");
        output.closeOutStreams();
    }

    const string r1 = readFile(stem + ".1.fq");
    const string r2 = readFile(stem + ".2.fq");
    check(countOccurrences(r1, "\n+\n") == threadCount * recordsPerThread,
        "concurrent paired output should retain every R1 record");
    check(countOccurrences(r2, "\n+\n") == threadCount * recordsPerThread,
        "concurrent paired output should retain every R2 record");
    removeSeedOutputFiles(stem);
}

void testConcurrentMergedAndPairedOutputLimits() {
    const string stem = "tests/reg_mixed_merged_limit";
    removeSeedOutputFiles(stem);
    auto options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-maxReadsPerOutput"] = "17";
    Filters filters(&options);
    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        output.setReadLimit(73);
        vector<std::thread> workers;
        for (int t = 0; t < 12; ++t) {
            workers.emplace_back([&, t] {
                for (int i = 0; i < 20; ++i) {
                    const string id = "mixed_" + to_string(t) + "_" + to_string(i);
                    auto first = make_shared<DNAunique>(makeDNA("ACGT", id + "/1"), -1);
                    if (t % 2) {
                        first->attachMerge(make_shared<DNAunique>(makeDNA("ACGTACGT", id), -1));
                        output.saveForWrite_merge(first, id);
                    } else {
                        auto second = makeDNA("TGCA", id + "/2");
                        first->prepareWrite(33); second->prepareWrite(33);
                        output.writeForWrite(first, 1, 0, second, 2, 0);
                    }
                }
            });
        }
        for (auto& worker : workers) { worker.join(); }
        output.closeOutStreams();
        check(output.getReadsWritten() == 73, "concurrent paired and merged output must share the global admission limit");
    }
    size_t total = 0;
    for (int generation = 0; generation < 5; ++generation) {
        const string firstPath = applyFileIT(stem + ".1.fq", generation);
        const string secondPath = applyFileIT(stem + ".2.fq", generation);
        const string mergedPath = applyFileIT(stem + ".merg.fq", generation);
        const auto first = countOccurrences(readFile(firstPath), "\n+\n"), second = countOccurrences(readFile(secondPath), "\n+\n");
        const auto merged = countOccurrences(readFile(mergedPath), "\n+\n");
        check(first == second, "mixed output rotation must preserve both mates");
        check(first + merged == (generation < 4 ? 17 : 5),
            "merged and paired records must share each generation's quota without losing records");
        total += first + merged;
    }
    check(total == 73, "mixed output files must contain all admitted records exactly once");
    removeSeedOutputFiles(stem);
}

void testConcurrentFastaQualityOutputOrder() {
    const string sequencePath = "tests/reg_parallel_output.fna";
    const string qualityPath = "tests/reg_parallel_output.qual";
    std::remove(sequencePath.c_str());
    std::remove(qualityPath.c_str());
    OptContainer options = makeFastaOutputOptions(sequencePath, qualityPath, 1);
    Filters filters(&options);

    constexpr int threadCount = 6;
    constexpr int recordsPerThread = 80;
    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, threadCount);
        vector<thread> writers;
        for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
            writers.emplace_back([&, threadIndex] {
                for (int recordIndex = 0; recordIndex < recordsPerThread; ++recordIndex) {
                    const string id = "fasta_qual_" + to_string(threadIndex) + "_" +
                        to_string(recordIndex);
                    output.writeForWrite(makeDNA(string(256, 'A'), id), 1, 0,
                        nullptr, -1, 100);
                }
            });
        }
        for (thread& writer : writers) { writer.join(); }
        output.closeOutStreams();
    }

    const vector<string> sequenceHeaders = fastaHeaders(sequencePath);
    const vector<string> qualityHeaders = fastaHeaders(qualityPath);
    check(sequenceHeaders.size() == threadCount * recordsPerThread,
        "concurrent FASTA/QUAL output should retain every record");
    check(sequenceHeaders == qualityHeaders,
        "concurrent FASTA and QUAL files must preserve identical record order");
    std::remove(sequencePath.c_str());
    std::remove(qualityPath.c_str());
}

void testConcurrentPairedFastaOutputOrder(bool writeQual = false, int recordsPerFile = 0) {
    const string r1Path = "tests/reg_parallel_output.1.fna";
    const string r2Path = "tests/reg_parallel_output.2.fna";
    const string q1Path = "tests/reg_parallel_output.1.qual";
    const string q2Path = "tests/reg_parallel_output.2.qual";
    OptContainer options = makeFastaOutputOptions(r1Path + "," + r2Path,
        writeQual ? q1Path + "," + q2Path : "", 2);
    options["-maxReadsPerOutput"] = to_string(recordsPerFile);
    Filters filters(&options);

    constexpr int threadCount = 12;
    constexpr int recordsPerThread = 80;
    const int expectedCount = recordsPerFile ? 499 : threadCount * recordsPerThread;
    const int generations = recordsPerFile ? (expectedCount + recordsPerFile - 1) / recordsPerFile : 1;
    auto cleanup = [&] {
        for (int generation = 0; generation <= generations; ++generation) {
            for (const string& path : {r1Path, r2Path, q1Path, q2Path}) {
                std::remove(applyFileIT(path, generation).c_str());
                std::remove(applyFileIT(path, generation, SingletonFileDescr).c_str());
                std::remove(applyFileIT(path + SingletonFileDescr, generation).c_str());
            }
        }
    };
    cleanup();
    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, threadCount);
        if (recordsPerFile) { output.setReadLimit(expectedCount); }
        vector<thread> writers;
        for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
            writers.emplace_back([&, threadIndex] {
                for (int recordIndex = 0; recordIndex < recordsPerThread; ++recordIndex) {
                    const string id = "paired_fasta_" + to_string(threadIndex) + "_" +
                        to_string(recordIndex);
                    output.writeForWrite(makeDNA(string(256, 'A'), id + "/1"), 1, 0,
                        makeDNA(string(256, 'T'), id + "/2"), 2, 0);
                }
            });
        }
        for (thread& writer : writers) { writer.join(); }
        check(output.getReadsWritten() == static_cast<uint64_t>(expectedCount),
            "concurrent FASTA output must account for exactly the admitted pairs");
        output.closeOutStreams();
    }

    size_t total = 0;
    vector<string> allHeaders;
    for (int generation = 0; generation < generations; ++generation) {
        vector<string> r1Headers = fastaHeaders(applyFileIT(r1Path, generation));
        vector<string> r2Headers = fastaHeaders(applyFileIT(r2Path, generation));
        total += r1Headers.size();
        allHeaders.insert(allHeaders.end(), r1Headers.begin(), r1Headers.end());
        check(!r1Headers.empty() && r2Headers.size() == r1Headers.size(),
            "each FASTA output generation must contain complete pairs");
        if (recordsPerFile) {
            check(r1Headers.size() <= static_cast<size_t>(recordsPerFile),
                "concurrent FASTA output must respect its per-file limit");
        }
        if (writeQual) {
            check(r1Headers == fastaHeaders(applyFileIT(q1Path, generation)) &&
                r2Headers == fastaHeaders(applyFileIT(q2Path, generation)),
                "both mates must preserve identical FASTA/QUAL order in each output generation");
        }
        for (string& header : r1Headers) {
            if (header.size() >= 2) { header.erase(header.size() - 2); }
        }
        for (string& header : r2Headers) {
            if (header.size() >= 2) { header.erase(header.size() - 2); }
        }
        check(r1Headers == r2Headers,
            "concurrent paired FASTA files must preserve identical pair order");
    }
    sort(allHeaders.begin(), allHeaders.end());
    check(total == static_cast<size_t>(expectedCount) &&
        adjacent_find(allHeaders.begin(), allHeaders.end()) == allHeaders.end(),
        "concurrent paired FASTA output must have no lost or duplicated pairs");
    check(!testFileExists(applyFileIT(r1Path, generations)),
        "the FASTA global limit must not create an empty extra generation");
    cleanup();
}

void testPairedFastaUsesOneQualityTier(bool fastq = false) {
    const string extension = fastq ? ".fq" : ".fna";
    const string primaryR1 = "tests/reg_fasta_tier.primary.1" + extension;
    const string primaryR2 = "tests/reg_fasta_tier.primary.2" + extension;
    const string secondaryR1 = "tests/reg_fasta_tier.secondary.1" + extension;
    const string secondaryR2 = "tests/reg_fasta_tier.secondary.2" + extension;
    const vector<string> paths = {primaryR1, primaryR2, secondaryR1, secondaryR2};
    for (int suppressPrimary : {0, 1}) {
        OptContainer options = fastq ? makeSeedOutputOptions(primaryR1, primaryR2) :
            makeFastaOutputOptions(primaryR1 + "," + primaryR2, "", 2);
        options[fastq ? "-o_fastq2" : "-o_fna2"] = secondaryR1 + "," + secondaryR2;
        options["-merge_pairs_seed"] = "0";
        options["-suppressOutput"] = to_string(suppressPrimary);
        Filters filters(&options);
        {
            OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
            for (int reversed = 0; reversed < 2; ++reversed) {
                const string id = "tiered_" + to_string(reversed);
                auto r1 = makeDNA("ACGT", id + "/1");
                auto r2 = makeDNA("TGCA", id + "/2");
                (reversed ? r2 : r1)->setPassed(true);
                (reversed ? r1 : r2)->setYellowQual(true);
                const int greenStream = suppressPrimary ? 100 : 0;
                output.writeForWrite(r1, 1, reversed ? 1 : greenStream,
                    r2, 2, reversed ? greenStream : 1);
            }
            // A failed mate must still be excluded, even with secondary output.
            output.writeForWrite(makeDNA("ACGT", "failed/1"), 1, 100,
                makeDNA("TGCA", "failed/2"), 2, 1);
            check(output.getReadsWritten() == 2,
                "mixed-tier passing mates must each reserve one complete pair, including with primary output suppressed");
            output.closeOutStreams();
        }

        check(readFile(primaryR1).empty() && readFile(primaryR2).empty(),
            "mates with different quality tiers must not be split into primary files");
        const string marker = fastq ? "@" : ">";
        for (int mate = 0; mate < 2; ++mate) {
            const string records = readFile(mate ? secondaryR2 : secondaryR1);
            check(countOccurrences(records, marker + "tiered_0/" + to_string(mate + 1) + "\n") == 1 &&
                countOccurrences(records, marker + "tiered_1/" + to_string(mate + 1) + "\n") == 1 &&
                records.find("failed/") == string::npos,
                "both passing mates must use their common less-stringent tier without admitting a failed mate");
        }
        for (const string& path : paths) {
            std::remove(path.c_str());
            std::remove((path + SingletonFileDescr).c_str());
            std::remove(applyFileIT(path, 0, SingletonFileDescr).c_str());
        }
    }
}

void testCombinedFileAndGlobalOutputLimits() {
    const string stem = "tests/reg_combined_output_limit";
    removeSeedOutputFiles(stem);
    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-merge_pairs_seed"] = "0";
    options["-maxReadsPerOutput"] = "3";
    Filters filters(&options);
    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        output.setReadLimit(3);
        for (int index = 0; index < 4; ++index) {
            output.writeForWrite(makeDNA("ACGT", "limited_" + to_string(index) + "/1"), 1, 0,
                makeDNA("TGCA", "limited_" + to_string(index) + "/2"), 2, 0);
        }
        check(output.getReadsWritten() == 3,
            "combined per-file and global limits must admit exactly the global allowance");
        output.closeOutStreams();
    }
    check(countOccurrences(readFile(stem + ".1.fq"), "\n+\n") == 3 &&
        countOccurrences(readFile(stem + ".2.fq"), "\n+\n") == 3,
        "combined output limits must retain three complete synchronized pairs");
    check(!testFileExists(applyFileIT(stem + ".1.fq", 1)) &&
        !testFileExists(applyFileIT(stem + ".2.fq", 1)),
        "a global-limit hit at the shard boundary must not create an empty generation");
    removeSeedOutputFiles(stem);
}

void testSubfilterMergePreservesOutputGeneration() {
    const string stem = "tests/reg_subfilter_rotation";
    for (int generation = 0; generation < 3; ++generation) {
        std::remove(applyFileIT(stem + ".1.fq", generation).c_str());
        std::remove(applyFileIT(stem + ".2.fq", generation).c_str());
        std::remove(applyFileIT(stem + ".1.fq", generation, SingletonFileDescr).c_str());
        std::remove(applyFileIT(stem + ".2.fq", generation, SingletonFileDescr).c_str());
        std::remove(applyFileIT(stem + ".merg.fq", generation).c_str());
    }

    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-merge_pairs_seed"] = "0";
    options["-maxReadsPerOutput"] = "3";
    Filters filters(&options);
    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, 2);
        for (int index = 0; index < 4; ++index) {
            const string id = "subfilter_rotation_" + to_string(index);
            output.writeForWrite(makeDNA("ACGT", id + "/1"), 1, 0,
                makeDNA("TGCA", id + "/2"), 2, 0);
        }
        check(filters.getFileIncrementor() == 1,
            "the fourth record should advance to output generation one");
        output.mergeSubFilters();
        check(filters.getFileIncrementor() == 1,
            "merging worker statistics must not restore a stale output generation");
        output.closeOutStreams();
    }

    for (int generation = 0; generation < 3; ++generation) {
        std::remove(applyFileIT(stem + ".1.fq", generation).c_str());
        std::remove(applyFileIT(stem + ".2.fq", generation).c_str());
        std::remove(applyFileIT(stem + ".1.fq", generation, SingletonFileDescr).c_str());
        std::remove(applyFileIT(stem + ".2.fq", generation, SingletonFileDescr).c_str());
        std::remove(applyFileIT(stem + ".merg.fq", generation).c_str());
    }
}

void testConcurrentOutputRotation() {
    const string stem = "tests/reg_parallel_rotation";
    constexpr int threadCount = 2;
    constexpr int recordsPerThread = 23;
    constexpr int recordsPerFile = 7;
    constexpr int expectedFiles =
        (threadCount * recordsPerThread + recordsPerFile - 1) / recordsPerFile;

    auto removeGeneration = [&](int generation) {
        std::remove(applyFileIT(stem + ".1.fq", generation).c_str());
        std::remove(applyFileIT(stem + ".2.fq", generation).c_str());
        std::remove(applyFileIT(stem + ".1.fq", generation, SingletonFileDescr).c_str());
        std::remove(applyFileIT(stem + ".2.fq", generation, SingletonFileDescr).c_str());
        std::remove(applyFileIT(stem + ".merg.fq", generation).c_str());
    };
    for (int generation = 0; generation < expectedFiles + 2; ++generation) {
        removeGeneration(generation);
    }

    OptContainer options = makeSeedOutputOptions(stem + ".1.fq", stem + ".2.fq");
    options["-merge_pairs_seed"] = "0";
    options["-maxReadsPerOutput"] = to_string(recordsPerFile);
    Filters filters(&options);
    {
        OutputStreamer output(&filters, &options, ios::out, nullptr, threadCount);
        vector<thread> writers;
        for (int threadIndex = 0; threadIndex < threadCount; ++threadIndex) {
            writers.emplace_back([&, threadIndex] {
                for (int recordIndex = 0; recordIndex < recordsPerThread; ++recordIndex) {
                    const string id = "rotated_" + to_string(threadIndex) + "_" +
                        to_string(recordIndex);
                    auto r1 = make_shared<DNA>("ACGT", id + "/1");
                    auto r2 = make_shared<DNA>("TGCA", id + "/2");
                    r1->setQual(vector<qual_score>(4, 30));
                    r2->setQual(vector<qual_score>(4, 30));
                    r1->seal();
                    r2->seal();
                    output.writeForWrite(r1, 1, 0, r2, 2, 0);
                }
            });
        }
        for (thread& writer : writers) { writer.join(); }
        check(output.getReadsWritten() == threadCount * recordsPerThread,
            "file rotation must not lose or double-count concurrent pairs");
        output.closeOutStreams();
    }

    size_t totalR1 = 0;
    size_t totalR2 = 0;
    for (int generation = 0; generation < expectedFiles; ++generation) {
        const size_t r1Count = countOccurrences(
            readFile(applyFileIT(stem + ".1.fq", generation)), "\n+\n");
        const size_t r2Count = countOccurrences(
            readFile(applyFileIT(stem + ".2.fq", generation)), "\n+\n");
        check(r1Count == r2Count,
            "every rotated output generation must keep R1/R2 synchronized");
        check(r1Count <= recordsPerFile,
            "a rotated output generation must respect maxReadsPerOutput");
        totalR1 += r1Count;
        totalR2 += r2Count;
    }
    check(totalR1 == threadCount * recordsPerThread && totalR2 == totalR1,
        "all concurrent pairs must survive output rotation exactly once");

    for (int generation = 0; generation < expectedFiles + 2; ++generation) {
        removeGeneration(generation);
    }
}

string derepTestSequence(int group) {
    string sequence(60, 'A');
    for (int i = 0; i < 8; ++i) {
        sequence[i] = "ACGT"[group % 4];
        group /= 4;
    }
    return sequence;
}

void removeDerepTestOutputs(const string& base) {
    for (const string& suffix : {".fq", ".fq.rest", ".map", ".1.hq.fq", ".2.hq.fq"}) {
        std::remove((base + suffix).c_str());
    }
}

OptContainer derepTestOptions(const string& base) {
    auto options = makeSeedOutputOptions("", "");
    options["-paired"] = "1";
    options["-o_dereplicate"] = base + ".fq";
    options["-derep_format"] = "fq";
    options["-derepSrchLen"] = "-1";
    options["-derep_reserve"] = "1"; // Exercise table growth in several shards.
    return options;
}

void testDerepPrefixStartupSelection() {
    const string base = "tests/reg_prefix_startup";
    const string optionFile = base + ".options";
    auto select = [&](const string& settings, bool expected, int paired = 1,
                      const string& search = "", const string& mode = "") {
        auto options = makeSeedOutputOptions("", "");
        options["-paired"] = to_string(paired);
        options["-o_dereplicate"] = base + ".fq";
        options["-derep_reserve"] = "1";
        if (!settings.empty()) {
            writeFile(optionFile, settings);
            options["-options"] = optionFile;
        }
        if (!search.empty()) { options["-derepSrchLen"] = search; }
        if (!mode.empty()) { options["-derepPrefix"] = mode; }
        Filters filters(&options);
        Dereplicate dereplicator(&options, &filters);
        check(dereplicator.usesPrefixDereplication() == expected,
            "startup prefix selection for " + settings + " search=" + search + " mode=" + mode);
    };
    select("", true);
    select("minSeqLength\t40\nmaxSeqLength\t100\n", true);
    select("minSeqLength\t40\nTruncateSequenceLength\t40\n", false);
    select("TruncateSequenceLength\t40\nminSeqLength\t60\n", false, 2);
    select("minSeqLength\t40\nTruncateSequenceLength\t40\n", false, 2, "-1");
    select("minSeqLength\t40\nTruncateSequenceLength\t60\n", true, 2);
    select("minSeqLength\t40\nmaxSeqLength\t40\n", false, 2);
    select("minSeqLength\t0.8\nTruncateSequenceLength\t40\n", true);
    select("minSeqLength\t40\n", false, 1, "20");
    select("minSeqLength\t40\nderepSrchLen\t20\n", false);
    select("minSeqLength\t40\nderepSrchLen\t20\n", true, 1, "-1");
    select("minSeqLength\t40\nderepPrefix\t0\n", false);
    select("minSeqLength\t40\nderepPrefix\t0\n", true, 1, "", "1");
    select("minSeqLength\t40\nTruncateSequenceLength\t40\n", true, 1, "", "1");
    select("", false, 1, "-1", "0");
    std::remove(optionFile.c_str());
    removeDerepTestOutputs(base);
}

void testDerepPrefixGroupsAndSnapshots() {
    const string prefix = derepTestSequence(93).substr(0, 40);
    const string a = prefix + "AA", b = prefix + "CC";
    struct Read { string sequence; int quality, sample; bool pass; };
    const vector<Read> reads = {
        {prefix, 20, 1, true}, {a, 40, 0, true}, {b, 40, 0, true},
        {a, 20, 1, true}, {prefix.substr(0, 35), 10, 1, false},
        {a + "TT", 20, 1, false}, {string(45, 'T'), 10, 1, false},
        {"T" + prefix, 30, 0, true}
    };
    for (int order = 0; order < 3; ++order) {
        const string base = "tests/reg_prefix_groups_" + to_string(order);
        removeDerepTestOutputs(base);
        auto options = derepTestOptions(base);
        Filters filters(&options);
        filters.SampleID = {"sample_a", "sample_b"};
        filters.collectStatistics[0]->BarcodeDetected.assign(2, -1);
        filters.collectStatistics[0]->BarcodeDetectedFail.assign(2, 0);
        Dereplicate dereplicator(&options, &filters);
        dereplicator.BCnamesAdding(&filters);
        auto add = [&](const Read& read, size_t id) {
            auto dna = makeDNA(read.sequence + "GG", "prefix_" + to_string(id));
            dna->setAllQual(read.quality);
            dna->cutSeqPseudo(static_cast<int>(read.sequence.size()));
            dna->setBCnumber(read.sample, 0);
            dna->setPassed(read.pass);
            dereplicator.addDNA(dna, nullptr);
        };
        // Failures alone never seed output, including an empty sequence.
        auto empty = makeDNA("A", "empty");
        empty->cutSeqPseudo(0); empty->setBCnumber(0, 0); empty->setPassed(false);
        dereplicator.addDNA(empty, nullptr);
        add(reads[6], 6);
        dereplicator.writeDereplDNA(&filters, "");
        check(readFile(base + ".fq").empty(), "failed-only prefixes must not seed output groups");
        dereplicator.reset();
        const int copies = order == 2 ? 300 : 1;
        if (order == 2) {
            vector<thread> workers;
            for (int t = 0; t < 4; ++t) {
                workers.emplace_back([&, t] {
                    for (int copy = 0; copy < copies; ++copy) {
                        for (size_t i = t; i < reads.size(); i += 4) { add(reads[i], i); }
                    }
                });
            }
            for (auto& worker : workers) { worker.join(); }
        } else {
            for (size_t i = 0; i < reads.size(); ++i) {
                const size_t index = order == 0 ? i : reads.size() - i - 1;
                add(reads[index], index);
            }
        }
        const auto report = dereplicator.writeDereplDNA(&filters, "");
        const string first = readFile(base + ".fq");
        ifstream output(base + ".fq");
        string header, sequence, plus, quality;
        int groups = 0;
        while (getline(output, header)) {
            getline(output, sequence); getline(output, plus); getline(output, quality);
            if (sequence == a) {
                check(header.find(";size=" + to_string(5 * copies) + ";") != string::npos,
                    "prefix aggregation must count every compatible read exactly once");
                check(quality == string(35, char(33 + 22)) + string(5, char(33 + 25)) + string(2, char(33 + 27)),
                    "prefix qualities must use position-specific coverage and all subgroup evidence");
            } else {
                check(sequence == b || sequence == "T" + prefix,
                    "longer conflicts stay separate and matches must start at position zero");
                check(header.find(";size=" + to_string(copies) + ";") != string::npos,
                    "an ambiguous short read must be assigned only once");
            }
            ++groups;
        }
        check(groups == 3, "prefix bridge must not join conflicting longer sequences");
        const string map = readFile(base + ".map");
        const size_t rowStart = map.find(";size=" + to_string(5 * copies) + ";");
        const string row = rowStart == string::npos ? "" : map.substr(rowStart, map.find('\n', rowStart) - rowStart);
        check(row.find("\t0:" + to_string(copies)) != string::npos &&
              row.find("\t1:" + to_string(4 * copies)) != string::npos,
            "prefix output must preserve sample counts across lengths and quality tiers");
        check(report.find("non-passing reads recovered") != string::npos,
            "snapshot must report deferred prefix recovery");
        dereplicator.writeDereplDNA(&filters, "");
        check(readFile(base + ".fq") == first, "repeated prefix snapshots must not mutate counts or qualities");
        // Extend an existing anchor after output; snapshots must rebuild membership.
        add({a + "GG", 30, 0, true}, 9);
        dereplicator.writeDereplDNA(&filters, "");
        const string next = readFile(base + ".fq");
        check(next.find(a + "GG\n") != string::npos && next.find("\n" + a + "\n") == string::npos,
            "a later longer passing prefix must replace the earlier output anchor");
        dereplicator.finishMap();
        removeDerepTestOutputs(base);
    }
}

void testConcurrentDereplicationCountsAndQualities() {
    const string base = "tests/reg_concurrent_derep";
    removeDerepTestOutputs(base);
    auto options = derepTestOptions(base);
    Filters filters(&options);
    constexpr int groups = 257, threads = 8, perThread = groups * 3 + 19;
    filters.SampleID = {"sample_a", "sample_b"};
    filters.collectStatistics[0]->BarcodeDetected.assign(2, perThread * threads / 2);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(2, 0);
    {
        Dereplicate dereplicator(&options, &filters);
        dereplicator.BCnamesAdding(&filters);
        vector<thread> workers;
        vector<exception_ptr> errors(threads);
        for (int t = 0; t < threads; ++t) {
            workers.emplace_back([&, t] {
                try {
                    for (int i = 0; i < perThread; ++i) {
                        const int group = i % groups;
                        auto dna = make_shared<DNA>(derepTestSequence(group), "group_" + to_string(group));
                        dna->setQual(vector<qual_score>(60, 20 + 2 * t));
                        if (!dna->seal()) { throw runtime_error("invalid derep fixture"); }
                        dna->setBCnumber(t % 2, 0);
                        dna->setPassed(true);
                        dereplicator.addDNA(dna, nullptr);
                    }
                } catch (...) { errors[t] = current_exception(); }
            });
        }
        for (auto& worker : workers) { worker.join(); }
        for (const auto& error : errors) { check(!error, "concurrent dereplication must not throw"); }
        // Producers have exited with nonempty tails smaller than the batch threshold.
        dereplicator.writeDereplDNA(&filters, "");
        ifstream output(base + ".fq");
        string header, sequence, plus, quality;
        int records = 0;
        vector<bool> seen(groups, false);
        while (getline(output, header)) {
            getline(output, sequence); getline(output, plus); getline(output, quality);
            int group = -1;
            for (int g = 0; g < groups; ++g) { if (sequence == derepTestSequence(g)) { group = g; break; } }
            check(group >= 0, "dereplication must preserve the full sequence key");
            if (group < 0) { continue; }
            check(!seen[group], "concurrent dereplication must emit each key exactly once");
            seen[group] = true;
            const int count = threads * (3 + (group < 19));
            check(header.find(";size=" + to_string(count) + ";") != string::npos,
                "concurrent dereplication must retain all counts, including exited-worker tails");
            check(plus == "+" && quality == string(60, static_cast<char>(33 + 27)),
                "concurrent dereplication must average every per-base quality contribution exactly");
            ++records;
        }
        check(records == groups, "concurrent dereplication must retain all sequence groups");
        const string map = readFile(base + ".map");
        for (int group = 0; group < groups; ++group) {
            const int sampleCount = (threads / 2) * (3 + (group < 19));
            const string prefix = "group_" + to_string(group) + ";size=" + to_string(2 * sampleCount) + ";";
            const size_t start = map.find(prefix + "\t");
            const string row = start == string::npos ? "" : map.substr(start, map.find('\n', start) - start);
            check(row.find("\t0:" + to_string(sampleCount)) != string::npos &&
                row.find("\t1:" + to_string(sampleCount)) != string::npos,
                "concurrent dereplication must preserve exact counts per sample");
        }
        dereplicator.reset();
        auto fresh = makeDNA(derepTestSequence(groups), "after_reset");
        fresh->setBCnumber(0, 0); fresh->setPassed(true);
        dereplicator.addDNA(fresh, nullptr);
        dereplicator.writeDereplDNA(&filters, "");
        check(readFile(base + ".fq").find("@after_reset;size=1;") == 0,
            "reset must clear every shard and permit new pending work");
        dereplicator.finishMap();
    }
    removeDerepTestOutputs(base);
}

void testDereplicationBatchOwnerReuse() {
    const string base = "tests/reg_derep_owner";
    removeDerepTestOutputs(base);
    auto options = derepTestOptions(base);
    Filters filters(&options);
    filters.SampleID = {"sample"};
    filters.collectStatistics[0]->BarcodeDetected.assign(1, 1);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
    alignas(Dereplicate) unsigned char storage[sizeof(Dereplicate)];
    for (int iteration = 0; iteration < 3; ++iteration) {
        auto* dereplicator = new (storage) Dereplicate(&options, &filters);
        dereplicator->BCnamesAdding(&filters);
        auto dna = makeDNA("ACGTACGT", "owner_" + to_string(iteration));
        dna->setBCnumber(0, 0); dna->setPassed(true);
        dereplicator->addDNA(dna, nullptr);
        // Leave the first owner's tail pending, then reuse its exact address.
        if (iteration > 0) {
            dereplicator->writeDereplDNA(&filters, "");
            check(readFile(base + ".fq").find("@owner_" + to_string(iteration) + ";size=1;") == 0,
                "thread caches must neither inherit old pending reads nor lose registration after owner address reuse");
        }
        dereplicator->~Dereplicate();
    }
    removeDerepTestOutputs(base);
}

void testDereplicationSnapshotDuringProduction() {
    const string base = "tests/reg_derep_snapshot";
    removeDerepTestOutputs(base);
    auto options = derepTestOptions(base);
    Filters filters(&options);
    constexpr int threads = 4, perThread = 4000, groups = 8;
    filters.SampleID = {"sample"};
    filters.collectStatistics[0]->BarcodeDetected.assign(1, threads * perThread);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
    Dereplicate dereplicator(&options, &filters);
    dereplicator.BCnamesAdding(&filters);
    atomic<int> ready{0};
    atomic<bool> go{false};
    vector<thread> workers;
    vector<exception_ptr> errors(threads);
    for (int t = 0; t < threads; ++t) {
        workers.emplace_back([&, t] {
            ++ready;
            while (!go.load()) { this_thread::yield(); }
            try {
                for (int i = 0; i < perThread; ++i) {
                    auto dna = make_shared<DNA>(derepTestSequence(i % groups), "snapshot_" + to_string(i % groups));
                    dna->setQual(vector<qual_score>(60, 30));
                    if (!dna->seal()) { throw runtime_error("invalid snapshot fixture"); }
                    dna->setBCnumber(0, 0); dna->setPassed(true);
                    dna->merge_seed_pos_ = (i / groups) % 3 - 1;
                    dereplicator.addDNA(dna, nullptr);
                }
            } catch (...) { errors[t] = current_exception(); }
        });
    }
    while (ready.load() != threads) { this_thread::yield(); }
    go.store(true);
    for (int i = 0; i < 20; ++i) {
        dereplicator.writeDereplDNA(&filters, "");
        this_thread::yield();
    }
    for (auto& worker : workers) { worker.join(); }
    for (const auto& error : errors) { check(!error, "snapshot must let every producer finish"); }
    dereplicator.writeDereplDNA(&filters, "");
    const string output = readFile(base + ".fq");
    check(countOccurrences(output, ";size=" + to_string(threads * perThread / groups) + ";") == groups,
        "snapshots during production must resume workers without dropping or duplicating updates");
    removeDerepTestOutputs(base);
}

void testDereplicationResumesAfterOutputFailure() {
    const string directory = "tests/reg_derep_missing_directory";
    const string base = directory + "/out";
    std::filesystem::remove_all(directory);
    auto options = derepTestOptions(base);
    Filters filters(&options);
    filters.SampleID = {"sample"};
    filters.collectStatistics[0]->BarcodeDetected.assign(1, 2);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
    Dereplicate dereplicator(&options, &filters);
    dereplicator.BCnamesAdding(&filters);
    for (int i = 0; i < 2; ++i) {
        auto dna = makeDNA("ACGTACGT", "resume");
        dna->setBCnumber(0, 0); dna->setPassed(true);
        dereplicator.addDNA(dna, nullptr);
        if (i == 0) {
            checkThrows([&] { dereplicator.writeDereplDNA(&filters, ""); },
                "dereplication must report an inaccessible output directory");
            std::filesystem::create_directory(directory);
        }
    }
    dereplicator.writeDereplDNA(&filters, "");
    check(readFile(base + ".fq").find("@resume;size=2;") == 0,
        "failed snapshots must resume producers and retain pending counts");
    removeDerepTestOutputs(base);
    std::filesystem::remove(directory);
}

} // namespace

int main(int argc, char** argv) {
    ini_DNAconstants();
    OptContainer globalOptions;
    setGlobalOptContainer(&globalOptions);

    try {
        if (argc == 2 && string(argv[1]) == "--derep") {
            testDerepPrefixStartupSelection();
            testDerepPrefixGroupsAndSnapshots();
            testConcurrentDereplicationCountsAndQualities();
            testDereplicationBatchOwnerReuse();
            testDerepSnapshotsRetainAllSubgroups();
            testMergedDerepQualityEvidence();
            testDereplicationSnapshotDuringProduction();
            testDereplicationResumesAfterOutputFailure();
        } else if (argc == 2 && string(argv[1]) == "--merge") {
            testMergeIdentityCutoff();
            testMergedQualityLookupMatchesLegacyCalculation();
            testMergedQualityBoundsForUnusualInputs();
            testFailedSeedSearchClearsMergeMetadata();
            testMergePreservesInputsAndLogicalTrim();
            testMergeOriginalOverlapAndBoundaries();
            testRejectedMergeClearsCachedConsensus();
            testMergedFilteringAndCachedLength();
            testSeedExtensionDoesNotRefilterProcessedPairs();
            testMergedSeedKeepsFullLength();
            testUnmergeablePairUsesOnlyPairedOutput();
        } else if (argc == 2 && string(argv[1]) == "--sync") {
            testPairedReadSynchronization();
            testDirectPairedDNAInput();
            testFastaQualitySynchronization();
            testMalformedBatchedFastqInput();
            testThreadedPipelineDrainsAfterWorkerFailure();
            testAllStreamResetResetsLocalReadLimit();
            testConcurrentFastaHeaderOrder();
            testNoBarcodeOutputArityAndSingleEnd();
            testConcurrentPairedOutput();
            testConcurrentFastaQualityOutputOrder();
            testConcurrentPairedFastaOutputOrder();
            testConcurrentPairedFastaOutputOrder(true);
            testConcurrentPairedFastaOutputOrder(true, 37);
            testPairedFastaUsesOneQualityTier();
            testPairedFastaUsesOneQualityTier(true);
            testCombinedFileAndGlobalOutputLimits();
            testConcurrentOutputRotation();
        } else {
        testLazyHeadersAndPhysicalTrimming();
        testUniformCoverageAndQualityHistory();
        testFastaQualityParsing();
        testMedianQualityDoesNotMutateRead();
        testQualityCachesFollowMutations();
        testQualityMathBoundsAndFractionalWindow();
        testDerepRepresentativeReplacementPreservesQualitySums();
        testDerepQualityEvidenceAcrossLengthChanges();
        testEnabledMergeStatistics();
		testDerepSetBestAbsorbsQualityEvidence();
        testDerepSnapshotsRetainAllSubgroups();
        testMergedDerepQualityEvidence();
		testOutputFinalizationFailureSurfaces();
		testRepresentativeScoringDirections();
        testSequenceBounds();
        testReversePrimerMatching();
        testPairedReadSynchronization();
		testOnlyPairSelectionUsesSingleInputTopology();
        testDirectPairedDNAInput();
        testFastaQualitySynchronization();
        testInputFormatMismatchIsRejected();
		testFilterOptionValidationDoesNotTruncateOutput();
        testDirectTripleInputRoutingAndReuse();
        testMalformedDirectFastqInput();
        testMalformedBatchedFastqInput();
        testThreadedPipelineDrainsAfterWorkerFailure();
        testAllStreamResetResetsLocalReadLimit();
        testConcurrentFastaHeaderOrder();
        testMappingRows();
        testHighDepthDereplicationCountsAndQualities();
        testDereplicationUsesLogicalTrimmedSequence();
        testDerepPrefixStartupSelection();
        testDerepPrefixGroupsAndSnapshots();
        testConcurrentDereplicationCountsAndQualities();
        testDereplicationBatchOwnerReuse();
        testDereplicationSnapshotDuringProduction();
        testDereplicationResumesAfterOutputFailure();
        testStatistics();
        testMergeIdentityCutoff();
        testMergedQualityLookupMatchesLegacyCalculation();
        testMergedQualityBoundsForUnusualInputs();
        testFailedSeedSearchClearsMergeMetadata();
        testMergePreservesInputsAndLogicalTrim();
        testMergeOriginalOverlapAndBoundaries();
        testRejectedMergeClearsCachedConsensus();
        testMergedFilteringAndCachedLength();
        testNoBarcodeOutputArityAndSingleEnd();
        testFallbackSeedUsesOnlySingletonOutput();
		testSeedExtensionDoesNotRefilterProcessedPairs();
            testMergedSeedKeepsFullLength();
        testUnmergeablePairUsesOnlyPairedOutput();
        testMergedSeedWriterReportsActualOutput();
        testMergedOutputRotationAndLimit();
        testConcurrentMergedAndPairedOutputLimits();
        testConcurrentPairedOutput();
        testConcurrentFastaQualityOutputOrder();
        testConcurrentPairedFastaOutputOrder();
        testConcurrentPairedFastaOutputOrder(true);
        testConcurrentPairedFastaOutputOrder(true, 37);
        testPairedFastaUsesOneQualityTier();
        testPairedFastaUsesOneQualityTier(true);
        testCombinedFileAndGlobalOutputLimits();
        testSubfilterMergePreservesOutputGeneration();
        testConcurrentOutputRotation();
        }
    }
    catch (const exception& error) {
        cerr << "Unexpected regression-test exception: " << error.what() << '\n';
        ++failures;
    }

    setGlobalOptContainer(nullptr);
    if (failures != 0) {
        cerr << failures << " regression test(s) failed\n";
        return 1;
    }
    cout << "All targeted regression tests passed\n";
    return 0;
}
