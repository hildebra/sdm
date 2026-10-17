// SPDX-License-Identifier: GPL-3.0-or-later
// Unit tests for the 2026-09 audit fixes and memory changes (sdm 3.55-3.60).
// Each test prints PASS/FAIL; the program exits non-zero if any check fails.
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "Common.h"
#include "DNA.h"
#include "DNAconsts.h"
#include "Filters.h"
#include "InputStream.h"
#include "OutputStreamer.h"
#include "PartitionWorkPool.h"
#include "ReadMerger.h"
#include "ThreadPool.h"
#include "Statistics.h"
#include "containers.h"

using namespace std;

namespace {

vector<string> failures;
string currentTest;

void check(bool condition, const string& message) {
    if (!condition) { failures.push_back(currentTest + ": " + message); }
}

void testHeaderStem() {
    // Pre-refactor semantics: split on any ':' or ' '.
    check(header_stem("M00123:45:000000000-ABCDE:1:1101:1000:2000") == "M00123:000000000-ABCDE:1:1101:1000",
        "7-token Illumina header");
    check(header_stem("inst:run:flow:lane:tile:x:y 1:N:0:ACGT") == "inst:lane:tile:x:y", "11-token Illumina header");
    check(header_stem("read7/1") == "read7", "slash suffix");
    bool threw = false;
    try { (void)header_stem("badheader"); } catch (const runtime_error&) { threw = true; }
    check(threw, "unparsable header throws");
}

void testSpaceX() {
    check(spaceX(-5).empty(), "negative width gives empty string (was a ~4 GB allocation)");
    check(spaceX(0).empty(), "zero width");
    check(spaceX(3) == "   ", "positive width");
}

void testIsDigits() {
    check(is_digits("0123"), "digits");
    check(!is_digits("12a"), "letter");
    check(!is_digits(string("1\xE9")), "high-bit byte is not a digit (no UB)");
}

void testGetFileNoPath() {
    string a = "dir/sub/b.fq", b = "C:\\x\\y.fq", c = "a/b\\c.fq", d = "plain.fq";
    check(getFileNoPath(a) == "b.fq", "POSIX path");
    check(getFileNoPath(b) == "y.fq", "Windows path");
    check(getFileNoPath(c) == "c.fq", "mixed separators");
    check(getFileNoPath(d) == "plain.fq", "no directory");
}

void testPrimerVariantTracker() {
    using Dir = PrimerVariantTracker::PrimerDirection;
    PrimerVariantTracker tracker;
    check(tracker.record("frag1", 0, Dir::Forward, "ACGT", "ACGT", false), "first forward record");
    check(!tracker.record("frag1", 0, Dir::Forward, "ACGT", "ACGT", false), "same fragment/direction counted once");
    check(tracker.record("frag1", 0, Dir::Reverse, "TTTT", "TTTA", true), "reverse direction counted separately");
    check(tracker.record("frag2", 0, Dir::Forward, "ACGT", "ACGA", true), "next fragment");
    for (int i = 0; i < 100000; ++i) { tracker.record("r" + to_string(i), 0, Dir::Forward, "ACGT", "ACGT", false); }
    check(tracker.totalForwardMatches == 100002 && tracker.totalReverseMatches == 1, "totals");
    check(tracker.forwardVariants.size() == 1 && tracker.reverseVariants.size() == 1, "variant tables");
    tracker.reset();
    check(tracker.totalForwardMatches == 0 && tracker.record("frag1", 0, Dir::Forward, "ACGT", "ACGT", false),
        "reset clears totals and the last-fragment memory");
}

void testDerepSubgroupMap() {
    DerepSubgroupMap m;
    auto a = make_shared<DNAunique>(), b = make_shared<DNAunique>(), c = make_shared<DNAunique>();
    m[5] = a; m[-1] = b; m[3] = c;
    vector<int> keys;
    for (const auto& entry : m) { keys.push_back(entry.first); }
    check(keys == vector<int>({-1, 3, 5}), "ascending key order, like std::map");
    check(m.find(3) != m.end() && m.find(3)->second == c, "find existing");
    check(m.find(4) == m.end(), "find missing");
    check(m[5] == a && m.size() == 3, "operator[] on existing key does not insert");
}

void testCompactOptionalString() {
    CompactOptionalString s;
    check(s.empty() && s.str().empty(), "default empty");
    s = string("ACGT-TGCA");
    CompactOptionalString copy = s;
    s = string();
    check(s.empty() && copy.str() == "ACGT-TGCA", "copies are independent");
    const string& asRef = copy;
    check(asRef == "ACGT-TGCA", "implicit const string& conversion");
}

void testMateHeaderRoundTrip() {
    check(DerepRead::mateHeaderFrom("x/1") == "x/2", "trailing /1 becomes /2");
    check(DerepRead::mateHeaderFrom("x") == "x", "CASAVA-style header (truncated at space) unchanged");
    check(DerepRead::mateHeaderFrom("x 1:N:0:ACGT") == "x 2:N:0:ACGT", "Illumina 1:N becomes 2:N");
    check(DerepRead::mateHeaderFrom("x/2") == "x/2", "no /1 suffix unchanged");

    DNAunique unique(make_shared<DNA>("ACGTACGT", "frag1/1 extra"), -1);
    unique.saveMem();
    auto mate = make_shared<DerepRead>(make_shared<DNA>("TTTTACGT", "frag1/2 extra"));
    auto merged = make_shared<DerepRead>(make_shared<DNA>("ACGTACGTTT", "frag1/1"));
    unique.attachPair(mate);
    unique.attachMerge(merged);
    unique.dropDerivableMateHeaders();
    check(mate->storedHeader().empty(), "derivable R2 header dropped");
    check(merged->storedHeader().empty(), "merged header equal to R1 dropped");
    mate->restoreMateHeader(unique.storedHeader());
    check(mate->storedHeader() == "frag1/2", "R2 header restored exactly");

    DNAunique other(make_shared<DNA>("ACGT", "fragA/1"), -1);
    other.saveMem();
    auto oddMate = make_shared<DerepRead>(make_shared<DNA>("TTTT", "fragB/2"));
    other.attachPair(oddMate);
    other.dropDerivableMateHeaders();
    check(oddMate->storedHeader() == "fragB/2", "non-derivable R2 header kept");
}

void testDerepReadIsSmaller() {
    check(sizeof(DerepRead) + 64 <= sizeof(DNAunique), "mate/merged objects omit DNAunique-only state");
}

// Exposes a protected member for the test only.
struct TestDNA : DNA {
    using DNA::DNA;
    using DNA::removeNonCanonicalNTs;
    using DNA::newId;
};

void testCutSeqWithShortQualities() {
    TestDNA d("ACGTACGTAC", "r");  // no qualities
    d.cutSeq(3);
    check(d.getSequence() == "ACG" && d.length() == 3, "3' cut without qualities");
    d.removeNonCanonicalNTs();
    TestDNA n("ACNNT", "n");
    n.removeNonCanonicalNTs();
    check(n.getSequence() == "ACAAT", "non-canonical bases replaced without qualities");
}

void testNegativeHomopolymerDisables() {
    DNA d("AAAAAAAAAAAACCCC", "h");
    int ambiguous = 0;
    check(d.scanSequenceChecks(-1, ambiguous), "negative maxHomonucleotide disables the check");
    check(!d.scanSequenceChecks(5, ambiguous), "positive threshold still rejects");
    check(d.HomoNTTrim(-1) == 0, "negative trimHomonucleotide disables trimming");
}

void testChangeHeadPverKeepsOriginal() {
    DNA d("ACGT", "r1 1:N:0:1");
    d.changeHeadPver(1);
    d.changeHeadPver(1);
    check(d.getOldId() == "r1 1:N:0:1", "original read ID not overwritten");
    check(d.getId() == "r1/1", "converted header");

    DNA labeled("ACGT", "r2 1:N:0:1");
    labeled.setNewID("S1___r2 orig_bc=AC");
    labeled.changeHeadPver(1);
    check(labeled.getId() == "S1___r2/1 orig_bc=AC", "sample label kept when converting to /1");
}

void testSetMateNumber() {
    DNA slash("ACGT", "frag/2");
    slash.setMateNumber(1);
    check(slash.getId() == "frag/1", "slash tag rewritten");
    DNA illumina("ACGT", "M1:9 2:N:0:ACGT+TT");
    illumina.setMateNumber(1);
    check(illumina.getId() == "M1:9 1:N:0:ACGT+TT", "Illumina tag rewritten, comment kept");
    DNA untagged("ACGT", "frag x");
    untagged.setMateNumber(2);
    check(untagged.getId() == "frag x", "header without mate tag unchanged");
    TestDNA labeled("ACGT", "frag/2");
    labeled.setNewID("S1___frag/2 orig_bc=AC");
    labeled.setMateNumber(1);
    check(labeled.newId() == "S1___frag/1 orig_bc=AC" && labeled.getOldId() == "frag/1", "relabeled header rewritten");
}

void testMatchSeqHighBitPrimer() {
    DNA d("ACGTACGTACGT", "r");
    const int hit = d.matchSeq(string("AC\xC3G"), 1, 6, 0);
    check(hit >= -1, "high-bit primer byte does not index DNA_IUPAC out of bounds");
    check(d.matchSeq("GTAC", 0, 8, 0) == 2, "normal primer match unchanged");
}

void testMergeStatsCountBeyondFloat() {
    mergeStats stats;
    const uint64_t n = (uint64_t(1) << 24) + 10;  // float ++ stops changing at 2^24
    for (uint64_t i = 0; i < n; ++i) { stats.logDistri(1, 0, 10, true); }
    stats.logDistri(1, 0, 10, false);
    check(stats.matchFwd.at(1) == n, "match counter keeps counting past 2^24");
    stats.prepStats();
    check(stats.percFwd.at(1) > 0.f && stats.percFwd.at(1) < 1e-6f, "mismatch fraction computed in double");
}

void testGoldenAxeMergeKeepsAllColumns() {
    auto a = make_shared<GAstats>(), b = make_shared<GAstats>();
    for (GAstats* g : {a.get(), b.get()}) {
        g->GAperBC.assign(2, 0); g->CNTperBC.assign(2, 0); g->GALENperBC.assign(2, 0); g->rdLENperBC.assign(2, 0);
        g->GAfailsPerBC.assign(2, 0); g->inCorPrimPerBC.assign(2, 0); g->missGAsPerBC.assign(2, 0);
        g->GALDISperBC.assign(2, {}); g->GALAMPLperBC.assign(2, {});
        g->SmplIDBC.assign(2, ""); g->BC1.assign(2, ""); g->BC2.assign(2, "");
    }
    b->inCorPrimPerBC = {3, 4};
    b->missGAsPerBC = {1, 2};
    a->addStats(b);
    check(a->inCorPrimPerBC == vector<uint64_t>({3, 4}), "IncorBCs column merged");
    check(a->missGAsPerBC == vector<uint64_t>({1, 2}), "MissedAmplics column merged");
}

void testDiffDNAClearKeepsPartialPairs() {
    auto ref = make_shared<const string>("ACGTACGT");
    auto mate = make_shared<const string>("TTTTACGT");
    diffDNA d(ref, false, false, mate, true);
    diffDNA shared = d;  // shared sequence data takes the copy-on-clear branch
    d.clear();
    bool ok = true;
    try { d.addPair("ACGTACGT", "", 1, 0); } catch (const exception&) { ok = false; }
    check(ok, "partial-pair mode survives clear() of shared data");
}

// --- I/O audit fixes (docs/io-audit.md) ---

string scratchPath(const string& name) {
    return (std::filesystem::temp_directory_path() / ("sdm_audit_" + name)).string();
}

void writeScratch(const string& path, const string& content) {
    ofstream out(path, ios::binary | ios::trunc);
    if (!out) { throw runtime_error("cannot create " + path); }
    out << content;
}

string fastqRecord(const string& id, const string& quals) {
    static const char* bases = "ACGT";
    string seq;
    for (size_t i = 0; i < quals.size(); ++i) { seq += bases[(i * 7 + id.size()) % 4]; }
    return "@" + id + "\n" + seq + "\n+\n" + quals + "\n";
}

// Offset detected from the probe block (the reader's path in read_sequences).
qual_score detectFromBlock(const string& path, int records) {
    InputStreamer input(false, 0, "0", "1", 1);
    int pairs = -1;
    if (!input.setupFastq("", path, pairs, "", false, false)) { throw runtime_error("cannot open " + path); }
    multi_tmp_lines probe(records);
    input.getDNAlines(&probe, records, false, true);
    input.finalizeFastqVersion();
    return input.fastQscore();
}

void testFastqOffsetAutoDetect() {
    // Finding 1: a first record with all bases >= Q26 and some at Q93 ('~')
    // used to select Phred+64 for the whole file.
    const string hifi = scratchPath("hifi.fq");
    writeScratch(hifi, fastqRecord("m1/ccs", string(10, '~') + string(10, 'I')) +
        fastqRecord("m2/ccs", string(20, '5')));
    check(detectFromBlock(hifi, 4) == 33, "HiFi-like qualities ('~' present) must select Phred+33");
    {
        // Immediate conversion path: the first record alone is conclusive here.
        InputStreamer input(false, 0, "0", "1", 1);
        int pairs = -1;
        input.setupFastq("", hifi, pairs, "", false, false);
        shared_ptr<DNA> read = input.getDNA(0);
        check(read != nullptr && read->getQual().size() == 20 && read->getQual()[0] == 93,
            "getDNA converts a '~' base as Q93 under Phred+33");
    }

    // Legacy Illumina 1.5 ('B'..'h') stays Phred+64 once the probe block is exhausted.
    const string legacy = scratchPath("legacy64.fq");
    writeScratch(legacy, fastqRecord("l1", string(8, 'h') + string(8, 'B')) +
        fastqRecord("l2", string(16, 'h')) + fastqRecord("l3", string(16, 'e')));
    check(detectFromBlock(legacy, 4) == 64, "Illumina 1.5 range (66..104) still selects Phred+64");

    // Undecided first record (94 '^' fits both), a later '!' settles on 33; the
    // old single-record rule returned 64 here.
    const string later = scratchPath("later33.fq");
    writeScratch(later, fastqRecord("u1", string(16, '^')) + fastqRecord("u2", string(15, '^') + "!"));
    check(detectFromBlock(later, 4) == 33, "evidence from later probe records selects Phred+33");

    // Modern Illumina header stays exempt.
    const string illumina = scratchPath("illumina.fq");
    writeScratch(illumina, fastqRecord("M0:1:F:1:1:1:1 1:N:0:1", string(16, '^')));
    check(detectFromBlock(illumina, 4) == 33, "' 1:N:' header selects Phred+33");

    for (const string& path : {hifi, legacy, later, illumina}) { std::remove(path.c_str()); }
}

// --- Dereplication audit fixes (docs/derep-audit.md) ---

string lcgSequence(size_t length, uint32_t seed) {
    static const char* bases = "ACGT";
    string s(length, 'A');
    uint32_t state = 0x9e3779b9u ^ seed;
    for (char& c : s) { state = state * 1664525u + 1013904223u; c = bases[(state >> 24) & 3]; }
    return s;
}

shared_ptr<DNA> derepRead(const string& sequence, const string& name, int quality, bool pass) {
    auto dna = make_shared<DNA>(sequence, name);
    dna->setQual(vector<qual_score>(sequence.size(), static_cast<qual_score>(quality)));
    check(dna->seal(), "seal test read");
    dna->setBCnumber(0, 0);
    dna->setPassed(pass);
    return dna;
}

string readScratch(const string& path) {
    ifstream in(path, ios::binary);
    return in ? string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>()) : string();
}

void testExactDerepRecoveryIsOrderIndependent() {
    // Finding 1: in exact-key mode a non-passing read was counted only if its
    // key already existed. Deferred reads are now attached at output when a
    // passing group exists (any order) and dropped otherwise (never in .rest).
    const string base = scratchPath("exact_order");
    const string shared = lcgSequence(80, 3), lonely = lcgSequence(80, 9);
    for (bool badFirst : {false, true}) {
        OptContainer opts;
        opts["-sample_sep"] = "___"; opts["-paired"] = "1";
        opts["-i_qual_offset"] = "33"; opts["-o_qual_offset"] = "33";
        opts["-5PR1cut"] = "0"; opts["-5PR2cut"] = "0";
        opts["-GoldenAxe"] = "0"; opts["-suppressOutput"] = "0"; opts["-threadIO"] = "0";
        opts["-merge_pairs_seed"] = "0"; opts["-count_chimeras"] = "F";
        opts["-o_dereplicate"] = base + ".fa"; opts["-derep_format"] = "fa";
        opts["-derepSrchLen"] = "60"; opts["-derep_reserve"] = "1";
        OptContainer* previousGlobal = getGlobalOptContainer();
        setGlobalOptContainer(&opts); // the record writers read runtime options from it
        Filters filters(&opts); filters.SampleID = {"S"};
        filters.collectStatistics[0]->BarcodeDetected.assign(1, -1);
        filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
        Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
        check(!derep.usesPrefixDereplication() && !derep.usesCoarseDereplication(), "fixture selects exact-key mode");
        auto passing = derepRead(shared, "good", 40, true), failed = derepRead(shared, "bad", 10, false);
        if (badFirst) { derep.addDNA(failed, nullptr); derep.addDNA(passing, nullptr); }
        else { derep.addDNA(passing, nullptr); derep.addDNA(failed, nullptr); }
        derep.addDNA(derepRead(lonely, "lonely", 10, false), nullptr);
        const string report = derep.writeDereplDNA(&filters, "");
        const string fa = readScratch(base + ".fa"), rest = readScratch(base + ".fa.rest"), map = readScratch(base + ".map");
        const string order = badFirst ? " (failed read first)" : " (passing read first)";
        check(fa.find(">good;size=2;") != string::npos, "the deferred failed read is recovered into its passing group" + order);
        check(fa.find(lonely) == string::npos && rest.find(lonely) == string::npos && map.find("lonely") == string::npos,
            "a deferred read without a passing group is dropped, not written to .rest or the map" + order);
        check(map.find("good;size=2;\t0:2") != string::npos, "the map carries the recovered sample count" + order);
        // A failed copy arriving after its passing key is counted directly; only
        // a copy arriving first is deferred and recovered at output.
        check(report.find(string("Exact-key matching: ") + (badFirst ? "1" : "0") + " non-passing reads recovered") != string::npos,
            "the report counts the recovered read" + order);
        derep.finishMap();
        setGlobalOptContainer(previousGlobal);
    }
    for (const char* suffix : {".fa", ".fa.rest", ".map", ".1.hq.fq", ".2.hq.fq"}) { std::remove((base + suffix).c_str()); }
}

void testMidQualityCountedOrPublished() {
    // A mid-quality record set is either counted by the dereplicator (key has a
    // passing dereplicate) or written to the mid-quality output, never both.
    const string base = scratchPath("mid_tier");
    const string optionsPath = base + ".options.txt";
    writeScratch(optionsPath, "minSeqLength\t10\nmaxSeqLength\t1000\nminAvgQuality\t30\n*minAvgQuality\t5\n"
        "maxAmbiguousNT\t-1\nmaxHomonucleotide\t0\ntrimHomonucleotide\t0\nTrimWindowThreshhold\t0\n"
        "QualWindowThreshhold\t0\nfastqVersion\t1\n");
    const string shared = lcgSequence(80, 21), alone = lcgSequence(80, 22);
    OptContainer opts;
    opts["-options"] = optionsPath;
    opts["-sample_sep"] = "___"; opts["-paired"] = "1";
    opts["-i_qual_offset"] = "33"; opts["-o_qual_offset"] = "33";
    opts["-5PR1cut"] = "0"; opts["-5PR2cut"] = "0";
    opts["-GoldenAxe"] = "0"; opts["-suppressOutput"] = "0"; opts["-threadIO"] = "0";
    opts["-merge_pairs_seed"] = "0"; opts["-count_chimeras"] = "F";
    opts["-o_fna"] = base + ".main.fna"; opts["-o_fna2"] = base + ".add.fna";
    opts["-o_qual"] = ""; opts["-o_fastq"] = ""; opts["-o_fastq2"] = ""; opts["-o_fastq_noBC"] = "";
    opts["-o_dereplicate"] = base + ".fa"; opts["-derep_format"] = "fa";
    opts["-derepSrchLen"] = "60"; opts["-derep_reserve"] = "1";
    OptContainer* previousGlobal = getGlobalOptContainer();
    setGlobalOptContainer(&opts);
    {
        Filters filters(&opts); filters.SampleID = {"S"};
        filters.collectStatistics[0]->BarcodeDetected.assign(1, -1);
        filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
        check(filters.doFilterAtAll(), "fixture options enable filtering and a second tier");
        auto derep = make_shared<Dereplicate>(&opts, &filters); derep->BCnamesAdding(&filters);
        OutputStreamer output(&filters, &opts, ios::out, nullptr, 0);
        output.attachDereplicator(derep);
        auto midRead = [&](const string& seq, const string& name) {
            auto d = derepRead(seq, name, 20, false); d->setYellowQual(true); return d;
        };
        // Same sequence of calls as process_DNA: stats/prepare, dereplicate, write unless withheld.
        auto process = [&](const shared_ptr<DNA>& d) {
            int cstream = 100;
            output.saveForWrite(d, 1, 0, cstream, false);
            const bool withheld = output.dereplicateDNA(d, nullptr, 0, 1, 2);
            if (!withheld) { output.writeForWrite(d, 1, cstream, nullptr, -1, 100, 0); }
            return withheld;
        };
        check(process(midRead(shared, "mid_before")), "a mid-quality read is withheld for the dereplicator");
        check(!process(derepRead(shared, "green", 40, true)), "a passing read is written by the caller");
        check(process(midRead(shared, "mid_after")), "a later mid-quality copy is withheld too");
        check(process(midRead(alone, "mid_alone")), "a mid-quality read without any passing copy is withheld");
        output.flushDereplication();
        const string report = derep->writeDereplDNA(&filters, "");
        output.closeOutStreams(true);
        const string main = readScratch(base + ".main.fna"), add = readScratch(base + ".add.fna"), fa = readScratch(base + ".fa");
        check(fa.find(">green;size=2;") != string::npos, "the mid-quality copy after the passing read is counted");
        check(main.find(">green") != string::npos && main.find("mid_") == string::npos, "the primary output holds only the passing read");
        check(add.find(">mid_before") != string::npos && add.find(">mid_alone") != string::npos,
            "uncounted mid-quality reads are published to the mid-quality output");
        check(add.find(">mid_after") == string::npos, "a counted mid-quality read is not published");
        check(report.find("2 published to the mid-quality output") != string::npos, "the report counts published sets");
        derep->finishMap();
    }
    setGlobalOptContainer(previousGlobal);
    for (const char* suffix : {".fa", ".fa.rest", ".map", ".1.hq.fq", ".2.hq.fq", ".main.fna", ".add.fna", ".options.txt"}) {
        std::remove((base + suffix).c_str());
    }
}

void testRejectedMergeLeavesNoMergeMarks() {
    // Finding 2: ReadMerger::merge marked both mates as merged before the
    // consensus was accepted; only DNA::adoptMergeStats marks them now.
    const string amplicon = lcgSequence(150, 5);
    string second = amplicon.substr(50);
    reverseTS(second);
    auto r1 = derepRead(amplicon.substr(0, 100), "pair/1", 40, true);
    auto r2 = derepRead(second, "pair/2", 40, true);
    ReadMerger merger;
    check(merger.findSeedForMerge(r1, r2), "fixture mates overlap by 50 bases");
    auto merged = merger.merge(r1, r2);
    check(merged != nullptr && merged->getMergeLength() == static_cast<int>(merged->length()),
        "the consensus carries its own merge length");
    check(r1->getMergeLength() == -1 && r2->getMergeLength() == -1 && r1->getMergeErrors() == -1,
        "merging leaves the mates unmarked until the consensus is accepted");
    r1->adoptMergeStats(*merged);
    check(r1->getMergeLength() == merged->getMergeLength() && r1->getMergeErrors() == merged->getMergeErrors(),
        "an accepted consensus marks the mate with its length and overlap errors");
}

void testBarcodeMatchPrefersBestOffset() {
    // Open finding 3: the first offset within budget won, so a one-error match at
    // offset 0 shadowed the exact barcode at offset 1.
    const string barcode = "ACGTTG";
    auto read = derepRead("CACGTTGTTTTTTTTTTTTT", "bc", 30, true); // exact at offset 1
    int errors = -1;
    check(read->matchSeq_tot(barcode, 1, 5, errors) == 1 && errors == 0, "exact match at offset 1 wins");
    auto shadow = derepRead("ACGTTAACGTTGTTTTTTTT", "bc2", 30, true); // 1 error at 0, exact at 6
    check(shadow->matchSeq_tot(barcode, 1, 8, errors) == 6 && errors == 0, "a later exact match beats an earlier one-error match");
    check(shadow->matchSeq_tot(barcode, 1, 5, errors) == 0 && errors == 1, "within the window the best available offset is used");
    auto none = derepRead("TTTTTTTTTTTTTTTTTTTT", "bc3", 30, true);
    check(none->matchSeq_tot(barcode, 1, 8, errors) == -1 && errors == 2, "no match reports budget + 1 errors");
}

void testThreadPoolWorkerSubmitRunsInline() {
    // Every worker submitting into its own full queue used to block forever.
    ThreadPool& pool = ThreadPool::instance();
    const size_t workers = pool.size(), perWorker = 8 * workers + 8;
    std::atomic<size_t> ran{0};
    std::mutex innerMutex;
    vector<std::future<void>> inner;
    vector<std::future<void>> outer;
    for (size_t w = 0; w < workers; ++w) {
        outer.push_back(pool.submit([&] {
            for (size_t i = 0; i < perWorker; ++i) {
                auto task = pool.submit([&ran] { ran.fetch_add(1); });
                std::lock_guard<std::mutex> guard(innerMutex);
                inner.push_back(std::move(task));
            }
        }));
    }
    for (auto& task : outer) { task.get(); }
    for (auto& task : inner) { task.get(); }
    check(ran.load() == workers * perWorker, "every task submitted by a worker ran exactly once");
}

void testPartitionPoolAfterEnqueueFailure() {
    // A throwing afterEnqueue used to unwind while published jobs pointed into the frame.
    PartitionWorkPool pool(4, 1);
    std::atomic<int> executed{0};
    const std::array<PartitionWorkPool::Work, 4> work{{{0}, {1}, {2}, {3}}};
    bool rethrown = false;
    try { pool.run(work, [&](size_t) { executed.fetch_add(1); }, [] { throw std::runtime_error("afterEnqueue"); }); }
    catch (const std::runtime_error&) { rethrown = true; }
    check(rethrown && executed.load() == 4, "the batch completes before the afterEnqueue error is rethrown");
    pool.run(work, [&](size_t) { executed.fetch_add(1); });
    check(executed.load() == 8, "the pool stays usable afterwards");
}

void testApplyFileITKeepsGzipAndDirectories() {
    check(applyFileIT("out.fq.gz", 1) == "out.1.fq.gz", "inner extension gets the iteration");
    check(applyFileIT("out.gz", 1) == "out.1.gz", "a .gz file without inner extension stays compressed");
    check(applyFileIT("out.gz", 0, ".singl") == "out.singl.gz", "a suffix goes before .gz");
    check(applyFileIT("out.fq", 2) == "out.2.fq", "plain extension");
    check(applyFileIT("/data/run.2024/out", 1) == "/data/run.2024/out.1", "a dot in a directory is not an extension");
    check(applyFileIT("/data/run.2024/out.fq", 1) == "/data/run.2024/out.1.fq", "directory dots do not disturb the basename");
    check(applyFileIT("out", 0, ".singl") == "out.singl", "no extension at all");
}

void testShortOverlapMergesNeedExactMatch() {
    // Finding 8: 12-19 base overlaps with a mismatch were accepted from 5-mer seeds.
    const string amplicon = lcgSequence(150, 31);
    auto mates = [&](size_t overlap, bool mismatch) {
        string second = amplicon.substr(100 - overlap);
        if (mismatch) { char& base = second[overlap / 2]; base = base == 'A' ? 'C' : 'A'; }
        reverseTS(second);
        return make_pair(derepRead(amplicon.substr(0, 100), "m/1", 40, true), derepRead(second, "m/2", 40, true));
    };
    ReadMerger merger;
    check(ReadMerger::minOverlap() == 16, "default minimum overlap is 16 bases");
    auto [r1a, r2a] = mates(14, false);
    check(!merger.findSeedForMerge(r1a, r2a), "a 14 base overlap is below the minimum");
    auto [r1b, r2b] = mates(18, false);
    check(merger.findSeedForMerge(r1b, r2b), "an exact 18 base overlap merges");
    auto [r1c, r2c] = mates(18, true);
    check(!merger.findSeedForMerge(r1c, r2c), "an 18 base overlap with a mismatch does not merge");
    auto [r1d, r2d] = mates(40, true);
    check(merger.findSeedForMerge(r1d, r2d), "a 40 base overlap tolerates a mismatch");
    ReadMerger::setMinOverlap(20);
    check(!merger.findSeedForMerge(r1b, r2b), "MergeMinOverlap raises the minimum");
    ReadMerger::setMinOverlap(16);
}

void testBetterPrefersCleanOverlap() {
    // Finding 4: overlap errors of merged consensus reads now enter the score.
    const string seq = lcgSequence(120, 41), merged = lcgSequence(200, 41);
    auto candidate = derepRead(seq, "cand", 35, true), reference = derepRead(seq, "ref", 35, true);
    auto candidateMerge = derepRead(merged, "candM", 35, true), referenceMerge = derepRead(merged, "refM", 35, true);
    for (auto& d : {candidate, reference, candidateMerge, referenceMerge}) { d->setMergeLength(200); }
    candidate->setMergeErrors(0, 0); candidateMerge->setMergeErrors(0, 0);
    reference->setMergeErrors(6, 20); referenceMerge->setMergeErrors(6, 20);
    auto history = make_shared<DNAunique>(reference, 0);
    float everBest = -1.f;
    check(whoIsBetter(candidate, nullptr, candidateMerge, reference, nullptr, referenceMerge, everBest, false, history.get()),
        "a consensus without overlap errors beats an otherwise equal one with six");
    everBest = -1.f;
    check(!whoIsBetter(reference, nullptr, referenceMerge, candidate, nullptr, candidateMerge, everBest, false, history.get()),
        "a consensus with overlap errors does not replace a clean one");
}

#if defined(__linux__)
void testOutputFinalizeFailureIsFatal() {
    // Finding 2: a delayed write failure on a non-retained main output (and
    // demultiplexed files) was only printed by ~ofbufstream; closeOutStreams
    // must raise it.
    const string outPath = scratchPath("full_out.fq");
    std::remove(outPath.c_str());
    std::filesystem::create_symlink("/dev/full", outPath);

    OptContainer options;
    options["-sample_sep"] = "___";
    options["-paired"] = "1";
    options["-onlyPair"] = "";
    options["-i_qual_offset"] = "33";
    options["-o_qual_offset"] = "33";
    options["-5PR1cut"] = "0";
    options["-5PR2cut"] = "0";
    options["-logLvsQ"] = "";
    options["-GoldenAxe"] = "0";
    options["-o_fastq"] = outPath;
    options["-o_fna"] = "";
    options["-o_qual"] = "";
    options["-o_fastq2"] = "";
    options["-o_fastq_noBC"] = "";
    options["-suppressOutput"] = "0";
    options["-threadIO"] = "0";
    options["-uparseVer"] = "";
    options["-count_chimeras"] = "F";
    OptContainer* previousGlobal = getGlobalOptContainer();
    setGlobalOptContainer(&options);
    bool raised = false;
    string message;
    {
        Filters filters(&options);
        OutputStreamer output(&filters, &options, ios::out, nullptr, 0);
        vector<string> lines{"@r1", "ACGTACGTAC", "+", "IIIIIIIIII"};
        shared_ptr<DNA> read = str2DNA(lines, true, 33, 0);
        int cstream = 0;
        output.saveForWrite(read, 1, 0, cstream, true);
        try { output.closeOutStreams(true); }
        catch (const exception& e) { raised = true; message = e.what(); }
    }
    setGlobalOptContainer(previousGlobal);
    std::remove(outPath.c_str());
    check(raised, "closeOutStreams must raise the finalize failure of a non-retained output");
    check(message.find(outPath) != string::npos, "the failure names the output file: " + message);
}
#endif

}  // namespace

int main(int argc, char** argv) {
    ini_DNAconstants();
    const string only = argc > 1 ? argv[1] : "";  // optional: run a single named test
    const vector<pair<string, function<void()>>> tests = {
        {"header_stem", testHeaderStem},
        {"spaceX", testSpaceX},
        {"is_digits", testIsDigits},
        {"getFileNoPath", testGetFileNoPath},
        {"primer_variant_tracker", testPrimerVariantTracker},
        {"derep_subgroup_map", testDerepSubgroupMap},
        {"compact_optional_string", testCompactOptionalString},
        {"mate_header_roundtrip", testMateHeaderRoundTrip},
        {"derep_read_size", testDerepReadIsSmaller},
        {"cutseq_short_qualities", testCutSeqWithShortQualities},
        {"negative_homopolymer", testNegativeHomopolymerDisables},
        {"change_head_pver", testChangeHeadPverKeepsOriginal},
        {"set_mate_number", testSetMateNumber},
        {"matchseq_high_bit", testMatchSeqHighBitPrimer},
        {"merge_stats_counts", testMergeStatsCountBeyondFloat},
        {"goldenaxe_merge", testGoldenAxeMergeKeepsAllColumns},
        {"diffdna_partial_pairs", testDiffDNAClearKeepsPartialPairs},
        {"fastq_offset_autodetect", testFastqOffsetAutoDetect},
        {"exact_derep_recovery_order_independent", testExactDerepRecoveryIsOrderIndependent},
        {"rejected_merge_leaves_no_marks", testRejectedMergeLeavesNoMergeMarks},
        {"mid_quality_counted_or_published", testMidQualityCountedOrPublished},
        {"apply_file_it_naming", testApplyFileITKeepsGzipAndDirectories},
        {"barcode_match_best_offset", testBarcodeMatchPrefersBestOffset},
        {"threadpool_worker_submit_inline", testThreadPoolWorkerSubmitRunsInline},
        {"partition_pool_after_enqueue_failure", testPartitionPoolAfterEnqueueFailure},
        {"short_overlap_merges_exact", testShortOverlapMergesNeedExactMatch},
        {"better_prefers_clean_overlap", testBetterPrefersCleanOverlap},
#if defined(__linux__)
        {"output_finalize_failure_fatal", testOutputFinalizeFailureIsFatal},
#endif
    };
    size_t ran = 0;
    for (const auto& [name, test] : tests) {
        if (!only.empty() && name != only) { continue; }
        ++ran;
        currentTest = name;
        const size_t before = failures.size();
        try { test(); } catch (const exception& e) { failures.push_back(name + ": unexpected exception: " + e.what()); }
        cout << (failures.size() == before ? "PASS " : "FAIL ") << name << "\n";
    }
    for (const auto& failure : failures) { cerr << "  " << failure << "\n"; }
    if (ran == 0) {
        cerr << "No test named '" << only << "'\n";
        return 2;
    }
    if (!failures.empty()) {
        cerr << failures.size() << " audit-fix check(s) failed\n";
        return 1;
    }
    cout << "All " << ran << " audit-fix unit tests passed\n";
    return 0;
}
