// Tests for the streaming candidate index, nucleotide deltas, and SDM integration.
#include "CoarseDereplication.h"
#include "Filters.h"
#include "containers.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <thread>
#include <atomic>
#include <functional>

namespace {
void check(bool ok, const std::string& message) { if (!ok) { throw std::runtime_error(message); } }
template<class F> void throws(F f) {
    bool caught = false; try { f(); } catch (const std::exception&) { caught = true; }
    check(caught, "Expected invalid input to throw");
}
std::string randomSequence(std::mt19937_64& rng, size_t length) {
    std::string sequence(length, 'A'); for (char& c : sequence) { c = "ACGT"[rng() % 4]; } return sequence;
}
void mutate(char& base) { base = base == 'A' ? 'C' : 'A'; }
std::shared_ptr<DNA> readDNA(const std::string& sequence, int quality = 25, int sample = 0, bool pass = true) {
    auto dna = std::make_shared<DNA>(sequence, "read_q" + std::to_string(quality));
    dna->setQual(std::vector<qual_score>(sequence.size(), quality));
    check(dna->seal(), "Seal test read"); dna->setBCnumber(sample, 0); dna->setPassed(pass); return dna;
}
size_t add(CoarseDereplication& index, const std::string& sequence, int quality = 25, bool pass = true) {
    return index.add(sequence, sequence, readDNA(sequence, quality, 0, pass), nullptr, nullptr, 0, pass);
}
void testExactOutput() {
    std::mt19937_64 rng(2026091418);
    const auto a = randomSequence(rng, 200), mate = randomSequence(rng, 230);
    auto b = a, rejected = a, differentMate = mate;
    mutate(b[31]); mutate(rejected[32]);
    for (size_t i = 0; i < 20; ++i) { mutate(differentMate[i]); }
    for (bool qualities : {false, true}) {
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, true, false, qualities, true);
        auto insert = [&](const auto& sequence, const auto& r2, int q, int sample, bool passing = true) {
            return index.add(sequence, sequence, readDNA(sequence, q, sample, passing),
                readDNA(r2, q, sample, passing), nullptr, sample, passing);
        };
        check(insert(a, mate, 20, 0) == 0 && insert(b, mate, 30, 1) == 0,
            "Nearby R1 keys share a coarse reference");
        check(index.groups().size() == 1 && index.exactGroups().size() == 2,
            "Coarse search retains two distinct exact output parents");
        check(index.groups()[0]->representatives.size() == 0,
            "Exact output does not duplicate coarse representatives/counts/qualities");
        const auto searched = index.stats().comparisons;
        check(insert(a, differentMate, 40, 1) == 0 && index.stats().comparisons == searched,
            "Exact R1 reuse bypasses search despite a divergent R2");
        check(insert(rejected, mate, 60, 1, false) == CoarseDereplication::npos,
            "A nearby coarse anchor cannot admit a new failed exact R1 key");
        check(index.stats().comparisons == searched && index.stats().rejectedExactMisses == 1,
            "A rejected exact-key miss performs no coarse candidate comparisons");
        check(insert(a, differentMate, 60, 2, false) == 0,
            "An established exact R1 still counts failed repeats");
        auto& exactA = *index.exactGroups().find(a)->second.representatives;
        auto& exactB = *index.exactGroups().find(b)->second.representatives;
        auto snapshot = exactA.best(true, false);
        check(snapshot->totalSum() == 3 && exactB.totalSum() == 1,
            "Abundances remain local to the exact R1 dereplicate");
        check(snapshot->getDerepCounts().sampleCounts() == diffDNA::SampleCounts({{0,1},{1,1},{2,1}}),
            "Exact R1 sample provenance survives coarse grouping");
        check(snapshot->getSeqPseudoView() == a && snapshot->getPair()->getSeqPseudoView() == differentMate &&
            snapshot->getQual() == std::vector<qual_score>(200, 40),
            "The better passing pair replaces both mates only within its exact R1 key");
        auto differences = snapshot->getDiffDNA();
        check(&differences->reference() == index.groups()[0]->reference.get(),
            "Exact parents share the immutable compression reference");
        check(differences->uniqueCount() == (qualities ? 2 : 1), "Retained full-pair variants keep their own counts");
        insert(a, mate, 45, 0);
        check(snapshot->totalSum() == 3 && exactA.totalSum() == 4, "Output snapshots leave live exact counts independent");
        check(insert(rejected, mate, 40, 1) == 0 && index.exactGroups().size() == 3 &&
            index.exactGroups().at(rejected).representatives->totalSum() == 1,
            "Rejecting a failed exact-key miss does not prevent its later passing anchor");
        index.clear();
        check(index.stats().queries == 0 && index.stats().rejectedExactMisses == 0,
            "Reset clears admission shortcut statistics");
        check(index.groups().empty() && index.exactGroups().empty() && snapshot->totalSum() == 3,
            "Reset clears exact parents while existing snapshots stay valid");
    }
}
void testCompactStorage(const std::string& base) {
    std::mt19937_64 rng(2026091531);
    const auto first = randomSequence(rng, 245), second = randomSequence(rng, 235);
    for (bool fasta : {true, false}) {
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, false, false, false, true);
        std::map<std::string, DNAuniqSet> ordinary;
        std::map<std::string, std::map<std::pair<std::string, std::string>, uint64_t>> expectedVariants;
        for (int i = 0; i < 1000; ++i) {
            auto r1 = first, r2 = second;
            r1[20] = "ACGT"[i % 4]; r1[21] = "ACGT"[i / 4 % 4];
            r1[220] = "ACGT"[i / 16 % 4]; r2[205] = "ACGT"[i / 64 % 4];
            const auto key = r1.substr(0, 200);
            auto a = readDNA(r1, 25 + (i / 16 % 3) * 7, i % 3);
            auto b = readDNA(r2, 25 + (i / 16 % 3) * 7, i % 3);
            auto q = a->getQual();
            // Exercise incompressible selected qualities, and variable coverage.
            for (size_t j = 0; j < q.size(); ++j) { q[j] -= (j * 17 + i) % 5; }
            a->setQual(q); a->setHeader("compact_" + std::to_string(i)); b->setHeader(a->getId());
            a->cutSeq(200 + i % 30, -1, true); b->cutSeq(200 + i % 20, -1, true);
            if (i >= 16 && i % 11 == 0) { a->setPassed(false); }
            if (i % 17 == 0) { b->setPassed(false); }
            if (i % 37 == 0) { b->delself(); } // Physical empty mate is still admitted ordinarily.
            std::shared_ptr<DNA> merged;
            if (i % 3) {
                merged = readDNA(r1 + r2.substr(130), 30 + i % 10, i % 3);
                a->setMergeLength(merged->length()); a->merge_seed_pos_ = i % 7;
            }
            auto oa = std::make_shared<DNA>(*a), ob = std::make_shared<DNA>(*b);
            auto om = merged ? std::make_shared<DNA>(*merged) : nullptr;
            const int position = i % 7;
            ordinary[key].addObservation(oa, ob, om, position, i % 3, fasta, true);
            check(index.add(key, a->getSeqPseudoView(), a, b, merged, i % 3, a->isGreenQual(), 200,
                fasta, true, position) != CoarseDereplication::npos, "Established compact key admitted");
            ++expectedVariants[key][{a->getSequence(), b->getSequence()}];
        }
        check(index.groups().size() < index.exactGroups().size(), "Exact keys share compact references");
        for (const auto& [key, assignment] : index.exactGroups()) {
            auto& compact = *assignment.representatives;
            check(compact.isCompact() && compact.begin() == compact.end(), "No persistent full DNA objects per compact exact key");
            auto expected = ordinary[key].best(true, !fasta), actual = compact.best(true, !fasta);
            check(expected->getDerepMap() == actual->getDerepMap(), "Compact sample counts");
            check(expected->getWhoIsBetterMergeStats().merged_compares == actual->getWhoIsBetterMergeStats().merged_compares &&
                expected->getWhoIsBetterMergeStats().non_merged_compares == actual->getWhoIsBetterMergeStats().non_merged_compares &&
                expected->getWhoIsBetterMergeStats().accumulated_merge_length == actual->getWhoIsBetterMergeStats().accumulated_merge_length,
                "Compact seed comparison history");
            expected->Count2Head(true); actual->Count2Head(true);
            check(expected->getId() == actual->getId(), "Compact selected representative ID");
            for (const auto& pair : {std::make_pair(std::shared_ptr<DerepRead>(expected), std::shared_ptr<DerepRead>(actual)),
                    std::make_pair(expected->getPair(), actual->getPair()),
                    std::make_pair(expected->getMerge(), actual->getMerge())}) {
                check(bool(pair.first) == bool(pair.second), "Compact mate/merge presence");
                if (!pair.first) { continue; }
                check(pair.first->getSequence() == pair.second->getSequence() && pair.first->length() == pair.second->length() &&
                    pair.first->getQual() == pair.second->getQual(), "Lossless full backing sequence and selected qualities");
            }
            expected->prepareDerepQualities(33); actual->prepareDerepQualities(33);
            { std::ofstream a(base + ".expected"), b(base + ".observed"); expected->writeDerepFastQ(a); actual->writeDerepFastQ(b); }
            std::ifstream a(base + ".expected"), b(base + ".observed");
            check(std::string(std::istreambuf_iterator<char>(a), {}) == std::string(std::istreambuf_iterator<char>(b), {}),
                "Compact main FASTQ keeps per-base quality evidence");
            std::map<std::pair<std::string, std::string>, uint64_t> observed;
            for (const auto& variants : compact.compactVariants()) {
                check(!variants.storesQualities(), "No quality vector allocated per full-fragment variant");
                variants.forEachPair([&](auto r1, auto r2, uint64_t count, auto*, auto*) { observed[{std::string(r1), std::string(r2)}] += count; });
            }
            check(observed == expectedVariants[key], "All 1000 admitted fragment observations retained as nucleotide deltas/counts");
        }
    }
}

void testOrdinaryStorageUpdates(const std::string& base) {
    std::mt19937_64 rng(2026091522);
    const auto first = randomSequence(rng, 200), mate = randomSequence(rng, 230);
    for (bool fasta : {false, true}) for (bool reverse : {false, true}) {
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, false, false, false, true);
        DNAuniqSet ordinary;
        for (int n = 0; n < 4; ++n) {
            int i = reverse ? 3 - n : n;
            auto a = readDNA(first, 20 + 10 * i, i % 2), b = readDNA(mate, 25 + 10 * i, i % 2);
            a->setHeader("pair" + std::to_string(i)); b->setHeader(a->getId());
            const int offset = i == 3 ? 20 : 10;
            // An independent call through the original ordinary API is the oracle.
            auto found = ordinary.find(offset);
            if (found == ordinary.end()) { ordinary.addNewDNAuniq(a, b, nullptr, offset, i % 2, fasta); }
            else { found->second->matchedDNA(a, b, nullptr, i % 2, fasta, true); }
            index.add(first, first, a, b, nullptr, i % 2, true, 200, fasta, true, offset);
        }
        auto expected = ordinary.best(true, !fasta);
        auto observed = index.exactGroups().at(first).representatives->best(true, !fasta);
        expected->Count2Head(true); observed->Count2Head(true);
        check(expected->getId() == observed->getId() && expected->getId().starts_with("pair2;"),
            "Storage keeps ordinary abundance-based merge-offset representative selection: " + expected->getId() + " / " + observed->getId());
        check(expected->getQual() == observed->getQual() && expected->getPair()->getQual() == observed->getPair()->getQual(),
            "Storage keeps both observed representative quality strings");
        expected->prepareDerepQualities(33); observed->prepareDerepQualities(33);
        { std::ofstream a(base + ".expected"), b(base + ".observed"); expected->writeDerepFastQ(a); observed->writeDerepFastQ(b); }
        std::ifstream a(base + ".expected"), b(base + ".observed");
        check(std::string(std::istreambuf_iterator<char>(a), {}) == std::string(std::istreambuf_iterator<char>(b), {}) && observed->totalSum() == 4 &&
            expected->getDerepMap() == observed->getDerepMap(), "Storage keeps ordinary main FASTQ qualities and sample counts");
    }
}
void testStorageOrderEffects() {
    std::mt19937_64 rng(2026091525);
    const auto first = randomSequence(rng, 200), mate = randomSequence(rng, 230);
    for (bool reverse : {false, true}) {
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, false, false, false, true);
        DNAuniqSet ordinary;
        for (int n = 0; n < 2; ++n) {
            const int i = reverse ? 1 - n : n;
            auto a = readDNA(first, 40), b = readDNA(mate, 40);
            a->setHeader("tie" + std::to_string(i)); b->setHeader(a->getId());
            if (n == 0) { ordinary.addNewDNAuniq(a, b, nullptr, 0, 0, false); }
            else { ordinary[0]->matchedDNA(a, b, nullptr, 0, false); }
            index.add(first, first, a, b, nullptr, 0, true, 200, false, true);
        }
        auto a = ordinary.best(true, true), b = index.exactGroups().at(first).representatives->best(true, true);
        a->Count2Head(true); b->Count2Head(true);
        const auto expected = std::string(reverse ? "tie1" : "tie0") + ";size=2;";
        check(a->getId() == expected && b->getId() == expected,
            "Exact ties retain the first observation in both paths; reversing order changes the ID");
        index.clear();
        bool anchored = false; unsigned ordinaryCount = 0;
        for (int n = 0; n < 3; ++n) {
            const bool pass = reverse ? n == 0 : n == 2;
            if (pass || anchored) { ++ordinaryCount; anchored = true; }
            index.add(first, first, readDNA(first, pass ? 40 : 10, 0, pass), nullptr, nullptr, 0, pass);
        }
        check(ordinaryCount == (reverse ? 3 : 1) && index.exactGroups().at(first).representatives->totalSum() == ordinaryCount,
            "Fixed-key recovery preserves ordinary before/after-anchor admission dependence");
    }
}
void testMixedPairAdmission() {
    std::mt19937_64 rng(2026091502);
    const auto first = randomSequence(rng, 200), second = randomSequence(rng, 230);
    for (const char* identity : {"100", "97"}) {
        CoarseDereplication index(CoarseDereplication::parseIdentity(identity), 4, false, false, true, true);
        DNAuniqSet ordinary;
        auto observe = [&](int q, int sample, bool r1Pass, bool r2Pass) {
            auto r1 = readDNA(first, q, sample, r1Pass), r2 = readDNA(second, q, sample, r2Pass);
            if (ordinary.size() == 0) { ordinary.addNewDNAuniq(r1, r2, nullptr, 0, sample, true); }
            else { ordinary[0]->matchedDNA(r1, r2, nullptr, sample, true, false); }
            check(index.add(first, first, r1, r2, nullptr, sample, r1Pass) == 0,
                "Passing R1 establishes an exact parent even when R2 failed");
        };
        observe(20, 0, true, false);
        observe(40, 1, true, false);
        observe(50, 2, false, true);
        auto retained = index.exactGroups().at(first).representatives->best(true, false);
        check(retained->getQual() == ordinary[0]->getQual() &&
            retained->getPair()->getQual() == ordinary[0]->getPair()->getQual() &&
            retained->getQual() == std::vector<qual_score>(200, 40),
            "R1 eligibility and paired replacement agree with ordinary dereplication");
        check(retained->getDerepCounts().sampleCounts() == ordinary[0]->getDerepCounts().sampleCounts(),
            "Mixed-quality fragments contribute the same sample counts");
        auto fresh = first; mutate(fresh[31]);
        check(index.add(fresh, fresh, readDNA(fresh, 40, 0, false), readDNA(second), nullptr, 0, false) ==
            CoarseDereplication::npos, "Passing R2 alone cannot establish a new exact R1 parent");
    }
}
void testExactReuse() {
    std::mt19937_64 rng(2026091409);
    const auto anchor = randomSequence(rng, 200);
    auto variant = anchor; mutate(variant[31]);
    CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, true, true, true);
    auto insert = [&](const std::string& sequence, int q, int sample, bool passed = true) {
        return index.add(sequence, sequence, readDNA(sequence, q, sample, passed), nullptr, nullptr, sample, passed);
    };
    check(insert(anchor, 20, 0) == 0 && insert(variant, 10, 0) == 0, "First coarse match establishes exact assignment");
    const auto searched = index.stats();
    check(insert(variant, 40, 1) == 0 && insert(variant, 60, 2, false) == 0,
        "Passing and failed exact repeats reuse the same immutable anchor");
    check(index.stats().exactHits == searched.exactHits + 2 && index.stats().comparisons == searched.comparisons,
        "Exact repeats bypass candidate verification");
    auto differences = index.groups()[0]->representatives[0]->getDiffDNA();
    check(differences->readCount() == 4 && differences->uniqueCount() == 2,
        "Exact hits still record every exact variant and abundance");
    differences->forEachWithQualities([&](auto sequence, auto count, const auto* samples, const auto& quality) {
        if (sequence != variant) { return; }
        check(count == 3 && samples && *samples == diffDNA::SampleCounts({{0,1},{1,1},{2,1}}),
            "Exact matches retain sample provenance");
        check(quality.seed && quality.seed->qualities == std::vector<qual_score>(200, 40),
            "Exact hits improve the selected quality, but failed reads cannot replace it");
    });
    const auto beforeSnapshot = index.stats();
    (void)index.refine(false, -1);
    check(index.find(variant) == 0 && index.stats().comparisons == beforeSnapshot.comparisons,
        "Refinement snapshots do not invalidate live exact assignments");
    index.clear();
    check(index.stats().queries == 0 && index.stats().exactHits == 0 && index.find(variant) == CoarseDereplication::npos,
        "Clear discards exact IDs and statistics");
    check(insert(variant, 60, 2, false) == CoarseDereplication::npos && insert(variant, 30, 0) == 0,
        "An earlier miss cannot hide a later passing anchor");
    check(insert(variant, 60, 2, false) == 0 && index.stats().exactHits == 1,
        "An admitted key is indexed after creation");

    // One effective search prefix can stand for several full nucleotide variants.
    CoarseDereplication capped(CoarseDereplication::identityScale, 4, true, false, true);
    auto tail = anchor; mutate(tail[170]);
    const auto prefix = anchor.substr(0, 100);
    check(capped.add(prefix, anchor, readDNA(anchor), nullptr, nullptr, 0, true, 100) == 0 &&
        capped.add(prefix, tail, readDNA(tail, 40), nullptr, nullptr, 1, true, 100) == 0,
        "Identical search keys reuse one parent");
    check(capped.stats().exactHits == 1 && capped.groups()[0]->representatives[0]->getDiffDNA()->uniqueCount() == 2,
        "Full tails remain distinct exact variants on exact hits");

    CoarseDereplication paired(CoarseDereplication::identityScale, 4, false, false, true);
    const auto mateA = randomSequence(rng, 200), mateB = randomSequence(rng, 200);
    auto first = readDNA(anchor, 10), second = readDNA(tail, 50);
    check(paired.add(prefix, anchor, first, readDNA(mateA, 10), nullptr, 0, true, 100) == 0 &&
        paired.add(prefix, tail, second, readDNA(mateB, 55), nullptr, 1, true, 100) == 0,
        "Exact R1 dereplication ignores R2 identity");
    const auto representative = paired.groups()[0]->representatives[0];
    check(paired.stats().exactHits == 1 && representative->getSeqPseudoView() == tail &&
        representative->getQual() == std::vector<qual_score>(200, 50) &&
        representative->getPair()->getSeqPseudoView() == mateB &&
        representative->getPair()->getQual() == std::vector<qual_score>(200, 55),
        "A better pair with different full R1 and R2 replaces both representative mates");
    check(representative->getDiffDNA()->uniqueCount() == 2 && representative->totalSum() == 2,
        "Distinct full pairs remain separate variants under one R1 dereplicate");

    // Assignments must persist; eviction would recheck R2 and split an old R1.
    CoarseDereplication many(CoarseDereplication::identityScale, 4, false, false, true);
    std::vector<std::string> keys;
    for (size_t i = 0; i < 1100; ++i) {
        keys.push_back(randomSequence(rng, 100));
        check(many.add(keys.back(), keys.back(), readDNA(keys.back()), readDNA(mateA), nullptr, 0, true) == i,
            "Distinct R1 keys establish exact dereplicates");
    }
    const auto before = many.stats();
    check(many.add(keys.front(), keys.front(), readDNA(keys.front(), 50), readDNA(mateB, 55), nullptr, 1, true) == 0 &&
        many.stats().exactHits == before.exactHits + 1 && many.stats().comparisons == before.comparisons,
        "Older R1 dereplicates survive index growth and accept different R2");
    check(many.groups()[0]->representatives.totalSum() == 2, "Persistent exact reuse preserves abundance");
    auto probe = keys[0]; probe.push_back('A');
    check(many.find(probe, mateA) == 0, "A new R1 prefix can find a coarse candidate");
    check(many.add(probe, probe, readDNA(probe), readDNA(mateB), nullptr, 0, true) == keys.size(),
        "Read-only searches do not create authoritative R1 assignments");

}

void testSeedsAndThresholds() {
    uint64_t word = 0;
    check(!CoarseDereplication::seedWord("ACGTAAAAAACGTACGTACG", 0, 4, word), "Exclude a seed containing a homopolymer");
    check(!CoarseDereplication::seedWord(std::string(20, 'A'), 0, 20, word), "Never index a pure homopolymer");
    check(CoarseDereplication::seedWord("ACGTAAACGTACGTACGTAC", 0, 4, word), "Retain seeds with runs below cutoff");
    check(!CoarseDereplication::seedWord("ACGTACGTNCGTACGTACGT", 0, 4, word), "Exclude ambiguous seed");
    for (const auto& bad : {"", "0", "100.1", "-99", "nan", "98x", "99.", "1.0000001"}) {
        throws([&] { CoarseDereplication::parseIdentity(bad); });
    }
    for (const auto& bad : {"0", "1", "21", "4x", "4.0", "-4"}) {
        throws([&] { CoarseDereplication::parseHomopolymer(bad); });
    }
    CoarseDereplication index(CoarseDereplication::parseIdentity("97.5"), 4, false);
    check(index.mismatchBudget(1500) == 37 && index.mismatchBudget(40) == 1, "Exact decimal mismatch boundary");
    CoarseDereplication tiny(CoarseDereplication::parseIdentity("99.999999"), 4, false);
    check(tiny.mismatchBudget(99999999) == 0 && tiny.mismatchBudget(100000000) == 1, "Six-decimal boundary");
}
void testDiffs() {
    const std::string anchor = "ACGTACGTACGT";
    diffDNA delta(anchor);
    std::vector<std::string> sequences = {anchor, anchor, "TCGTACGTACGA", "ACGT", anchor + "NGA", "", "ACGTNCGTACGT", "TCGTACGTACGA", "ACGT"};
    std::map<std::string, uint64_t> expected;
    for (const auto& sequence : sequences) { delta.add(sequence); ++expected[sequence]; }
    auto verify = [&](const diffDNA& d) {
        std::map<std::string, uint64_t> observed;
        d.forEach([&](std::string_view sequence, uint64_t count) { observed[std::string(sequence)] += count; });
        check(observed == expected, "Nucleotide haplotypes, suffixes, shorter and empty reads round trip");
        check(d.readCount() == sequences.size() && d.referenceCount() == 2 && d.uniqueCount() == expected.size(), "Delta multiplicities and unique sequences");
    };
    verify(delta);
    std::ostringstream fasta; delta.writeFasta(fasta, "cluster;size=9;");
    std::istringstream fastaInput(fasta.str()); std::string header, sequence;
    std::map<std::string, uint64_t> fastaCounts;
    while (std::getline(fastaInput, header)) {
        check(std::getline(fastaInput, sequence).good(), "Subcluster FASTA sequence line");
        const auto size = header.find(";size=");
        check(header.starts_with(">cluster.sub") && size != std::string::npos && header.find(";size=", size + 1) == std::string::npos,
            "Each subcluster has one abundance attribute and a parent label");
        fastaCounts[sequence] += std::stoull(header.substr(size + 6));
    }
    check(fastaCounts == expected, "FASTA exports one exact subcluster per sequence with its own count");
    diffDNA repeated(anchor); repeated.add("TCGTACGTACGA", 5); const auto patchBytes = repeated.payloadBytes();
    repeated.add("TCGTACGTACGA", 7);
    check(repeated.uniqueCount() == 1 && repeated.readCount() == 12 && repeated.payloadBytes() == patchBytes,
        "Repeated non-reference reads increment a count without copying their patch");
    diffDNA large(anchor); large.add("TCGTACGTACGA", std::numeric_limits<uint64_t>::max());
    throws([&] { large.add(anchor); });
    check(large.readCount() == std::numeric_limits<uint64_t>::max(), "Overflow rejection preserves counts");
    // Version 1 stored repeated patches separately. Import coalesces them.
    std::string legacy = "SDMDIFF1";
    for (unsigned char byte : std::initializer_list<unsigned char>{1, 'c', 1, 'A', 3, 1, 8, 1, 1, 0, 'C', 1, 1, 0, 'C'}) { legacy.push_back(static_cast<char>(byte)); }
    std::istringstream legacyInput(legacy, std::ios::binary); std::string legacyLabel;
    auto legacyDiff = diffDNA::read(legacyInput, legacyLabel);
    check(legacyDiff->readCount() == 3 && legacyDiff->uniqueCount() == 2, "Version 1 imports exact repeats into counted subclusters");
    std::ostringstream out(std::ios::binary); delta.write(out, "cluster;size=9;"); delta.write(out, "second");
    std::istringstream in(out.str(), std::ios::binary); std::string label;
    auto decoded = diffDNA::read(in, label); check(label == "cluster;size=9;", "Cluster label survives serialization"); verify(*decoded);
    verify(*diffDNA::read(in, label)); check(label == "second" && !diffDNA::read(in, label), "Concatenated blocks and clean EOF");
    const auto bytes = out.str();
    for (size_t cut = 1; cut < bytes.size() / 2; ++cut) {
        throws([&] { std::istringstream broken(bytes.substr(0, cut), std::ios::binary); diffDNA::read(broken, label); });
    }
    diffDNA exact(anchor); for (size_t i = 0; i < 10000; ++i) { exact.add(anchor); }
    check(exact.payloadBytes() == 0 && exact.readCount() == 10000, "Identical reads require only a counter");
}
void testUnifiedCounts() {
    diffDNA counts;
    check(!counts.hasSequenceData() && !counts.tracksSamples() && counts.readCount() == 0,
        "Default diffDNA retains only the base sample counts");
    counts.addSample(7, 3); counts.addSample(2, 4); counts.addSample(7, 5);
    counts.addSample(-1, 99); counts.addSample(9, 0);
    check(counts.readCount() == 12 && counts.sampleCounts() == diffDNA::SampleCounts({{7, 8}, {2, 4}}),
        "Sparse sample counts coalesce repeated samples without counting unassigned reads");
    check(counts.uniqueCount() == 0 && counts.payloadBytes() == 0 && counts.referenceCount() == 0,
        "Exact count-only mode leaves all nucleotide entries empty");
    size_t visited = 0; counts.forEach([&](auto, auto) { ++visited; });
    check(visited == 0, "Count-only mode does not invent a nucleotide observation");
    throws([&] { counts.reference(); });
    throws([&] { std::ostringstream out; counts.write(out, "no_sequence"); });
    auto saved = counts;
    counts.add("ignored", diffDNA::SampleCounts{{2, 6}, {9, 2}});
    check(counts.readCount() == 20 && saved.readCount() == 12 && saved.sampleCounts()[1].second == 4,
        "A count-only copy remains independent after another sample update");
    counts.merge(saved); counts.merge(counts);
    check(counts.readCount() == 64 && counts.sampleCounts() == diffDNA::SampleCounts({{7, 32}, {2, 28}, {9, 4}}),
        "Counts-only merges, including self-merge, preserve each sample exactly");
    counts.clear(); check(counts.readCount() == 0 && counts.sampleCounts().empty(), "Reset clears base counts");
    counts.addSample(0, std::numeric_limits<uint64_t>::max());
    throws([&] { counts.addSample(1); }); throws([&] { counts.merge(saved); });
    check(counts.readCount() == std::numeric_limits<uint64_t>::max() && counts.sampleCounts().size() == 1,
        "Overflow across samples is rejected before changing the sparse counts");
    throws([&] { saved.add("ignored", diffDNA::SampleCounts{{-1, 1}}); });
    check(saved.readCount() == 12, "Invalid sample provenance does not change counts");

    auto unique = std::make_shared<DNAunique>(readDNA("ACGT", 30, 3), 3);
    unique->incrementSampleCounterBy(1, 5);
    check(!unique->getDiffDNA() && !unique->getDerepCounts().hasSequenceData() &&
        unique->getDerepCounts().readCount() == 6 && unique->totalSum() == 6,
        "Exact DNAunique uses its embedded diffDNA for every count without a variant allocation");
    check(&unique->getDerepMap() == &unique->getDerepCounts().sampleCounts(),
        "The legacy sparse-map API exposes the same storage as diffDNA");
    auto output = unique->cloneForOutput(); unique->incrementSampleCounter(1);
    check(output->totalSum() == 6 && unique->totalSum() == 7, "Exact output snapshots retain independent counts");
    output->transferOccurence(unique);
    check(output->totalSum() == 13 && output->getDerepMap() == read_occ({{3, 2}, {1, 11}}),
        "Legacy occurrence merging delegates to the unified count store");
    output->clearDerepCounts();
    check(output->getDerepCounts().readCount() == 0 && output->getDerepMap().empty() && output->totalSum() == 1,
        "Clearing removes counts while preserving the legacy uncounted-intermediate fallback");
}
void testUnifiedVariantCounts() {
    // Independently accumulate observations across samples and distinct lengths.
    // Merge sources with different immutable references, exercising rebasing.
    const std::vector<std::string> sequences = {"ACGT", "TCGT", "AC", "ACGTNN", ""};
    for (bool tracked : {false, true}) {
        diffDNA left("ACGT", tracked), right("TCGT", tracked);
        std::map<std::string, std::map<int, uint64_t>> expected;
        std::mt19937_64 rng(22091);
        for (size_t i = 0; i < 2000; ++i) {
            const auto& sequence = sequences[rng() % sequences.size()];
            const int sample = rng() % 11; const uint64_t copies = 1 + rng() % 20;
            (i % 2 ? left : right).add(sequence, copies, sample);
            expected[sequence][sample] += copies;
        }
        auto verify = [&](const diffDNA& value, uint64_t multiplier = 1) {
            std::map<int, uint64_t> aggregates, actualSamples;
            std::map<std::string, uint64_t> variants, actualVariants;
            uint64_t total = 0;
            for (const auto& variant : expected) {
                for (const auto& sample : variant.second) {
                    const auto count = sample.second * multiplier;
                    aggregates[sample.first] += count; variants[variant.first] += count; total += count;
                }
            }
            for (const auto& entry : value.sampleCounts()) { actualSamples[entry.first] += entry.second; }
            value.forEach([&](auto sequence, auto count) { actualVariants[std::string(sequence)] += count; });
            check(actualSamples == aggregates && actualVariants == variants && value.readCount() == total,
                "Merged aggregate sample counts and nucleotide counts equal an independent observation tally");
            if (value.tracksSamples()) {
                std::map<std::string, std::map<int, uint64_t>> actual;
                value.forEachWithSamples([&](auto sequence, auto count, const auto& samples) {
                    uint64_t sum = 0;
                    for (const auto& sample : samples) {
                        actual[std::string(sequence)][sample.first] += sample.second / multiplier;
                        sum += sample.second;
                    }
                    check(sum == count, "Each variant's multiplicity is the sum of its sparse sample counts");
                });
                check(actual == expected, "Merge preserves each variant's sample provenance");
            }
        };
        left.merge(right); verify(left);
        auto snapshot = left; left.merge(left); verify(left, 2); verify(snapshot);
        left.clear();
        check(left.hasSequenceData() && left.tracksSamples() == tracked && left.reference() == "ACGT" &&
            left.readCount() == 0 && left.uniqueCount() == 0 && left.sampleCounts().empty(),
            "Reset clears samples and variants together while retaining the mode and reference");
        verify(snapshot);
        // Read the aggregate before mutation to check invalidation of the tracked cache.
        (void)left.sampleCounts(); left.add("ACGT", diffDNA::SampleCounts{{7, 2}, {2, 3}});
        (void)left.sampleCounts(); left.add("TCGT", 4, 7);
        const auto aggregate = left.sampleCounts();
        const std::map<int, uint64_t> actual(aggregate.begin(), aggregate.end());
        check(actual == std::map<int, uint64_t>({{2, 3}, {7, 6}}) && left.readCount() == 9,
            "Multi-sample recording updates totals once and invalidates the aggregate view");
        throws([&] { left.addSample(0); });
        diffDNA anonymous; anonymous.addSample(0);
        throws([&] { left.merge(anonymous); });
        check(left.readCount() == 9, "Reject count-only observations when nucleotide variants must be retained");
        diffDNA countOnly; countOnly.merge(left);
        check(!countOnly.hasSequenceData() && countOnly.readCount() == 9 && countOnly.sampleCounts() == aggregate,
            "A count-only destination can absorb samples without allocating nucleotide data");
        std::ostringstream out; left.write(out, "counted"); std::istringstream in(out.str()); std::string label;
        auto restored = diffDNA::read(in, label);
        check(label == "counted" && restored->readCount() == 9 && restored->referenceCount() == 5 && restored->uniqueCount() == 2 &&
            restored->sampleCounts().empty() && !restored->tracksSamples(),
            "The existing nucleotide-only sidecar preserves multiplicities without inventing sample labels");
        if (tracked) {
            throws([&] { left.merge(*restored); });
            check(left.readCount() == 9, "Tracked merge rejects a sidecar that lacks sample provenance");
        }
    }
}
void testIndexAgainstExhaustive() {
    std::mt19937_64 rng(2026091219);
    // Include length-bin boundaries, homopolymer-only and ambiguous reads.
    std::vector<std::string> templates;
    for (size_t length : {19, 20, 99, 127, 128, 129, 199, 200, 201, 220, 250, 254, 255, 256, 1408, 1499, 1500, 1535, 1536, 1600}) {
        templates.push_back(randomSequence(rng, length));
    }
    templates.push_back(std::string(1500, 'A')); templates.push_back(std::string(1500, 'N'));
    for (const auto& identity : {"100", "99", "98", "97.5", "97", "90", "0.1"}) {
        for (unsigned homopolymer : {2, 4, 6}) {
            CoarseDereplication index(CoarseDereplication::parseIdentity(identity), homopolymer, false);
            std::vector<std::string> anchors;
            std::map<std::string, size_t> exactKeys;
            for (size_t r = 0; r < 700; ++r) {
                std::string sequence = templates[rng() % templates.size()];
                if (r % 5 == 0 && sequence.size() > 30) { sequence.resize(sequence.size() - rng() % 20); }
                for (char& base : sequence) { if (rng() % 100 < 2) { mutate(base); } }
                // A repeated exact key reuses its group. Otherwise cluster output
                // assigns to the nearest anchor within budget (lowest mismatch
                // rate); ties go to the earliest anchor.
                size_t expected = CoarseDereplication::npos, bestErrors = 0, bestLength = 1;
                const auto known = exactKeys.find(sequence);
                if (known != exactKeys.end()) { expected = known->second; }
                else for (size_t id = 0; id < anchors.size(); ++id) {
                    const size_t length = std::min(anchors[id].size(), sequence.size());
                    size_t errors = 0;
                    for (size_t i = 0; i < length; ++i) { errors += anchors[id][i] != sequence[i]; }
                    if (errors > index.mismatchBudget(length)) { continue; }
                    if (expected == CoarseDereplication::npos || errors * bestLength < bestErrors * length) {
                        expected = id; bestErrors = errors; bestLength = length;
                        if (errors == 0) { break; }
                    }
                }
                check(index.find(sequence) == expected, "Seed lookup disagrees with exhaustive nearest-anchor match at " + std::string(identity));
                if (expected == CoarseDereplication::npos) { expected = anchors.size(); anchors.push_back(sequence); }
                exactKeys.emplace(sequence, expected);
                check(add(index, sequence) == expected, "Streaming assignment disagrees with exhaustive clustering");
            }
            check(index.groups().size() == anchors.size(), "Cluster count agrees with exhaustive clustering");
        }
    }
    CoarseDereplication low(CoarseDereplication::parseIdentity("99"), 4, false);
    add(low, std::string(1500, 'A'));
    std::string query(1500, 'A'); for (size_t i = 0; i < 15; ++i) { query[i * 50] = 'C'; }
    check(low.find(query) == 0 && low.stats().fallbackQueries > 0, "Homopolymer fallback retains a boundary match");
    query.back() = 'C'; check(low.find(query) == CoarseDereplication::npos, "Fallback rejects over-budget mismatch");
    check(low.stats().indexedSeeds == 0, "Homopolymer anchor has no search seeds");
}
void testAmpliconSeedPruning() {
    std::mt19937_64 rng(2026091402);
    auto sequence = [&](size_t length) {
        auto result = randomSequence(rng, length);
        // Avoid homopolymer exclusions so this isolates length-bound pruning.
        for (size_t i = 1; i < length; ++i) {
            if (result[i] == result[i - 1]) { mutate(result[i]); }
        }
        return result;
    };
    const auto identity = CoarseDereplication::parseIdentity("97");
    for (size_t length = 200; length <= 255; ++length) {
        CoarseDereplication index(identity, 4, false);
        std::vector<std::string> anchors;
        for (size_t i = 0; i < 32; ++i) {
            anchors.push_back(sequence(length));
            check(add(index, anchors.back()) == i, "Unrelated amplicons seed separate clusters");
        }
        const auto before = index.stats();
        check(index.find(anchors.back()) == anchors.size() - 1, "Find the last amplicon anchor");
        check(index.stats().comparisons - before.comparisons < 4,
            "200-255 bp exact queries prune unrelated anchors instead of scanning the length bin");
        auto boundary = anchors.back();
        for (size_t i = 0; i < index.mismatchBudget(length); ++i) { mutate(boundary[i * 20]); }
        check(index.find(boundary) == anchors.size() - 1, "Seed pruning retains a 97 percent boundary match");
        mutate(boundary[index.mismatchBudget(length) * 20]);
        check(index.find(boundary) == CoarseDereplication::npos, "Reject a query one substitution outside 97 percent");
        check(index.stats().fallbackQueries == 0, "Seed-rich 200-255 bp queries do not require fallback");
    }
    // Bounds must expand in both directions as differently sized anchors arrive.
    for (bool longerFirst : {false, true}) {
        CoarseDereplication index(identity, 4, false);
        const auto shorter = sequence(200), longer = sequence(255);
        add(index, longerFirst ? longer : shorter);
        add(index, longerFirst ? shorter : longer);
        auto longQuery = longer;
        for (size_t i = 0; i < 7; ++i) { mutate(longQuery[i * 20]); }
        check(index.find(longQuery) == (longerFirst ? 0 : 1), "Use the longest actual overlap to bound mismatches");
        auto shortQuery = shorter + sequence(55);
        for (size_t i = 0; i < 6; ++i) { mutate(shortQuery[i * 20]); }
        check(index.find(shortQuery) == (longerFirst ? 1 : 0), "Short anchors remain reachable despite unconstrained query tails");
        check(index.stats().fallbackQueries == 0, "Mixed 200/255 bp anchors can use eligible shared-prefix seeds");
        index.clear();
        add(index, longer);
        check(index.find(longer) == 0 && index.stats().fallbackQueries == 0, "Cleared indexes rebuild populated length bounds");
    }
    // Bounds describe the immutable search key, not the full retained sequence.
    CoarseDereplication capped(identity, 4, true);
    const auto full = sequence(255);
    auto dna = readDNA(full);
    capped.add(std::string_view(full).substr(0, 200), full, dna, nullptr, nullptr, 0, true);
    auto query = full;
    for (size_t i = 200; i < query.size(); ++i) { mutate(query[i]); }
    for (size_t i = 0; i < 6; ++i) { mutate(query[i * 20]); }
    check(capped.find(query) == 0 && capped.stats().fallbackQueries == 0,
        "Search-key bounds retain prefix matches when stored tails differ");
}
void testRepresentativeAndFixedAnchor() {
    std::mt19937_64 rng(55); const auto a = randomSequence(rng, 1500);
    auto b = a, c = a;
    for (size_t i = 0; i < 15; ++i) { mutate(b[i * 20]); mutate(c[i * 20]); }
    for (size_t i = 0; i < 15; ++i) { mutate(c[i * 20 + 1]); }
    CoarseDereplication index(CoarseDereplication::parseIdentity("99"), 4, true);
    add(index, a, 10); check(add(index, b, 40) == 0, "Inclusive 99 percent boundary");
    check(add(index, c, 40) == 1, "Representative replacement cannot chain clusters");
    auto representative = index.groups()[0]->representatives.best(true, false);
    check(representative->getSeqPseudoView() == b, "Use existing quality selection for the representative");
    check(representative->quality_sum_per_base_.empty(), "Coarse mode retains no quality sums");
    check(representative->getDiffDNA()->reference() == a, "Deltas retain the immutable anchor after replacement");
    std::map<std::string, uint64_t> observed;
    representative->getDiffDNA()->forEach([&](std::string_view s, uint64_t n) { observed[std::string(s)] += n; });
    check(observed.size() == 2 && observed[a] == 1 && observed[b] == 1, "Both original observations survive quality replacement");
    check(add(index, a, 10, false) == 0 && index.groups()[0]->representatives.totalSum() == 3,
        "Matching failed reads contribute counts without replacing representative");
    check(representative->getDiffDNA()->readCount() == 2 &&
        index.groups()[0]->representatives[0]->getDiffDNA()->readCount() == 3,
        "A retained output snapshot keeps independent diff counts when processing resumes");
    CoarseDereplication disabled(CoarseDereplication::parseIdentity("99"), 4, false);
    check(add(disabled, a, 25, false) == CoarseDereplication::npos && disabled.groups().empty(), "Failed reads cannot seed clusters");
    add(disabled, a); check(!disabled.groups()[0]->representatives[0]->getDiffDNA(), "No per-read delta allocation by default");
}
OptContainer options(const std::string& base) {
    OptContainer opts;
    opts["-sample_sep"] = "___"; opts["-paired"] = "1";
    opts["-i_qual_offset"] = "33"; opts["-o_qual_offset"] = "33";
    opts["-5PR1cut"] = "0"; opts["-5PR2cut"] = "0";
    opts["-GoldenAxe"] = "0"; opts["-suppressOutput"] = "0"; opts["-threadIO"] = "0";
    opts["-merge_pairs_seed"] = "0"; opts["-count_chimeras"] = "F";
    opts["-o_dereplicate"] = base + ".fq"; opts["-derep_format"] = "fq";
    opts["-derepCoarseClusters"] = "1";
    opts["-derepIdentity"] = "99"; opts["-derepSrchLen"] = "-1"; opts["-derep_reserve"] = "1";
    return opts;
}
std::string fileText(const std::string& path) {
    std::ifstream in(path, std::ios::binary); check(bool(in), "Open " + path);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
void testRefinement() {
    // Sample provenance is opt-in, coalesces exact repeats, and shares snapshots safely.
    diffDNA samples("ACGT", true); samples.add("ACGT", 3, 0); samples.add("TCGT", 2, 1);
    diffDNA snapshot = samples; samples.add("TCGT", 4, 2);
    uint64_t oldCount = 0;
    snapshot.forEachWithSamples([&](auto, uint64_t n, const auto& counts) {
        uint64_t sum = 0; for (auto sample : counts) { sum += sample.second; }
        check(sum == n, "Snapshot sample counts equal variant count"); oldCount += sum;
    });
    check(oldCount == 5 && samples.readCount() == 9, "Sample provenance uses copy on write");
    throws([&] { samples.add("ACGT"); });
    throws([&] { diffDNA("ACGT").forEachWithSamples([](auto, auto, const auto&) {}); });

    std::mt19937_64 rng(2026091301);
    for (size_t length : {200, 1500}) {
        const auto a = randomSequence(rng, length);
        auto b = a, move = a, tie = a;
        const size_t changes = length == 200 ? 8 : 60;
        for (size_t i = 0; i < changes; ++i) { mutate(b[i]); }
        for (size_t i = 0; i < changes / 2 + 1; ++i) { mutate(move[i]); }
        for (size_t i = 0; i < changes / 2; ++i) { mutate(tie[i]); }
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, true, true);
        // move (closer to b) arrives while a is the only anchor; exact-key reuse
        // keeps its later copies there until final reassignment.
        check(add(index, a, 40) == 0, "Initial anchor");
        for (int i = 0; i < 5; ++i) {
            check(index.add(move, move, readDNA(move, 20, i % 2), nullptr, nullptr, i % 2, true) == 0, "Streaming admission to the only anchor; exact reuse afterwards");
        }
        check(add(index, b, 40) == 1, "Separate later anchor");
        add(index, tie, 20);
        auto refined = index.refine(false, -1);
        check(refined.movedReads == 5 && refined.movedVariants == 1, "Move the complete exact subcluster once");
        check(refined.representatives[0]->totalSum() == 2 && refined.representatives[1]->totalSum() == 6, "Reassigned cluster abundances");
        const auto& target = refined.representatives[1]->getDerepMap();
        check(target.size() == 2 && target[0] == read_occ_entry(0, 4) && target[1] == read_occ_entry(1, 2), "Move every sample's counts with the variant");
        check(refined.representatives[1]->getDiffDNA()->reference() == b, "Rebase exported differences onto the destination reference");
        uint64_t tied = 0;
        refined.representatives[0]->getDiffDNA()->forEach([&](auto seq, auto n) { if (seq == tie) { tied += n; } });
        check(tied == 1, "Equal identities keep the original assignment");
        check(index.groups()[0]->representatives[0]->totalSum() == 7, "Refinement does not mutate live counts");
        check(index.refine(false, -1).movedReads == 5, "Repeated final passes do not compound moves");
        add(index, move, 20);
        check(index.refine(false, -1).movedReads == 6 && refined.representatives[1]->totalSum() == 6, "Resume admission after a retained snapshot");
    }
    {
        auto a = randomSequence(rng, 200), b = a, x = a;
        for (size_t i = 150; i < 158; ++i) { mutate(b[i]); }
        for (size_t i = 150; i < 155; ++i) { mutate(x[i]); }
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, false, true);
        auto insertMerged = [&](const auto& sequence, int q) {
            return index.add(sequence, sequence, readDNA(sequence.substr(0, 100), q), nullptr, readDNA(sequence, q), 0, true);
        };
        // x is admitted while a is the only anchor; it is closer to b over the full merged length.
        insertMerged(a, 40); insertMerged(x, 20); insertMerged(b, 40);
        check(index.refine(true, -1).movedReads == 1, "Refinement searches selected merged representatives");
        check(index.refine(true, 100).movedReads == 0, "A search cap excludes differences outside the compared prefix");
        check(index.refine(false, -1).movedReads == 0, "R1 comparison uses the shorter representative length");
    }
    {
        const auto a = randomSequence(rng, 200);
        auto b = a; for (size_t p : {0, 1, 150, 151, 152, 153, 154, 155}) { mutate(b[p]); }
        auto full = b; full[1] = a[1]; // Tied across the logical 50 bp, closer to B across the hidden tail.
        auto shortRead = readDNA(full, 20); shortRead->cutSeqPseudo(50);
        CoarseDereplication index(CoarseDereplication::parseIdentity("97"), 4, true, true, true);
        check(add(index, a, 50) == 0 && add(index, b, 50) == 1, "Separate retained-quality search anchors");
        check(index.add(shortRead->getSeqPseudoView(), shortRead->getSeqPseudoView(), shortRead,
            nullptr, nullptr, 0, true) == 0, "Short logical observation enters the first qualifying anchor");
        auto refined = index.refine(false, -1);
        check(refined.movedReads == 0, "Full HQ tails do not change logical reassignment identity or ties");
        bool found = false;
        refined.representatives[0]->getDiffDNA()->forEachWithQualities([&](auto sequence, auto, const auto*, const auto& q) {
            if (sequence == full) { found = true; check(q.seed->qualities.size() == 200 && q.seed->stats.logical == 50,
                "Reassignment preserves full quality payload and separate logical length"); }
        });
        check(found, "Reassignment exports the original full observation");
    }
    // Independent exhaustive nearest-representative oracle: bin boundaries,
    // ambiguity/homopolymer fallback, caps, quality replacement, and failed reads.
    for (const auto& identity : {"100", "99", "98", "97", "90"}) {
        for (int cap : {-1, 180}) {
            CoarseDereplication index(CoarseDereplication::parseIdentity(identity), 4, true, true);
            std::vector<std::string> templates;
            for (size_t length : {19, 127, 128, 129, 200, 255, 256, 1500}) { templates.push_back(randomSequence(rng, length)); }
            templates.push_back(std::string(200, 'A')); templates.push_back(std::string(200, 'N'));
            auto key = [cap](std::string_view seq) { return seq.substr(0, cap > 0 ? std::min(seq.size(), size_t(cap)) : seq.size()); };
            struct Observation { std::string sequence; size_t original; int sample; };
            std::vector<Observation> observations;
            for (int i = 0; i < 400; ++i) {
                auto sequence = templates[rng() % templates.size()];
                for (char& base : sequence) { if (rng() % 100 < 2) { mutate(base); } }
                const int sample = i % 3; const bool passed = i % 7 != 0;
                auto dna = readDNA(sequence, 10 + rng() % 31, sample, passed);
                const auto id = index.add(key(sequence), sequence, dna, nullptr, nullptr, sample, passed);
                if (id != CoarseDereplication::npos) { observations.push_back({sequence, id, sample}); }
            }
            using Variants = std::map<std::string, uint64_t>;
            std::vector<Variants> expected(index.groups().size());
            std::vector<std::map<int,uint64_t>> expectedSamples(index.groups().size());
            uint64_t moved = 0;
            for (const auto& observation : observations) {
                auto query = key(observation.sequence);
                auto distance = [&](size_t id) {
                    auto rep = key(index.groups()[id]->representatives[0]->getSeqPseudoView());
                    uint64_t length = std::min(query.size(), rep.size()), errors = 0;
                    for (size_t i = 0; i < length; ++i) { errors += query[i] != rep[i]; }
                    return std::pair<uint64_t,uint64_t>(errors, length);
                };
                size_t target = observation.original; auto best = distance(target);
                for (size_t i = 0; i < expected.size(); ++i) {
                    const auto candidate = distance(i);
                    if (candidate.first <= index.mismatchBudget(candidate.second) && candidate.first * best.second < best.first * candidate.second) {
                        target = i; best = candidate;
                    }
                }
                expected[target][observation.sequence]++; expectedSamples[target][observation.sample]++;
                moved += target != observation.original;
            }
            auto result = index.refine(false, cap);
            check(result.movedReads == moved, "Nearest index agrees with independent exhaustive reassignment");
            for (size_t i = 0; i < expected.size(); ++i) {
                const auto& dna = result.representatives[i];
                if (expected[i].empty()) { check(!dna, "Empty output group omitted"); continue; }
                Variants observed;
                dna->getDiffDNA()->forEach([&](auto seq, auto count) { observed[std::string(seq)] += count; });
                check(observed == expected[i], "Every variant reaches its nearest qualifying representative");
                std::map<int,uint64_t> counted;
                for (auto sample : dna->getDerepMap()) { counted[sample.first] += sample.second; }
                check(counted == expectedSamples[i], "Exhaustive sample abundance conservation");
            }
        }
    }
}
void testRefinementOutput(const std::string& base) {
    std::mt19937_64 rng(2026091302); auto a = randomSequence(rng, 200), b = a, x = a;
    for (size_t i = 0; i < 8; ++i) { mutate(b[i]); }
    for (size_t i = 0; i < 5; ++i) { mutate(x[i]); }
    for (bool reassign : {false, true}) {
        for (bool store : {false, true}) {
            const auto stem = base + std::to_string(reassign) + std::to_string(store);
            auto opts = options(stem); opts["-derepIdentity"] = "97";
            opts["-derepReassign"] = reassign ? "1" : "0"; opts["-derepStoreDiffs"] = store ? "1" : "0";
            opts["-derepSubclusterFasta"] = store ? "1" : "0";
            Filters filters(&opts); filters.SampleID = {"A", "B"};
            filters.collectStatistics[0]->BarcodeDetected.assign(2, -1);
            filters.collectStatistics[0]->BarcodeDetectedFail.assign(2, 0);
            Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
            auto insert = [&](const auto& seq, const auto& name, int q, int sample) {
                auto dna = readDNA(seq, q, sample); dna->setHeader(name); derep.addDNA(dna, nullptr);
            };
            // The variant (closer to B) is admitted while A is the only anchor; only
            // explicit reassignment moves it to B at output.
            insert(a, "cluster_A", 40, 0);
            for (int i = 0; i < 3; ++i) { insert(x, "variant", 20, 0); }
            insert(b, "cluster_B", 40, 1);
            const auto report = derep.writeDereplDNA(&filters, "");
            auto first = fileText(stem + ".fq"); auto map = fileText(stem + ".map");
            const auto winner = reassign ? "cluster_B" : "cluster_A";
            check(first.find(std::string(winner) + ";size=4;") != std::string::npos, "Output parent reflects reassigned abundance");
            check((report.find("Final reassignment: 3 reads") != std::string::npos) == reassign, "Reassignment report and default-off behavior");
            if (reassign) { check(map.find("cluster_B;size=4;\t0:3\t1:1") != std::string::npos || map.find("cluster_B;size=4;\t1:1\t0:3") != std::string::npos, "Output sample map follows the moved variant"); }
            if (store) {
                std::ifstream input(stem + ".diff", std::ios::binary); std::string label;
                bool found = false;
                while (auto delta = diffDNA::read(input, label)) {
                    uint64_t copies = 0; delta->forEach([&](auto seq, auto n) { if (seq == x) { copies += n; } });
                    if (copies) { found = true; check(copies == 3 && label.starts_with(winner), "Delta variant belongs to its output parent"); }
                }
                check(found && fileText(stem + ".subclusters.fna").find(";size=3;\n" + x) != std::string::npos, "Exact FASTA variants retain own abundance");
            } else { check(!std::filesystem::exists(stem + ".diff"), "Internal reassignment storage does not enable diff output"); }
            derep.writeDereplDNA(&filters, "");
            check(fileText(stem + ".fq") == first && fileText(stem + ".map") == map, "Repeated refined snapshots are stable");
            insert(x, "variant", 20, 1); derep.writeDereplDNA(&filters, "");
            check(fileText(stem + ".fq").find(std::string(winner) + ";size=5;") != std::string::npos, "Refined output includes resumed input once");
        }
    }
}
void testIntegration(const std::string& base) {
    std::mt19937_64 rng(35); const auto a = randomSequence(rng, 1500); auto b = a; mutate(b[45]);
    for (bool store : {false, true}) {
        const auto stem = base + (store ? "_diff" : "_plain");
        auto opts = options(stem); opts["-derepStoreDiffs"] = store ? "1" : "0"; opts["-derepSubclusterFasta"] = store ? "1" : "0";
        Filters filters(&opts); filters.SampleID = {"sample_a", "sample_b"};
        filters.collectStatistics[0]->BarcodeDetected.assign(2, -1);
        filters.collectStatistics[0]->BarcodeDetectedFail.assign(2, 0);
        Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
        check(derep.usesCoarseDereplication() && !derep.usesPrefixDereplication(), "Startup selects coarse mode");
        derep.addDNA(readDNA(a, 10, 0), nullptr);
        derep.writeDereplDNA(&filters, ""); // Establish anchor before parallel arrivals.
        std::vector<std::thread> workers;
        for (int t = 0; t < 4; ++t) {
            workers.emplace_back([&, t] { for (int r = 0; r < 300; ++r) { derep.addDNA(readDNA(b, 40, t % 2), nullptr); } });
        }
        for (auto& worker : workers) { worker.join(); }
        derep.writeDereplDNA(&filters, "");
        const auto first = fileText(stem + ".fq");
        check(first.find(";size=1201;") != std::string::npos && first.find(b) != std::string::npos, "Count every assigned read and select high quality sequence");
        check(first.find(std::string(1500, 'I')) != std::string::npos, "Representative FASTQ qualities are preserved, not averaged");
        const auto map = fileText(stem + ".map");
        check(map.find("\t0:601") != std::string::npos && map.find("\t1:600") != std::string::npos, "All sample counts survive concurrent batches");
        derep.writeDereplDNA(&filters, ""); check(fileText(stem + ".fq") == first, "Snapshots do not change live clusters");
        if (store) {
            std::ifstream in(stem + ".diff", std::ios::binary); std::string label;
            auto diff = diffDNA::read(in, label);
            check(diff && diff->readCount() == 1201 && diff->reference() == a && diff->uniqueCount() == 2, "Delta sidecar counts exact repeats and retains anchor");
            check(label.find(";size=1201;") != std::string::npos, "Sidecar labels match dereplicate headers");
            uint64_t aCount = 0, bCount = 0;
            diff->forEach([&](std::string_view s, uint64_t n) { if (s == a) { aCount += n; } else if (s == b) { bCount += n; } });
            check(aCount == 1 && bCount == 1200 && !diffDNA::read(in, label), "Sidecar restores all nucleotide observations");
            const auto subclusters = fileText(stem + ".subclusters.fna");
            check(subclusters.find(";size=1;\n" + a) != std::string::npos &&
                subclusters.find(";size=1200;\n" + b) != std::string::npos && subclusters.find(";size=1201;") == std::string::npos,
                "Optional subcluster FASTA exports separate exact counts, not the coarse abundance");
        } else { check(!std::filesystem::exists(stem + ".diff"), "Diff storage is opt-in"); }
        derep.reset(); derep.addDNA(readDNA(b, 30), nullptr); derep.writeDereplDNA(&filters, "");
        check(fileText(stem + ".fq").find(";size=1;") != std::string::npos, "Reset releases clusters and counts");
        derep.finishMap(); check(fileText(stem + ".map").starts_with("#SMPLS"), "Map finalization still works");
    }
    {
        auto invalid = options(base + "_invalid"); invalid["-derepSubclusterFasta"] = "1";
        Filters filters(&invalid);
        throws([&] { Dereplicate noDiffs(&invalid, &filters); });
        invalid["-derepStoreDiffs"] = "1"; invalid["-o_dereplicate"] = base + "_collision.diff";
        throws([&] { Dereplicate collision(&invalid, &filters); });
        invalid["-derepCoarseClusters"] = "yes"; throws([&] { Dereplicate badClusters(&invalid, &filters); });
        invalid["-derepCoarseClusters"] = "1";
        invalid["-derepReassign"] = "yes"; throws([&] { Dereplicate badReassign(&invalid, &filters); });
    }
    // Delta output follows per-run naming and includes groups below min copies.
    {
        const auto stem = base + "_runs";
        auto opts = options(stem); opts["-derepStoreDiffs"] = "1";
        opts["-derepReassign"] = "1"; opts["-derepPerSR"] = "1"; opts["-min_derep_copies"] = "2";
        Filters filters(&opts); filters.SampleID = {"sample"};
        filters.collectStatistics[0]->BarcodeDetected.assign(1, -1);
        filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
        Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
        const auto other = randomSequence(rng, 1500);
        derep.addDNA(readDNA(a), nullptr); derep.addDNA(readDNA(a), nullptr);
        derep.addDNA(readDNA(other), nullptr);
        const auto diffPath = stem + ".runA.diff";
        std::filesystem::create_directory(diffPath);
        throws([&] { derep.writeDereplDNA(&filters, "run A"); });
        std::filesystem::remove(diffPath);
        derep.writeDereplDNA(&filters, "run A");
        check(fileText(stem + ".runA.fq").find(";size=2;") != std::string::npos,
            "Per-run output retains accepted cluster");
        check(fileText(stem + ".fq.rest").find(";size=1;") != std::string::npos,
            "Below-cutoff cluster goes to rest output");
        std::ifstream input(diffPath, std::ios::binary); std::string label; uint64_t count = 0; size_t groups = 0;
        while (auto delta = diffDNA::read(input, label)) { count += delta->readCount(); ++groups; }
        check(count == 3 && groups == 2, "All groups survive output failure and appear in the delta sidecar");
        const auto firstRun = fileText(diffPath);
        derep.reset(); derep.addDNA(readDNA(b), nullptr); derep.writeDereplDNA(&filters, "runB");
        check(fileText(diffPath) == firstRun && std::filesystem::exists(stem + ".runB.diff"),
            "Separate sequencing runs retain separate delta sidecars");
    }
    // File options and CLI precedence, plus truncation of clustering keys only.
    const auto stem = base + "_options";
    { std::ofstream settings(stem + ".txt"); settings << "minSeqLength\t100\nderepIdentity\t98.5\r\nderepStoreDiffs\t1\nderepReassign\t1\r\nderepSeedHomopolymer\t6\nTruncateSequenceLength\t100\n"; }
    auto opts = options(stem); opts.data().erase("-derepIdentity"); opts["-options"] = stem + ".txt";
    Filters filters(&opts); check(opts.get("-derepIdentity") == "98.5", "Option-file identity is trimmed and loaded");
    check(opts.get("-derepReassign") == "1", "Option-file reassignment is trimmed and loaded");
    filters.SampleID = {"sample"}; filters.collectStatistics[0]->BarcodeDetected.assign(1, -1);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
    Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
    auto tail = a; for (size_t i = 100; i < tail.size(); ++i) { mutate(tail[i]); }
    derep.addDNA(readDNA(a), nullptr); derep.addDNA(readDNA(tail), nullptr); derep.writeDereplDNA(&filters, "");
    check(fileText(stem + ".fq").find(";size=2;") != std::string::npos, "Truncation caps coarse matching even when lengths are fixed");
    std::ifstream input(stem + ".diff", std::ios::binary); std::string label;
    auto differences = diffDNA::read(input, label); uint64_t tailCount = 0;
    differences->forEach([&](std::string_view s, uint64_t n) { if (s == tail) { tailCount += n; } });
    check(tailCount == 1, "A search cap does not lose differences in the processed sequence suffix");
    auto explicitOpts = options(base + "_override"); explicitOpts["-options"] = stem + ".txt";
    explicitOpts["-derepReassign"] = "0"; explicitOpts["-derepIdentity"] = "100"; explicitOpts["-derepStoreDiffs"] = "0";
    Filters overrideFilters(&explicitOpts); Dereplicate exact(&explicitOpts, &overrideFilters);
    check(!exact.usesCoarseDereplication() && !exact.usesPrefixDereplication(), "CLI precedence and default fixed-length fast path");
}
using FastqVariants = std::map<std::string, std::pair<uint64_t, std::string>>;
FastqVariants parseVariantFastq(const std::string& text) {
    std::istringstream in(text); std::string header, sequence, plus, quality; FastqVariants records;
    while (std::getline(in, header)) {
        check(bool(std::getline(in, sequence)) && bool(std::getline(in, plus)) && bool(std::getline(in, quality)), "Complete FASTQ record");
        check(header.starts_with("@") && plus == "+" && sequence.size() == quality.size(), "Valid reconstructed FASTQ lengths");
        const auto pos = header.find(";size="); check(pos != std::string::npos && header.find(";size=", pos + 1) == std::string::npos, "One exact-variant size annotation");
        check(!records.contains(sequence), "One FASTQ record per exact variant in this fixture");
        records[sequence] = {std::stoull(header.substr(pos + 6)), quality};
    }
    return records;
}
void testRetainedQualities() {
    std::mt19937_64 rng(2026091401); std::vector<std::string> sequences;
    for (size_t length : {20, 21, 50, 127, 128, 200, 1500}) { sequences.push_back(randomSequence(rng, length)); }
    sequences.push_back(sequences[0] + "NN"); sequences.push_back(sequences[0].substr(0, 15));
    for (bool tracked : {false, true}) {
        diffDNA variants(sequences[0], tracked, true);
        std::map<std::string, std::shared_ptr<DNAunique>> expected;
        std::map<std::string, uint64_t> counts; std::map<int,uint64_t> samples;
        for (size_t i = 0; i < 1200; ++i) {
            const auto& sequence = sequences[i % sequences.size()]; const int sample = i % 4;
            const size_t repetition = i / sequences.size();
            const bool pass = repetition % 7 != 0;
            // Revisit each variant with different eligibility, hidden tails and profiles.
            auto read = readDNA(sequence + std::string(repetition % 3, 'N'), 10 + rng() % 36, sample, pass);
            auto profile = read->getQual();
            for (auto& value : profile) { value = std::clamp<int>(int(value) + int(rng() % 9) - 4, 2, 50); }
            read->setQual(std::move(profile));
            read->cutSeqPseudo(sequence.size());
            if (i % 3) { read->setFwdPrimCut(); }
            if (i % 5) { read->setRevPrimCut(); }
            const auto& fullSequence = read->getSequence();
            auto& winner = expected[fullSequence];
            if (!winner) { winner = std::make_shared<DNAunique>(read, -1); }
            else if (pass && !winner->isGreenQual()) { winner->takeOverDNA(read, nullptr, nullptr); }
            else { winner->matchedDNA(read, nullptr, nullptr, -1, true, false); }
            variants.addQuality(*read, 1, sample); ++counts[fullSequence]; ++samples[sample];
        }
        variants.forEachWithQualities([&](auto sequence, auto count, const auto* occurrence, const auto& quality) {
            const auto key = std::string(sequence); const auto& winner = expected.at(key);
            check(count == counts[key] && quality.seed->qualities == std::vector<qual_score>(winner->getQual().begin(), winner->getQual().begin() + sequence.size()), "Selected qualities agree with uncompressed DNAunique selection");
            check(quality.seed->stats.ambiguous == winner->numNonCanonicalDNA(true) && quality.seed->stats.backing == winner->mem_length(), "Retain scoring statistics from hidden backing tails");
            check(quality.history.non_merged_compares == winner->getWhoIsBetterMergeStats().non_merged_compares, "Comparison history survives quality replacement");
            check(bool(occurrence) == tracked, "Sample provenance remains optional within quality variants");
        });
        std::map<int,uint64_t> actualSamples(variants.sampleCounts().begin(), variants.sampleCounts().end());
        check(actualSamples == samples && variants.readCount() == 1200, "Quality selection counts every observation once");
        std::ostringstream out; variants.writeFastq(out, "parent;size=1200; description", 64);
        const auto restored = parseVariantFastq(out.str());
        for (const auto& entry : expected) {
            const auto& q = entry.second->getQual(); std::string encoded(entry.first.size(), ' ');
            for (size_t i=0; i<encoded.size(); ++i) { encoded[i] = static_cast<char>(std::min(126, int(q[i]) + 64)); }
            check(restored.at(entry.first) == std::make_pair(counts[entry.first], encoded), "FASTQ reconstructs winning qualities at requested offset");
        }
        auto snapshot = variants; auto incoming = readDNA(sequences[0], 50); incoming->setFwdPrimCut(); incoming->setRevPrimCut();
        variants.addQuality(*incoming, 1, 0); std::ostringstream unchanged; snapshot.writeFastq(unchanged, "parent;size=1200; description", 64);
        check(unchanged.str() == out.str() && snapshot.readCount() == 1200, "Retained snapshots preserve previous quality vectors and counts");
        diffDNA merged(sequences[1], tracked, true); merged.merge(snapshot);
        std::ostringstream rebased; merged.writeFastq(rebased, "other");
        check(merged.readCount() == 1200 && parseVariantFastq(rebased.str()).size() == counts.size(), "Rebasing merges carry every variant quality payload");
        merged.merge(merged); check(merged.readCount() == 2400 && merged.uniqueCount() == counts.size(), "Quality-aware self merge coalesces exact variants");
        std::ostringstream binary; snapshot.write(binary, "parent"); std::istringstream input(binary.str()); std::string label;
        auto nucleotideOnly = diffDNA::read(input, label);
        check(nucleotideOnly->readCount() == 1200 && !nucleotideOnly->storesQualities(), "Existing SDMDIFF2 remains nucleotide-only");
        throws([&] { std::ostringstream f; nucleotideOnly->writeFastq(f, "missing"); });
        throws([&] { merged.merge(*nucleotideOnly); });
        merged.clear(); check(merged.storesQualities() && merged.uniqueCount() == 0 && merged.readCount() == 0, "Clearing removes qualities and counts together");
        auto absent = readDNA(sequences[0]); absent->setQual({}); throws([&] { merged.addQuality(*absent, 1, 0); });
        auto paired = readDNA(sequences[0]); paired->setMergeLength(100); throws([&] { merged.addQuality(*paired, 1, 0); });
        throws([&] { merged.add(sequences[0], 1, 0); });
        check(merged.readCount() == 0 && merged.uniqueCount() == 0, "Invalid quality observations do not create variants");
    }
    diffDNA limited(sequences[0], false, true); auto read = readDNA(sequences[0]);
    limited.addQuality(*read, std::numeric_limits<uint64_t>::max(), 0);
    throws([&] { limited.addQuality(*read, 1, 0); });
    check(limited.readCount() == std::numeric_limits<uint64_t>::max(), "Retained-quality count overflow is rejected");
}
void testQualityOutput(const std::string& base) {
    std::mt19937_64 rng(2026091402); auto a = randomSequence(rng, 200), b = a, x = a;
    for (size_t i=0;i<8;++i) { mutate(b[i]); } for(size_t i=0;i<5;++i) { mutate(x[i]); }
    for (bool reassign : {false,true}) { for (bool store : {false,true}) {
        const auto stem=base+std::to_string(reassign)+std::to_string(store);
        auto opts=options(stem);opts["-derepIdentity"]="97";opts["-derepStoreQuals"]="1";
        opts["-derepStoreDiffs"]=store?"1":"0";opts["-derepReassign"]=reassign?"1":"0";opts["-min_derep_copies"]="3";
        Filters filters(&opts);filters.SampleID={"A","B"};filters.collectStatistics[0]->BarcodeDetected.assign(2,-1);filters.collectStatistics[0]->BarcodeDetectedFail.assign(2,0);
        Dereplicate derep(&opts,&filters);derep.BCnamesAdding(&filters);
        // x is 5 substitutions from a and 3 from b. Admitted while only a exists, its exact
        // key stays in a's cluster (exact reuse) until final reassignment moves it to b.
        derep.addDNA(readDNA(a,50,0),nullptr);
        derep.addDNA(readDNA(x,20,0),nullptr);derep.addDNA(readDNA(x,40,1),nullptr);
        derep.addDNA(readDNA(b,50,1),nullptr);derep.addDNA(readDNA(x,60,0,false),nullptr);
        const auto report=derep.writeDereplDNA(&filters,"");const auto before=fileText(stem+".1.hq.fq");const auto observed=parseVariantFastq(before);
        check(observed.size()==3 && observed.at(x)==std::make_pair(uint64_t(3),std::string(200,char(40+33))),"Failed duplicate cannot replace passing selected qualities");
        check(observed.at(a).first==1 && observed.at(b).first==1,"FASTQ includes groups below the parent minimum-copy cutoff");
        check((report.find("Final reassignment: 3 reads")!=std::string::npos)==reassign,"Quality retention and reassignment agree on moved abundance");
        check(std::filesystem::exists(stem+".diff")==store,"Quality storage does not implicitly write nucleotide sidecars");
        check(!std::filesystem::exists(stem + ".subclusters.fq"), "Retained variants use only the normal HQ FASTQ path");
        derep.writeDereplDNA(&filters,"");check(fileText(stem+".1.hq.fq")==before,"Repeated quality exports are stable");
        derep.addDNA(readDNA(x,45,0),nullptr);derep.writeDereplDNA(&filters,"");
        check(parseVariantFastq(fileText(stem+".1.hq.fq")).at(x)==std::make_pair(uint64_t(4),std::string(200,char(45+33))),"Resumed input updates selected quality and abundance once");
    }}
    auto opts=options(base+"_flags");opts["-derepStoreQuals"]="yes";Filters filters(&opts);
    throws([&]{Dereplicate invalid(&opts,&filters);});opts["-derepStoreQuals"]="1";filters.setPaired()=2;
    { Dereplicate paired(&opts,&filters); throws([&]{paired.addDNA(readDNA(a),nullptr);}); } filters.setPaired()=1;
    opts["-i_fna"]="without_qualities.fa";throws([&]{Dereplicate noQual(&opts,&filters);});opts["-i_fna"]="";
    opts["-derepIdentity"]="100";opts["-derepPerSR"]="1";Dereplicate exact(&opts,&filters);exact.BCnamesAdding(&filters);
    check(exact.usesCoarseDereplication(),"Quality retention selects streaming mode at 100 percent");
    throws([&]{exact.addDNA(readDNA(a),readDNA(b));});
    exact.addDNA(readDNA(a),nullptr);exact.writeDereplDNA(&filters,"run1");
    check(parseVariantFastq(fileText(base+"_flags.1.hq.fq")).at(a).first==1,"Per-run retained variants use the standard cumulative HQ filename");
    exact.reset(); exact.addDNA(readDNA(b, 40), nullptr); exact.writeDereplDNA(&filters, "run2");
    const auto cumulative = parseVariantFastq(fileText(base + "_flags.1.hq.fq"));
    check(cumulative.size() == 2 && cumulative.at(a).first == 1 && cumulative.at(b).first == 1,
        "Per-run retained HQ output appends both runs using the normal filename");
    check(std::filesystem::exists(base + "_flags.run1.fq") && std::filesystem::exists(base + "_flags.run2.fq") &&
        !std::filesystem::exists(base + "_flags.run1.1.hq.fq") &&
        !std::filesystem::exists(base + "_flags.run1.subclusters.fq"),
        "Only main per-run outputs get run suffixes; no duplicate variant FASTQ is emitted");
    const auto config=base+"_options.txt";{std::ofstream f(config);f<<"derepStoreQuals\t1\n";}
    auto loaded=options(base+"_loaded");loaded["-options"]=config;Filters fromFile(&loaded);check(loaded.get("-derepStoreQuals")=="1","Options file enables selected qualities");
    auto overridden=options(base+"_override");overridden["-options"]=config;overridden["-derepStoreQuals"]="0";Filters fromCli(&overridden);check(overridden.get("-derepStoreQuals")=="0","CLI disables options-file quality retention");
}
using PairKey = std::pair<std::string, std::string>;
using PairCounts = std::map<PairKey, uint64_t>;
struct PairFastqRecord { PairKey sequence, quality; uint64_t count; std::string name; };
std::vector<PairFastqRecord> pairedFastq(const std::string& first, const std::string& second) {
    std::istringstream a(first), b(second); std::vector<PairFastqRecord> result;
    std::string h1,h2,s1,s2,p1,p2,q1,q2;
    while (std::getline(a,h1)) {
        check(bool(std::getline(b,h2)) && bool(std::getline(a,s1)) && bool(std::getline(b,s2)) &&
            bool(std::getline(a,p1)) && bool(std::getline(b,p2)) && bool(std::getline(a,q1)) && bool(std::getline(b,q2)), "Complete paired FASTQ records");
        const auto pos=h1.rfind(";size=");
        check(h1==h2 && h1.starts_with("@") && pos!=std::string::npos && h1.find(";size=")==pos &&
            p1=="+" && p2=="+" && s1.size()==q1.size() && s2.size()==q2.size(), "Paired FASTQ headers, qualities, and abundances stay synchronized");
        result.push_back({{s1,s2},{q1,q2},std::stoull(h1.substr(pos+6)),h1});
    }
    check(!std::getline(b,h2), "No unmatched R2 records"); return result;
}
void testCoarsePartitions(const std::string& base) {
    constexpr size_t threads = 8, variants = 128, perThread = 1024;
    std::mt19937_64 rng(2026091403);
    const auto r1 = randomSequence(rng, 200);
    std::vector<std::string> firsts, mates;
    for (size_t i = 0; i < variants; ++i) { firsts.push_back(randomSequence(rng, 200)); mates.push_back(randomSequence(rng, 250)); }
    auto opts = options(base);
    opts["-threads"] = std::to_string(threads); opts["-derepIdentity"] = "97";
    opts["-derepStoreQuals"] = "1"; opts["-derepStoreDiffs"] = "1";
    // Partitions exist for exact-parent output; explicit cluster output uses the global index.
    opts["-derepCoarseClusters"] = "0";
    Filters filters(&opts); filters.setPaired() = 2; filters.SampleID = {"A", "B", "C"};
    filters.collectStatistics[0]->BarcodeDetected.assign(3, -1);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(3, 0);
    Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
    check(derep.coarsePartitionCount() == 64, "Multiple workers select independent coarse indexes by default");
    std::atomic<bool> start{false};
    std::vector<std::thread> workers;
    std::vector<std::exception_ptr> errors(threads);
    for (size_t t = 0; t < threads; ++t) {
        workers.emplace_back([&, t] {
            while (!start.load()) { std::this_thread::yield(); }
            try {
                for (size_t i = 0; i < perThread; ++i) {
                    auto a = readDNA(firsts[i % variants], 10 + t * 5, t % 3);
                    auto b = readDNA(mates[i % variants], 11 + t * 5, t % 3);
                    a->setHeader("fragment_" + std::to_string(t) + "_" + std::to_string(i));
                    derep.addDNA(a, b);
                }
            } catch (...) { errors[t] = std::current_exception(); }
        });
    }
    start.store(true);
    for (int i = 0; i < 4; ++i) { derep.writeDereplDNA(&filters, ""); std::this_thread::yield(); }
    for (auto& worker : workers) { worker.join(); }
    for (const auto& error : errors) { check(!error, "Partition snapshots let every producer finish"); }
    const auto report = derep.writeDereplDNA(&filters, "");
    check(report.find("Coarse work pool:") != std::string::npos && report.find("Coarse exact reuse:") != std::string::npos,
        "Partitioned admission reports work sharing and exact reuse");
    const auto active = report.find("; active ");
    check(active != std::string::npos && std::stoul(report.substr(active + 9)) >= 32,
        "Different exact R1 keys spread over many coarse partitions");
    const auto first = fileText(base + ".1.hq.fq"), second = fileText(base + ".2.hq.fq");
    const auto records = pairedFastq(first, second);
    PairCounts counts;
    check(records.size() == variants, "Exact pairs occur once across partitions");
    for (const auto& record : records) {
        counts[record.sequence] += record.count;
        check(record.quality == PairKey(std::string(200, char(33 + 45)), std::string(250, char(33 + 46))),
            "Independent indexes retain both qualities from the best fragment");
    }
    for (size_t i = 0; i < variants; ++i) {
        check(counts.at({firsts[i], mates[i]}) == threads * perThread / variants, "All duplicate fragments reach their owning partition");
    }
    uint64_t binaryCount = 0;
    std::ifstream binary(base + ".diff", std::ios::binary); std::string label;
    while (auto differences = diffDNA::read(binary, label)) {
        check(differences->storesPairs(), "Every partition writes linked paired differences");
        binaryCount += differences->readCount();
    }
    check(binaryCount == threads * perThread, "Partitioned binary output counts every fragment once");
    std::vector<uint64_t> sampleCounts(3);
    std::istringstream map(fileText(base + ".map")); std::string line;
    while (std::getline(map, line)) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream fields(line); std::string field;
        std::getline(fields, field, '\t');
        while (std::getline(fields, field, '\t')) {
            const auto colon = field.find(':');
            if (colon != std::string::npos) { sampleCounts.at(std::stoul(field.substr(0, colon))) += std::stoull(field.substr(colon + 1)); }
        }
    }
    check(sampleCounts == std::vector<uint64_t>({3 * perThread, 3 * perThread, 2 * perThread}),
        "Partitioned sample maps conserve each sample's fragment count");
    derep.writeDereplDNA(&filters, "");
    check(fileText(base + ".1.hq.fq") == first && fileText(base + ".2.hq.fq") == second,
        "Repeated partition snapshots do not alter live data");
    // A failed exact pair uses the same owner and cannot replace its passing witness.
    derep.addDNA(readDNA(firsts[0], 60, 0, false), readDNA(mates[0], 60, 0, false));
    derep.writeDereplDNA(&filters, "");
    uint64_t resumed = 0;
    for (const auto& record : pairedFastq(fileText(base + ".1.hq.fq"), fileText(base + ".2.hq.fq"))) {
        resumed += record.count;
        check(record.quality == PairKey(std::string(200, char(33 + 45)), std::string(250, char(33 + 46))),
            "Failed pairs cannot replace passing quality witnesses after a snapshot");
    }
    check(resumed == threads * perThread + 1, "Resuming partitioned admission adds a fragment exactly once");
    derep.reset();
    derep.addDNA(readDNA(firsts[0]), readDNA(mates[0]));
    derep.writeDereplDNA(&filters, "");
    auto reset = pairedFastq(fileText(base + ".1.hq.fq"), fileText(base + ".2.hq.fq"));
    check(reset.size() == 1 && reset[0].count == 1, "Reset clears every coarse partition");
    derep.finishMap();
    derep.addDNA(readDNA(firsts[1]), readDNA(mates[1]));
    derep.writeDereplDNA(&filters, "");
    reset = pairedFastq(fileText(base + ".1.hq.fq"), fileText(base + ".2.hq.fq"));
    check(reset.size() == 1 && reset[0].sequence.second == mates[1] && reset[0].count == 1,
        "Map finalization clears every coarse partition before reuse");

    // Even mutually compatible coarse groups remain independent by default.
    derep.reset();
    for (size_t i = 0; i < variants; ++i) {
        auto mate = mates[0]; mutate(mate[i]);
        auto nearFirst = r1; mutate(nearFirst[i]);
        auto dna = readDNA(nearFirst); dna->setHeader("near_" + std::to_string(i));
        derep.addDNA(dna, readDNA(mate));
    }
    derep.writeDereplDNA(&filters, "");
    size_t independentGroups = 0; uint64_t independentCount = 0;
    std::ifstream separate(base + ".diff", std::ios::binary);
    while (auto differences = diffDNA::read(separate, label)) { ++independentGroups; independentCount += differences->readCount(); }
    check(independentGroups >= 32 && independentCount == variants,
        "Default partition output performs no consolidation of compatible coarse groups");
    check(pairedFastq(fileText(base + ".1.hq.fq"), fileText(base + ".2.hq.fq")).size() == variants,
        "Independent coarse groups retain all exact variants");

    auto config = options(base + "_mode"); config["-threads"] = "8";
    config["-derepStoreQuals"] = "1"; config["-derepIdentity"] = "97";
    Filters single(&config);
    { Dereplicate clusters(&config, &single); check(clusters.coarsePartitionCount() == 1, "Explicit coarse-cluster output retains the global index: partitions cannot form clusters"); }
    config["-derepCoarseClusters"] = "0";
    { Dereplicate parallel(&config, &single); check(parallel.coarsePartitionCount() == 64, "Exact-parent output with multiple workers uses independent partitions"); }
    config["-derepGlobal"] = "1";
    { Dereplicate global(&config, &single); check(global.coarsePartitionCount() == 1, "Explicit global grouping retains one index"); }
    config["-derepGlobal"] = "0"; config["-derepReassign"] = "1";
    { Dereplicate global(&config, &single); check(global.coarsePartitionCount() == 1, "Explicit reassignment retains global representative search"); }
    config["-derepReassign"] = "0"; config["-derepIdentity"] = "100";
    { Dereplicate exact(&config, &single); check(exact.coarsePartitionCount() == 1, "100 percent prefix matching stays global"); }
    config["-derepIdentity"] = "97"; config["-threads"] = "1";
    { Dereplicate serial(&config, &single); check(serial.coarsePartitionCount() == 1, "One worker retains the original global order"); }
    config["-derepGlobal"] = "yes";
    throws([&] { Dereplicate invalid(&config, &single); });
    const auto settings = base + "_global.txt";
    { std::ofstream file(settings); file << "derepGlobal\t1\r\n"; }
    auto fromFile = options(base + "_file"); fromFile["-options"] = settings; fromFile["-threads"] = "8";
    Filters loaded(&fromFile);
    { Dereplicate global(&fromFile, &loaded); check(global.coarsePartitionCount() == 1, "Option file enables global coarse grouping"); }
    auto overridden = options(base + "_override"); overridden["-options"] = settings;
    overridden["-threads"] = "8"; overridden["-derepGlobal"] = "0"; overridden["-derepCoarseClusters"] = "0";
    Filters explicitParallel(&overridden);
    { Dereplicate parallel(&overridden, &explicitParallel); check(parallel.coarsePartitionCount() == 64, "CLI can override the option file to select parallel partitions"); }
}
void testPairedVariants() {
    std::mt19937_64 rng(2026091501);
    const auto r1=randomSequence(rng,200), r2=randomSequence(rng,250);
    std::vector<PairKey> variants{{r1,r2},{r1,r2},{r1,r2},{r1.substr(0,128),r2.substr(0,127)},
        {r1+"NN",r2+"AA"},{r1.substr(0,20),r2.substr(0,21)},{randomSequence(rng,1500),randomSequence(rng,1450)}};
    mutate(variants[1].second[33]); mutate(variants[2].first[71]);
    variants.emplace_back(r1.substr(0,199), r1.back()+r2); // Same concatenation as reference, different mate boundary.
    auto object=[&](bool tracked,bool qualities,const std::string& first,const std::string& second) {
        return diffDNA(std::make_shared<const std::string>(first),tracked,qualities,std::make_shared<const std::string>(second));
    };
    for (bool tracked : {false,true}) {
        auto delta=object(tracked,true,r1,r2); PairCounts counts;
        std::map<PairKey,std::shared_ptr<DNAunique>> expected;
        for (size_t i=0;i<800;++i) {
            const auto& pair=variants[i%variants.size()]; const auto repetition=i/variants.size();
            auto a=readDNA(pair.first+std::string(repetition%3,'N'),10+rng()%35,i%3,true);
            auto b=readDNA(pair.second+std::string(repetition%2,'N'),10+rng()%35,i%3,repetition%9!=0);
            a->cutSeqPseudo(pair.first.size()); b->cutSeqPseudo(pair.second.size());
            if(i%5) {a->setFwdPrimCut();} if(i%7) {b->setRevPrimCut();}
            const PairKey fullPair{a->getSequence(),b->getSequence()};
            auto& winner=expected[fullPair];
            if (!winner) { winner=std::make_shared<DNAunique>(a,-1); winner->takeOverDNA(a,b,nullptr); }
            else if(a->isGreenQual() && b->isGreenQual()) {
                if(!winner->isGreenQual() || !winner->getPair()->isGreenQual()) {winner->takeOverDNA(a,b,nullptr);}
                else {winner->matchedCoarseDNA(a,b,nullptr,-1);}
            }
            delta.addPairQuality(*a,*b,1,i%3); ++counts[fullPair];
        }
        PairCounts restored;
        delta.forEachPair([&](auto first,auto second,auto count,const auto* samples,const auto* quality){
            const PairKey key{first,second}; restored[key]+=count; const auto& winner=expected.at(key);
            const auto& q1=winner->getQual(); const auto& q2=winner->getPair()->getQual();
            check(quality && quality->seed && quality->seed->mate &&
                quality->seed->qualities==std::vector<qual_score>(q1.begin(),q1.begin()+first.size()) &&
                quality->seed->mate->qualities==std::vector<qual_score>(q2.begin(),q2.begin()+second.size()), "Both qualities match one uncompressed selected fragment");
            check(quality->seed->stats.ambiguous==winner->numNonCanonicalDNA(true) &&
                quality->seed->mate->stats.ambiguous==winner->getPair()->numNonCanonicalDNA(true), "Preserve both mates' scoring summaries");
            check(quality->history.non_merged_compares==winner->getWhoIsBetterMergeStats().non_merged_compares, "Paired comparison history agrees with full DNA selection");
            check(bool(samples)==tracked,"One optional shared sample vector per pair");
        });
        check(restored==counts && delta.readCount()==800 && delta.uniqueCount()==counts.size(),"Exact pair identity includes both mate boundaries and counts fragments once");
        uint64_t sampleTotal=0; for(const auto& s:delta.sampleCounts()){sampleTotal+=s.second;} check(sampleTotal==800,"Shared pair sample totals");
        std::ostringstream f1,f2;delta.writeFastq(f1,f2,"parent;size=800; description",64);
        check(pairedFastq(f1.str(),f2.str()).size()==counts.size(),"Export every paired variant once");
        auto snapshot=delta;auto better1=readDNA(r1,60),better2=readDNA(r2,60);better1->setFwdPrimCut();
        delta.addPairQuality(*better1,*better2,1,0);
        std::ostringstream old1,old2;snapshot.writeFastq(old1,old2,"parent;size=800; description",64);
        check(f1.str()==old1.str() && f2.str()==old2.str(),"Paired snapshots preserve both old quality vectors");
        auto rebased=object(tracked,true,variants[2].first,variants[1].second);rebased.merge(snapshot);rebased.merge(rebased);
        check(rebased.readCount()==1600 && rebased.uniqueCount()==counts.size(),"Rebase both mates and self-merge counts once");
        std::ostringstream binary;snapshot.write(binary,"pair");check(binary.str().starts_with("SDMDIFF3"),"Pairs use explicit linked binary format");
        std::istringstream input(binary.str());std::string label;auto loaded=diffDNA::read(input,label);PairCounts imported;
        loaded->forEachPair([&](auto first,auto second,auto n,const auto*,const auto* q){imported[PairKey(first,second)]+=n;check(q==nullptr,"Binary is nucleotide-only");});
        check(imported==counts && loaded->storesPairs() && !loaded->storesQualities(),"Paired binary reconstructs exact linked variants");
        for(size_t length : {size_t(0),size_t(9),binary.str().size()-1}) {
            if(!length)continue;std::istringstream truncated(binary.str().substr(0,length));throws([&]{diffDNA::read(truncated,label);});
        }
        throws([&]{snapshot.forEach([](auto,auto){});});
        throws([&]{std::ostringstream out;snapshot.writeFastq(out,"single");});
        throws([&]{delta.merge(*loaded);});
        auto single=diffDNA(r1,false,true);single.addQuality(*better1,1,0);throws([&]{delta.merge(single);});
        delta.clear();check(delta.storesPairs() && delta.storesQualities() && delta.readCount()==0,"Clear keeps pair topology without stale variants");
        throws([&]{delta.addQuality(*better1,1,0);});
        better2->setQual({});throws([&]{delta.addPairQuality(*better1,*better2,1,0);});
        check(delta.readCount()==0 && delta.uniqueCount()==0,"Missing mate qualities do not partially count a pair");
    }
}
void testPairedCandidateOracle() {
    std::mt19937_64 rng(2026091502);
    std::vector<PairKey> templates;
    for(size_t n:{20,127,128,200,255,256}) {templates.emplace_back(randomSequence(rng,n),randomSequence(rng,n+31));}
    templates.emplace_back(std::string(200,'A'),std::string(200,'N'));
    for (const auto* identity : {"97","99","100"}) {
        CoarseDereplication index(CoarseDereplication::parseIdentity(identity),4,false,false,true);
        std::vector<PairKey> anchors;
        std::map<std::string, size_t> exactR1;
        for(size_t i=0;i<400;++i) {
            auto pair=templates[rng()%templates.size()];
            for(size_t n=0;n<rng()%9;++n){mutate(pair.first[rng()%pair.first.size()]);}
            for(size_t n=0;n<rng()%9;++n){mutate(pair.second[rng()%pair.second.size()]);}
            // Each mate must pass its own budget; the nearest anchor by combined
            // mismatch rate wins (cluster output), ties go to the earliest anchor.
            auto distance=[&](const std::string& x,const std::string& y,uint64_t& errors,uint64_t& overlap){size_t length=std::min(x.size(),y.size());uint64_t e=0;for(size_t p=0;p<length;++p)e+=x[p]!=y[p];if(e>index.mismatchBudget(length))return false;errors+=e;overlap+=length;return true;};
            const auto exact = exactR1.find(pair.first);
            size_t expected = exact == exactR1.end() ? anchors.size() : exact->second;
            if (exact == exactR1.end()) {
                uint64_t bestErrors=0,bestOverlap=1;
                for(size_t a=0;a<anchors.size();++a){uint64_t e=0,o=0;if(!distance(pair.first,anchors[a].first,e,o)||!distance(pair.second,anchors[a].second,e,o))continue;
                    if(expected==anchors.size()||e*bestOverlap<bestErrors*o){expected=a;bestErrors=e;bestOverlap=o;if(e==0)break;}}
            }
            exactR1.emplace(pair.first, expected);
            if(expected==anchors.size())anchors.push_back(pair);
            auto a=readDNA(pair.first),b=readDNA(pair.second);
            check(index.add(pair.first,pair.first,a,b,nullptr,0,true)==expected,"Exact R1 reuse followed by coarse two-mate search agrees with the oracle");
        }
        check(index.groups().size()==anchors.size(),"Exact and approximate paired cluster counts agree with oracle");
        throws([&]{add(index,templates[0].first);});
    }
}
void testPairedOutput(const std::string& base) {
    std::mt19937_64 rng(2026091503);const auto first=randomSequence(rng,200);auto a=randomSequence(rng,200),b=a,x=a;
    for(size_t i=0;i<8;++i)mutate(b[i]);for(size_t i=0;i<5;++i)mutate(x[i]);
    auto firstB = first, firstX = first;
    for (size_t i=0;i<8;++i) mutate(firstB[i]);
    for (size_t i=0;i<5;++i) mutate(firstX[i]);
    for(bool reassign:{false,true})for(bool binary:{false,true}) {
        const auto stem=base+std::to_string(reassign)+std::to_string(binary);
        auto opts=options(stem);opts["-derepStoreQuals"]="1";opts["-derepIdentity"]="97";opts["-derepStoreDiffs"]=binary?"1":"0";
        opts["-derepSubclusterFasta"]=binary?"1":"0";opts["-derepReassign"]=reassign?"1":"0";opts["-min_derep_copies"]="3";
        Filters filters(&opts);filters.setPaired()=2;filters.SampleID={"A","B"};filters.collectStatistics[0]->BarcodeDetected.assign(2,-1);filters.collectStatistics[0]->BarcodeDetectedFail.assign(2,0);
        Dereplicate derep(&opts,&filters);derep.BCnamesAdding(&filters);
        auto insert=[&](const std::string& second,int q1,int q2,int sample,bool pass=true){auto r1=readDNA(second==b?firstB:second==x?firstX:first,q1,sample),r2=readDNA(second,q2,sample,pass);r1->setHeader(second==b?"B":"A");derep.addDNA(r1,r2);};
        // x joins A's cluster while it is the only anchor (exact reuse keeps it there);
        // explicit reassignment later moves it to the closer B.
        insert(a,50,50,0);insert(x,20,25,0);insert(x,40,45,1);insert(b,50,50,1);insert(x,60,60,0,false);
        const auto report=derep.writeDereplDNA(&filters,"");const auto f1=fileText(stem+".1.hq.fq"),f2=fileText(stem+".2.hq.fq");
        const auto records=pairedFastq(f1,f2);check(records.size()==3,"Three full paired variants remain distinct");
        uint64_t counted=0;for(const auto& r:records){counted+=r.count;if(r.sequence.second==x){check(r.count==3&&r.quality==PairKey(std::string(200,char(73)),std::string(200,char(78))),"Failed R2 cannot replace either selected quality vector");check(r.name.starts_with(reassign?"@B.sub":"@A.sub"),"Explicit paired reassignment selects the closer representative");}}
        check(counted==5 && ((report.find("Final reassignment: 3 reads")!=std::string::npos)==reassign),"Fragment counts and paired final reassignment agree");
        check(std::filesystem::exists(stem+".diff")==binary,"Paired binary export stays optional");
        if(binary){check(fileText(stem+".diff").starts_with("SDMDIFF3"),"Linked pair sidecar selected");check(std::filesystem::exists(stem+".subclusters.1.fna")&&std::filesystem::exists(stem+".subclusters.2.fna"),"Both mate FASTA sidecars exist");}
        derep.writeDereplDNA(&filters,"");check(fileText(stem+".1.hq.fq")==f1&&fileText(stem+".2.hq.fq")==f2,"Repeated paired output is stable");
        insert(x,50,55,1);derep.writeDereplDNA(&filters,"");uint64_t after=0;for(const auto&r:pairedFastq(fileText(stem+".1.hq.fq"),fileText(stem+".2.hq.fq")))after+=r.count;check(after==6,"Resumed pair increments one fragment");
    }
    // Exact R1 reuse takes precedence even when R2 differs within its cap.
    CoarseDereplication capped(CoarseDereplication::parseIdentity("100"),4,false,false,true);
    auto shortRead=readDNA(first.substr(0,40));auto secondRead=readDNA(a);auto withinCap=a,beyondCap=a;mutate(withinCap[45]);mutate(beyondCap[70]);
    check(capped.add(shortRead->getSeqPseudoView(),shortRead->getSeqPseudoView(),shortRead,secondRead,nullptr,0,true,50)==0,"First pair anchors both search prefixes");
    check(capped.add(shortRead->getSeqPseudoView(),shortRead->getSeqPseudoView(),shortRead,readDNA(withinCap),nullptr,0,true,50)==0,"Exact R1 remains one dereplicate despite an R2 mismatch within the cap");
    check(capped.add(shortRead->getSeqPseudoView(),shortRead->getSeqPseudoView(),shortRead,readDNA(beyondCap),nullptr,0,true,50)==0,"R2 differences beyond search cap stay exact variants in the same parent");
    auto trimmed=readDNA(a);trimmed->cutSeqPseudo(0);
    check(capped.add(shortRead->getSeqPseudoView(),shortRead->getSeqPseudoView(),shortRead,trimmed,nullptr,0,true,50)==0,"A logically hidden mate retains its full HQ sequence and fragment count");
    auto removed=readDNA(a);removed->cutSeq(0,-1);
    check(capped.add(shortRead->getSeqPseudoView(),shortRead->getSeqPseudoView(),shortRead,removed,nullptr,0,true,50)==CoarseDereplication::npos,"A physically removed mate cannot supply an HQ pair");
    auto opts=options(base+"_caps");opts["-derepStoreQuals"]="1";opts["-derepIdentity"]="100";opts["-derepSrchLen"]="50";opts["-derepPerSR"]="1";
    Filters filters(&opts);filters.setPaired()=2;filters.SampleID={"S"};Dereplicate derep(&opts,&filters);derep.BCnamesAdding(&filters);
    auto short1=first.substr(0,40),different2=a;mutate(different2[45]);
    derep.addDNA(readDNA(short1),readDNA(a));derep.addDNA(readDNA(short1),readDNA(different2));
    derep.writeDereplDNA(&filters,"run1");
    check(pairedFastq(fileText(base+"_caps.1.hq.fq"),fileText(base+"_caps.2.hq.fq")).size()==2,"Apply search cap to R2 independently of a shorter R1");
    throws([&]{derep.addDNA(readDNA(first),nullptr);});
    check(derep.addDNA(readDNA(first),trimmed),"Logical trimming does not discard a retained full mate");
    check(!derep.addDNA(readDNA(first),removed),"Do not buffer a fragment whose mate was physically removed");
    for(const auto* flag:{"-merge_pairs_derep","-merge_pairs_filter","-merge_pairs_demulti"}){opts[flag]="1";throws([&]{Dereplicate merged(&opts,&filters);});opts[flag]="0";}
}
}
#ifndef SDM_COARSE_TEST_NO_MAIN
// Global-index commits from snapshot probes must equal serial admission: the
// same anchors in the same order, the same exact-key assignments, counts and
// retained variants, including when a probe's snapshot is older than earlier
// commits (anchors created in between are found by the delta search).
void testSnapshotCommits() {
    struct Config { uint64_t identity; bool diffs, qualities, exactOutput, paired; };
    const Config configs[] = {
        {97000000, false, true, true, true},    // retained pairs, exact R1 output
        {97000000, true, false, false, false},  // coarse clusters with deltas
        {99000000, false, false, true, false},  // storage-only exact keys
        {100000000, false, true, true, true},   // 100% with retained pairs
    };
    for (const auto& config : configs) {
        std::mt19937_64 rng(2026092401 + config.identity);
        std::vector<std::string> templates;
        for (int t = 0; t < 30; ++t) { templates.push_back(randomSequence(rng, 200)); }
        struct Fragment { std::string r1, r2; int quality, sample; bool pass; };
        std::vector<Fragment> fragments;
        for (int i = 0; i < 1500; ++i) {
            Fragment fragment{templates[rng() % templates.size()], templates[rng() % 7], 20 + int(rng() % 20),
                int(rng() % 4), rng() % 10 != 0};
            for (int m = int(rng() % 6); m > 0; --m) { mutate(fragment.r1[rng() % fragment.r1.size()]); }
            for (int m = int(rng() % 3); m > 0; --m) { mutate(fragment.r2[rng() % fragment.r2.size()]); }
            fragments.push_back(std::move(fragment));
        }
        auto makeIndex = [&] {
            return CoarseDereplication(config.identity, 4, config.diffs, false, config.qualities, config.exactOutput);
        };
        auto reads = [&](const Fragment& fragment) {
            auto r1 = readDNA(fragment.r1, fragment.quality, fragment.sample, fragment.pass);
            auto r2 = config.paired ? readDNA(fragment.r2, fragment.quality, fragment.sample, fragment.pass) : nullptr;
            return std::make_pair(r1, r2);
        };
        auto serial = makeIndex();
        std::vector<size_t> serialIds;
        for (const auto& fragment : fragments) {
            auto [r1, r2] = reads(fragment);
            serialIds.push_back(serial.add(fragment.r1, fragment.r1, r1, r2, nullptr, fragment.sample, fragment.pass));
        }
        auto snapshot = makeIndex();
        std::vector<size_t> snapshotIds;
        CoarseDereplication::SearchScratch scratch;
        const size_t batch = 16;
        // Probe two batches before committing the first, so the second batch's
        // probes miss anchors that the first batch creates.
        std::vector<std::pair<std::shared_ptr<DNA>, std::shared_ptr<DNA>>> pendingReads;
        std::vector<CoarseDereplication::Probe> probes;
        auto probeRange = [&](size_t begin, size_t end) {
            for (size_t i = begin; i < end; ++i) {
                auto pair = reads(fragments[i]);
                probes.push_back(snapshot.probe(fragments[i].r1, pair.first, pair.second, fragments[i].pass, -1, scratch));
                pendingReads.push_back(std::move(pair));
            }
        };
        std::vector<CoarseDereplication::Update> updates;
        size_t committed = 0;
        auto commitRange = [&](size_t end) {
            updates.clear();
            for (size_t i = committed; i < end; ++i) {
                const auto& fragment = fragments[i];
                snapshotIds.push_back(snapshot.add(fragment.r1, fragment.r1, pendingReads[i].first, pendingReads[i].second,
                    nullptr, fragment.sample, fragment.pass, -1, true, false, 0, &probes[i], &updates));
            }
            for (const auto& update : updates) { snapshot.apply(update); }
            committed = end;
        };
        size_t probed = 0;
        while (committed < fragments.size()) {
            while (probed < fragments.size() && probed < committed + 2 * batch) {
                const size_t end = std::min(fragments.size(), probed + batch);
                probeRange(probed, end);
                probed = end;
            }
            commitRange(std::min(fragments.size(), committed + batch));
        }
        check(snapshotIds == serialIds, "Snapshot commits assign every fragment to the serial cluster");
        check(snapshot.groups().size() == serial.groups().size(), "Snapshot commits create the serial anchors");
        for (size_t g = 0; g < serial.groups().size(); ++g) {
            check(snapshot.groups()[g]->key() == serial.groups()[g]->key(), "Snapshot anchors keep creation order");
            check(snapshot.groups()[g]->representatives.totalSum() == serial.groups()[g]->representatives.totalSum(),
                "Snapshot cluster counts match serial admission");
        }
        check(snapshot.exactGroups().size() == serial.exactGroups().size(), "Snapshot exact keys match serial admission");
        for (const auto& [key, assignment] : serial.exactGroups()) {
            const auto& other = snapshot.exactGroups().at(key);
            check(other.group == assignment.group, "Snapshot exact keys keep their serial cluster");
            if (assignment.representatives) {
                check(other.representatives && other.representatives->totalSum() == assignment.representatives->totalSum(),
                    "Snapshot exact-key counts match serial admission");
            }
        }
        check(snapshot.stats().exactHits == serial.stats().exactHits, "Exact reuse is counted as in serial admission");
    }
}

// A large group's variant FASTQ is written from one sorted snapshot, decoded
// and sorted in parallel tasks and formatted in ranges: the concatenated
// ranges must equal the serial export, for single-end and paired variants.
void testParallelVariantSnapshot() {
    std::mt19937_64 rng(2026092402);
    const auto r1 = randomSequence(rng, 150), r2 = randomSequence(rng, 160);
    diffDNA paired(std::make_shared<const std::string>(r1), true, true, std::make_shared<const std::string>(r2));
    diffDNA single(std::make_shared<const std::string>(r1), false, true);
    for (size_t i = 0; i < 10000; ++i) {
        auto a = r1, b = r2;
        mutate(a[rng() % a.size()]); mutate(a[rng() % a.size()]); mutate(b[rng() % b.size()]);
        auto x = readDNA(a, 10 + rng() % 30, i % 3), y = readDNA(b, 10 + rng() % 30, i % 3);
        paired.addPairQuality(*x, *y, 1, i % 3);
        single.addQuality(*x, 1, i % 3);
    }
    check(paired.uniqueCount() > 4096 && single.uniqueCount() > 4096, "Enough variants to split decoding and sorting");
    size_t calls = 0;
    const diffDNA::ParallelFor parallel = [&calls](size_t count, const std::function<void(size_t)>& task) {
        ++calls;
        std::vector<std::thread> threads;
        for (size_t i = 0; i < count; ++i) { threads.emplace_back(task, i); }
        for (auto& thread : threads) { thread.join(); }
    };
    const std::string cluster = "parent;size=10000; description";
    std::ostringstream serial1, serial2, serialSingle;
    paired.writeFastq(serial1, serial2, cluster, 33);
    single.writeFastq(serialSingle, cluster, 33);
    const auto pairedSnapshot = paired.sortedFragments(&parallel);
    const auto singleSnapshot = single.sortedFragments(&parallel);
    check(calls > 2 && pairedSnapshot.size() == paired.uniqueCount() && singleSnapshot.size() == single.uniqueCount(),
        "Parallel snapshots keep every variant");
    std::ostringstream ranged1, ranged2, rangedSingle;
    for (size_t begin = 0; begin < pairedSnapshot.size(); begin += 777) {
        std::ostringstream part1, part2;
        paired.writeFastq(part1, part2, cluster, pairedSnapshot, begin, std::min(pairedSnapshot.size(), begin + 777), 33);
        ranged1 << part1.str(); ranged2 << part2.str();
    }
    for (size_t begin = 0; begin < singleSnapshot.size(); begin += 1000) {
        single.writeFastq(rangedSingle, cluster, singleSnapshot, begin, std::min(singleSnapshot.size(), begin + 1000), 33);
    }
    check(ranged1.str() == serial1.str() && ranged2.str() == serial2.str(), "Ranged paired export equals the serial export");
    check(rangedSingle.str() == serialSingle.str(), "Ranged single-end export equals the serial export");
}
// Several workers sealing numbered blocks out of order through the global
// index (concurrent probes, non-blocking ordered commits, update lanes), then
// the parallel output (large groups in ranges), must write exactly what one
// worker adding the reads in input order writes.
void testParallelGlobalDerep(const std::string& base) {
    std::mt19937_64 rng(2026092401);
    struct Pair { std::string first, second; int q1, q2, sample; bool pass; };
    std::vector<Pair> pairs;
    const auto r1 = randomSequence(rng, 200), r2 = randomSequence(rng, 200);
    for (size_t i = 0; i < 9000; ++i) { // one group with thousands of exact variants
        Pair pair{r1, r2, int(20 + rng() % 20), int(20 + rng() % 20), int(i % 2), rng() % 10 != 0};
        if (i % 3) { mutate(pair.first[10 + rng() % 180]); }
        if (i % 4) { mutate(pair.second[10 + rng() % 180]); }
        pairs.push_back(pair);
    }
    for (size_t group = 0; group < 400; ++group) { // many small groups
        const auto s1 = randomSequence(rng, 200), s2 = randomSequence(rng, 200);
        for (size_t copy = 0, copies = 1 + rng() % 6; copy < copies; ++copy) {
            Pair pair{s1, s2, int(20 + rng() % 20), int(20 + rng() % 20), int(copy % 2), rng() % 5 != 0};
            if (copy % 2) { mutate(pair.first[10 + rng() % 180]); }
            pairs.push_back(pair);
        }
    }
    std::shuffle(pairs.begin(), pairs.end(), rng);
    auto run = [&](int threads, const std::string& stem) {
        auto opts = options(stem);
        opts["-derepStoreQuals"] = "1"; opts["-derepIdentity"] = "97"; opts["-derepGlobal"] = "1";
        opts["-min_derep_copies"] = "2"; opts["-threads"] = std::to_string(threads);
        Filters filters(&opts); filters.setPaired() = 2; filters.SampleID = {"A", "B"};
        filters.collectStatistics[0]->BarcodeDetected.assign(2, -1); filters.collectStatistics[0]->BarcodeDetectedFail.assign(2, 0);
        Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
        constexpr size_t kBlock = 37;
        const size_t blocks = (pairs.size() + kBlock - 1) / kBlock;
        std::vector<uint64_t> numbers(blocks);
        for (auto& number : numbers) { number = derep.reserveBlock(); } // the reader numbers blocks in input order
        std::vector<size_t> claims(blocks);
        for (size_t i = 0; i < blocks; ++i) { claims[i] = i; }
        if (threads > 1) { std::shuffle(claims.begin(), claims.end(), std::mt19937_64(7)); }
        std::atomic<size_t> next{0};
        auto worker = [&](unsigned seed) {
            std::mt19937_64 jitter(seed);
            for (size_t claim; (claim = next.fetch_add(1)) < blocks;) {
                const size_t block = claims[claim];
                derep.beginBlock(numbers[block]);
                for (size_t i = block * kBlock; i < std::min(pairs.size(), (block + 1) * kBlock); ++i) {
                    const auto& pair = pairs[i];
                    derep.addDNA(readDNA(pair.first, pair.q1, pair.sample), readDNA(pair.second, pair.q2, pair.sample, pair.pass));
                }
                if (jitter() % 3 == 0) { std::this_thread::yield(); }
                derep.endBlock();
            }
        };
        std::vector<std::thread> workers;
        for (int t = 0; t < threads; ++t) { workers.emplace_back(worker, unsigned(t + 1)); }
        for (auto& thread : workers) { thread.join(); }
        derep.writeDereplDNA(&filters, "");
        std::vector<std::string> files;
        for (const char* suffix : {".fq", ".fq.rest", ".map", ".1.hq.fq", ".2.hq.fq"}) { files.push_back(fileText(stem + suffix)); }
        return files;
    };
    const auto serial = run(1, base + "_t1");
    check(pairedFastq(serial[3], serial[4]).size() > 4096, "One group exports enough variants for parallel ranges");
    for (int threads : {2, 4, 8}) {
        check(run(threads, base + "_t" + std::to_string(threads)) == serial, "Global-index output is independent of the worker count");
    }
}
// More output groups than the parallel output sort's threshold (2 x 16384; about
// five sorted runs here, so the pairwise merge carries an odd run), with many equal
// counts: the written order must be count descending, then key ascending, exactly
// as a serial std::sort orders them.
void testParallelOutputSort(const std::string& base) {
    std::mt19937_64 rng(2026092403);
    auto opts = options(base);
    opts["-derepIdentity"] = "100"; opts["-derepCoarseClusters"] = "0"; opts["-min_derep_copies"] = "1";
    Filters filters(&opts); filters.SampleID = {"A"};
    filters.collectStatistics[0]->BarcodeDetected.assign(1, -1); filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
    Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
    std::map<std::string, uint64_t> counts;
    for (size_t group = 0; group < 90000; ++group) {
        const auto sequence = randomSequence(rng, 60);
        const uint64_t copies = 1 + rng() % 4;
        for (uint64_t copy = 0; copy < copies; ++copy) { derep.addDNA(readDNA(sequence), nullptr); }
        counts[sequence] += copies;
    }
    derep.writeDereplDNA(&filters, "");
    std::vector<std::pair<std::string, uint64_t>> expected(counts.begin(), counts.end());
    std::sort(expected.begin(), expected.end(), [](const auto& a, const auto& b) {
        return a.second != b.second ? a.second > b.second : a.first < b.first;
    });
    std::istringstream records(fileText(base + ".fq"));
    std::vector<std::string> written;
    for (std::string header, sequence, plus, quality; std::getline(records, header) && std::getline(records, sequence) &&
            std::getline(records, plus) && std::getline(records, quality);) {
        written.push_back(sequence);
    }
    check(written.size() == expected.size(), "Every exact group is written once");
    for (size_t i = 0; i < written.size(); ++i) {
        check(written[i] == expected[i].first, "Output order is count descending, then key ascending");
    }
}
int main() {
    ini_DNAconstants();
    OptContainer globalOptions; setGlobalOptContainer(&globalOptions);
    const auto directory = std::filesystem::temp_directory_path() / ("sdm-coarse-test-" + std::to_string(std::random_device{}()));
    std::filesystem::create_directories(directory);
    try {
        testExactOutput(); testCompactStorage((directory / "compact-storage").string()); testOrdinaryStorageUpdates((directory / "ordinary-storage").string()); testStorageOrderEffects(); testMixedPairAdmission(); testExactReuse(); testRetainedQualities(); testQualityOutput((directory / "quality").string());
        testCoarsePartitions((directory / "partitions").string());
        testPairedVariants(); testPairedCandidateOracle(); testPairedOutput((directory / "pairs").string());
        testSeedsAndThresholds(); testDiffs(); testUnifiedCounts(); testUnifiedVariantCounts(); testIndexAgainstExhaustive(); testAmpliconSeedPruning(); testRepresentativeAndFixedAnchor(); testRefinement();
        testSnapshotCommits();
        testParallelVariantSnapshot(); testParallelGlobalDerep((directory / "parallel-global").string());
        testParallelOutputSort((directory / "parallel-sort").string());
        testIntegration((directory / "run").string());
        testRefinementOutput((directory / "refine").string());
        std::filesystem::remove_all(directory);
        std::cout << "Coarse dereplication and diffDNA tests passed\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << "\nFixtures: " << directory << '\n'; return 1;
    }
}

#endif
