// Focused search audit: independent edit-distance oracle plus legacy boundaries.
#include "Filters.h"
#include <random>
#include <set>
#include <stdexcept>

namespace {
void require(bool ok, const char* message) {
    if (!ok) { throw std::runtime_error(message); }
}
shared_ptr<DNA> read(const string& sequence) {
    auto d = make_shared<DNA>(sequence, "search_test");
    d->setQual(vector<qual_score>(sequence.size(), 30));
    require(d->seal(), "fixture must seal");
    return d;
}
class SearchFilters : public Filters {
public:
    using Filters::Filters;
    using Filters::cutPrimer;
    using Filters::cutPrimerRev;
    void enableDualBarcodes() { bDoBarcode2 = true; }
};
int distance(const string& a, const string& b) {
    vector<vector<int>> matrix(a.size() + 1, vector<int>(b.size() + 1));
    for (size_t i = 0; i <= a.size(); ++i) { matrix[i][0] = i; }
    for (size_t j = 0; j <= b.size(); ++j) { matrix[0][j] = j; }
    for (size_t i = 1; i <= a.size(); ++i) {
        for (size_t j = 1; j <= b.size(); ++j) {
            matrix[i][j] = min({matrix[i-1][j]+1, matrix[i][j-1]+1,
                matrix[i-1][j-1]+(a[i-1] != b[j-1])});
        }
    }
    return matrix.back().back();
}
void editOracle() {
    std::mt19937 rng(351);
    auto dna = [&](int n) {
        string s;
        while (n-- > 0) { s += "ACGT"[rng() % 4]; }
        return s;
    };
    for (int trial = 0; trial < 2000; ++trial) {
        const string barcode = dna(4 + rng() % 11);
        const string sequence = trial % 2 ? dna(rng() % 30) : dna(rng() % 10) + barcode + dna(rng() % 10);
        const int budget = rng() % 4, window = 1 + rng() % 24;
        int best = budget + 1;
        std::set<int> stops;
        for (int start = 0; start < min(window, (int)sequence.size()); ++start) {
            for (int stop = start + 1; stop <= (int)sequence.size(); ++stop) {
                const int cost = distance(barcode, sequence.substr(start, stop - start));
                if (cost > budget) { continue; }
                if (cost < best) { best = cost; stops.clear(); }
                if (cost == best) { stops.insert(stop); }
            }
        }
        const auto hit = DNA::matchSeq_tot(sequence, barcode, budget, window);
        require((hit.start < 0) == stops.empty(), "edit oracle: hit presence");
        require(hit.errors == best, "edit oracle: best edit distance");
        require(hit.ambiguous == (stops.size() > 1), "edit oracle: trim ambiguity");
        if (hit.start >= 0) {
            require(stops.count(hit.stop) && hit.start < window &&
                distance(barcode, sequence.substr(hit.start, hit.stop - hit.start)) == best,
                "edit oracle: selected coordinates");
        }
    }
}
void primerBounds() {
    auto d = read("TTACGTTT");
    require(d->matchSeq("ACGT", 0, 3, 0) == 2, "inclusive last forward start");
    require(d->matchSeq("ACGT", 0, 2, 0) == -1, "exclusive forward limit");
    require(d->matchSeq("", 0, 8, 0) == -1, "empty forward primer must not match");
    require(d->matchSeqRev("", 0, -1) == -1, "empty reverse primer must not match");
    require(d->matchSeqRev("ACGT", 0, -1, -8) == 2, "negative rear offset is clamped safely");
    require(d->matchSeqRev("TT", 0, -1, 100) == -1, "rear offset beyond read must not match");
    auto partial = read("GGGA");
    require(partial->matchSeqRev("ACGT", 0, -1, 1) == -1, "explicit offset cannot admit a partial primer");
    auto trimmed = read("TTTTACGT");
    trimmed->cutSeqPseudo(4);
    int errors = 0;
    require(trimmed->matchSeq_tot("ACGT", 0, 20, errors) == -1, "barcode search respects logical length");
    require(read("ACGT")->matchSeqRev("ACGT", 0, -1) == 0, "full read reverse primer at position zero");
    require(read("AGGT")->matchSeq("ARGT", 0, 1, 0) == 0, "IUPAC primer matching remains supported");
    OptContainer opts;
    SearchFilters filters(&opts);
    filters.PrimerL = {"ACGT"}; filters.PrimerL_RC = {"ACGT"};
    filters.PrimerR = {"TGCA"}; filters.PrimerR_RC = {"TGCA"};
    bool rc = false;
    require(filters.cutPrimer(d, -1, rc, 0), "absent per-sample forward primer is safe");
    require(filters.cutPrimerRev(d, -1, rc), "absent per-sample reverse primer is safe");
    require(!filters.checkIfPrimerHits(d, -1) && !filters.checkIfRevPrimerHits(d, 20) &&
        !filters.findPrimer(d, -1, false, 0), "public primer searches validate sample indices");
    filters.removePrimer(false);
    auto retained = read(string(80, 'A') + "TGCA" + "TTTT");
    rc = true;
    require(filters.cutPrimerRev(retained, 0, rc) && retained->getRevPrimDetect() &&
        retained->getSequence() == string(80, 'A') + "TGCA",
        "retained reverse primer is detected and trailing adapter is removed");
}
void legacyBarcodes() {
    OptContainer opts;
    opts["-sample_sep"] = "__";
    SearchFilters filters(&opts);
    filters.allResize(2);
    filters.Barcode = {"ACGT", "TGCAACGT"};
    filters.BarcodePreStats();
    filters.setTagErrs(1);
    string observed;
    int errors = 0;
    auto shortRead = read("ACGTAA");
    require(filters.findTag(shortRead, observed, errors, true, 0, false, false) == 0,
        "short read must still match a shorter barcode in a variable-length map");
    filters.Barcode = {"ACGTCAGT", "TGCATCGA"};
    filters.BarcodePreStats();
    auto rear = read(string(50, 'G') + "ACGTAAGT" + "TT");
    require(filters.findTag(rear, observed, errors, true, 0, false, true) == 0,
        "approximate rear barcode search must examine the rear");
    require(errors == 1 && observed == "ACGTAAGT", "barcode errors and observed sequence propagate to caller");
    auto front = read("ACGTAAGT" + string(50, 'G'));
    require(filters.findTag(front, observed, errors, true, 0, false, true) == -1,
        "rear search must not rescue a front-only approximate barcode");
    filters.removeTag(false);
    filters.SampleID = {"sampleA", "sampleB"};
    auto kept = read("ACGTCAGTAAAA");
    require(filters.detectCutBC(kept, true) == 0 && kept->getNewIDshort().find("sampleA__") == 0,
        "retaining barcode bases must still assign the output header");
    require(kept->getSequence() == "ACGTCAGTAAAA", "keepBarcodeSeq retains sequence");
    filters.enableDualBarcodes();
    filters.Barcode = {"ACGTCAGTGCTA", "TGCATCGACAGT"};
    filters.Barcode2 = {"ATCG", "GCTA"};
    filters.BarcodePreStats(); filters.setTagErrs(0); filters.removeTag(true);
    require(filters.detectCutBC(read("ATCGTT"), false) == 0,
        "short mate 2 barcode must not be checked against mate 1's minimum length");
}
}
int main() {
    ini_DNAconstants();
    try {
        editOracle();
        primerBounds();
        legacyBarcodes();
    } catch (const std::exception& e) {
        cerr << "Search regression failed: " << e.what() << '\n';
        return 1;
    }
    cout << "Search audit and 2,000 edit-distance oracle cases passed\n";
}
