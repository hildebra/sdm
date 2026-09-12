#include "HtsInput.h"
#include "Common.h"
#include <array>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <zlib.h>
#ifdef SDM_WITH_HTSLIB
#include <htslib/hts.h>
#include <htslib/sam.h>
#endif

namespace {
void requireLocalFile(const std::string& path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) {
        throw std::runtime_error("HTS input requires a readable, regular local file: " + path);
    }
}
}

std::string HtsInput::detectFormat(const std::string& filename) {
    // Recognize extensions even in a portable build, to give an actionable error.
    std::string lower = filename;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (lower.ends_with(".gz")) lower.resize(lower.size() - 3);
    for (const std::string ext : {"sam", "bam", "cram"}) {
        if (lower.ends_with("." + ext)) return "-i_" + ext;
    }
    // gzread transparently reads plain files too. Probe only a small prefix;
    // BAM's magic is inside its BGZF stream, CRAM's is uncompressed.
    struct stat info;
    if (stat(filename.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) return "";
    gzFile input = gzopen(filename.c_str(), "rb");
    if (!input) return "";
    char prefix[4096];
    const int count = gzread(input, prefix, sizeof(prefix));
    gzclose(input);
    if (count <= 0) return "";
    const std::string text(prefix, count);
    if (text.starts_with(std::string("BAM\1", 4))) return "-i_bam";
    if (text.starts_with("CRAM")) return "-i_cram";
    if (text.starts_with("@HD\t") || text.starts_with("@SQ\t") ||
        text.starts_with("@RG\t") || text.starts_with("@PG\t") || text.starts_with("@CO\t"))
        return "-i_sam";
    // Headerless SAM: unlike a FASTQ record, the first line has >= 11 fields.
    const auto end = text.find('\n');
    const auto line = text.substr(0, end);
    if (!line.empty() && line[0] != '@' && line[0] != '>' &&
        std::count(line.begin(), line.end(), '\t') >= 10) return "-i_sam";
    return "";
}

#ifdef SDM_WITH_HTSLIB
namespace {
using FilePtr = std::unique_ptr<samFile, decltype(&hts_close)>;
using HeaderPtr = std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)>;
using RecordPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

// Compare digit runs numerically without integer overflow. Both natural and
// lexicographic query-name sorts are accepted, consistently over the whole file.
int naturalCompare(const std::string& a, const std::string& b) {
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (std::isdigit(static_cast<unsigned char>(a[i])) &&
            std::isdigit(static_cast<unsigned char>(b[j]))) {
            size_t ae = i, be = j;
            while (ae < a.size() && std::isdigit(static_cast<unsigned char>(a[ae]))) ++ae;
            while (be < b.size() && std::isdigit(static_cast<unsigned char>(b[be]))) ++be;
            size_t az = i, bz = j;
            while (az < ae && a[az] == '0') ++az;
            while (bz < be && b[bz] == '0') ++bz;
            if (ae - az != be - bz) return ae - az < be - bz ? -1 : 1;
            const int cmp = a.compare(az, ae - az, b, bz, be - bz);
            if (cmp) return cmp;
            i = ae; j = be;
        } else {
            if (a[i] != b[j]) return a[i] < b[j] ? -1 : 1;
            ++i; ++j;
        }
    }
    return i == a.size() ? (j == b.size() ? 0 : -1) : 1;
}

std::string stringTag(const bam1_t* record, const char* tag) {
    const uint8_t* data = bam_aux_get(record, tag);
    if (!data) return "";
    if (*data != 'Z') throw std::runtime_error(std::string(tag) + " tag must have type Z");
    const char* value = bam_aux2Z(data);
    if (!value) throw std::runtime_error(std::string("Invalid ") + tag + " tag");
    return value;
}
}

struct HtsInput::Impl {
    std::string path, reference, onlyPair;
    bool tags;
    int requestedPairs, sourcePairs = 1, outputPairs = 1;
    FilePtr file{nullptr, &hts_close};
    HeaderPtr header{nullptr, &sam_hdr_destroy};
    RecordPtr next{bam_init1(), &bam_destroy1};
    bool haveNext = false, natural = true, lexical = true;
    std::string previousName;
    std::array<std::vector<std::string>, 2> pending;
    int expectedPosition = 0;

    Impl(const std::string& p, int pairs, const std::string& only,
         const std::string& ref, bool bc)
        : path(p), reference(ref), onlyPair(only), tags(bc), requestedPairs(pairs) {
        if (!next) throw std::bad_alloc();
        if (pairs != -1 && pairs != 1 && pairs != 2)
            throw std::runtime_error("HTS input supports -paired 1 or 2, not a three-file layout");
        if (!only.empty() && only != "1" && only != "2")
            throw std::runtime_error("HTS -onlyPair must be 1 or 2");
        open();
    }
    [[noreturn]] void fail(const std::string& message) const {
        throw std::runtime_error(message + " [" + path + "]");
    }
    void advance() {
        int result;
        do {
            result = sam_read1(file.get(), header.get(), next.get());
        } while (result >= 0 && (next->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)));
        if (result < -1) fail("Failed to decode HTS record (check file integrity and CRAM reference)");
        haveNext = result >= 0;
    }
    void open() {
        requireLocalFile(path);
        file.reset(sam_open(path.c_str(), "r"));
        if (!file) fail("Cannot open HTS input");
        const auto format = hts_get_format(file.get())->format;
        if ((format == bam || format == cram) && hts_check_EOF(file.get()) == 0)
            fail("Missing HTS EOF marker; input may be truncated");
        if (format != sam && format != bam && format != cram) fail("Expected SAM, BAM or CRAM input");
        if (format == cram && !reference.empty() &&
            hts_set_fai_filename(file.get(), reference.c_str()) < 0)
            fail("Cannot load -cramRef " + reference);
        header.reset(sam_hdr_read(file.get()));
        if (!header) fail("Cannot read HTS header");
        pending = {}; expectedPosition = 0; previousName.clear(); natural = lexical = true;
        advance();
        sourcePairs = haveNext ? ((next->core.flag & BAM_FPAIRED) ? 2 : 1)
                               : (requestedPairs == 2 ? 2 : 1);
        if (haveNext && onlyPair.empty() && requestedPairs != -1 && requestedPairs != sourcePairs)
            fail("-paired does not match HTS records; use -paired 2 for paired input");
        if (!onlyPair.empty() && sourcePairs != 2) fail("-onlyPair requires paired HTS input");
        outputPairs = onlyPair.empty() ? sourcePairs : 1;
        if (sourcePairs == 2) {
            kstring_t order{0, 0, nullptr};
            sam_hdr_find_tag_hd(header.get(), "SO", &order);
            const std::string sortOrder = order.s ? order.s : "";
            free(order.s);
            if (sortOrder == "coordinate")
                fail("Paired HTS input must be query-name sorted; run samtools sort -n or -N");
        }
    }
    std::vector<std::string> convert(const bam1_t* record) {
        const int length = record->core.l_qseq;
        if (length <= 0) fail("Missing SEQ for " + std::string(bam_get_qname(record)));
        std::vector<std::string> lines(4);
        lines[0] = "@" + std::string(bam_get_qname(record));
        if (sourcePairs == 2) lines[0] += (record->core.flag & BAM_FREAD1) ? "/1" : "/2";
        lines[1].resize(length); lines[2] = "+"; lines[3].resize(length);
        const uint8_t* seq = bam_get_seq(record);
        const uint8_t* qual = bam_get_qual(record);
        for (int i = 0; i < length; ++i) {
            const char base = seq_nt16_str[bam_seqi(seq, i)];
            if (base == '=') fail("Reference-equality SEQ bases are unsupported; provide explicit bases");
            lines[1][i] = base;
            if (qual[i] > 93) fail("Missing or out-of-range QUAL for " + lines[0]);
            lines[3][i] = static_cast<char>(qual[i] + 33);
        }
        if (record->core.flag & BAM_FREVERSE) {
            reverseTS(lines[1]);
            std::reverse(lines[3].begin(), lines[3].end());
        }
        if (tags) {
            auto bc = stringTag(record, "BC");
            std::transform(bc.begin(), bc.end(), bc.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            if (!bc.empty() && (bc.front() == '-' || bc.back() == '-' ||
                std::count(bc.begin(), bc.end(), '-') > 1 ||
                bc.find_first_not_of("ACGTN-") != std::string::npos))
                fail("BC tag must contain one barcode or two barcodes separated by '-'");
            lines.push_back(std::move(bc));
        }
        return lines;
    }
    bool loadSet() {
        if (!haveNext) return false;
        const std::string name = bam_get_qname(next.get());
        if (name.empty() || name == "*") fail("HTS reads require a query name");
        if (sourcePairs == 2 && !previousName.empty()) {
            lexical = lexical && previousName < name;
            natural = natural && naturalCompare(previousName, name) < 0;
            if (!lexical && !natural) fail("Paired HTS records are not query-name sorted near " + name);
        }
        std::array<std::vector<std::string>, 2> reads;
        do {
            const auto flag = next->core.flag;
            if (bool(flag & BAM_FPAIRED) != (sourcePairs == 2))
                fail("Mixed singleton and paired HTS records near " + name);
            int slot = 0;
            if (sourcePairs == 2) {
                if (bool(flag & BAM_FREAD1) == bool(flag & BAM_FREAD2))
                    fail("Paired HTS record requires exactly one READ1/READ2 flag: " + name);
                slot = (flag & BAM_FREAD2) ? 1 : 0;
            }
            if (!reads[slot].empty()) fail("Duplicate primary record for " + name);
            reads[slot] = convert(next.get());
            advance();
        } while (haveNext && name == bam_get_qname(next.get()));
        if (sourcePairs == 2 && (reads[0].empty() || reads[1].empty()))
            fail("Missing adjacent mate for " + name + "; paired input must contain complete query-name sorted pairs");
        if (tags) {
            std::string bc = reads[0][4];
            if (sourcePairs == 2) {
                if (bc.empty()) bc = reads[1][4];
                if (!reads[1][4].empty() && bc != reads[1][4]) fail("Conflicting BC tags for " + name);
            }
            if (bc.empty()) fail("Missing BC:Z tag for " + name + " with -htsBarcode BC");
            reads[0][4] = bc;
            if (sourcePairs == 2) reads[1][4] = bc;
        }
        previousName = name;
        pending = std::move(reads);
        if (!onlyPair.empty()) {
            if (onlyPair == "2") pending[0] = std::move(pending[1]);
            pending[1].clear();
        }
        return true;
    }
};

HtsInput::HtsInput(const std::string& f, int p, const std::string& only,
                   const std::string& ref, bool tags)
    : impl_(std::make_unique<Impl>(f, p, only, ref, tags)) {}
HtsInput::~HtsInput() = default;
std::string HtsInput::version() { return hts_version(); }
int HtsInput::pairCount() const { return impl_->outputPairs; }
const std::string& HtsInput::filename() const { return impl_->path; }
bool HtsInput::hasData() const {
    return impl_->haveNext || !impl_->pending[0].empty() || !impl_->pending[1].empty();
}
void HtsInput::reset() { impl_->open(); }
bool HtsInput::read(std::vector<std::string>& record, int position) {
    record.clear();
    if (position < 0 || position >= impl_->outputPairs) return false;
    if (!hasData()) return false;
    if (position != impl_->expectedPosition)
        impl_->fail("Read HTS mates in R1/R2 order");
    if (position == 0 && !impl_->loadSet()) return false;
    record = std::move(impl_->pending[position]);
    impl_->pending[position].clear();
    impl_->expectedPosition = (position + 1) % impl_->outputPairs;
    return true;
}
#else
struct HtsInput::Impl {};
HtsInput::HtsInput(const std::string&, int, const std::string&, const std::string&, bool) {
    throw std::runtime_error("SAM/BAM/CRAM input requires HTSlib support; rebuild sdm with make WITH_HTS=1");
}
HtsInput::~HtsInput() = default;
std::string HtsInput::version() { return "disabled"; }
bool HtsInput::read(std::vector<std::string>&, int) { return false; }
bool HtsInput::hasData() const { return false; }
int HtsInput::pairCount() const { return 1; }
void HtsInput::reset() {}
const std::string& HtsInput::filename() const { static const std::string empty; return empty; }
#endif
