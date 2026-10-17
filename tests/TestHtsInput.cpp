#include "InputStream.h"
#include "HtsInput.h"
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <filesystem>
#include <functional>
#include <fstream>
#include <stdexcept>
#include <unistd.h>

namespace fs = std::filesystem;
namespace {
void check(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
void rejects(const std::function<void()>& fn, const std::string& message) {
    try { fn(); } catch (const std::runtime_error& error) {
        check(std::string(error.what()).find(message) != std::string::npos,
              "Wrong error: " + std::string(error.what()) + "; expected " + message);
        return;
    }
    throw std::runtime_error("Expected rejection: " + message);
}
void write(const fs::path& p, const std::string& data) { std::ofstream(p) << data; }
void convert(const std::string& in, const std::string& out, const std::string& format,
             const std::string& reference = "") {
    std::unique_ptr<samFile, decltype(&hts_close)> src(sam_open(in.c_str(), "r"), &hts_close);
    check(bool(src), "Open fixture source");
    std::unique_ptr<sam_hdr_t, decltype(&sam_hdr_destroy)> hdr(sam_hdr_read(src.get()), &sam_hdr_destroy);
    const char* mode = format == "bam" ? "wb" : "wc";
    std::unique_ptr<samFile, decltype(&hts_close)> dst(sam_open(out.c_str(), mode), &hts_close);
    check(bool(dst) && bool(hdr), "Open fixture output/header");
    if (format != "bam") {
        if (reference.empty()) check(hts_set_opt(dst.get(), CRAM_OPT_NO_REF, 1) == 0, "no_ref");
        else check(hts_set_fai_filename(dst.get(), reference.c_str()) == 0, "Set reference");
        if (format == "embedded") check(hts_set_opt(dst.get(), CRAM_OPT_EMBED_REF, 1) == 0, "embed_ref");
    }
    check(sam_hdr_write(dst.get(), hdr.get()) >= 0, "Write fixture header");
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> rec(bam_init1(), &bam_destroy1);
    int result;
    while ((result = sam_read1(src.get(), hdr.get(), rec.get())) >= 0)
        check(sam_write1(dst.get(), hdr.get(), rec.get()) >= 0, "Write fixture record");
    check(result == -1, "Read fixture EOF");
    check(sam_close(dst.release()) == 0, "Close fixture output");
}
std::string header(const std::string& order = "queryname") {
    return "@HD\tVN:1.6\tSO:" + order + "\n@SQ\tSN:chr1\tLN:1000\n";
}
std::string rec(const std::string& name, int flags, const std::string& sequence = "ACGA",
                const std::string& quality = "IIII", const std::string& tags = "") {
    return name + "\t" + std::to_string(flags) + "\tchr1\t1\t60\t4M\t*\t0\t0\t" + sequence + "\t" + quality + tags + "\n";
}
void consume(const fs::path& path, int pairs = -1, bool tags = false, const std::string& reference = "") {
    HtsInput input(path.string(), pairs, "", reference, tags);
    std::vector<std::string> record;
    while (input.read(record, 0)) {
        if (input.pairCount() == 2) check(input.read(record, 1), "Missing test R2");
    }
}
}
int main(int argc, char** argv) {
    try {
        if (argc >= 5 && std::string(argv[1]) == "--convert") {
            convert(argv[2], argv[3], argv[4], argc == 6 ? argv[5] : "");
            return 0;
        }
        ini_DNAconstants();
        char pattern[] = "/tmp/sdm-hts-test-XXXXXX";
        const char* temporary = mkdtemp(pattern);
        check(temporary != nullptr, "Create temporary directory");
        const fs::path root(temporary);
        struct Cleanup { fs::path p; ~Cleanup() { fs::remove_all(p); } } cleanup{root};
        const auto sam = root / "reads.sam";
        write(sam, header() + rec("r1", 145, "TCGT", "$#\"!", "\tBC:Z:ACGT-TGCA") +
              rec("r1", 65, "ACGA", "!\"#$") + rec("r1", 321) + rec("r1", 2113) +
              rec("r2", 65) + rec("r2", 129, "ACGA", "IIII", "\tBC:Z:ACGT-TGCA"));
        const auto ref = root / "reference.fa";
        write(ref, ">chr1\n" + std::string(1000, 'A') + "\n");
        check(fai_build(ref.c_str()) == 0, "Build reference index");
        for (const std::string format : {"sam", "bam", "cram", "external", "embedded"}) {
            const auto path = format == "sam" ? sam : root / (format + (format == "bam" ? ".bam" : ".cram"));
            if (format != "sam") convert(sam, path, format,
                (format == "external" || format == "embedded") ? ref.string() : "");
            check(detectSeqFmt(path) == (format == "sam" ? "-i_sam" : format == "bam" ? "-i_bam" : "-i_cram"), "Format detection");
            // Exercise the normal InputStreamer APIs, including fixed-quality
            // decoding even when the caller explicitly requests Phred+64.
            InputStreamer input(false, 64, "0", "1", 2);
            input.configureHts(format == "external" ? ref.string() : "", true);
            int pairs = -1;
            check(input.setupFastq("", path, pairs, ""), "Open input");
            check(pairs == 2 && input.fastQscore() == 33, "Pair inference / fixed quality offset");
            auto first = input.getDNApairs();
            check(first[0] && first[1] && first[0]->getSequence() == "ACGA" && first[1]->getSequence() == "ACGA", "Reverse complement restoration");
            check(first[0]->getQual() == first[1]->getQual() && first[0]->getQual()[0] == 0 && first[0]->getQual()[3] == 3, "Quality orientation / offset");
            check(first[0]->sampleBarcode() == "ACGT-TGCA" && first[1]->sampleBarcode() == "ACGT-TGCA", "Mate barcode inheritance");
            input.allStreamReset();
            input.setMaxRdsRead(2);
            multi_tmp_lines batch(2);
            check(!input.getDNAlines(&batch, 2, false, true) && batch.size() == 1, "Limit must preserve complete pairs");
            auto decoded = str2DNA(std::move(batch.tmp[0][1]), true, 0, 1);
            check(decoded->getSequence() == "ACGA" && decoded->getQual()[3] == 3, "Batch conversion");
            input.allStreamReset();
            input.setMaxRdsRead(0);
            check(input.getDNA(0)->getId() == "r1/1" && input.getDNA(1)->getId() == "r1/2", "Direct paired API");
            check(input.getDNApairs()[0]->getId() == "r2/1", "Skip secondary/supplementary");
            check(!input.getDNApairs()[0] && !input.checkInFileStatus(), "Clean EOF");
            input.allStreamReset();
            check(input.checkInFileStatus(), "Reset EOF");
            HtsInput only(path, -1, "2", format == "external" ? ref.string() : "", false);
            std::vector<std::string> record;
            check(only.pairCount() == 1 && only.read(record, 0) && record[0] == "@r1/2", "onlyPair 2");
        }
        fs::copy_file(root / "bam.bam", root / "extensionless");
        check(detectSeqFmt(root / "extensionless") == "-i_bam", "BAM magic detection");
        fs::copy_file(root / "cram.cram", root / "no-extension");
        check(detectSeqFmt(root / "no-extension") == "-i_cram", "CRAM magic detection");
        const auto bad = root / "bad.sam";
        auto invalid = [&](const std::string& data, const std::string& error, bool tags = false) {
            write(bad, data); rejects([&] { consume(bad, -1, tags); }, error);
        };
        invalid(header("coordinate") + rec("r1", 65) + rec("r1", 129), "query-name sorted");
        invalid(header() + rec("r1", 65), "Missing adjacent mate");
        invalid(header() + rec("r1", 65) + rec("r2", 129), "Missing adjacent mate");
        invalid(header() + rec("r1", 65) + rec("r1", 65), "Duplicate primary");
        invalid(header() + rec("r1", 193), "exactly one READ1/READ2");
        invalid(header() + rec("r1", 65) + rec("r1", 129) + rec("r2", 0), "Mixed singleton");
        invalid(header() + rec("r2", 65) + rec("r2", 129) + rec("r1", 65) + rec("r1", 129), "not query-name sorted");
        invalid(header() + rec("r1", 0, "ACGA", "*"), "Missing or out-of-range QUAL");
        invalid(header() + "r1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n", "Missing SEQ");
        invalid(header() + rec("r1", 0), "Missing BC:Z", true);
        invalid(header() + rec("r1", 0, "ACGA", "IIII", "\tBC:i:7"), "type Z", true);
        invalid(header() + rec("r1", 65, "ACGA", "IIII", "\tBC:Z:AAAA") +
                rec("r1", 129, "ACGA", "IIII", "\tBC:Z:CCCC"), "Conflicting BC", true);
        invalid(header() + "broken\n", "Failed to decode");
        invalid(header() + rec("r1", 65) + rec("r1", 129) + rec("r01", 65) + rec("r01", 129) +
                rec("r1", 65) + rec("r1", 129), "not query-name sorted");
        for (const auto& names : {std::vector<std::string>{"r1", "r2", "r10"},
                                 std::vector<std::string>{"r1", "r10", "r2"}}) {
            std::string data = header();
            for (const auto& name : names) data += rec(name, 65) + rec(name, 129);
            write(bad, data); consume(bad, 2);
        }
        // SAM has no compulsory header for unmapped singleton records.
        write(root / "headerless", "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGA\tIIII\n");
        check(detectSeqFmt(root / "headerless") == "-i_sam", "Headerless SAM detection");
        consume(root / "headerless", 1);
        // Singleton input may be coordinate sorted; it must still be single-end.
        write(bad, header("coordinate") + rec("z", 0) + rec("a", 16));
        consume(bad, 1);
        rejects([&] { consume(sam, 1); }, "-paired does not match");
        // A corrupted binary stream must not silently become EOF.
        const auto truncated = root / "truncated.bam";
        fs::copy_file(root / "bam.bam", truncated);
        fs::resize_file(truncated, fs::file_size(truncated) / 2);
        rejects([&] { consume(truncated); }, "HTS");
        const auto wrong = root / "wrong.fa";
        write(wrong, ">chr1\n" + std::string(1000, 'C') + "\n");
        check(fai_build(wrong.c_str()) == 0, "Build wrong reference");
        rejects([&] { consume(root / "external.cram", 2, false, wrong); }, "decode");
        // Hide the UR reference and disable remote/cache fallback for this test.
        fs::rename(ref, root / "hidden.fa");
        setenv("REF_PATH", "/nonexistent/sdm-hts-reference/%s", 1);
        setenv("REF_CACHE", "/nonexistent/sdm-hts-cache/%s", 1);
        rejects([&] { consume(root / "external.cram"); }, "decode");
        consume(root / "embedded.cram");
        consume(root / "cram.cram");
        std::cout << "All HTS input tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "HTS test failure: " << error.what() << '\n';
        return 1;
    }
}
