// Production dereplication benchmark; reuse the integration-test fixtures.
// Substitution simulation with the earlier conserved/variable 16S-like model.
#define SDM_COARSE_TEST_NO_MAIN
#include "TestCoarseDerep.cpp"
#include <chrono>
#include <iomanip>
#include "CoarseSimulation.h"
#if defined(__linux__)
#include <sys/resource.h>
#endif

int main(int argc, char** argv) {
    if (argc < 2 || argc > 10) {
        std::cerr << "Usage: TestCoarseBenchmark <exact|identity-percent> [reads=20000] [diffs=0|1] [length-jitter=0..200] [length=1500] [quality=20] [reassign=0|1] [retain-qualities=0|1] [fa|fq]\n";
        return 2;
    }
    ini_DNAconstants(); OptContainer globals; setGlobalOptContainer(&globals);
    const std::string identity = argv[1];
    const size_t reads = argc > 2 ? std::stoull(argv[2]) : 20000;
    const bool store = argc > 3 && std::string(argv[3]) == "1";
    const size_t jitter = argc > 4 ? std::stoull(argv[4]) : 0;
    const size_t meanLength = argc > 5 ? std::stoull(argv[5]) : 1500;
    const unsigned quality = argc > 6 ? std::stoul(argv[6]) : 20;
    const bool reassign = argc > 7 && std::string(argv[7]) == "1";
    const bool retainQualities = argc > 8 && std::string(argv[8]) == "1";
    const std::string format = argc > 9 ? argv[9] : "fa";
    if (format != "fa" && format != "fq") { return 2; }
    // Supplying the final argument uses variable quality profiles in BOTH modes.
    const bool variableQualities = argc > 8;
    if (!reads || reads > 1000000 || jitter > 200 || meanLength <= jitter || meanLength > 10000 || quality > 60) { return 2; }
    const auto directory = std::filesystem::temp_directory_path() / ("sdm-coarse-bench-" + std::to_string(std::random_device{}()));
    std::filesystem::create_directories(directory);
    const std::string stem = (directory / "clusters").string();
    auto opts = options(stem); opts["-derepIdentity"] = identity == "exact" ? "100" : identity;
    opts["-derepPrefix"] = "0"; opts["-derepStoreDiffs"] = store ? "1" : "0";
    opts["-derepReassign"] = reassign ? "1" : "0";
    opts["-derepStoreQuals"] = retainQualities ? "1" : "0";
    // Default to FASTA; FASTQ also exercises exact per-base quality sums.
    opts["-derep_format"] = format; opts["-o_dereplicate"] = stem + "." + format;
    CoarseSimulation simulation(meanLength, quality, jitter);
    Filters filters(&opts); filters.SampleID = {"sample"};
    filters.collectStatistics[0]->BarcodeDetected.assign(1, -1);
    filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
    Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
    const auto start = std::chrono::steady_clock::now();
    for (size_t r = 0; r < reads; ++r) {
        auto sequence = simulation.next();
        const int observedQuality = variableQualities ? std::clamp(int(quality) + int(r % 3) * 6 - 6, 2, 60) : quality;
        auto dna = readDNA(sequence, observedQuality);
        dna->setHeader("read_" + std::to_string(r));
        derep.addDNA(dna, nullptr);
    }
    const auto ingested = std::chrono::steady_clock::now();
    const auto report = derep.writeDereplDNA(&filters, "");
    uint64_t movedReads = 0;
    const auto movedPosition = report.find("Final reassignment: ");
    if (movedPosition != std::string::npos) {
        auto count = report.substr(movedPosition + 20);
        count.erase(std::remove(count.begin(), count.end(), ','), count.end());
        movedReads = std::stoull(count);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    size_t clusters = 0; uint64_t counted = 0;
    std::ifstream fasta(stem + "." + format); std::string line;
    while (std::getline(fasta, line)) {
        if (!line.empty() && line.front() == (format == "fq" ? '@' : '>')) {
            ++clusters; const auto pos = line.rfind(";size="); check(pos != std::string::npos, "Missing benchmark abundance");
            counted += std::stoull(line.substr(pos + 6));
            if (format == "fq") {
                for (int field = 0; field < 3; ++field) { check(bool(std::getline(fasta, line)), "Complete benchmark FASTQ record"); }
            }
        }
    }
    check(counted == reads, "Benchmark must count every read");
    size_t exactVariants = 0; uint64_t qualityBytes = 0;
    if (retainQualities) {
        std::ifstream fastq(stem + ".1.hq.fq");
        std::string sequence, plus, qualities; uint64_t variantCount = 0;
        while (std::getline(fastq, line)) {
            check(bool(std::getline(fastq, sequence)) && bool(std::getline(fastq, plus)) && bool(std::getline(fastq, qualities)), "Complete benchmark FASTQ");
            const auto pos = line.rfind(";size=");
            check(line.starts_with("@") && pos != std::string::npos && plus == "+" && sequence.size() == qualities.size(), "Valid benchmark FASTQ");
            variantCount += std::stoull(line.substr(pos + 6)); ++exactVariants;
        }
        check(variantCount == reads, "Retained qualities preserve every read count");
        qualityBytes = std::filesystem::file_size(stem + ".1.hq.fq");
    }
    long peakKiB = 0;
#if defined(__linux__)
    rusage usage{}; getrusage(RUSAGE_SELF, &usage); peakKiB = usage.ru_maxrss;
    // ru_maxrss can include the launcher's pre-exec high-water mark. VmHWM
    // describes this executable's address space, without that inherited floor.
    std::ifstream status("/proc/self/status"); std::string statusLine;
    while (std::getline(status, statusLine)) {
        if (statusLine.starts_with("VmHWM:")) { peakKiB = std::stol(statusLine.substr(6)); break; }
    }
#endif
    // Byte-for-byte output comparison across implementations, including HQ mates
    // and sample maps. File names are stable; the temporary directory is ignored.
    uint64_t outputChecksum = 1469598103934665603ULL;
    std::vector<std::filesystem::path> outputs;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) { outputs.push_back(entry.path()); }
    }
    std::sort(outputs.begin(), outputs.end());
    for (const auto& path : outputs) {
        std::ifstream input(path, std::ios::binary);
        char buffer[8192];
        while (input.read(buffer, sizeof(buffer)) || input.gcount()) {
            for (std::streamsize i = 0; i < input.gcount(); ++i) {
                outputChecksum = (outputChecksum ^ static_cast<unsigned char>(buffer[i])) * 1099511628211ULL;
            }
        }
    }
    uint64_t diffBytes = store ? std::filesystem::file_size(stem + ".diff") : 0;
    std::cout << "identity,reads,diffs,length_jitter,length,quality,reassign,moved_reads,clusters,seconds,peak_rss_kib,diff_file_bytes,simulated_error_fraction,retain_qualities,variable_qualities,exact_subclusters,subcluster_fastq_bytes,format,output_checksum,ingest_seconds,write_seconds\n"
        << identity << ',' << reads << ',' << store << ',' << jitter << ',' << meanLength << ',' << quality << ',' << reassign << ',' << movedReads << ',' << clusters << ',' << std::setprecision(8)
        << seconds << ',' << peakKiB << ',' << diffBytes << ',' << double(simulation.errors) / simulation.bases << ',' << retainQualities << ',' << variableQualities << ',' << exactVariants << ',' << qualityBytes << ',' << format << ',' << outputChecksum << ',' << std::chrono::duration<double>(ingested - start).count() << ',' << seconds - std::chrono::duration<double>(ingested - start).count() << '\n';
    std::filesystem::remove_all(directory);
}
