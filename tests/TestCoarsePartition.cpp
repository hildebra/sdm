// Investigation only: real coarse indexes, owned partitions, centroid-only merge.
// Not a production dereplication backend; preloaded reads are retained for validation.
#define SDM_COARSE_TEST_NO_MAIN
#include "TestCoarseDerep.cpp"
#include "CoarseSimulation.h"
#include <atomic>
#include <chrono>
#include <iomanip>
#include <numeric>

using Clock = std::chrono::steady_clock;
static double elapsed(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }
static uint64_t mix(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}
static std::vector<uint64_t> seeds(std::string_view sequence) {
    std::vector<uint64_t> result;
    for (size_t pos = 0; pos + 20 <= sequence.size(); pos += 20) {
        uint64_t word;
        if (CoarseDereplication::seedWord(sequence, pos, 4, word)) { result.push_back(word | ((pos / 20) << 40)); }
    }
    return result;
}
int main(int argc, char** argv) {
    if (argc != 8) {
        std::cerr << "Usage: TestCoarsePartition <exact|exact_diff|global|roundrobin|fullhash|min20|min8|rare20|supported20> threads buckets reads length quality identity\n";
        return 2;
    }
    ini_DNAconstants(); OptContainer globals; setGlobalOptContainer(&globals);
    const std::string mode = argv[1], identity = argv[7];
    const size_t threads = std::stoull(argv[2]), buckets = std::stoull(argv[3]), n = std::stoull(argv[4]);
    const size_t length = std::stoull(argv[5]); const unsigned quality = std::stoul(argv[6]);
    check(threads && threads <= 64 && buckets && buckets <= 512 && n && n <= 1000000 && length >= 20 && length <= 10000,
        "Invalid benchmark dimensions");
    const uint64_t threshold = CoarseDereplication::parseIdentity(identity);
    CoarseSimulation simulation(length, quality);
    std::vector<std::string> input; input.reserve(n);
    std::vector<size_t> source; source.reserve(n);
    for (size_t i = 0; i < n; ++i) { input.push_back(simulation.next()); source.push_back(simulation.lastTemplate); }
    double routeSeconds = 0, localSeconds = 0, mergeSeconds = 0, waitSeconds = 0;
    size_t localGroups = 0, finalGroups = 0, violations = 0;
    double maxBucketFraction = 0, pairTogether = 0;
    std::vector<size_t> workerReads(threads, 0);
    std::vector<size_t> route(n, 0), assignment(n, 0);
    const auto allStart = Clock::now();
    if (mode == "exact" || mode == "exact_diff") {
        check(identity == "100", "Production exact cases require identity 100");
        const auto dir = std::filesystem::temp_directory_path() / ("sdm-partition-" + std::to_string(std::random_device{}()));
        std::filesystem::create_directories(dir); const auto stem = (dir / "clusters").string();
        auto opts = options(stem); opts["-derepIdentity"] = "100"; opts["-derepPrefix"] = "0";
        opts["-derepStoreDiffs"] = mode == "exact_diff" ? "1" : "0";
        opts["-derep_format"] = "fa"; opts["-o_dereplicate"] = stem + ".fa";
        Filters filters(&opts); filters.SampleID = {"sample"};
        filters.collectStatistics[0]->BarcodeDetected.assign(1, -1);
        filters.collectStatistics[0]->BarcodeDetectedFail.assign(1, 0);
        Dereplicate derep(&opts, &filters); derep.BCnamesAdding(&filters);
        std::atomic<size_t> next{0}; std::vector<std::thread> workers;
        const auto localStart = Clock::now();
        for (size_t t = 0; t < threads; ++t) {
            workers.emplace_back([&, t] {
                for (size_t begin; (begin = next.fetch_add(256)) < n;) {
                    const size_t end = std::min(n, begin + 256);
                    for (size_t i = begin; i < end; ++i) {
                        auto dna = readDNA(input[i], quality); dna->setHeader("read_" + std::to_string(i));
                        derep.addDNA(dna, nullptr); ++workerReads[t];
                    }
                }
            });
        }
        for (auto& worker : workers) { worker.join(); } localSeconds = elapsed(localStart);
        const auto mergeStart = Clock::now(); derep.writeDereplDNA(&filters, ""); mergeSeconds = elapsed(mergeStart);
        std::ifstream output(stem + ".fa"); std::string line; uint64_t counted = 0;
        while (std::getline(output, line)) {
            if (line.starts_with('>')) { ++finalGroups; counted += std::stoull(line.substr(line.rfind(";size=") + 6)); }
        }
        check(counted == n, "Production exact benchmark lost observations"); localGroups = finalGroups;
        std::filesystem::remove_all(dir);
    } else {
        const auto routeStart = Clock::now();
        robin_hood::unordered_map<uint64_t, size_t> frequencies;
        if (mode == "rare20" || mode == "supported20") {
            for (size_t i = 0; i < std::min(n, size_t(4096)); ++i) { for (auto seed : seeds(input[i])) { ++frequencies[seed]; } }
        }
        std::vector<std::vector<size_t>> ids(buckets);
        std::vector<std::vector<size_t>> templateBuckets(50, std::vector<size_t>(buckets));
        for (size_t i = 0; i < n; ++i) {
            uint64_t signature = std::numeric_limits<uint64_t>::max();
            size_t rarest = std::numeric_limits<size_t>::max();
            if (mode == "global") { route[i] = 0; }
            else if (mode == "roundrobin") { route[i] = i % buckets; }
            else if (mode == "fullhash") { route[i] = robin_hood::hash<std::string>{}(input[i]) % buckets; }
            else {
                for (uint64_t seed : seeds(input[i])) {
                    const uint64_t rank = mix(mode == "min8" ? ((seed & ((1ULL << 40) - 1)) >> 24) : seed);
                    if (mode == "rare20" || mode == "supported20") {
                        const auto found = frequencies.find(seed); const size_t freq = found == frequencies.end() ? 0 : found->second;
                        if (mode == "supported20" && freq < 3) { continue; }
                        if (freq < rarest || (freq == rarest && rank < signature)) { rarest = freq; signature = rank; }
                    } else { signature = std::min(signature, rank); }
                }
                // Reads with no eligible seed share one fallback bucket.
                route[i] = mix(signature) % buckets;
            }
            ids[route[i]].push_back(i); ++templateBuckets[source[i]][route[i]];
        }
        uint64_t same = 0, pairs = 0;
        for (const auto& counts : templateBuckets) {
            uint64_t total = 0;
            for (auto count : counts) { total += count; same += uint64_t(count) * (count ? count - 1 : 0); }
            pairs += total * (total ? total - 1 : 0);
        }
        pairTogether = pairs ? double(same) / pairs : 0;
        for (const auto& bucket : ids) { maxBucketFraction = std::max(maxBucketFraction, double(bucket.size()) / n); }
        routeSeconds = elapsed(routeStart);
        std::vector<std::unique_ptr<CoarseDereplication>> indexes;
        for (size_t b = 0; b < buckets; ++b) { indexes.push_back(std::make_unique<CoarseDereplication>(threshold, 4, false)); }
        std::vector<std::vector<size_t>> firstRead(buckets);
        std::mutex globalMutex; std::atomic<size_t> next{0};
        std::vector<double> waits(threads, 0); std::vector<std::thread> workers;
        const auto localStart = Clock::now();
        for (size_t t = 0; t < threads; ++t) {
            workers.emplace_back([&, t] {
                auto insert = [&](size_t i, size_t bucket, const std::shared_ptr<DNA>& dna) {
                    const auto before = indexes[bucket]->groups().size();
                    assignment[i] = indexes[bucket]->add(input[i], input[i], dna, nullptr, nullptr, 0, true);
                    if (indexes[bucket]->groups().size() != before) { firstRead[bucket].push_back(i); }
                    ++workerReads[t];
                };
                if (mode == "global") {
                    for (size_t begin; (begin = next.fetch_add(256)) < n;) {
                        const size_t end = std::min(n, begin + 256);
                        std::vector<std::shared_ptr<DNA>> reads; reads.reserve(end - begin);
                        for (size_t i = begin; i < end; ++i) { reads.push_back(readDNA(input[i], quality)); }
                        const auto startWait = Clock::now(); std::lock_guard<std::mutex> lock(globalMutex);
                        waits[t] += elapsed(startWait);
                        for (size_t i = begin; i < end; ++i) { insert(i, 0, reads[i - begin]); }
                    }
                } else {
                    for (size_t b; (b = next.fetch_add(1)) < buckets;) {
                        for (auto i : ids[b]) { insert(i, b, readDNA(input[i], quality)); }
                    }
                }
            });
        }
        for (auto& worker : workers) { worker.join(); }
        localSeconds = elapsed(localStart); waitSeconds = std::accumulate(waits.begin(), waits.end(), 0.);
        for (const auto& index : indexes) { localGroups += index->groups().size(); }
        const auto mergeStart = Clock::now();
        struct Child { size_t first, bucket, id; };
        std::vector<Child> children;
        std::vector<std::vector<size_t>> parent(buckets);
        for (size_t b = 0; b < buckets; ++b) {
            parent[b].resize(indexes[b]->groups().size());
            for (size_t id = 0; id < firstRead[b].size(); ++id) { children.push_back({firstRead[b][id], b, id}); }
        }
        std::sort(children.begin(), children.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        CoarseDereplication merged(threshold, 4, false);
        for (const auto& child : children) {
            auto anchor = indexes[child.bucket]->groups()[child.id]->key();
            parent[child.bucket][child.id] = merged.add(anchor, anchor, readDNA(std::string(anchor), quality), nullptr, nullptr, 0, true);
        }
        finalGroups = merged.groups().size(); mergeSeconds = elapsed(mergeStart);
        if (mode == "global") { mergeSeconds = 0; } // No merge is required by the shared baseline.
        // Validation is outside the reported route/local/merge timings.
        std::vector<uint64_t> counts(finalGroups, 0);
        for (size_t i = 0; i < n; ++i) {
            const size_t target = parent[route[i]][assignment[i]]; ++counts[target];
            const auto anchor = merged.groups()[target]->key(); size_t errors = 0;
            const size_t overlap = std::min(anchor.size(), input[i].size());
            for (size_t p = 0; p < overlap; ++p) { errors += anchor[p] != input[i][p]; }
            violations += errors > merged.mismatchBudget(overlap);
        }
        check(std::accumulate(counts.begin(), counts.end(), uint64_t(0)) == n, "Partition accounting lost observations");
        if (mode == "global") { check(violations == 0, "Shared baseline violated its anchor threshold"); }
    }
    long peakKiB = 0; std::ifstream status("/proc/self/status"); std::string line;
    while (std::getline(status, line)) { if (line.starts_with("VmHWM:")) { peakKiB = std::stol(line.substr(6)); } }
    const double maxWorkerFraction = double(*std::max_element(workerReads.begin(), workerReads.end())) / n;
    std::cout << "mode,threads,buckets,reads,length,quality,identity,route_seconds,local_seconds,merge_or_output_seconds,total_seconds,mutex_wait_thread_seconds,local_clusters,final_clusters,outside_final_anchor,max_bucket_fraction,max_worker_fraction,same_template_pair_bucket_fraction,peak_rss_kib\n"
        << mode << ',' << threads << ',' << buckets << ',' << n << ',' << length << ',' << quality << ',' << identity << ','
        << std::setprecision(8) << routeSeconds << ',' << localSeconds << ',' << mergeSeconds << ',' << routeSeconds + localSeconds + mergeSeconds << ','
        << waitSeconds << ',' << localGroups << ',' << finalGroups << ',' << violations << ',' << maxBucketFraction << ',' << maxWorkerFraction << ',' << pairTogether << ',' << peakKiB << '\n';
}
