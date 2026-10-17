// Isolated representation experiment; production retains signed byte qualities.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct PackedSix {
    explicit PackedSix(size_t count) : bytes((count * 6 + 7) / 8 + 1, 0) {}
    unsigned get(size_t i) const {
        const size_t bit = i * 6, word = bit / 8;
        return ((unsigned(bytes[word]) | (unsigned(bytes[word + 1]) << 8)) >> (bit % 8)) & 63;
    }
    void set(size_t i, unsigned value) {
        const size_t bit = i * 6, word = bit / 8; const unsigned shift = bit % 8;
        unsigned pair = unsigned(bytes[word]) | (unsigned(bytes[word + 1]) << 8);
        pair = (pair & ~(63U << shift)) | (value << shift);
        bytes[word] = pair & 255; bytes[word + 1] = pair >> 8;
    }
    std::vector<uint8_t> bytes;
};
template<class Byte> struct Plain {
    explicit Plain(size_t count) : bytes(count) {}
    unsigned get(size_t i) const { return bytes[i]; }
    void set(size_t i, unsigned value) { bytes[i] = value; }
    std::vector<Byte> bytes;
};
long peakRSS() {
    std::ifstream input("/proc/self/status"); std::string line;
    while (std::getline(input, line)) { if (line.starts_with("VmHWM:")) { return std::stol(line.substr(6)); } }
    return 0;
}
template<class Storage> void run(const std::string& mode, size_t count, size_t passes) {
    const auto start = std::chrono::steady_clock::now();
    Storage qualities(count);
    for (size_t i = 0; i < count; ++i) { qualities.set(i, (i * 17 + i / 11) % 61); }
    const auto initialized = std::chrono::steady_clock::now();
    uint64_t sum = 0;
    for (size_t pass = 0; pass < passes; ++pass) {
        // Representative quality scans and mutation, in 150-base reads.
        for (size_t i = 0; i < count; ++i) { sum += qualities.get(i); }
        for (size_t start = 0; start + 150 <= count; start += 150) {
            for (size_t i = 0; i < 75; ++i) {
                const auto left = qualities.get(start + i), right = qualities.get(start + 149 - i);
                qualities.set(start + i, right); qualities.set(start + 149 - i, left);
            }
        }
    }
    const auto finished = std::chrono::steady_clock::now();
    // Include final content so all work remains observable and comparable.
    for (size_t i = 0; i < count; ++i) { sum += qualities.get(i) * (i % 31); }
    std::cout << "mode,bases,passes,storage_bytes,initialization_seconds,scan_reverse_seconds,peak_rss_kib,checksum\n"
        << mode << ',' << count << ',' << passes << ',' << qualities.bytes.size() << ',' << std::setprecision(9)
        << std::chrono::duration<double>(initialized-start).count() << ','
        << std::chrono::duration<double>(finished-initialized).count() << ',' << peakRSS() << ',' << sum << '\n';
}
int main(int argc, char** argv) {
    if (argc != 4) { return 2; }
    const std::string mode=argv[1];const size_t count=std::stoull(argv[2]),passes=std::stoull(argv[3]);
    if (mode=="signed_byte") { run<Plain<int8_t>>(mode,count,passes); }
    else if (mode=="unsigned_byte") { run<Plain<uint8_t>>(mode,count,passes); }
    else if (mode=="packed_six") { run<PackedSix>(mode,count,passes); }
    else { return 2; }
}
