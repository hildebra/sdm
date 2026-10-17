#include "Common.h"
#include "OutputFileStream.h"
#include "ThreadPool.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <zlib.h>

// DNA_trans comes from DNAconsts.o (linked by the Makefile test rule).

namespace {

constexpr std::array<uint8_t, 28> kBgzfEof = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

uint16_t loadLe16(const uint8_t* data) {
    return static_cast<uint16_t>(data[0]) |
        static_cast<uint16_t>(static_cast<uint16_t>(data[1]) << 8);
}

uint32_t loadLe32(const uint8_t* data) {
    return static_cast<uint32_t>(data[0]) |
        (static_cast<uint32_t>(data[1]) << 8) |
        (static_cast<uint32_t>(data[2]) << 16) |
        (static_cast<uint32_t>(data[3]) << 24);
}

uint64_t loadLe64(const uint8_t* data) {
    uint64_t value = 0;
    for (size_t i = 0; i < 8; ++i) {
        value |= static_cast<uint64_t>(data[i]) << (i * 8);
    }
    return value;
}

std::vector<uint8_t> readBytes(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    require(input.is_open(), "could not open " + path.string());
    input.seekg(0, std::ios::end);
    const std::streamoff size = input.tellg();
    require(size >= 0, "could not size " + path.string());
    input.seekg(0, std::ios::beg);
    std::vector<uint8_t> result(static_cast<size_t>(size));
    if (!result.empty()) {
        input.read(reinterpret_cast<char*>(result.data()), size);
        require(input.good(), "could not read " + path.string());
    }
    return result;
}

std::string decompressGzip(const std::filesystem::path& path) {
    gzFile input = gzopen(path.string().c_str(), "rb");
    require(input != nullptr, "zlib could not open " + path.string());
    std::string result;
    std::array<char, 32768> buffer{};
    for (;;) {
        const int count = gzread(input, buffer.data(), buffer.size());
        if (count > 0) {
            result.append(buffer.data(), static_cast<size_t>(count));
            continue;
        }
        if (count == 0) {
            break;
        }
        int errorNumber = Z_OK;
        const char* error = gzerror(input, &errorNumber);
        gzclose(input);
        throw std::runtime_error(
            "zlib decompression failed: " + std::string(error ? error : "unknown"));
    }
    require(gzclose(input) == Z_OK, "zlib close failed for " + path.string());
    return result;
}

void writeOutput(
    const std::filesystem::path& path,
    const std::string& data,
    std::ios_base::openmode mode,
    bool parallel = true) {
    auto output = openOutputFile(path.string(), std::ios::out | mode, 4096, parallel);
    output->write(data.data(), static_cast<std::streamsize>(data.size()));
    finalizeOutputFile(output, path.string());
    output.reset();
}

void validateBgzfAndGzi(
    const std::filesystem::path& path,
    uint64_t expectedUncompressedSize) {
    const std::vector<uint8_t> compressed = readBytes(path);
    require(!compressed.empty(), "BGZF output must not be empty");

    std::vector<std::pair<uint64_t, uint64_t>> expectedIndex;
    uint64_t compressedOffset = 0;
    uint64_t uncompressedOffset = 0;
    size_t dataBlocks = 0;
    bool sawEof = false;
    while (compressedOffset < compressed.size()) {
        require(compressedOffset + 18 <= compressed.size(), "truncated BGZF header");
        const uint8_t* block = compressed.data() + compressedOffset;
        require(block[0] == 0x1f && block[1] == 0x8b && block[2] == 0x08,
            "invalid gzip header in BGZF block");
        require((block[3] & 0x04) != 0 && block[12] == 0x42 && block[13] == 0x43,
            "missing BGZF BC extra field");
        const uint64_t blockSize = static_cast<uint64_t>(loadLe16(block + 16)) + 1;
        require(blockSize >= 26 && blockSize <= 65536, "invalid BGZF BSIZE");
        require(compressedOffset + blockSize <= compressed.size(), "truncated BGZF block");

        const bool isEof = blockSize == kBgzfEof.size() &&
            std::equal(kBgzfEof.begin(), kBgzfEof.end(), block);
        if (isEof) {
            sawEof = true;
            require(compressedOffset + blockSize == compressed.size(),
                "BGZF EOF marker must be the final block");
            compressedOffset += blockSize;
            break;
        }

        if (uncompressedOffset != 0) {
            expectedIndex.emplace_back(compressedOffset, uncompressedOffset);
        }
        const uint32_t blockUncompressed = loadLe32(
            block + static_cast<size_t>(blockSize) - 4);
        require(blockUncompressed > 0, "data BGZF block has zero ISIZE");
        uncompressedOffset += blockUncompressed;
        compressedOffset += blockSize;
        ++dataBlocks;
    }

    require(sawEof, "missing canonical BGZF EOF marker");
    require(compressedOffset == compressed.size(), "unparsed bytes follow BGZF data");
    require(uncompressedOffset == expectedUncompressedSize,
        "BGZF ISIZE totals do not match source data");
    require(dataBlocks > 1 || expectedUncompressedSize == 0,
        "large BGZF fixture did not produce multiple blocks");

    const std::vector<uint8_t> index = readBytes(path.string() + ".gzi");
    require(index.size() >= 8, "truncated .gzi header");
    const uint64_t count = loadLe64(index.data());
    require(count == expectedIndex.size(), "unexpected .gzi entry count");
    require(index.size() == 8 + count * 16, "unexpected .gzi file size");
    for (size_t i = 0; i < expectedIndex.size(); ++i) {
        const uint8_t* pair = index.data() + 8 + i * 16;
        require(loadLe64(pair) == expectedIndex[i].first,
            "incorrect compressed offset in .gzi");
        require(loadLe64(pair + 8) == expectedIndex[i].second,
            "incorrect uncompressed offset in .gzi");
    }
}

std::string makeFixture(size_t size) {
    std::string data(size, '\0');
    uint32_t state = 0x12345678u;
    for (char& value : data) {
        state = state * 1664525u + 1013904223u;
        value = static_cast<char>(state >> 24);
    }
    return data;
}

void expectFinalizeFailure(const std::filesystem::path& path) {
    bool failed = false;
    try {
        auto output = openOutputFile(path.string(), std::ios::out | std::ios::trunc, 64, true);
        *output << std::string(4096, 'X');
        finalizeOutputFile(output, path.string());
    }
    catch (const std::exception&) {
        failed = true;
    }
    require(failed, "delayed write failure was not reported for " + path.string());
}

} // namespace

int main() {
    try {
        ThreadPool::configure_gzip(4);
        const std::filesystem::path directory =
            std::filesystem::temp_directory_path() / "sdm_bgzf_output_tests";
        std::filesystem::create_directories(directory);

        // The index is opt-in (-o_gzi 1): by default none is written and a stale one is removed.
        require(!bgzfIndexOutput(), ".gzi output must be off by default");
        const std::filesystem::path noIndexPath = directory / "noindex.gz";
        {
            std::ofstream stale(noIndexPath.string() + ".gzi", std::ios::binary | std::ios::trunc);
            stale << "stale";
        }
        writeOutput(noIndexPath, makeFixture(100000), std::ios::trunc, true);
        require(!std::filesystem::exists(noIndexPath.string() + ".gzi"),
            "a .gzi was written (or a stale one kept) although the index is disabled");
        setBgzfIndexOutput(true);

        require(isGZfile("reads.gz"), ".gz suffix should select BGZF output");
        require(isGZfile("reads.GZIP"), ".gzip suffix should be case-insensitive");
        require(!isGZfile("reads.gz.tmp"), "gzip suffix matching must be anchored");

        const std::string fixture = makeFixture(420000);
        const std::filesystem::path parallelPath = directory / "parallel.GZIP";
        writeOutput(parallelPath, fixture, std::ios::trunc, true);
        require(decompressGzip(parallelPath) == fixture,
            "parallel BGZF output did not round-trip");
        validateBgzfAndGzi(parallelPath, fixture.size());

        const std::filesystem::path appendPath = directory / "append.gz";
        const std::string first = fixture.substr(0, 135000);
        const std::string second = fixture.substr(135000);
        writeOutput(appendPath, first, std::ios::trunc, true);
        writeOutput(appendPath, second, std::ios::app, true);
        require(decompressGzip(appendPath) == first + second,
            "appended BGZF output did not round-trip");
        validateBgzfAndGzi(appendPath, fixture.size());

        const std::filesystem::path sequentialPath = directory / "sequential.gz";
        writeOutput(sequentialPath, fixture, std::ios::trunc, false);
        require(decompressGzip(sequentialPath) == fixture,
            "sequential BGZF output did not round-trip");
        validateBgzfAndGzi(sequentialPath, fixture.size());

        const std::filesystem::path emptyPath = directory / "empty.gz";
        writeOutput(emptyPath, "", std::ios::trunc, true);
        require(readBytes(emptyPath).size() == kBgzfEof.size(),
            "empty BGZF output should contain only the canonical EOF marker");
        require(decompressGzip(emptyPath).empty(), "empty BGZF output was not readable");
        validateBgzfAndGzi(emptyPath, 0);

#if defined(__linux__)
        expectFinalizeFailure("/dev/full");
        const std::filesystem::path fullGzipPath = directory / "full.gz";
        std::error_code removeError;
        std::filesystem::remove(fullGzipPath, removeError);
        std::filesystem::create_symlink("/dev/full", fullGzipPath);
        expectFinalizeFailure(fullGzipPath);
#endif

        std::cout << "BGZF output tests passed\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "BGZF output test failed: " << error.what() << '\n';
        return 1;
    }
}
