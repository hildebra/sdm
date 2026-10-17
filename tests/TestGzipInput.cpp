#include "Common.h"
#include "InputFileStream.h"
#include "InputStream.h"
#include "OutputFileStream.h"
#include "ThreadPool.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

// DNA_trans comes from DNAconsts.o (linked by the Makefile test rule).

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

class CerrCapture {
public:
    CerrCapture() : previous_(std::cerr.rdbuf(buffer_.rdbuf())) {}
    ~CerrCapture() { std::cerr.rdbuf(previous_); }
    std::string str() const { return buffer_.str(); }

private:
    std::ostringstream buffer_;
    std::streambuf* previous_;
};

std::string makeFixture(size_t size) {
    std::string data(size, '\0');
    uint32_t state = 0x9e3779b9u;
    for (char& value : data) {
        state = state * 1664525u + 1013904223u;
        value = static_cast<char>(state >> 24);
    }
    return data;
}

std::vector<uint8_t> readBytes(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    require(input.is_open(), "could not open " + path.string());
    input.seekg(0, std::ios::end);
    const std::streamoff size = input.tellg();
    require(size >= 0, "could not size " + path.string());
    input.seekg(0, std::ios::beg);
    std::vector<uint8_t> bytes(static_cast<size_t>(size));
    if (!bytes.empty()) {
        input.read(reinterpret_cast<char*>(bytes.data()), size);
        require(input.gcount() == size, "could not read " + path.string());
    }
    return bytes;
}

void writeBytes(const fs::path& path, const std::vector<uint8_t>& bytes) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(output.is_open(), "could not create " + path.string());
    if (!bytes.empty()) {
        output.write(
            reinterpret_cast<const char*>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
    }
    require(output.good(), "could not write " + path.string());
}

void writePlain(const fs::path& path, const std::string& data) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    require(output.is_open(), "could not create " + path.string());
    output.write(data.data(), static_cast<std::streamsize>(data.size()));
    require(output.good(), "could not write " + path.string());
}

void writeGzip(const fs::path& path, const std::string& data) {
    gzFile output = gzopen(path.string().c_str(), "wb");
    require(output != nullptr, "gzopen failed for " + path.string());
    const int written = gzwrite(
        output, data.data(), static_cast<unsigned>(data.size()));
    const int closeResult = gzclose(output);
    require(written == static_cast<int>(data.size()), "gzwrite failed");
    require(closeResult == Z_OK, "gzclose failed");
}

void writeBgzf(const fs::path& path, const std::string& data) {
    auto output = openOutputFile(
        path.string(), std::ios::out | std::ios::trunc, 4096, true);
    output->write(data.data(), static_cast<std::streamsize>(data.size()));
    finalizeOutputFile(output, path.string());
}

std::string readAll(std::unique_ptr<std::istream> input) {
    std::string result;
    std::array<char, 32768> buffer{};
    for (;;) {
        input->read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize count = input->gcount();
        if (count > 0) {
            result.append(buffer.data(), static_cast<size_t>(count));
        }
        if (input->eof()) {
            break;
        }
        if (!input->good()) {
            throw std::runtime_error("input stream failed without reaching EOF");
        }
    }
    return result;
}

void testOrdinaryGzip(const fs::path& directory) {
    const std::string data = makeFixture(150000);
    const fs::path path = directory / "ordinary.gzip";
    writeGzip(path, data);
    require(readAll(openInputFile(path.string())) == data,
        "ordinary gzip did not round-trip");
    require(!inputFileUsesParallelBgzf(*openInputFile(path.string())),
        "ordinary gzip without an index should stay sequential");
}

void testTruncatedGzipWarns(const fs::path& directory) {
    const std::string data = makeFixture(190000);
    const fs::path path = directory / "truncated.gz";
    writeGzip(path, data);
    std::vector<uint8_t> bytes = readBytes(path);
    require(bytes.size() > 8, "gzip fixture is too small");
    bytes.resize(bytes.size() - 8);
    writeBytes(path, bytes);

    CerrCapture capture;
    const std::string recovered = readAll(openInputFile(path.string()));
    require(recovered == data, "truncated gzip did not return recovered data");
    require(capture.str().find("truncated gzip input") != std::string::npos,
        "truncated gzip did not produce a warning");
}

void testCorruptGzipStillFails(const fs::path& directory) {
    const std::string data = makeFixture(80000);
    const fs::path path = directory / "bad-crc.gz";
    writeGzip(path, data);
    std::vector<uint8_t> bytes = readBytes(path);
    require(bytes.size() > 8, "gzip fixture is too small");
    bytes[bytes.size() - 8] ^= 0x80;
    writeBytes(path, bytes);

    bool failed = false;
    try {
        (void)readAll(openInputFile(path.string()));
    }
    catch (const std::exception&) {
        failed = true;
    }
    require(failed, "CRC-corrupt gzip should remain a hard data error");
}

void testIndexedParallelBgzf(const fs::path& directory) {
    const std::string data = makeFixture(620000);
    const fs::path path = directory / "indexed.fastq.gz";
    writeBgzf(path, data);
    require(fs::exists(path.string() + ".gzi"), "BGZF writer did not create .gzi");

    auto parallel = openInputFile(path.string(), std::ios::in, 4096, true);
    require(inputFileUsesParallelBgzf(*parallel),
        "valid indexed BGZF did not select the parallel reader");
    require(readAll(std::move(parallel)) == data,
        "parallel indexed BGZF did not preserve block order");

    auto sequential = openInputFile(path.string(), std::ios::in, 4096, false);
    require(!inputFileUsesParallelBgzf(*sequential),
        "parallel=false should select sequential gzip");
    require(readAll(std::move(sequential)) == data,
        "sequential BGZF fallback did not round-trip");

    writeBytes(path.string() + ".gzi", std::vector<uint8_t>{
        1, 0, 0, 0, 0, 0, 0, 0
    });
    CerrCapture capture;
    auto fallback = openInputFile(path.string(), std::ios::in, 4096, true);
    require(!inputFileUsesParallelBgzf(*fallback),
        "invalid .gzi should fall back to sequential gzip");
    require(readAll(std::move(fallback)) == data,
        "invalid-index fallback did not round-trip");
    require(capture.str().find("falling back to sequential gzip") != std::string::npos,
        "invalid .gzi fallback did not produce a warning");
}

void testMissingBgzfEofWarns(const fs::path& directory) {
    const std::string data = makeFixture(310000);
    const fs::path path = directory / "missing-eof.gz";
    writeBgzf(path, data);
    std::vector<uint8_t> bytes = readBytes(path);
    require(bytes.size() > 28, "BGZF fixture is too small");
    bytes.resize(bytes.size() - 28);
    writeBytes(path, bytes);

    CerrCapture capture;
    auto input = openInputFile(path.string(), std::ios::in, 4096, true);
    require(inputFileUsesParallelBgzf(*input),
        "complete indexed blocks should remain parallel without BGZF EOF");
    require(readAll(std::move(input)) == data,
        "BGZF without EOF did not return recovered data");
    require(capture.str().find("EOF marker is missing") != std::string::npos,
        "missing BGZF EOF did not produce a warning");
}

void testExactBufferEof(const fs::path& directory) {
    const fs::path path = directory / "exact-buffer.txt";
    const std::string data = "0123456789abcdef";
    writePlain(path, data);
    ifbufstream input(path.string(), data.size(), false, false);
    std::string line;
    require(input.getlines(line), "exact-buffer input lost its final line");
    require(line == data, "exact-buffer input changed its final line");
    require(input.eof(), "exact-buffer input did not mark EOF after refill");
    require(!input.getlines(line), "exact-buffer input returned a phantom line");
}

void testGzipSuffixDetection(const fs::path& directory) {
    const fs::path path = directory / "reads.fastq.gzip";
    writeGzip(path, "@r1\nACGT\n+\nIIII\n");
    require(detectSeqFmt(path.string()) == "-i_fastq",
        ".fastq.gzip suffix was not recognized as FASTQ");
}

} // namespace

int main() {
    try {
        ThreadPool::configure_gzip(4);
        setBgzfIndexOutput(true); // the indexed-input tests need the writer's .gzi
        const fs::path directory =
            fs::temp_directory_path() / "sdm_gzip_input_tests";
        fs::create_directories(directory);

        testOrdinaryGzip(directory);
        testTruncatedGzipWarns(directory);
        testCorruptGzipStillFails(directory);
        testIndexedParallelBgzf(directory);
        testMissingBgzfEofWarns(directory);
        testExactBufferEof(directory);
        testGzipSuffixDetection(directory);

        std::cout << "Gzip input tests passed\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "Gzip input test failed: " << error.what() << '\n';
        return 1;
    }
}
