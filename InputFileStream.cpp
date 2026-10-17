#include "InputFileStream.h"

#include "ThreadPool.h"
#include "ThreadedGzStream.h"

#include <isa-l.h>
#include <sys/stat.h>
#include <zstd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace {

// The threads a compressed input file gets (ThreadedGzStreambuf::set_threads) on top of
// its own inflating thread, from the gzip worker capacity (-threads): workers for BGZF
// blocks, zstd frames and gzip members, which wait (not spin) while sdm takes the bytes
// slower than they decompress them (on a 1.7 GB FASTQ, 1-8 workers: BGZF 1.2-3.5 GB/s,
// zstd frames 1.5-5.2 GB/s, peak RSS 9-35 MB). Single-member gzip stays on its file's own
// thread (ParallelGzip off, unless SDM_INFLATE_THREADS asks for it): one ISA-L thread
// inflated 1.8 GB/s of it, more than sdm takes from a file, and pieces on 2-4 threads
// added only 4-25% for 150-300 MB more per file (.claude/archive/gz-port-2026-10).
struct CompressedInputThreads {
    unsigned inflate = 1;
    unsigned gzip = 1;
};

CompressedInputThreads compressedInputThreads(bool requested) {
    if (!requested) {
        return {};
    }
    const size_t hw = std::thread::hardware_concurrency();
    const size_t threads = ThreadPool::configured_gzip_thread_count();
    if ((hw > 0 && hw < 2) || threads < 2) {
        return {};
    }
    CompressedInputThreads result;
    result.inflate = static_cast<unsigned>(std::clamp<size_t>(threads / 2, 2, 8));
    return result;
}

void warnInput(const std::string& filename, const std::string& message) {
    static std::mutex warningMutex;
    std::lock_guard<std::mutex> lock(warningMutex);
    std::cerr << "Warning: " << message << ": " << filename << '\n';
}

// Whether a read error of ThreadedGzStreambuf only ends the file early, everything read
// before it being the file's (as sdm 3.65 read such files, with a warning): the file is
// truncated (a member, block or frame cut off, a BGZF file without its end-of-file
// block), or what follows a whole gzip member or BGZF block is no member or block (a
// damaged header, other data). Corrupt compressed data (a wrong CRC-32 or length,
// invalid deflate or zstd data) remains an error: bytes already read may be wrong.
bool endsEarly(const std::string& error) {
    return error.find("(truncated file?)") != std::string::npos ||
        error.find("is no gzip member") != std::string::npos ||
        error.find("no BGZF block at byte") != std::string::npos;
}

// ThreadedGzStreambuf, whose read error is a warning where the file only ends early
// (endsEarly: reading goes on with what was read before the damage) and is thrown where
// reading reaches it otherwise, rather than read as the end of the file.
class CompressedInputStreambuf final : public sdm::ThreadedGzStreambuf {
public:
    explicit CompressedInputStreambuf(std::string filename)
        : filename_(std::move(filename)) {
    }

protected:
    int_type underflow() override {
        const int_type next = sdm::ThreadedGzStreambuf::underflow();
        if (traits_type::eq_int_type(next, traits_type::eof()) && read_failed()) {
            const std::string error = read_error_message();
            if (!endsEarly(error)) {
                throw std::runtime_error("Could not read input file " + filename_ + ": " + error);
            }
            if (!warned_) {
                warned_ = true;
                warnInput(filename_,
                    "damaged compressed input, continuing with the data read before the damage (" + error + ")");
            }
        }
        return next;
    }

private:
    std::string filename_;
    bool warned_ = false;
};

class CompressedInputStream final : public std::istream {
public:
    CompressedInputStream(const std::string& filename, CompressedInputThreads threads)
        : std::istream(nullptr), buffer_(filename) {
        rdbuf(&buffer_); // also clears the badbit of the null buffer
        buffer_.set_threads(threads.inflate, threads.gzip);
        if (!buffer_.open(filename.c_str())) {
            throw std::runtime_error(
                "Could not open input file: " + filename + " (" + buffer_.open_error() + ")");
        }
        exceptions(std::ios::badbit);
    }

    ~CompressedInputStream() override {
        exceptions(std::ios::goodbit);
        buffer_.close();
    }

    const CompressedInputStreambuf& buffer() const { return buffer_; }

private:
    CompressedInputStreambuf buffer_;
};

// Whether a regular file begins as gzip (BGZF or not) or zstd (a frame or a skippable
// frame): only those are given to ThreadedGzStreambuf, plain files to std::ifstream.
bool startsCompressed(const std::string& filename) {
    std::ifstream input(filename, std::ios::binary);
    unsigned char head[4] = {};
    input.read(reinterpret_cast<char*>(head), sizeof(head));
    const std::streamsize got = input.gcount();
    if (got >= 2 && head[0] == 0x1f && head[1] == 0x8b) {
        return true;
    }
    return got == 4 &&
        ((head[0] == 0x28 && head[1] == 0xb5 && head[2] == 0x2f && head[3] == 0xfd) ||
         ((head[0] & 0xf0) == 0x50 && head[1] == 0x2a && head[2] == 0x4d && head[3] == 0x18));
}

} // namespace

std::unique_ptr<std::istream> openInputFile(
    const std::string& filename,
    std::ios_base::openmode mode,
    size_t bufferSize,
    bool parallelGzip) {
    (void)bufferSize; // the decompressing streams keep blocks of their own (ThreadedGzStreambuf::kBlockSize)
    struct stat info {};
    const bool regular = ::stat(filename.c_str(), &info) == 0 && S_ISREG(info.st_mode);
    // A pipe (FIFO, /dev/stdin, process substitution) cannot be peeked at: ThreadedGzStreambuf
    // tells its format from the bytes as it reads them, on one thread.
    if (!regular || startsCompressed(filename)) {
        return std::make_unique<CompressedInputStream>(
            filename, compressedInputThreads(parallelGzip));
    }
    auto input = std::make_unique<std::ifstream>(
        filename, mode | std::ios::in | std::ios::binary);
    if (!input->is_open()) {
        throw std::runtime_error("Could not open input file: " + filename);
    }
    return input;
}

bool inputFileUsesParallelBgzf(const std::istream& stream) {
    const auto* compressed = dynamic_cast<const CompressedInputStream*>(&stream);
    return compressed && compressed->buffer().opened_as_bgzf() &&
        compressed->buffer().inflate_threads() > 1;
}

std::string compressionLibraries() {
    return "gzip/BGZF: ISA-L " + std::to_string(ISAL_MAJOR_VERSION) + "." +
        std::to_string(ISAL_MINOR_VERSION) + "." + std::to_string(ISAL_PATCH_VERSION) +
        " (.gz output: BGZF); zstd: libzstd " + ZSTD_versionString();
}
