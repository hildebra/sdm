#include "OutputFileStream.h"

#include "Common.h"
#include "ThreadPool.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>
#include <limits>
#include <stdexcept>
#include <streambuf>
#include <utility>
#include <vector>

#ifdef _gzipread
#include <zlib.h>
#endif

namespace {

class FinalizableOutput {
public:
    virtual ~FinalizableOutput() = default;
    virtual void finalize() = 0;
};

#ifdef _gzipread

constexpr size_t kBgzfBlockSize = 0xff00;
constexpr size_t kBgzfMaxBlockSize = 0x10000;
constexpr size_t kBgzfHeaderSize = 18;
constexpr size_t kBgzfFooterSize = 8;

bool shouldUseParallelGzip(bool requested) {
    if (!requested) {
        return false;
    }
    const size_t hw = std::thread::hardware_concurrency();
    if (hw > 0 && hw < 2) {
        return false;
    }
    return ThreadPool::configured_gzip_thread_count() >= 2;
}

constexpr std::array<uint8_t, 18> kBgzfHeader = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00,
    0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x00, 0x00
};

constexpr std::array<uint8_t, 28> kBgzfEof = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

void storeLe16(uint8_t* dst, uint16_t value) {
    dst[0] = static_cast<uint8_t>(value);
    dst[1] = static_cast<uint8_t>(value >> 8);
}

void storeLe32(uint8_t* dst, uint32_t value) {
    dst[0] = static_cast<uint8_t>(value);
    dst[1] = static_cast<uint8_t>(value >> 8);
    dst[2] = static_cast<uint8_t>(value >> 16);
    dst[3] = static_cast<uint8_t>(value >> 24);
}

void storeLe64(uint8_t* dst, uint64_t value) {
    for (size_t i = 0; i < 8; ++i) {
        dst[i] = static_cast<uint8_t>(value >> (i * 8));
    }
}

uint16_t loadLe16(const uint8_t* src) {
    return static_cast<uint16_t>(src[0]) |
        static_cast<uint16_t>(static_cast<uint16_t>(src[1]) << 8);
}

uint32_t loadLe32(const uint8_t* src) {
    return static_cast<uint32_t>(src[0]) |
        (static_cast<uint32_t>(src[1]) << 8) |
        (static_cast<uint32_t>(src[2]) << 16) |
        (static_cast<uint32_t>(src[3]) << 24);
}

struct CompressedBlock {
    std::vector<uint8_t> bytes;
    uint32_t uncompressedSize = 0;
};

CompressedBlock makeStoredBgzfBlock(const std::vector<uint8_t>& input) {
    if (input.size() > std::numeric_limits<uint16_t>::max()) {
        throw std::runtime_error("BGZF stored block input exceeds 65535 bytes");
    }

    const size_t totalSize = kBgzfHeaderSize + 5 + input.size() + kBgzfFooterSize;
    if (totalSize > kBgzfMaxBlockSize) {
        throw std::runtime_error("BGZF stored block exceeds 65536 bytes");
    }

    CompressedBlock result;
    result.bytes.resize(totalSize);
    result.uncompressedSize = static_cast<uint32_t>(input.size());
    std::copy(kBgzfHeader.begin(), kBgzfHeader.end(), result.bytes.begin());
    storeLe16(result.bytes.data() + 16, static_cast<uint16_t>(totalSize - 1));

    uint8_t* body = result.bytes.data() + kBgzfHeaderSize;
    body[0] = 0x01; // BFINAL=1, BTYPE=00.
    const uint16_t length = static_cast<uint16_t>(input.size());
    storeLe16(body + 1, length);
    storeLe16(body + 3, static_cast<uint16_t>(~length));
    if (!input.empty()) {
        std::memcpy(body + 5, input.data(), input.size());
    }

    const uint32_t checksum = static_cast<uint32_t>(
        crc32(crc32(0L, Z_NULL, 0), input.data(), static_cast<uInt>(input.size())));
    storeLe32(result.bytes.data() + totalSize - 8, checksum);
    storeLe32(result.bytes.data() + totalSize - 4, result.uncompressedSize);
    return result;
}

CompressedBlock compressBgzfBlock(std::vector<uint8_t> input) {
    if (input.empty() || input.size() > kBgzfBlockSize) {
        throw std::runtime_error("invalid BGZF input block size");
    }

    CompressedBlock result;
    result.bytes.resize(kBgzfMaxBlockSize);
    result.uncompressedSize = static_cast<uint32_t>(input.size());

    z_stream stream{};
    stream.next_in = input.data();
    stream.avail_in = static_cast<uInt>(input.size());
    stream.next_out = result.bytes.data() + kBgzfHeaderSize;
    stream.avail_out = static_cast<uInt>(
        kBgzfMaxBlockSize - kBgzfHeaderSize - kBgzfFooterSize);

    const int initResult = deflateInit2(
        &stream, Z_BEST_SPEED, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY);
    if (initResult != Z_OK) {
        throw std::runtime_error(
            "deflateInit2 failed with code " + std::to_string(initResult));
    }

    const int deflateResult = deflate(&stream, Z_FINISH);
    const size_t compressedBodySize = static_cast<size_t>(stream.total_out);
    const bool useStoredBlock = stream.avail_out == 0;
    const int endResult = deflateEnd(&stream);

    if (useStoredBlock) {
        return makeStoredBgzfBlock(input);
    }
    if (deflateResult != Z_STREAM_END) {
        throw std::runtime_error(
            "deflate failed with code " + std::to_string(deflateResult));
    }
    if (endResult != Z_OK) {
        throw std::runtime_error(
            "deflateEnd failed with code " + std::to_string(endResult));
    }

    const size_t totalSize = kBgzfHeaderSize + compressedBodySize + kBgzfFooterSize;
    if (totalSize > kBgzfMaxBlockSize) {
        return makeStoredBgzfBlock(input);
    }

    std::copy(kBgzfHeader.begin(), kBgzfHeader.end(), result.bytes.begin());
    storeLe16(result.bytes.data() + 16, static_cast<uint16_t>(totalSize - 1));
    const uint32_t checksum = static_cast<uint32_t>(
        crc32(crc32(0L, Z_NULL, 0), input.data(), static_cast<uInt>(input.size())));
    storeLe32(result.bytes.data() + totalSize - 8, checksum);
    storeLe32(result.bytes.data() + totalSize - 4, result.uncompressedSize);
    result.bytes.resize(totalSize);
    return result;
}

class ParallelBgzfOstream final : public std::ostream,
                                  private std::streambuf,
                                  public FinalizableOutput {
public:
    ParallelBgzfOstream(
        const std::string& filename,
        std::ios_base::openmode mode,
        size_t requestedBufferSize,
        bool parallel)
        : std::ostream(this),
          filename_(filename),
          parallel_(parallel),
          maxPending_(parallel
              ? std::min<size_t>(8, std::max<size_t>(2,
                    ThreadPool::configured_gzip_thread_count() * 2))
              : 0) {
        (void)requestedBufferSize;
        buffer_.resize(kBgzfBlockSize);
        setp(reinterpret_cast<char*>(buffer_.data()),
            reinterpret_cast<char*>(buffer_.data() + buffer_.size()));

        prepareAppend(mode);
        file_.open(filename_, mode | std::ios::binary);
        if (!file_.is_open()) {
            throw std::runtime_error("Error opening output file: " + filename_);
        }
        file_.exceptions(std::ios::badbit | std::ios::failbit);
    }

    ~ParallelBgzfOstream() override {
        try {
            finalize();
        }
        catch (...) {
        }
    }

    void finalize() override {
        if (closed_) {
            rethrowFailure();
            return;
        }

        try {
            if (!failure_) {
                queueCurrentBuffer();
            }
            drainPending();
            if (!failure_) {
                file_.write(
                    reinterpret_cast<const char*>(kBgzfEof.data()),
                    static_cast<std::streamsize>(kBgzfEof.size()));
                file_.flush();
            }
        }
        catch (...) {
            rememberFailure(std::current_exception());
        }

        try {
            if (file_.is_open()) {
                file_.close();
            }
        }
        catch (...) {
            rememberFailure(std::current_exception());
        }

        if (!failure_) {
            try {
                writeGzi();
            }
            catch (...) {
                rememberFailure(std::current_exception());
            }
        }
        closed_ = true;
        rethrowFailure();
    }

protected:
    std::streambuf::int_type overflow(
        std::streambuf::int_type ch) override {
        using traits = std::streambuf::traits_type;
        if (!queueCurrentBuffer()) {
            return traits::eof();
        }
        if (!traits::eq_int_type(ch, traits::eof())) {
            *pptr() = traits::to_char_type(ch);
            pbump(1);
        }
        return traits::not_eof(ch);
    }

    std::streamsize xsputn(const char* data, std::streamsize count) override {
        if (count <= 0 || failure_ || closed_) {
            return 0;
        }

        std::streamsize written = 0;
        while (written < count) {
            if (pptr() == epptr() && !queueCurrentBuffer()) {
                break;
            }
            const std::streamsize available = epptr() - pptr();
            const std::streamsize amount = std::min(available, count - written);
            std::memcpy(pptr(), data + written, static_cast<size_t>(amount));
            pbump(static_cast<int>(amount));
            written += amount;
            if (pptr() == epptr() && !queueCurrentBuffer()) {
                break;
            }
        }
        return written;
    }

    int sync() override {
        if (closed_) {
            return failure_ ? -1 : 0;
        }
        if (!queueCurrentBuffer()) {
            return -1;
        }
        drainPending();
        return failure_ ? -1 : 0;
    }

private:
    struct IndexEntry {
        uint64_t compressedOffset;
        uint64_t uncompressedOffset;
    };

    bool queueCurrentBuffer() {
        if (failure_ || closed_) {
            return false;
        }
        const size_t size = static_cast<size_t>(pptr() - pbase());
        if (size == 0) {
            return true;
        }

        std::vector<uint8_t> input(size);
        std::memcpy(input.data(), buffer_.data(), size);
        setp(reinterpret_cast<char*>(buffer_.data()),
            reinterpret_cast<char*>(buffer_.data() + buffer_.size()));

        try {
            if (!parallel_) {
                writeBlock(compressBgzfBlock(std::move(input)));
                return !failure_;
            }
            pending_.emplace_back(ThreadPool::gzip_instance().submit(
                [block = std::move(input)]() mutable {
                    return compressBgzfBlock(std::move(block));
                }));
            if (pending_.size() >= maxPending_) {
                writeNextPending();
            }
        }
        catch (...) {
            rememberFailure(std::current_exception());
        }
        return !failure_;
    }

    void writeNextPending() {
        if (pending_.empty()) {
            return;
        }
        std::future<CompressedBlock> future = std::move(pending_.front());
        pending_.pop_front();
        try {
            CompressedBlock block = future.get();
            if (!failure_) {
                writeBlock(std::move(block));
            }
        }
        catch (...) {
            rememberFailure(std::current_exception());
        }
    }

    void drainPending() {
        while (!pending_.empty()) {
            writeNextPending();
        }
    }

    void writeBlock(CompressedBlock block) {
        if (block.bytes.empty() || block.bytes.size() > kBgzfMaxBlockSize ||
            block.uncompressedSize == 0) {
            throw std::runtime_error("invalid compressed BGZF block");
        }
        if (uncompressedOffset_ != 0) {
            index_.push_back({compressedOffset_, uncompressedOffset_});
        }
        file_.write(
            reinterpret_cast<const char*>(block.bytes.data()),
            static_cast<std::streamsize>(block.bytes.size()));
        compressedOffset_ += block.bytes.size();
        uncompressedOffset_ += block.uncompressedSize;
    }

    void prepareAppend(std::ios_base::openmode mode) {
        if ((mode & std::ios::app) == 0 || !std::filesystem::exists(filename_)) {
            return;
        }

        std::error_code sizeError;
        const uintmax_t fileSize = std::filesystem::file_size(filename_, sizeError);
        if (sizeError || fileSize == 0) {
            if (sizeError) {
                throw std::runtime_error(
                    "Could not inspect BGZF file for append: " + filename_);
            }
            return;
        }

        std::ifstream input(filename_, std::ios::binary);
        if (!input) {
            throw std::runtime_error("Could not inspect BGZF file for append: " + filename_);
        }

        uint64_t compressed = 0;
        uint64_t uncompressed = 0;
        uint64_t dataEnd = 0;
        bool foundEof = false;
        while (compressed < fileSize) {
            std::array<uint8_t, kBgzfHeaderSize> header{};
            input.read(reinterpret_cast<char*>(header.data()), header.size());
            if (input.gcount() != static_cast<std::streamsize>(header.size()) ||
                !std::equal(kBgzfHeader.begin(), kBgzfHeader.begin() + 16, header.begin())) {
                throw std::runtime_error(
                    "Cannot append: existing file is not standard BGZF: " + filename_);
            }
            const uint64_t blockSize = static_cast<uint64_t>(loadLe16(header.data() + 16)) + 1;
            if (blockSize < kBgzfHeaderSize + kBgzfFooterSize ||
                blockSize > kBgzfMaxBlockSize || compressed + blockSize > fileSize) {
                throw std::runtime_error(
                    "Cannot append: invalid BGZF block size in " + filename_);
            }

            std::vector<uint8_t> rest(static_cast<size_t>(blockSize - header.size()));
            input.read(reinterpret_cast<char*>(rest.data()), rest.size());
            if (input.gcount() != static_cast<std::streamsize>(rest.size())) {
                throw std::runtime_error(
                    "Cannot append: truncated BGZF block in " + filename_);
            }
            const uint32_t blockUncompressed = loadLe32(rest.data() + rest.size() - 4);
            const bool isCanonicalEof = blockSize == kBgzfEof.size() &&
                std::equal(
                    kBgzfEof.begin(),
                    kBgzfEof.begin() + static_cast<std::ptrdiff_t>(header.size()),
                    header.begin()) &&
                std::equal(
                    kBgzfEof.begin() + static_cast<std::ptrdiff_t>(header.size()),
                    kBgzfEof.end(), rest.begin());
            if (isCanonicalEof) {
                foundEof = true;
                dataEnd = compressed;
                compressed += blockSize;
                if (compressed != fileSize) {
                    throw std::runtime_error(
                        "Cannot append: data follows BGZF EOF marker in " + filename_);
                }
                break;
            }

            if (uncompressed != 0) {
                index_.push_back({compressed, uncompressed});
            }
            compressed += blockSize;
            uncompressed += blockUncompressed;
            dataEnd = compressed;
        }

        compressedOffset_ = dataEnd;
        uncompressedOffset_ = uncompressed;
        if (foundEof) {
            std::error_code resizeError;
            std::filesystem::resize_file(filename_, dataEnd, resizeError);
            if (resizeError) {
                throw std::runtime_error(
                    "Could not remove BGZF EOF marker before append: " + filename_);
            }
        }
    }

    void writeGzi() {
        const std::string indexFilename = filename_ + ".gzi";
        std::ofstream indexFile(indexFilename, std::ios::binary | std::ios::trunc);
        if (!indexFile.is_open()) {
            throw std::runtime_error("Could not open gzip index: " + indexFilename);
        }
        indexFile.exceptions(std::ios::badbit | std::ios::failbit);
        std::array<uint8_t, 8> encoded{};
        storeLe64(encoded.data(), static_cast<uint64_t>(index_.size()));
        indexFile.write(reinterpret_cast<const char*>(encoded.data()), encoded.size());
        for (const IndexEntry& entry : index_) {
            std::array<uint8_t, 16> pair{};
            storeLe64(pair.data(), entry.compressedOffset);
            storeLe64(pair.data() + 8, entry.uncompressedOffset);
            indexFile.write(reinterpret_cast<const char*>(pair.data()), pair.size());
        }
        indexFile.close();
    }

    void rememberFailure(std::exception_ptr failure) {
        if (!failure_) {
            failure_ = std::move(failure);
        }
    }

    void rethrowFailure() const {
        if (failure_) {
            std::rethrow_exception(failure_);
        }
    }

    std::string filename_;
    std::ofstream file_;
    std::vector<uint8_t> buffer_;
    std::deque<std::future<CompressedBlock>> pending_;
    std::vector<IndexEntry> index_;
    uint64_t compressedOffset_ = 0;
    uint64_t uncompressedOffset_ = 0;
    bool parallel_ = true;
    size_t maxPending_ = 2;
    bool closed_ = false;
    std::exception_ptr failure_;
};

#endif

} // namespace

std::unique_ptr<std::ostream> openOutputFile(
    const std::string& filename,
    std::ios_base::openmode mode,
    size_t bufferSize,
    bool parallelGzip) {
    if (isGZfile(filename)) {
#ifdef _gzipread
        const bool useParallelGzip = shouldUseParallelGzip(parallelGzip);
        return std::make_unique<ParallelBgzfOstream>(
            filename, mode, bufferSize, useParallelGzip);
#else
        (void)bufferSize;
        (void)parallelGzip;
        throw std::runtime_error("gzip output requested but support is not compiled: " + filename);
#endif
    }

    auto output = std::make_unique<std::ofstream>(filename, mode | std::ios::binary);
    if (!output->is_open()) {
        throw std::runtime_error("Error opening output file: " + filename);
    }
    return output;
}

void finalizeOutputFile(
    std::unique_ptr<std::ostream>& stream,
    const std::string& filename) {
    if (!stream) {
        return;
    }
    finalizeOutputFile(*stream, filename);
}

void finalizeOutputFile(std::ostream& stream, const std::string& filename) {
    if (auto* finalizable = dynamic_cast<FinalizableOutput*>(&stream)) {
        finalizable->finalize();
        return;
    }

    stream.flush();
    bool failed = !stream.good();
    if (auto* file = dynamic_cast<std::ofstream*>(&stream)) {
        file->close();
        failed = failed || file->fail();
    }
    if (failed) {
        throw std::runtime_error("Failed to finalize output file: " + filename);
    }
}

