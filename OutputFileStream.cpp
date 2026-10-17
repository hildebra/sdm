#include "OutputFileStream.h"

#include "Bgzf.h"
#include "Common.h"
#include "ThreadPool.h"

#include <zstd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <streambuf>
#include <thread>
#include <utility>
#include <vector>

namespace {

std::atomic<bool> g_writeBgzfIndex{false};

class FinalizableOutput {
public:
    virtual ~FinalizableOutput() = default;
    virtual void finalize() = 0;
};

namespace bgzf = sdm::bgzf;

constexpr size_t kBgzfBlockSize = bgzf::kBlockInput;
constexpr size_t kBgzfMaxBlockSize = bgzf::kMaxBlock;
constexpr size_t kBgzfHeaderSize = bgzf::kHeaderBytes;
constexpr size_t kBgzfFooterSize = bgzf::kFooterBytes;

// zstd output: independent frames of up to kZstdFrameInput content bytes each, every
// one with its content size and checksum (as pzstd writes them), so that readers can
// decompress them on several threads (ThreadedGzStreambuf) and any zstd reads the file.
constexpr size_t kZstdFrameInput = size_t{1} << 20;
constexpr int kZstdLevel = 3; // zstd's default

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
static_assert(sizeof(bgzf::kEof) == kBgzfEof.size(), "the BGZF end-of-file block");

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

    const uint32_t checksum = bgzf::Crc32(input.data(), input.size());
    storeLe32(result.bytes.data() + totalSize - 8, checksum);
    storeLe32(result.bytes.data() + totalSize - 4, result.uncompressedSize);
    return result;
}

// One BGZF block, deflated with ISA-L at level 1 (bgzf::Deflater, one per thread): its
// bytes depend on the input alone, not on the thread that compresses it, so the file
// is the same however many threads write it.
CompressedBlock compressBgzfBlock(std::vector<uint8_t> input) {
    if (input.empty() || input.size() > kBgzfBlockSize) {
        throw std::runtime_error("invalid BGZF input block size");
    }

    thread_local bgzf::Deflater deflater;
    CompressedBlock result;
    result.bytes.resize(kBgzfMaxBlockSize);
    result.uncompressedSize = static_cast<uint32_t>(input.size());
    const size_t compressedBodySize = deflater.Compress(
        reinterpret_cast<const char*>(input.data()), input.size(),
        result.bytes.data() + kBgzfHeaderSize,
        kBgzfMaxBlockSize - kBgzfHeaderSize - kBgzfFooterSize);
    if (compressedBodySize == 0) { // does not fit: store the block
        return makeStoredBgzfBlock(input);
    }

    const size_t totalSize = kBgzfHeaderSize + compressedBodySize + kBgzfFooterSize;
    std::copy(kBgzfHeader.begin(), kBgzfHeader.end(), result.bytes.begin());
    storeLe16(result.bytes.data() + 16, static_cast<uint16_t>(totalSize - 1));
    const uint32_t checksum = bgzf::Crc32(input.data(), input.size());
    storeLe32(result.bytes.data() + totalSize - 8, checksum);
    storeLe32(result.bytes.data() + totalSize - 4, result.uncompressedSize);
    result.bytes.resize(totalSize);
    return result;
}

// One zstd frame of the input, with its content size and checksum (one compression
// context per thread).
CompressedBlock compressZstdFrame(std::vector<uint8_t> input) {
    if (input.size() > kZstdFrameInput) {
        throw std::runtime_error("invalid zstd input frame size");
    }
    struct Context {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        Context() {
            if (cctx &&
                (ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, kZstdLevel)) ||
                 ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1)) ||
                 ZSTD_isError(ZSTD_CCtx_setParameter(cctx, ZSTD_c_contentSizeFlag, 1)))) {
                ZSTD_freeCCtx(cctx);
                cctx = nullptr;
            }
        }
        ~Context() { ZSTD_freeCCtx(cctx); }
    };
    thread_local Context context;
    if (!context.cctx) {
        throw std::runtime_error("cannot allocate a zstd compression context");
    }
    CompressedBlock result;
    result.uncompressedSize = static_cast<uint32_t>(input.size());
    result.bytes.resize(ZSTD_compressBound(input.size()));
    const size_t written = ZSTD_compress2(
        context.cctx, result.bytes.data(), result.bytes.size(), input.data(), input.size());
    if (ZSTD_isError(written)) {
        throw std::runtime_error(
            std::string("zstd compression failed: ") + ZSTD_getErrorName(written));
    }
    result.bytes.resize(written);
    return result;
}

// Compressed output in independent blocks, compressed in order on one thread or on the
// gzip thread pool and written in order: BGZF (gzip) or zstd frames.
class ParallelBlockOstream final : private std::streambuf,
                                   public std::ostream,
                                   public FinalizableOutput {
public:
    enum class Format { Bgzf, Zstd };

    ParallelBlockOstream(
        const std::string& filename,
        std::ios_base::openmode mode,
        Format format,
        bool parallel)
        : std::ostream(static_cast<std::streambuf*>(this)),
          filename_(filename),
          format_(format),
          parallel_(parallel),
          // Blocks in flight per stream: BGZF's are 64 KB, zstd's 1 MiB frames (one per pool thread, at most 8).
          maxPending_(!parallel ? 0
              : format == Format::Bgzf
                  ? std::min<size_t>(64, std::max<size_t>(2, ThreadPool::configured_gzip_thread_count() * 2))
                  : std::min<size_t>(8, std::max<size_t>(2, ThreadPool::configured_gzip_thread_count()))) {
        buffer_.resize(format_ == Format::Bgzf ? kBgzfBlockSize : kZstdFrameInput);
        setp(reinterpret_cast<char*>(buffer_.data()),
            reinterpret_cast<char*>(buffer_.data() + buffer_.size()));

        if (format_ == Format::Bgzf) {
            prepareAppend(mode);
        }
        else {
            std::error_code sizeError;
            const bool hasContent = (mode & std::ios::app) != 0 &&
                std::filesystem::file_size(filename_, sizeError) > 0 && !sizeError;
            emptyFrameIfNothing_ = !hasContent;
        }
        file_.open(filename_, mode | std::ios::binary);
        if (!file_.is_open()) {
            throw std::runtime_error("Error opening output file: " + filename_);
        }
        file_.exceptions(std::ios::badbit | std::ios::failbit);
    }

    ~ParallelBlockOstream() override {
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
                if (format_ == Format::Bgzf) {
                    file_.write(
                        reinterpret_cast<const char*>(kBgzfEof.data()),
                        static_cast<std::streamsize>(kBgzfEof.size()));
                }
                else if (emptyFrameIfNothing_ && compressedOffset_ == 0) {
                    // No content: one empty frame, so that the file is zstd for every reader.
                    const CompressedBlock empty = compressZstdFrame({});
                    file_.write(
                        reinterpret_cast<const char*>(empty.bytes.data()),
                        static_cast<std::streamsize>(empty.bytes.size()));
                }
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

        if (!failure_ && format_ == Format::Bgzf) {
            try {
                if (bgzfIndexOutput()) {
                    writeGzi();
                }
                else {
                    removeGzi();
                }
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

    static CompressedBlock compress(Format format, std::vector<uint8_t> input) {
        return format == Format::Bgzf
            ? compressBgzfBlock(std::move(input))
            : compressZstdFrame(std::move(input));
    }

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
                writeBlock(compress(format_, std::move(input)));
                return !failure_;
            }
            pending_.emplace_back(ThreadPool::gzip_instance().submit(
                [format = format_, block = std::move(input)]() mutable {
                    return compress(format, std::move(block));
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
        if (block.bytes.empty() || block.uncompressedSize == 0 ||
            (format_ == Format::Bgzf && block.bytes.size() > kBgzfMaxBlockSize)) {
            throw std::runtime_error("invalid compressed output block");
        }
        if (format_ == Format::Bgzf && uncompressedOffset_ != 0) {
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

        // Only block headers and ISIZE trailers are needed to rebuild the index.
        // Seek over the compressed payloads (unbuffered) instead of reading the
        // whole existing file, which made repeated appends quadratic.
        std::ifstream input;
        input.rdbuf()->pubsetbuf(nullptr, 0);
        input.open(filename_, std::ios::binary);
        if (!input) {
            throw std::runtime_error("Could not inspect BGZF file for append: " + filename_);
        }

        uint64_t compressed = 0;
        uint64_t uncompressed = 0;
        uint64_t dataEnd = 0;
        bool foundEof = false;
        while (compressed < fileSize) {
            std::array<uint8_t, kBgzfHeaderSize> header{};
            input.seekg(static_cast<std::streamoff>(compressed));
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

            // The EOF candidate is compared completely; other blocks only need ISIZE.
            const bool eofCandidate = blockSize == kBgzfEof.size();
            std::vector<uint8_t> rest(eofCandidate ? static_cast<size_t>(blockSize - header.size()) : 4);
            if (!eofCandidate) {
                input.seekg(static_cast<std::streamoff>(compressed + blockSize - rest.size()));
            }
            input.read(reinterpret_cast<char*>(rest.data()), rest.size());
            if (input.gcount() != static_cast<std::streamsize>(rest.size())) {
                throw std::runtime_error(
                    "Cannot append: truncated BGZF block in " + filename_);
            }
            const uint32_t blockUncompressed = loadLe32(rest.data() + rest.size() - 4);
            const bool isCanonicalEof = eofCandidate &&
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

    // An index from an earlier run would describe different blocks.
    void removeGzi() {
        std::error_code ignored;
        std::filesystem::remove(filename_ + ".gzi", ignored);
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
    Format format_ = Format::Bgzf;
    bool emptyFrameIfNothing_ = false; // zstd: the file holds no frame unless this stream writes one
    bool parallel_ = true;
    size_t maxPending_ = 2;
    bool closed_ = false;
    std::exception_ptr failure_;
};

} // namespace

void setBgzfIndexOutput(bool enabled) {
    g_writeBgzfIndex.store(enabled, std::memory_order_relaxed);
}

bool bgzfIndexOutput() {
    return g_writeBgzfIndex.load(std::memory_order_relaxed);
}

std::unique_ptr<std::ostream> openOutputFile(
    const std::string& filename,
    std::ios_base::openmode mode,
    size_t bufferSize,
    bool parallelGzip) {
    (void)bufferSize;
    if (isGZfile(filename) || isZstdFile(filename)) {
        return std::make_unique<ParallelBlockOstream>(
            filename, mode,
            isGZfile(filename) ? ParallelBlockOstream::Format::Bgzf : ParallelBlockOstream::Format::Zstd,
            shouldUseParallelGzip(parallelGzip));
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
