#include "InputFileStream.h"

#include "Common.h"
#include "ThreadPool.h"

#include <algorithm>
#include <array>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef _gzipread
#include <zlib.h>
#endif
#if defined(_gzipread) && defined(SDM_WITH_LIBDEFLATE)
#include <libdeflate.h>
#endif

namespace {

constexpr size_t kDefaultReadBufferSize = 64 * 1024;
constexpr size_t kMaxBgzfBlockSize = 64 * 1024;

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

#ifdef _gzipread

constexpr std::array<uint8_t, 28> kBgzfEof = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

void warnInput(const std::string& filename, const std::string& message) {
    static std::mutex warningMutex;
    std::lock_guard<std::mutex> lock(warningMutex);
    std::cerr << "Warning: " << message << ": " << filename << '\n';
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

class GzipStreambuf final : public std::streambuf {
public:
    GzipStreambuf(const std::string& filename, size_t requestedBufferSize)
        : filename_(filename),
          buffer_(std::max<size_t>(
              4096, requestedBufferSize == 0
                  ? kDefaultReadBufferSize
                  : requestedBufferSize)) {
        file_ = gzopen(filename_.c_str(), "rb");
        if (!file_) {
            throw std::runtime_error("Could not open gzip input: " + filename_);
        }
        const unsigned zlibBufferSize = static_cast<unsigned>(
            std::min<size_t>(buffer_.size(), std::numeric_limits<unsigned>::max()));
        (void)gzbuffer(file_, zlibBufferSize);
        setg(buffer_.data(), buffer_.data(), buffer_.data());
    }

    ~GzipStreambuf() override {
        close();
    }

protected:
    int_type underflow() override {
        if (gptr() != egptr()) {
            return traits_type::to_int_type(*gptr());
        }
        if (!file_) {
            return traits_type::eof();
        }

        const unsigned request = static_cast<unsigned>(
            std::min<size_t>(buffer_.size(), static_cast<size_t>(std::numeric_limits<int>::max())));
        const int count = gzread(file_, buffer_.data(), request);
        if (count > 0) {
            setg(buffer_.data(), buffer_.data(), buffer_.data() + count);
            return traits_type::to_int_type(*gptr());
        }

        int errorCode = Z_OK;
        const char* errorText = gzerror(file_, &errorCode);
        if (errorCode == Z_BUF_ERROR) {
            warnTruncated(errorText);
            return traits_type::eof();
        }
        if (count < 0 || (errorCode != Z_OK && errorCode != Z_STREAM_END)) {
            throw std::runtime_error(
                "Gzip decompression failed for " + filename_ + ": " +
                (errorText ? std::string(errorText) : std::string("unknown zlib error")));
        }
        return traits_type::eof();
    }

private:
    void warnTruncated(const char* detail = nullptr) {
        if (truncationWarned_) {
            return;
        }
        truncationWarned_ = true;
        std::string message = "truncated gzip input; continuing with recovered data";
        if (detail && *detail) {
            message += " (" + std::string(detail) + ")";
        }
        warnInput(filename_, message);
    }

    void close() noexcept {
        if (!file_) {
            return;
        }
        const int result = gzclose(file_);
        file_ = nullptr;
        if (result == Z_BUF_ERROR) {
            warnTruncated();
        }
        else if (result != Z_OK) {
            warnInput(filename_, "gzip close reported an error");
        }
    }

    std::string filename_;
    gzFile file_ = nullptr;
    std::vector<char> buffer_;
    bool truncationWarned_ = false;
};

class GzipInputStream final : public std::istream {
public:
    GzipInputStream(const std::string& filename, size_t bufferSize)
        : std::istream(new GzipStreambuf(filename, bufferSize)) {
        exceptions(std::ios::badbit);
    }

    ~GzipInputStream() override {
        exceptions(std::ios::goodbit);
        std::streambuf* buffer = rdbuf();
        rdbuf(nullptr);
        delete buffer;
    }
};

class IndexUnavailable final : public std::runtime_error {
public:
    explicit IndexUnavailable(const std::string& message)
        : std::runtime_error(message) {
    }
};

struct GziEntry {
    uint64_t compressedOffset = 0;
    uint64_t uncompressedOffset = 0;
};

struct BgzfBlock {
    uint64_t compressedOffset = 0;
    uint32_t compressedSize = 0;
    uint64_t uncompressedOffset = 0;
    uint32_t uncompressedSize = 0;
};

std::vector<uint8_t> readAt(
    std::ifstream& input,
    uint64_t offset,
    size_t size,
    const std::string& filename) {
    if (offset > static_cast<uint64_t>(std::numeric_limits<std::streamoff>::max())) {
        throw IndexUnavailable("offset exceeds stream limits");
    }
    std::vector<uint8_t> data(size);
    input.clear();
    input.seekg(static_cast<std::streamoff>(offset), std::ios::beg);
    if (!input.good()) {
        throw IndexUnavailable("could not seek in " + filename);
    }
    if (size != 0) {
        input.read(
            reinterpret_cast<char*>(data.data()),
            static_cast<std::streamsize>(size));
        if (input.gcount() != static_cast<std::streamsize>(size)) {
            throw IndexUnavailable("truncated BGZF block in " + filename);
        }
    }
    return data;
}

std::vector<GziEntry> readGzi(const std::string& indexFilename) {
    std::error_code sizeError;
    const uintmax_t rawSize = std::filesystem::file_size(indexFilename, sizeError);
    if (sizeError || rawSize < 8 || (rawSize - 8) % 16 != 0) {
        throw IndexUnavailable("invalid .gzi file size");
    }
    if (rawSize > static_cast<uintmax_t>(std::numeric_limits<size_t>::max())) {
        throw IndexUnavailable(".gzi file is too large");
    }

    std::ifstream input(indexFilename, std::ios::binary);
    if (!input.is_open()) {
        throw IndexUnavailable("could not open .gzi index");
    }
    std::vector<uint8_t> bytes(static_cast<size_t>(rawSize));
    input.read(
        reinterpret_cast<char*>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    if (input.gcount() != static_cast<std::streamsize>(bytes.size())) {
        throw IndexUnavailable("could not read complete .gzi index");
    }

    const uint64_t count = loadLe64(bytes.data());
    const uint64_t encodedCount = static_cast<uint64_t>((bytes.size() - 8) / 16);
    if (count != encodedCount) {
        throw IndexUnavailable(".gzi entry count does not match its file size");
    }

    std::vector<GziEntry> entries;
    entries.reserve(static_cast<size_t>(count) + 1);
    entries.push_back({0, 0});
    for (uint64_t i = 0; i < count; ++i) {
        const uint8_t* pair = bytes.data() + 8 + static_cast<size_t>(i) * 16;
        GziEntry entry{loadLe64(pair), loadLe64(pair + 8)};
        if (entry.compressedOffset <= entries.back().compressedOffset ||
            entry.uncompressedOffset <= entries.back().uncompressedOffset) {
            throw IndexUnavailable(".gzi offsets are not strictly increasing");
        }
        entries.push_back(entry);
    }
    return entries;
}

struct ParsedBgzfBlock {
    uint32_t compressedSize = 0;
    uint32_t uncompressedSize = 0;
    bool canonicalEof = false;
};

ParsedBgzfBlock parseBgzfBlock(
    std::ifstream& input,
    const std::string& filename,
    uint64_t offset,
    uint64_t fileSize) {
    if (offset > fileSize || fileSize - offset < 18) {
        throw IndexUnavailable("truncated BGZF header");
    }

    const std::vector<uint8_t> base = readAt(input, offset, 12, filename);
    if (base[0] != 0x1f || base[1] != 0x8b || base[2] != Z_DEFLATED ||
        (base[3] & 0x04) == 0 || (base[3] & 0xe0) != 0) {
        throw IndexUnavailable("indexed input is not BGZF");
    }

    const uint16_t extraLength = loadLe16(base.data() + 10);
    if (extraLength < 6 || fileSize - offset < 12u + extraLength + 8u) {
        throw IndexUnavailable("invalid BGZF extra field");
    }
    const std::vector<uint8_t> extra = readAt(input, offset + 12, extraLength, filename);

    bool foundBlockSize = false;
    uint32_t compressedSize = 0;
    size_t position = 0;
    while (position + 4 <= extra.size()) {
        const uint16_t subfieldLength = loadLe16(extra.data() + position + 2);
        const size_t next = position + 4 + subfieldLength;
        if (next > extra.size()) {
            throw IndexUnavailable("malformed BGZF extra subfield");
        }
        if (extra[position] == 'B' && extra[position + 1] == 'C' &&
            subfieldLength == 2) {
            compressedSize =
                static_cast<uint32_t>(loadLe16(extra.data() + position + 4)) + 1;
            foundBlockSize = true;
        }
        position = next;
    }
    if (!foundBlockSize || compressedSize < 26 ||
        compressedSize > kMaxBgzfBlockSize ||
        compressedSize > fileSize - offset) {
        throw IndexUnavailable("invalid BGZF BSIZE");
    }

    const std::vector<uint8_t> footer =
        readAt(input, offset + compressedSize - 4, 4, filename);
    ParsedBgzfBlock result;
    result.compressedSize = compressedSize;
    result.uncompressedSize = loadLe32(footer.data());
    if (compressedSize == kBgzfEof.size()) {
        const std::vector<uint8_t> candidate =
            readAt(input, offset, kBgzfEof.size(), filename);
        result.canonicalEof =
            std::equal(kBgzfEof.begin(), kBgzfEof.end(), candidate.begin());
    }
    return result;
}

std::vector<BgzfBlock> loadIndexedBgzf(
    const std::string& filename,
    const std::string& indexFilename) {
    const std::vector<GziEntry> entries = readGzi(indexFilename);

    std::error_code sizeError;
    const uintmax_t rawFileSize = std::filesystem::file_size(filename, sizeError);
    if (sizeError || rawFileSize > std::numeric_limits<uint64_t>::max()) {
        throw IndexUnavailable("could not determine compressed file size");
    }
    const uint64_t fileSize = static_cast<uint64_t>(rawFileSize);
    if (fileSize == 0) {
        throw IndexUnavailable("indexed gzip file is empty");
    }

    std::ifstream input(filename, std::ios::binary);
    if (!input.is_open()) {
        throw IndexUnavailable("could not open indexed gzip input");
    }

    const ParsedBgzfBlock first =
        parseBgzfBlock(input, filename, entries.front().compressedOffset, fileSize);
    if (first.canonicalEof) {
        if (entries.size() != 1 || first.compressedSize != fileSize) {
            throw IndexUnavailable("invalid empty BGZF index");
        }
        return {};
    }

    std::vector<BgzfBlock> blocks;
    blocks.reserve(entries.size());
    uint64_t dataEnd = 0;
    for (size_t i = 0; i < entries.size(); ++i) {
        const GziEntry& entry = entries[i];
        const ParsedBgzfBlock parsed =
            parseBgzfBlock(input, filename, entry.compressedOffset, fileSize);
        if (parsed.canonicalEof) {
            throw IndexUnavailable(".gzi points to the BGZF EOF marker");
        }

        const uint64_t blockEnd =
            entry.compressedOffset + parsed.compressedSize;
        if (i + 1 < entries.size()) {
            const GziEntry& next = entries[i + 1];
            if (blockEnd != next.compressedOffset) {
                throw IndexUnavailable(".gzi compressed offsets do not match BGZF blocks");
            }
            if (entry.uncompressedOffset + parsed.uncompressedSize !=
                next.uncompressedOffset) {
                throw IndexUnavailable(".gzi uncompressed offsets do not match BGZF blocks");
            }
        }
        blocks.push_back({
            entry.compressedOffset,
            parsed.compressedSize,
            entry.uncompressedOffset,
            parsed.uncompressedSize
        });
        dataEnd = blockEnd;
    }

    if (dataEnd == fileSize) {
        warnInput(
            filename,
            "BGZF EOF marker is missing; input may be truncated, continuing with indexed data");
        return blocks;
    }
    if (fileSize - dataEnd != kBgzfEof.size()) {
        throw IndexUnavailable("data follows indexed BGZF blocks");
    }
    const std::vector<uint8_t> eof =
        readAt(input, dataEnd, kBgzfEof.size(), filename);
    if (!std::equal(kBgzfEof.begin(), kBgzfEof.end(), eof.begin())) {
        throw IndexUnavailable("indexed BGZF data is not followed by the canonical EOF marker");
    }
    return blocks;
}

std::vector<char> inflateBgzfBlock(
    std::vector<uint8_t> compressed,
    BgzfBlock block,
    const std::string& filename) {
    std::vector<char> output(static_cast<size_t>(block.uncompressedSize) + 1);

#ifdef SDM_WITH_LIBDEFLATE
    // One decompressor per thread. It checks the gzip header, CRC-32 and ISIZE.
    struct Decompressor {
        libdeflate_decompressor* handle = libdeflate_alloc_decompressor();
        ~Decompressor() { libdeflate_free_decompressor(handle); }
    };
    thread_local Decompressor decompressor;
    if (!decompressor.handle) {
        throw std::runtime_error("libdeflate_alloc_decompressor failed");
    }
    size_t consumed = 0, produced = 0;
    const libdeflate_result status = libdeflate_gzip_decompress_ex(decompressor.handle,
        compressed.data(), compressed.size(), output.data(), output.size(), &consumed, &produced);
    if (status != LIBDEFLATE_SUCCESS || produced != block.uncompressedSize || consumed != compressed.size()) {
        throw std::runtime_error(
            "Invalid BGZF block in " + filename + " at compressed offset " +
            std::to_string(block.compressedOffset));
    }
#else
    z_stream stream{};
    stream.next_in = compressed.data();
    stream.avail_in = static_cast<uInt>(compressed.size());
    stream.next_out = reinterpret_cast<Bytef*>(output.data());
    stream.avail_out = static_cast<uInt>(output.size());

    const int initResult = inflateInit2(&stream, 15 + 16);
    if (initResult != Z_OK) {
        throw std::runtime_error(
            "inflateInit2 failed for " + filename + " at BGZF offset " +
            std::to_string(block.compressedOffset));
    }

    const int inflateResult = inflate(&stream, Z_FINISH);
    const size_t produced = static_cast<size_t>(stream.total_out);
    const uInt remainingInput = stream.avail_in;
    const int endResult = inflateEnd(&stream);
    if (inflateResult != Z_STREAM_END || endResult != Z_OK ||
        produced != block.uncompressedSize || remainingInput != 0) {
        throw std::runtime_error(
            "Invalid BGZF block in " + filename + " at compressed offset " +
            std::to_string(block.compressedOffset));
    }
#endif
    output.resize(produced);
    return output;
}

class ParallelBgzfInputMarker {
public:
    virtual ~ParallelBgzfInputMarker() = default;
};

class ParallelBgzfStreambuf final : public std::streambuf {
public:
    ParallelBgzfStreambuf(
        std::string filename,
        std::vector<BgzfBlock> blocks)
        : filename_(std::move(filename)),
          input_(filename_, std::ios::binary),
          blocks_(std::move(blocks)),
          maxPending_(std::min<size_t>(
              8,
              std::max<size_t>(
                  2,
                  ThreadPool::configured_gzip_thread_count() * 2))) {
        if (!input_.is_open()) {
            throw std::runtime_error("Could not open indexed BGZF input: " + filename_);
        }
        setg(nullptr, nullptr, nullptr);
    }

protected:
    int_type underflow() override {
        if (gptr() && gptr() != egptr()) {
            return traits_type::to_int_type(*gptr());
        }

        for (;;) {
            fillPending();
            if (pending_.empty()) {
                return traits_type::eof();
            }

            std::future<std::vector<char>> future = std::move(pending_.front());
            pending_.pop_front();
            current_ = future.get();
            fillPending();
            if (current_.empty()) {
                continue;
            }
            setg(current_.data(), current_.data(), current_.data() + current_.size());
            return traits_type::to_int_type(*gptr());
        }
    }

private:
    std::vector<uint8_t> readCompressedBlock(const BgzfBlock& block) {
        if (block.compressedOffset >
            static_cast<uint64_t>(std::numeric_limits<std::streamoff>::max())) {
            throw std::runtime_error("BGZF offset exceeds stream limits: " + filename_);
        }
        std::vector<uint8_t> compressed(block.compressedSize);
        input_.clear();
        input_.seekg(
            static_cast<std::streamoff>(block.compressedOffset),
            std::ios::beg);
        input_.read(
            reinterpret_cast<char*>(compressed.data()),
            static_cast<std::streamsize>(compressed.size()));
        if (input_.gcount() != static_cast<std::streamsize>(compressed.size())) {
            throw std::runtime_error(
                "Could not read indexed BGZF block from " + filename_);
        }
        return compressed;
    }

    void fillPending() {
        while (nextBlock_ < blocks_.size() && pending_.size() < maxPending_) {
            const BgzfBlock block = blocks_[nextBlock_++];
            std::vector<uint8_t> compressed = readCompressedBlock(block);
            pending_.emplace_back(ThreadPool::gzip_instance().submit(
                [compressed = std::move(compressed), block, filename = filename_]() mutable {
                    return inflateBgzfBlock(
                        std::move(compressed), block, filename);
                }));
        }
    }

    std::string filename_;
    std::ifstream input_;
    std::vector<BgzfBlock> blocks_;
    std::deque<std::future<std::vector<char>>> pending_;
    std::vector<char> current_;
    size_t nextBlock_ = 0;
    size_t maxPending_ = 2;
};

class ParallelBgzfInputStream final : public std::istream,
                                      public ParallelBgzfInputMarker {
public:
    ParallelBgzfInputStream(
        const std::string& filename,
        std::vector<BgzfBlock> blocks)
        : std::istream(new ParallelBgzfStreambuf(
              filename, std::move(blocks))) {
        exceptions(std::ios::badbit);
    }

    ~ParallelBgzfInputStream() override {
        exceptions(std::ios::goodbit);
        std::streambuf* buffer = rdbuf();
        rdbuf(nullptr);
        delete buffer;
    }
};

// Plain gzip streams cannot be inflated in parallel, but each input file (R1,
// R2, MID) can be inflated on its own thread, ahead of and concurrently with
// record parsing on the reader thread. Decompression errors are rethrown on
// the reading thread at the position where they occurred.
class ReadaheadStreambuf final : public std::streambuf {
public:
    ReadaheadStreambuf(std::unique_ptr<std::istream> source, size_t chunkSize, size_t depth)
        : source_(std::move(source)), chunkSize_(chunkSize), depth_(depth) {
        thread_ = std::thread([this] { run(); });
    }
    ReadaheadStreambuf(const ReadaheadStreambuf&) = delete;
    ReadaheadStreambuf& operator=(const ReadaheadStreambuf&) = delete;
    ~ReadaheadStreambuf() override {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        space_.notify_all();
        if (thread_.joinable()) {
            thread_.join();
        }
    }

protected:
    int_type underflow() override {
        if (gptr() != egptr()) {
            return traits_type::to_int_type(*gptr());
        }
        std::unique_lock<std::mutex> lock(mutex_);
        if (!current_.empty()) {
            spare_.push_back(std::move(current_));
            current_.clear();
        }
        data_.wait(lock, [this] { return !ready_.empty() || done_; });
        if (ready_.empty()) {
            setg(nullptr, nullptr, nullptr);
            if (error_) {
                std::exception_ptr error = error_;
                error_ = nullptr;
                lock.unlock();
                std::rethrow_exception(error);
            }
            return traits_type::eof();
        }
        current_ = std::move(ready_.front());
        ready_.pop_front();
        // Wake the producer once half the queue is free, not for every chunk: each
        // wake-up of a sleeping thread is a futex call on this, the reading, thread.
        const bool wakeProducer = ready_.size() <= depth_ / 2;
        lock.unlock();
        if (wakeProducer) { space_.notify_one(); }
        setg(current_.data(), current_.data(), current_.data() + current_.size());
        return traits_type::to_int_type(*gptr());
    }

private:
    void run() {
        try {
            for (;;) {
                std::vector<char> chunk;
                {
                    std::unique_lock<std::mutex> lock(mutex_);
                    space_.wait(lock, [this] { return stop_ || ready_.size() < depth_; });
                    if (stop_) {
                        return;
                    }
                    if (!spare_.empty()) {
                        chunk = std::move(spare_.back());
                        spare_.pop_back();
                    }
                }
                chunk.resize(chunkSize_);
                source_->read(chunk.data(), static_cast<std::streamsize>(chunk.size()));
                chunk.resize(static_cast<size_t>(source_->gcount()));
                const bool ended = chunk.size() < chunkSize_ || !*source_;
                bool wakeConsumer;
                {
                    std::lock_guard<std::mutex> lock(mutex_);
                    if (!chunk.empty()) {
                        ready_.push_back(std::move(chunk));
                    }
                    done_ = ended;
                    // The consumer only waits on an empty queue.
                    wakeConsumer = ended || ready_.size() == 1;
                }
                if (wakeConsumer) { data_.notify_one(); }
                if (ended) {
                    return;
                }
            }
        }
        catch (...) {
            std::lock_guard<std::mutex> lock(mutex_);
            error_ = std::current_exception();
            done_ = true;
            data_.notify_one();
        }
    }

    std::unique_ptr<std::istream> source_;
    const size_t chunkSize_;
    const size_t depth_;
    std::mutex mutex_;
    std::condition_variable data_, space_;
    std::deque<std::vector<char>> ready_;
    std::vector<std::vector<char>> spare_;
    std::vector<char> current_;
    bool done_ = false;
    bool stop_ = false;
    std::exception_ptr error_;
    std::thread thread_;
};

class ReadaheadInputStream final : public std::istream {
public:
    explicit ReadaheadInputStream(std::unique_ptr<std::istream> source)
        : std::istream(new ReadaheadStreambuf(std::move(source), 256 * 1024, 8)) {
        exceptions(std::ios::badbit);
    }

    ~ReadaheadInputStream() override {
        exceptions(std::ios::goodbit);
        std::streambuf* buffer = rdbuf();
        rdbuf(nullptr);
        delete buffer;
    }
};

#endif

} // namespace

std::unique_ptr<std::istream> openInputFile(
    const std::string& filename,
    std::ios_base::openmode mode,
    size_t bufferSize,
    bool parallelGzip) {
    if (!isGZfile(filename)) {
        auto input = std::make_unique<std::ifstream>(
            filename, mode | std::ios::in | std::ios::binary);
        if (!input->is_open()) {
            throw std::runtime_error("Could not open input file: " + filename);
        }
        return input;
    }

#ifdef _gzipread
    const std::string indexFilename = filename + ".gzi";
    std::error_code indexError;
    const bool hasIndex = std::filesystem::exists(indexFilename, indexError);
    const bool useParallelGzip = shouldUseParallelGzip(parallelGzip);
    if (useParallelGzip && !indexError && hasIndex) {
        try {
            return std::make_unique<ParallelBgzfInputStream>(
                filename,
                loadIndexedBgzf(filename, indexFilename));
        }
        catch (const IndexUnavailable& error) {
            warnInput(
                filename,
                "cannot use .gzi index (" + std::string(error.what()) +
                    "); falling back to sequential gzip");
        }
    }

    auto sequential = std::make_unique<GzipInputStream>(filename, bufferSize);
    if (useParallelGzip) {
        return std::make_unique<ReadaheadInputStream>(std::move(sequential));
    }
    return sequential;
#else
    (void)mode;
    (void)bufferSize;
    (void)parallelGzip;
    throw std::runtime_error("gzip input requested but support is not compiled: " + filename);
#endif
}

bool inputFileUsesParallelBgzf(const std::istream& stream) {
#ifdef _gzipread
    return dynamic_cast<const ParallelBgzfInputMarker*>(&stream) != nullptr;
#else
    (void)stream;
    return false;
#endif
}
