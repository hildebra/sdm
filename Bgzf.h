// Bgzf.h - gzip as sdm writes it, and reads it where it can, with ISA-L (igzip).
//
// Ported from protal (github.com/hildebra/protalFastHit, branch audit-fixes, 956e18b, src/IO/Bgzf.h);
// the measurements quoted here and in ThreadedGzStream.h / ParallelGzip.* are protal's (its docs/claude).
//
// sdm writes every gzip output (.gz, .gzip; OutputFileStream.cpp) as BGZF (SAM specification, section
// 4.1): gzip members of at most 64 KB, each with its size in a "BC" extra field, then an empty member
// as end-of-file marker. zcat, gzip and zlib read it as any multi-member gzip file; htslib reads
// its blocks independently. Every block is compressed and decompressed whole, with ISA-L's
// stateless (de)compression. Other gzip files, e.g. the single-member .fq.gz of sequencers, are read
// with ISA-L's streaming inflate (ThreadedGzStream.h).
//
// ISA-L (protal docs/claude/2026-10-06-performance-profiling, the ISA-L section): it inflates about 1.6-1.8x
// as fast as zlib-ng (libdeflate was 10-20% faster on whole blocks but cannot stream); at level 1 it
// deflates about 7x as fast as libdeflate's level 6, into files about 16% larger. Its levels 1 and 2 wrote
// the same bytes through its SSE4.2, AVX and AVX2 code (AVX-512 not checked), level 3 not (its AVX2
// path matches differently). Level 1 writes the same bytes on any thread since 2026-10-08 (Deflater).
#pragma once

#include <isa-l/crc.h>
#include <isa-l/igzip_lib.h>
#include <sys/mman.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <memory>
#include <new>
#include <string>
#include <vector>

#if defined(__linux__) && !defined(MAP_FIXED_NOREPLACE)
#define MAP_FIXED_NOREPLACE 0x100000  // Linux 4.17; an older kernel takes the address as a hint (MapStream checks)
#endif

namespace sdm::bgzf {
    inline constexpr size_t kBlockInput = 0xff00;      // input bytes per block, as htslib
    inline constexpr size_t kMaxBlock = 0x10000;       // a block's size must fit its 16-bit field
    inline constexpr size_t kHeaderBytes = 18, kFooterBytes = 8;
    inline constexpr uint32_t kLevel = 1;              // ISA-L's: level 2 compresses no better, 3 differs by CPU
    inline constexpr uint32_t kLevelBuffer = ISAL_DEF_LVL1_DEFAULT;  // fixed: the output depends on it
    inline constexpr unsigned char kHeader[16] = { 0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0 };
    inline constexpr unsigned char kEof[28] = { 0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0,
                                                0x1b, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    inline void PutLE16(unsigned char* p, uint32_t v) {
        p[0] = static_cast<unsigned char>(v & 0xff);
        p[1] = static_cast<unsigned char>(v >> 8 & 0xff);
    }

    inline void PutLE32(unsigned char* p, uint32_t v) {
        for (int i = 0; i < 4; i++) p[i] = static_cast<unsigned char>(v >> (8 * i) & 0xff);
    }

    inline uint32_t LE32(unsigned char const* p) {
        return uint32_t(p[0]) | uint32_t(p[1]) << 8 | uint32_t(p[2]) << 16 | uint32_t(p[3]) << 24;
    }

    // Where a Deflater's ISA-L stream is put. ISA-L's level 1 (the SSE4.2, AVX and AVX2 versions of
    // igzip_icf_body_h1_gr_bt.asm, 2.31.0 to 2.32.1 and its master of 2026-10) hashes the input's third
    // byte, at the start of a stream without history (every stateless call), from a register it does not
    // load there (its level 0 does), which holds the address of the isal_zstream: bits 16-47 of it chose a
    // hash bucket, and with it now and then a match (in the block that showed it, the bucket of its first
    // four bytes: their next occurrence lost its match). So a block's bytes depended on the thread (and the
    // run) that compressed it, though they always held the same content (2026-10-08: one block of a
    // simulated sample in 50 unit suite runs; protal docs/claude/2026-10-08-bgzf-thread-bytes). A Deflater's
    // stream is therefore mapped where the CRC-32C of those bits (ISA-L's hash) is kStreamHash under the
    // level's hash mask: every stream hashes that byte into the same bucket, and a block gets the same
    // bytes on any thread. Where no such place is mapped, the stream is on the heap (the output is valid,
    // but may differ by thread).
    namespace detail {
        static_assert(kLevel == 1, "the hash mask is level 1's");
        inline constexpr uint32_t kStreamHashMask = IGZIP_LVL1_HASH_SIZE - 1;  // a block's mask is a part of it
        inline constexpr uint32_t kStreamHash = 0;
        inline constexpr uintptr_t kStreamStep = uintptr_t{1} << 16;  // the bits below do not count

        inline constexpr auto kCrc32cTable = [] {
            std::array<uint32_t, 256> table{};
            for (uint32_t i = 0; i < 256; i++) {
                uint32_t crc = i;
                for (int bit = 0; bit < 8; bit++) crc = (crc >> 1) ^ (0x82f63b78u & (0u - (crc & 1u)));
                table[i] = crc;
            }
            return table;
        }();

        // The CRC-32C of a 32-bit word, from 0 and without the final complement: the x86 crc32
        // instruction on a 32-bit operand, ISA-L's hash.
        inline uint32_t Crc32cWord(uint32_t word) {
            uint32_t crc = 0;
            for (int byte = 0; byte < 4; byte++) crc = (crc >> 8) ^ kCrc32cTable[(crc ^ (word >> (8 * byte))) & 0xff];
            return crc;
        }

        inline bool IsStreamPlace(uintptr_t address) {
            return (Crc32cWord(static_cast<uint32_t>(address >> 16)) & kStreamHashMask) == kStreamHash;
        }

        // `size` bytes of zeroes at an address that IsStreamPlace (one in 8192 places 64 KB apart), mapped
        // below where the kernel maps now (the first free such place: one freed is used again); null if
        // none was found within 64 GB, or the system has no MAP_FIXED_NOREPLACE (off Linux).
        inline void* MapStream(size_t size) {
#ifndef MAP_FIXED_NOREPLACE
            (void)size;
            return nullptr;
#else
            void* const probe = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
            if (probe == MAP_FAILED) return nullptr;
            uintptr_t at = reinterpret_cast<uintptr_t>(probe) & ~(kStreamStep - 1);
            munmap(probe, size);
            int failed = 0;
            for (int step = 0; step < (1 << 20) && failed < 4 && at > (uintptr_t{1} << 32); step++, at -= kStreamStep) {
                if (!IsStreamPlace(at)) continue;
                void* const p = mmap(reinterpret_cast<void*>(at), size, PROT_READ | PROT_WRITE,
                                     MAP_PRIVATE | MAP_ANONYMOUS | MAP_FIXED_NOREPLACE, -1, 0);
                if (p == reinterpret_cast<void*>(at)) return p;
                if (p == MAP_FAILED && errno == EEXIST) continue;  // in use: further down
                if (p != MAP_FAILED) munmap(p, size);  // mapped elsewhere: the flag was taken for a hint
                failed++;
            }
            return nullptr;
#endif
        }
    }

    // An ISA-L compressor at kLevel, for one thread: raw deflate of a whole block. Its state and level
    // buffer (~80 and ~280 KB) are kept for the thread's blocks; the state is mapped by MapStream.
    class Deflater {
    public:
        Deflater() : m_level_buffer(kLevelBuffer) {
            void* const place = detail::MapStream(sizeof(isal_zstream));
            m_mapped = place != nullptr;
            m_stream = m_mapped ? new (place) isal_zstream() : new isal_zstream();
        }
        ~Deflater() {
            if (m_mapped) munmap(m_stream, sizeof(isal_zstream));
            else delete m_stream;
        }
        Deflater(Deflater const&) = delete;
        Deflater& operator=(Deflater const&) = delete;

        // Raw deflate of [in, in + size) into out; the compressed size, or 0 if it does not fit. A block
        // that does not shrink is stored (ISA-L does it), so kBlockInput bytes always fit the BGZF limit.
        size_t Compress(char const* in, size_t size, unsigned char* out, size_t capacity) {
            isal_zstream& stream = *m_stream;
            isal_deflate_stateless_init(&stream);  // raw deflate (IGZIP_DEFLATE), the whole input at once
            stream.level = kLevel;
            stream.level_buf = m_level_buffer.data();
            stream.level_buf_size = static_cast<uint32_t>(m_level_buffer.size());
            stream.end_of_stream = 1;
            stream.flush = NO_FLUSH;
            stream.next_in = reinterpret_cast<uint8_t*>(const_cast<char*>(in));
            stream.avail_in = static_cast<uint32_t>(size);
            stream.next_out = out;
            stream.avail_out = static_cast<uint32_t>(std::min<size_t>(capacity, UINT32_MAX));
            return isal_deflate_stateless(&stream) == COMP_OK ? stream.total_out : 0;
        }

    private:
        isal_zstream* m_stream = nullptr;
        bool m_mapped = false;
        std::vector<uint8_t> m_level_buffer;
    };

    // The CRC-32 of gzip, as zlib's crc32(0, data, size).
    inline uint32_t Crc32(void const* data, size_t size) {
        return crc32_gzip_refl(0, static_cast<unsigned char const*>(data), size);
    }

    // Appends [data, data + size) to out as BGZF blocks, deflated at kLevel (a block that does not
    // shrink is stored). False only if ISA-L fails.
    inline bool Compress(char const* data, size_t size, std::string& out) {
        thread_local Deflater deflater;
        constexpr size_t kPayload = kMaxBlock - kHeaderBytes - kFooterBytes;
        for (size_t offset = 0; offset < size; offset += kBlockInput) {
            size_t const n = std::min(kBlockInput, size - offset);
            size_t const start = out.size();
            out.resize(start + kMaxBlock);
            auto* block = reinterpret_cast<unsigned char*>(out.data() + start);
            size_t const compressed = deflater.Compress(data + offset, n, block + kHeaderBytes, kPayload);
            if (compressed == 0) return false;
            size_t const total = kHeaderBytes + compressed + kFooterBytes;
            std::memcpy(block, kHeader, sizeof(kHeader));
            PutLE16(block + 16, static_cast<uint32_t>(total - 1));
            PutLE32(block + kHeaderBytes + compressed, Crc32(data + offset, n));
            PutLE32(block + kHeaderBytes + compressed + 4, static_cast<uint32_t>(n));
            out.resize(start + total);
        }
        return true;
    }

    // Whether the 18 bytes at h begin a BGZF block (gzip, FEXTRA with the BC subfield first).
    inline bool IsBlockHeader(unsigned char const* h) {
        return h[0] == 0x1f && h[1] == 0x8b && h[2] == 8 && (h[3] & 4) && h[10] == 6 && h[11] == 0 &&
               h[12] == 'B' && h[13] == 'C' && h[14] == 2 && h[15] == 0;
    }

    // The size of the block whose header is at h, header and footer included.
    inline size_t BlockSize(unsigned char const* h) {
        return (size_t(h[16]) | size_t(h[17]) << 8) + 1;
    }

    // Decompresses the whole block [block, block + size) into out (capacity bytes, at least its
    // content): its content size in n, or false and error for a corrupt block.
    inline bool DecompressBlock(unsigned char const* block, size_t size, char* out, size_t capacity, size_t& n, std::string& error) {
        thread_local std::unique_ptr<inflate_state> const state = std::make_unique<inflate_state>();  // ~85 KB
        if (size < kHeaderBytes + kFooterBytes || !IsBlockHeader(block) || BlockSize(block) != size) {
            error = "not a BGZF block";
            return false;
        }
        uint32_t const crc = LE32(block + size - 8), isize = LE32(block + size - 4);
        if (isize > capacity) {
            error = "a BGZF block holds more than 64 KB";
            return false;
        }
        // Raw deflate, the whole block at once: its content must be exactly isize bytes.
        isal_inflate_init(state.get());
        state->crc_flag = ISAL_DEFLATE;
        state->next_in = const_cast<uint8_t*>(block + kHeaderBytes);
        state->avail_in = static_cast<uint32_t>(size - kHeaderBytes - kFooterBytes);
        state->next_out = reinterpret_cast<uint8_t*>(out);
        state->avail_out = isize;
        if (isal_inflate_stateless(state.get()) != ISAL_DECOMP_OK || state->total_out != isize) {
            error = "a corrupt BGZF block";
            return false;
        }
        if (Crc32(out, isize) != crc) {
            error = "a BGZF block fails its CRC check";
            return false;
        }
        n = isize;
        return true;
    }

    // Whether the file begins with a BGZF block header.
    inline bool StartsAsBgzf(std::string const& path) {
        std::ifstream is(path, std::ios::binary);
        unsigned char b[kHeaderBytes] = {};
        is.read(reinterpret_cast<char*>(b), sizeof(b));
        return is.gcount() == static_cast<std::streamsize>(sizeof(b)) && IsBlockHeader(b);
    }

    // Whether the file ends with the BGZF end-of-file block.
    inline bool EndsWithEof(std::string const& path) {
        std::error_code ec;
        auto const size = std::filesystem::file_size(path, ec);
        if (ec || size < sizeof(kEof)) return false;
        std::ifstream is(path, std::ios::binary);
        is.seekg(static_cast<std::streamoff>(size - sizeof(kEof)));
        unsigned char b[sizeof(kEof)] = {};
        is.read(reinterpret_cast<char*>(b), sizeof(b));
        return is.gcount() == static_cast<std::streamsize>(sizeof(b)) && std::memcmp(b, kEof, sizeof(b)) == 0;
    }

    // Writes BGZF to a file as its content arrives, on one thread: the blocks start every kBlockInput bytes of the
    // content, however it arrives, so the file depends on the content alone, and the content is never on disk
    // uncompressed. Close() writes the rest and the end-of-file block; check Error() then.
    class Writer {
    public:
        explicit Writer(std::string path) : m_path(std::move(path)), m_out(std::fopen(m_path.c_str(), "wb")) {
            if (!m_out) m_error = "cannot write " + m_path + ": " + std::strerror(errno);
        }
        ~Writer() {
            if (m_out) std::fclose(m_out);
        }
        Writer(Writer const&) = delete;
        Writer& operator=(Writer const&) = delete;

        void Write(char const* data, size_t size) {
            if (!m_error.empty()) return;
            m_buffer.append(data, size);
            if (m_buffer.size() >= kFlush) Flush(false);
        }

        // Writes what is left and the end-of-file block, and closes the file; false (see Error) if anything failed.
        bool Close() {
            if (!m_out) return false;
            Flush(true);
            if (m_error.empty() && std::fwrite(kEof, 1, sizeof(kEof), m_out) != sizeof(kEof)) {
                m_error = "writing " + m_path + " failed: " + std::strerror(errno);
            }
            if (std::fclose(m_out) != 0 && m_error.empty()) m_error = "writing " + m_path + " failed: " + std::strerror(errno);
            m_out = nullptr;
            return m_error.empty();
        }

        std::string const& Error() const { return m_error; }

    private:
        static constexpr size_t kFlush = 64 * kBlockInput;  // ~4 MB of content, whole blocks, at a time

        // Compresses the buffer's whole blocks (all of it if `all`) and writes them.
        void Flush(bool all) {
            size_t const n = all ? m_buffer.size() : m_buffer.size() / kBlockInput * kBlockInput;
            if (n == 0 || !m_error.empty()) return;
            m_packed.clear();
            if (!Compress(m_buffer.data(), n, m_packed)) {
                m_error = "compressing for " + m_path + " failed";
                return;
            }
            if (std::fwrite(m_packed.data(), 1, m_packed.size(), m_out) != m_packed.size()) {
                m_error = "writing " + m_path + " failed: " + std::strerror(errno);
                return;
            }
            m_buffer.erase(0, n);
        }

        std::string m_path;
        std::FILE* m_out;
        std::string m_buffer, m_packed, m_error;
    };
}
