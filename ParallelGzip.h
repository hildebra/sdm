// ParallelGzip.h - one gzip file inflated on several threads with ISA-L, for ThreadedGzStreambuf.
//
// A gzip member (as gzip and most tools write .fastq.gz files: the whole file one member) cannot be cut
// at block boundaries as BGZF can: its deflate blocks are found only by decoding it from the start, and
// refer back up to 32 KB. Here the file is cut into pieces (kChunkBytes of it each), each beginning at a
// deflate block found by its header near the piece's start. Worker threads decode pieces ahead with ISA-L,
// twice, with the unknown 32 KB before the piece as two sets of labels: where the two outputs differ the
// byte is a copy from that window, and the labels say from where, so the bytes are filled in once the piece
// before is done. The reader's own thread decodes in order, with the true window, the pieces no worker took,
// and chains the windows. A piece is used only if the one before ends exactly where it begins, so a block
// header found by chance is never used; ISA-L's stock API is all it needs
// (protal docs/claude/2026-10-12-isal-parallel-inflate). Read gives the file's bytes in order, its members one
// after the other; each member's CRC-32 and length are checked, and that the file ends with a whole
// member (or zero bytes of padding). Small members (as tools that compress a file in chunks write) hold
// too few block starts for pieces: where they follow a member, Read stops, and the caller reads the rest by
// members (StoppedAt; protal docs/claude/2026-10-13-gzip-members-regression).
#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

namespace sdm {

    // Whether read files are asked to be read ahead of where they are read (ParallelGzip's mapping, ThreadedGzStreambuf's
    // reads): SDM_READ_AHEAD=0 in the environment turns it off, to compare runs reading their input cold.
    inline bool ReadAheadWanted() {
        static bool const wanted = [] {
            char const* env = std::getenv("SDM_READ_AHEAD");
            return !(env && std::strcmp(env, "0") == 0);
        }();
        return wanted;
    }

    class ParallelGzip {
    public:
        // Bytes of the compressed file per piece: on FASTQ ~4 MB of output, decoded twice ahead (1 MiB pieces were the
        // fastest of 1, 4 and 16 MiB in protal docs/claude/2026-10-12-isal-parallel-inflate, and held the least memory).
        static constexpr size_t kChunkBytes = size_t{1} << 20;

        // A reader of the regular file open as fd, from its start, on `threads` threads (the caller's and threads - 1
        // of its own; at least 2): the file is mapped, so its offset and the descriptor stay the caller's. Null if the
        // reader cannot be made (error then says why: not a regular file, no gzip header at its start). With
        // `max_member`, Read stops after a member that the next two members follow within max_member bytes each
        // (StoppedAt).
        static std::unique_ptr<ParallelGzip> Open(int fd, unsigned threads, std::string& error, size_t max_member = 0);

        ~ParallelGzip();
        ParallelGzip(ParallelGzip const&) = delete;
        ParallelGzip& operator=(ParallelGzip const&) = delete;

        // Up to size inflated bytes into out, the next ones of the file; 0 at its end, where it stopped (StoppedAt), or
        // on an error (Error() then says what: a corrupt member, a CRC-32 that does not match, a truncated file).
        size_t Read(char* out, size_t size);

        std::string const& Error() const { return m_error; }

        // Once Read gives 0: the byte where the member begins that Read stopped before (small members follow), every
        // member before it read and checked, and how many they were; kNotStopped at the file's end or an error.
        static constexpr uint64_t kNotStopped = UINT64_MAX;
        uint64_t StoppedAt() const;
        size_t MembersRead() const;
        // The members begun so far (the file's members, once it is read; the caller's thread's, while Read is not running).
        size_t MembersBegun() const;

        // How the file's ranges were decoded so far: pieces decoded ahead by the workers and used, ranges the reader's
        // thread decoded in order, ranges decoded on to their member's end (no piece start after them, or a corrupt or
        // too large piece), block starts found by their headers that were none, and pieces decoded ahead in vain.
        struct Counts {
            size_t ahead = 0, in_order = 0, streamed = 0, dropped_starts = 0;
            size_t unused = 0;  // pieces decoded ahead that were not used (a failed decode, another end)
        };
        Counts counts() const;

        // For tests: the bytes of the compressed file per piece of the readers opened from now on (0: kChunkBytes), so
        // that a small file has many, and how far ahead of the reader's thread the workers begin (pieces; 0: 1 per
        // thread and 2), so that they decode pieces the reader then takes even in a small file.
        static void SetPieceBytes(size_t bytes, size_t gap = 0) {
            PieceBytesSetting().store(bytes, std::memory_order_relaxed);
            GapSetting().store(gap, std::memory_order_relaxed);
        }

    private:
        ParallelGzip();
        static std::atomic<size_t>& PieceBytesSetting();
        static std::atomic<size_t>& GapSetting();

        struct Impl;
        std::unique_ptr<Impl> m_impl;
        std::string m_error;
    };

}  // namespace sdm
