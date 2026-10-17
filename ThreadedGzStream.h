#pragma once

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zstd.h>
#include <zstd_errors.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <istream>
#include <memory>
#include <mutex>
#include <new>
#include <streambuf>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "Bgzf.h"
#include "GzipMember.h"
#include "ParallelGzip.h"

namespace sdm {

    // Ported from protal (github.com/hildebra/protalFastHit, branch audit-fixes, 956e18b,
    // src/IO/ThreadedGzStream.h); sdm opens it through openInputFile (InputFileStream.cpp), which
    // gives each file its threads (set_threads) and turns a read error into an exception.
    //
    // Reads a gzip-compressed (or plain) file as igzstream does, but inflates it in a thread of its
    // own, up to kBlocks blocks ahead of the reader: sdm's reader thread (ifbufstream, which takes
    // the files of a pair and their MID file in step) only copies bytes that are already inflated,
    // and the files of a pair inflate in parallel. A BGZF file (Bgzf.h: bgzip's, sdm's own .gz
    // output) is inflated block by block, any other gzip file streamed, both with ISA-L (igzip,
    // Bgzf.h); a single-member .fq.gz inflates at about 1 GB/s of FASTQ on one thread, 1.6-1.8x
    // zlib-ng's speed (protal docs/claude/2026-10-06-performance-profiling). A zstd file (its frames,
    // one after the other, skippable ones too, any window zstd --long writes) is read with libzstd's
    // streaming decompression. The format is told from the first bytes, not from the name; a file of
    // none of them is read as it is.
    //
    // As with igzstream, a truncated or corrupt file reads as one that ends early; read_failed()
    // tells the two apart once reading has stopped. A BGZF file must end with its end-of-file
    // block, so a cut at a block boundary is found too (it may go on with gzip members of other
    // kinds, as cat a.bgzf.gz b.gz writes); in other gzip files, what follows a member must be
    // another member (or zero bytes of padding), so a damaged member header is an error rather
    // than the end; a zstd file must not end inside a frame, and holds zstd frames only. The file
    // is opened once and read from that descriptor only, so pipes and process substitution
    // (<(zcat ...)) work; only a regular file is checked for BGZF.
    // The buffer is read by one thread at a time (the reader lock serialises it); its inflating
    // thread is stopped by close().
    //
    // A BGZF file's blocks are independent deflate streams, so they can inflate on several threads
    // (set_threads, SetInflateThreads): one thread then reads the compressed blocks in
    // order and hands them out in batches of up to 1 MB of output, worker threads inflate the
    // batches with ISA-L, and the reader takes them in file order, as if one thread had inflated
    // them: the same bytes, and an error (a corrupt block, a cut file) where the file has it, with
    // everything before it read. A regular gzip file of several members whose second member begins
    // within kMaxGzipMember bytes (as tools that compress a file in chunks write it) is read the same
    // way, by members: the dispatcher cuts it where a member's header begins (GzipMember.h) and the
    // workers inflate whole members, each from its start. A header found by chance in a member's
    // data cuts a member wrongly; its batch then fails, and the file from that batch on is read on
    // one thread, as one thread reads it (damage fails there as on one thread). Every batch before
    // it ended exactly where the next began, so what was read before is the file's (protal docs/claude/
    // 2026-10-13-gzip-members-regression). One gzip member cannot be split this way (its deflate blocks are
    // found only by decoding it from the start, and refer back up to 32 KB): another gzip file, if
    // it is a regular file and gets more than one thread, is inflated in pieces on that many threads
    // with ISA-L (ParallelGzip.h: this file's thread decodes in order what the others have not decoded
    // ahead, and copies the output into the blocks). Where small members follow one of its members
    // (cat of a large gzip file and a file of many members), the pieces stop there and the rest is
    // read by members, this file's thread their dispatcher (SwitchToMembers). With one thread, or from
    // a pipe (which would have to be kept), ISA-L inflates it here on one thread. The bytes are the same either way; a cut file,
    // a corrupt member, a wrong CRC-32 or length, a damaged member header or data after the last member
    // fail either way, in other words (what was read before may differ: a corrupt member's bytes up to
    // where the damage shows).
    //
    // A zstd file's frames are independent too: a regular zstd file on several threads (as many as a BGZF file's)
    // is read frame by frame by the dispatcher (whole frames, their sizes from their headers, batched into output
    // blocks; a block grows to hold a larger frame) and the workers decompress them (ZSTD_decompressDCtx). pzstd and
    // sdm's own .zst output (a frame per 1 MiB of content) write many frames; zstd, also with -T, writes one, which is
    // read on one thread all the same. A frame whose header gives no size (a stream's: cat | zstd), or one above
    // kMaxZstdFrame, and the rest of the file after it, are decompressed as a stream on the dispatcher's thread, as
    // one thread does (the whole file on the file's own thread if its first frame is one). The bytes are the same
    // either way; a cut or corrupt file fails either way.
    //
    // Several files (open with a list: a sample's lanes) read as one stream, one after the other, each in
    // its own format: a file whose last line lacks its '\n' gets one, TakeLines ends a batch with its file
    // (so that the mate files of lanes stay in step, see TakeLines), and a truncated or corrupt file, or one
    // that cannot be opened when its turn comes, ends the stream there, named in read_error_message().
    class ThreadedGzStreambuf : public std::streambuf {
    public:
        // Called as each file of the stream is opened (the first by open, the others as the stream reaches them).
        using OnOpen = std::function<void(std::string const& path, ThreadedGzStreambuf const& buf)>;

        static constexpr size_t kBlockSize = size_t{1} << 20;  // inflated bytes per block
        static constexpr size_t kBlocks = 4;                   // with one inflating thread; 2 per worker + 2 with more
        static constexpr uint32_t kInputBuffer = 1u << 17;     // compressed bytes per read()
        static constexpr unsigned kMaxInflateThreads = 16;
        // A zstd frame decompressed whole by a worker (a regular file on several threads) holds at most this: one whose
        // header gives no size, or a larger one, is decompressed as a stream on the file's own thread with the rest of
        // the file (an output block grows to hold one frame of up to this).
        static constexpr size_t kMaxZstdFrame = size_t{64} << 20;
        static constexpr size_t kZstdFrameHeaderMax = 18;  // a zstd frame header's most bytes (zstd.h's static section)
        // The most threads for a gzip file that is not BGZF (ParallelGzip): a piece decoded ahead costs ~2.6
        // times one decoded in order (two passes, and its copies from the unknown window filled in: on FASTQ ~40% of
        // the bytes), so each thread adds less (on 4 cores 1.7, 2.0, 2.35, 2.4 GB/s of FASTQ on 1-4 threads,
        // protal docs/claude/2026-10-12-isal-parallel-inflate).
        static constexpr unsigned kMaxGzipThreads = 8;
        // A regular gzip file that is not BGZF but whose second member begins within this many bytes is read by members
        // on as many threads as a BGZF file (ReadGzipMembersJob). A later member longer than this (no header after it
        // within this many bytes), or holding more than kMaxZstdFrame, is inflated as a stream on the file's own thread
        // with the rest of the file.
        static constexpr size_t kMaxGzipMember = size_t{16} << 20;

        // Threads that inflate a BGZF file (or decompress a zstd file's frames) opened from now on (1: the file's own
        // thread does it; at most kMaxInflateThreads), unless the stream was given its own (set_threads).
        // SDM_INFLATE_THREADS in the environment overrides both, and GzipThreads (for comparing runs).
        static void SetInflateThreads(unsigned threads) {
            InflateThreadsSetting().store(std::clamp(threads, 1u, kMaxInflateThreads), std::memory_order_relaxed);
        }
        static unsigned InflateThreads() {
            if (unsigned const n = InflateThreadsOfTheEnvironment()) return n;
            return InflateThreadsSetting().load(std::memory_order_relaxed);
        }

        // Threads that inflate a gzip file that is not BGZF (a regular file) opened from now on: more than one inflate
        // it in pieces (ParallelGzip); 1: ISA-L on the file's own thread.
        static void SetGzipThreads(unsigned threads) {
            GzipThreadsSetting().store(std::clamp(threads, 1u, kMaxInflateThreads), std::memory_order_relaxed);
        }
        static unsigned GzipThreads() {
            if (unsigned const n = InflateThreadsOfTheEnvironment()) return n;
            return GzipThreadsSetting().load(std::memory_order_relaxed);
        }

        // This stream's threads for the files it opens from now on, instead of SetInflateThreads' and
        // SetGzipThreads' (0 each: theirs); SDM_INFLATE_THREADS still overrides them. sdm sets them per file
        // (openInputFile), as files open at once may want different ones.
        void set_threads(unsigned inflate, unsigned gzip) {
            m_inflate_threads = inflate == 0 ? 0 : std::clamp(inflate, 1u, kMaxInflateThreads);
            m_gzip_threads = gzip == 0 ? 0 : std::clamp(gzip, 1u, kMaxInflateThreads);
        }

        // For tests: kMaxGzipMember for the files opened from now on (0: no file is read by members).
        static void SetMaxGzipMemberBytes(size_t bytes) { MaxGzipMemberSetting().store(bytes, std::memory_order_relaxed); }

        ThreadedGzStreambuf() = default;
        ThreadedGzStreambuf(ThreadedGzStreambuf const&) = delete;
        ThreadedGzStreambuf& operator=(ThreadedGzStreambuf const&) = delete;
        ~ThreadedGzStreambuf() override { close(); }

        // False if the file cannot be read; open_error() then says why.
        bool open(char const* path) {
            return open(std::vector<std::string>{ path });
        }

        // Several files, read one after the other as one stream; each is opened once the one before it has been read.
        // False if the first cannot be read (open_error() says why); a later one that cannot be read ends the stream
        // there (read_failed()).
        bool open(std::vector<std::string> paths) {
            if (is_open() || paths.empty()) return false;
            m_paths = std::move(paths);
            m_path = 0;
            m_files_begun = 0;
            m_files_read.clear();
            m_last_byte = '\n';
            return OpenFile(m_paths.front().c_str());
        }

        // set_on_open(f): f is called as each file is opened from the next open() on.
        void set_on_open(OnOpen on_open) { m_on_open = std::move(on_open); }

        // The file being read (the last one opened, or the one that could not be opened).
        std::string const& current_path() const {
            static std::string const none;
            return m_path < m_paths.size() ? m_paths[m_path] : none;
        }

        // The files the stream reads (one, unless opened with a list).
        size_t files() const { return m_paths.size(); }

    private:
        bool OpenFile(char const* path) {
            if (is_open()) return false;
            m_open_error.clear();
            m_fd = ::open(path, O_RDONLY | O_CLOEXEC);
            if (m_fd < 0) {
                m_open_error = std::strerror(errno);
                return false;
            }
            struct stat st {};
            if (::fstat(m_fd, &st) != 0 || S_ISDIR(st.st_mode)) {
                m_open_error = S_ISDIR(st.st_mode) ? "it is a directory" : std::strerror(errno);
                ::close(m_fd);
                m_fd = -1;
                return false;
            }
            // A regular file is read in order: the kernel's read-ahead doubled (POSIX_FADV_SEQUENTIAL), and more asked for
            // as it is read (ReadUpTo).
            m_regular = S_ISREG(st.st_mode) && ReadAheadWanted();
            m_read_pos = m_read_ahead = 0;
#ifdef POSIX_FADV_SEQUENTIAL
            if (m_regular) ::posix_fadvise(m_fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
            // pread leaves the offset where it is; a pipe cannot be peeked at, and is read as gzip.
            unsigned char head[bgzf::kHeaderBytes] = {};
            ssize_t const peeked = S_ISREG(st.st_mode) ? ::pread(m_fd, head, sizeof(head), 0) : 0;
            m_bgzf = peeked == static_cast<ssize_t>(sizeof(head)) && bgzf::IsBlockHeader(head);
            m_offset = 0;
            m_members = 0;
            m_parallel.reset();
            m_parallel_error.clear();
            // A zstd file's frames are independent, as BGZF blocks are: on as many threads (pzstd and sdm's .zst output
            // write many frames; zstd, also with -T, one, read on one thread all the same). A file whose first frame does
            // not say its size (a stream's, `cat x | zstd`) is read as a stream from the start, on one thread.
            bool const zstd = !m_bgzf && peeked >= 4 && IsZstdFrame(head);
            bool const sized = zstd && ZSTD_getFrameContentSize(head, static_cast<size_t>(peeked)) <= kMaxZstdFrame;
            // A gzip file that is not BGZF whose second member begins near its start: by members, as BGZF.
            bool const gzip = !m_bgzf && !zstd && peeked >= 2 && head[0] == 0x1f && head[1] == 0x8b;
            m_max_member = MaxGzipMemberSetting().load(std::memory_order_relaxed);
            m_gzip_members = gzip && StreamInflateThreads() > 1 && m_max_member > 0 && SecondMemberWithin(static_cast<size_t>(st.st_size));
            m_members_read = 0;
            m_member_restarts = 0;
            unsigned const threads = m_bgzf || zstd || m_gzip_members ? StreamInflateThreads() : StreamGzipThreads();
            m_zstd_frames = sized && threads > 1;
            if (m_bgzf) {
                m_block_in.resize(bgzf::kMaxBlock);
                m_last_empty = false;
            } else if (m_zstd_frames || m_gzip_members) {
                m_block_in.assign(size_t{4} << 20, 0);
                m_in_pos = m_in_end = 0;
                m_in_eof = false;
                m_mode = m_zstd_frames ? Mode::Zstd : Mode::Look;
            } else if (gzip && threads > 1) {
                // A gzip file that is not BGZF, on several threads (ParallelGzip); on one if that cannot be had
                // (m_parallel_error says why). Where small members follow a member, the rest by members (Inflate).
                m_parallel = ParallelGzip::Open(m_fd, threads, m_parallel_error, StreamInflateThreads() > 1 ? m_max_member : 0);
            }
            if (!m_bgzf && !m_zstd_frames && !m_gzip_members && !m_parallel && !StartGzip(0, m_open_error)) {
                ::close(m_fd);
                m_fd = -1;
                return false;
            }
            m_opened_bgzf = m_bgzf;
            m_opened_zstd_frames = m_zstd_frames;
            m_opened_gzip_members = m_gzip_members;
            m_workers = m_bgzf || m_zstd_frames || m_gzip_members ? threads : 1;
            m_parallel_threads = m_parallel ? threads : 0;
            m_parallel_counts = {};
            m_switched_to_members = false;
            m_switched_at = m_restart_at = kNowhere;
            m_format = gzip || m_bgzf ? 'g' : zstd ? 'z' : peeked > 0 ? 'p' : 0;  // a pipe's: once read (FillStream)
            m_opened_threads = m_parallel_threads > 0 ? m_parallel_threads : m_workers;
            // Read in pieces, a file may go on by members (SwitchToMembers): blocks for as many workers.
            unsigned const block_workers = m_parallel && StreamInflateThreads() > 1 ? StreamInflateThreads() : m_workers;
            size_t const blocks = block_workers > 1 ? 2 * size_t{block_workers} + 2 : kBlocks;
            m_blocks.assign(blocks, std::vector<char>(kBlockSize));
            m_jobs.assign(blocks, Job{});
            m_free.clear();
            m_ready.clear();
            m_todo.clear();
            for (size_t i = 0; i < blocks; i++) m_free.push_back(i);
            m_current = kNone;
            m_stop = false;
            m_done = false;
            m_failed = false;
            m_error.clear();
            m_dispatch_error.clear();
            setg(nullptr, nullptr, nullptr);
            if (m_workers > 1) {
                m_thread = std::thread(&ThreadedGzStreambuf::Dispatch, this);
                for (unsigned w = 0; w < m_workers; w++) m_worker_threads.emplace_back(&ThreadedGzStreambuf::Work, this);
            } else {
                m_thread = std::thread(&ThreadedGzStreambuf::Inflate, this);
            }
            m_files_begun++;
            m_last_byte = '\n';
            if (m_on_open) m_on_open(path, *this);
            return true;
        }

        // After the end of a file read whole, the next file of the stream: false after the last, after a file that ended
        // in an error (nothing more is read), or if the next cannot be opened (an error, as for a damaged file).
        bool NextFile() {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (!m_error.empty() || m_failed) return false;
            }
            if (m_path + 1 >= m_paths.size()) return false;
            close();
            m_path++;
            if (OpenFile(m_paths[m_path].c_str())) return true;
            std::lock_guard<std::mutex> lock(m_mutex);
            m_error = "it cannot be read (" + m_open_error + ")";
            return false;
        }

    public:
        // The threads inflating the open file: BGZF blocks, or another gzip file's pieces (1: the file's own thread).
        unsigned inflate_threads() const { return m_parallel_threads > 0 ? m_parallel_threads : m_workers; }

        // Whether the open file began as BGZF (else gzip, zstd or uncompressed).
        bool opened_as_bgzf() const { return m_opened_bgzf; }

        // Whether the open file, zstd, began with its frames decompressed on several threads (inflate_threads()).
        bool opened_as_zstd_frames() const { return m_opened_zstd_frames; }

        // Whether the open file, gzip but not BGZF, is inflated in pieces on several threads (ParallelGzip).
        bool parallel_gzip() const { return m_parallel_threads > 0; }

        // Whether the open file, gzip of several members but not BGZF, began with its members inflated on several threads
        // (inflate_threads()), and how often reading it went back to one thread for the rest of the file (a batch of
        // members that did not inflate as the file was cut: damage, or a header found by chance in a member's data).
        bool opened_as_gzip_members() const { return m_opened_gzip_members; }
        size_t gzip_member_restarts() const { return m_member_restarts; }

        // Whether the open file, read in pieces, went on by members where small members followed a member (once read).
        bool switched_to_gzip_members() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_switched_to_members;
        }

        // How each file of the stream was read, in the order read, once closed (close(), or the stream going on to the
        // next file): its path and, for a run's log, e.g. "gzip of 22606 members, inflated by members on 4 threads".
        struct FileRead {
            std::string path, how;
            unsigned threads = 1;  // the threads that read it: inflating ones and, with workers, the file's own
        };
        std::vector<FileRead> const& files_read() const { return m_files_read; }

        // Why a gzip file that would have been inflated on several threads is inflated on one (empty if none).
        std::string const& parallel_gzip_error() const { return m_parallel_error; }

        // How the open file's pieces were decoded (ParallelGzip::Counts; none if not in pieces), once it is read.
        ParallelGzip::Counts parallel_gzip_counts() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_parallel ? m_parallel->counts() : m_parallel_counts;
        }

        bool is_open() const { return m_fd >= 0; }

        // Why the last open() failed.
        std::string const& open_error() const { return m_open_error; }

        void close() {
            if (!is_open()) return;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            if (m_thread.joinable()) m_thread.join();
            for (auto& worker : m_worker_threads) worker.join();
            m_worker_threads.clear();
            // The file's own thread and its workers (BGZF blocks, zstd frames, gzip members), or its pieces' threads (the
            // file's own among them; on by members, its workers instead), or the file's own thread alone.
            bool const workers = (m_opened_bgzf && m_opened_threads > 1) || m_opened_zstd_frames || m_opened_gzip_members;
            unsigned const threads = workers                    ? m_opened_threads + 1
                                     : m_switched_at != kNowhere ? std::max(m_opened_threads, m_switch_workers + 1)
                                                                 : m_opened_threads;
            m_files_read.push_back({ current_path(), HowRead(), threads });
            m_jobs.clear();
            m_parallel.reset();  // its threads end; its copy of the descriptor is closed
            m_parallel_threads = 0;
            m_inflate.reset();
            if (m_zstd) ZSTD_freeDCtx(m_zstd);
            m_zstd = nullptr;
            ::close(m_fd);
            m_fd = -1;
            setg(nullptr, nullptr, nullptr);
            m_blocks.clear();
        }

        // A truncated or corrupt file reads as one that ends early; this tells them apart.
        bool read_failed() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return !m_error.empty();
        }

        // Why reading failed; with several files, after the name of the file it failed in.
        std::string read_error_message() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_error.empty() || m_paths.size() < 2) return m_error;
            return current_path() + ": " + m_error;
        }

        // Appends up to `lines` whole lines to out, each with its '\n' (a last line without one gets
        // it, as getline reads such a line too), and returns how many. The lines are found with
        // memchr in the inflated block and copied at once, so a reader lock is held for little
        // more than a copy (protal's BufferedFastxReader::LoadBatch; sdm reads with read()).
        // Of several files, a batch ends with its file once it has a line: the batches of the mate files
        // of a sample's lanes end with the same pair if the lanes' mates match, and otherwise a lane's
        // extra reads show (protal's SeqReaderPE fails the sample) instead of pairing the next lane's reads wrongly.
        size_t TakeLines(size_t lines, std::string& out) {
            size_t taken = 0;
            bool partial = false;  // out ends inside a line
            bool ending = false;   // the file ended inside out's last line: the batch ends with that line
            while (taken < lines) {
                if (gptr() == egptr()) {
                    size_t const files = m_files_begun;
                    if (traits_type::eq_int_type(underflow(), traits_type::eof())) break;
                    if (m_files_begun != files) {  // the next file
                        if (partial) ending = true;    // its line gets the '\n' the last file lacked
                        else if (taken > 0) break;
                    }
                }
                char* const begin = gptr();
                char* const end = egptr();
                char* p = begin;
                while (taken < lines) {
                    auto* const newline = static_cast<char*>(std::memchr(p, '\n', static_cast<size_t>(end - p)));
                    if (!newline) {
                        p = end;
                        break;
                    }
                    p = newline + 1;
                    taken++;
                }
                out.append(begin, static_cast<size_t>(p - begin));
                partial = p[-1] != '\n';
                gbump(static_cast<int>(p - begin));
                if (ending && !partial) break;
            }
            if (partial && taken < lines) {
                out.push_back('\n');
                taken++;
            }
            return taken;
        }

    protected:
        // The current file's next block; at its end the next file's (NextFile), after a '\n' if the file's
        // last line lacked one.
        int_type underflow() override {
            if (gptr() < egptr()) return traits_type::to_int_type(*gptr());
            if (egptr() > eback()) m_last_byte = egptr()[-1];
            while (is_open()) {
                int_type const next = UnderflowFile();
                if (!traits_type::eq_int_type(next, traits_type::eof())) return next;
                bool const open_line = m_last_byte != '\n';
                if (!NextFile()) break;
                if (open_line) {
                    setg(m_newline, m_newline, m_newline + 1);
                    return traits_type::to_int_type('\n');
                }
            }
            return traits_type::eof();
        }

    private:
        int_type UnderflowFile() {
            if (!is_open()) return traits_type::eof();
            std::unique_lock<std::mutex> lock(m_mutex);
            if (m_current != kNone) {
                m_free.push_back(m_current);
                m_current = kNone;
                m_cv.notify_all();
            }
            // The blocks in file order: the next one once it is inflated (a worker may finish a later one first).
            while (!m_failed) {
                m_cv.wait(lock, [this] { return (!m_ready.empty() && m_jobs[m_ready.front()].complete) || (m_ready.empty() && m_done); });
                if (m_ready.empty()) {
                    // The end, or what stopped the dispatcher, after everything before it was read.
                    if (!m_dispatch_error.empty() && m_error.empty()) m_error = m_dispatch_error;
                    break;
                }
                size_t const index = m_ready.front();
                m_ready.pop_front();
                Job const& job = m_jobs[index];
                if (!job.error.empty() && job.gzip) {
                    // Members that did not inflate as the file was cut (damaged, or cut wrongly at a header found by
                    // chance): the file from the batch's first member on is read on one thread, which fails where one
                    // thread fails. Every batch before ended where this one begins, at a member's end.
                    uint64_t const from = job.offset;
                    lock.unlock();
                    std::string error;
                    bool const restarted = ReadOnOneThreadFrom(from, error);
                    lock.lock();
                    if (restarted) continue;
                    if (m_error.empty()) m_error = error;
                    m_failed = true;
                    break;
                }
                if (job.gzip) m_members_read += job.members;
                if (!job.error.empty()) {
                    // A corrupt block: what precedes it in the batch is read, then the file ends here.
                    if (m_error.empty()) m_error = job.error;
                    m_failed = true;
                    m_stop = true;
                    m_cv.notify_all();
                }
                if (job.size == 0) {
                    m_free.push_back(index);
                    m_cv.notify_all();
                    continue;
                }
                m_current = index;
                char* data = m_blocks[index].data();
                setg(data, data, data + job.size);
                return traits_type::to_int_type(*gptr());
            }
            setg(nullptr, nullptr, nullptr);
            return traits_type::eof();
        }

    private:
        static constexpr size_t kNone = SIZE_MAX;

        // An output block's batch: with several inflating threads, the compressed BGZF blocks (or zstd frames) to inflate
        // into it.
        struct Job {
            std::vector<unsigned char> in;               // whole BGZF blocks (zstd frames), one after the other
            std::vector<size_t> ends;                    // where each ends in `in`
            bool zstd = false;                           // zstd frames (Work: ZSTD_decompressDCtx), else BGZF blocks
            bool gzip = false;                           // gzip members (Work: InflateGzipMembers), as the dispatcher cut them
            bool last = false;                           // a gzip job to the file's end: zero bytes may follow its members
            size_t members = 0;                          // a gzip job: the members inflated
            size_t first = 0;                            // a zstd job: the number of its first frame in the file (1 on)
            uint64_t offset = 0;                         // the file offset of the first
            size_t size = 0;                             // inflated bytes in the output block
            bool complete = false;                       // inflated (or filled): the reader may take it
            std::string error;                           // the block that could not be inflated, if one
        };


        void Inflate() {
            bool end = false;
            while (!end) {
                size_t index;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return !m_free.empty() || m_stop; });
                    if (m_stop) break;
                    index = m_free.front();
                    m_free.pop_front();
                }
                char* data = m_blocks[index].data();
                size_t size = 0;
                std::string error;
                end = m_parallel ? FillParallel(data, size, error)
                    : m_bgzf     ? FillBgzf(data, size, error)
                                 : FillStream(data, size, error);
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    m_jobs[index].size = size;
                    m_jobs[index].error.clear();
                    m_jobs[index].complete = true;
                    if (size > 0) m_ready.push_back(index);
                    else m_free.push_back(index);
                    if (!error.empty() && m_error.empty()) m_error = error;
                }
                m_cv.notify_all();
                // Read in pieces up to small members: the rest by members, this thread their dispatcher.
                if (end && error.empty() && m_parallel && m_parallel->StoppedAt() != ParallelGzip::kNotStopped) {
                    std::string switch_error;
                    if (SwitchToMembers(switch_error)) {
                        Dispatch();  // which says when it is done
                        return;
                    }
                    std::lock_guard<std::mutex> lock(m_mutex);
                    if (m_error.empty()) m_error = switch_error;
                }
            }
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_done = true;
            }
            m_cv.notify_all();
        }

        // How the open file was read, its threads ended (close, for files_read): its format, its members (BGZF blocks,
        // zstd frames), on how many threads and how, and where it went on by members or on one thread.
        std::string HowRead() const {
            auto const count = [](size_t n, char const* one, char const* more) { return std::to_string(n) + ' ' + (n == 1 ? one : more); };
            std::string const threads = m_opened_threads > 1 ? " on " + std::to_string(m_opened_threads) + " threads" : " on one thread";
            std::string how;
            if (m_opened_bgzf) {
                how = "BGZF, " + count(m_members, "block", "blocks") + ", inflated" + threads;
            } else if (m_opened_zstd_frames) {
                how = "zstd, " + count(m_members, "frame", "frames") + ", decompressed" + threads;
            } else if (m_opened_gzip_members) {
                how = "gzip of " + count(m_members, "member", "members") + ", inflated by members" + threads;
            } else if (m_parallel_threads > 0) {
                how = "gzip of " + count(m_parallel ? m_parallel->MembersBegun() : m_members, "member", "members") +
                      ", inflated in pieces" + threads;
                if (m_switched_at != kNowhere) how += ", by members from byte " + std::to_string(m_switched_at) + " on";
            } else if (m_format == 'g') {
                how = "gzip of " + count(m_members, "member", "members") + ", inflated on one thread";
            } else if (m_format == 'z') {
                how = "zstd, " + count(m_members, "frame", "frames") + ", decompressed on one thread";
            } else {
                how = "uncompressed";
            }
            if (m_restart_at != kNowhere) {
                how += ", on one thread from byte " + std::to_string(m_restart_at) + " on (members that did not inflate as cut)";
            }
            if (!m_parallel_error.empty()) how += " (not in pieces: " + m_parallel_error + ")";
            return how;
        }

        // A gzip file read in pieces whose reader stopped before a member that small members follow (ParallelGzip::
        // StoppedAt), every member before it read and checked: the rest of the file read by members from there (as a
        // file of many members from its start), on as many workers as a BGZF file, this thread their dispatcher. The
        // pieces' threads end and the mapping goes first. False if the file cannot be read on from there (error set).
        bool SwitchToMembers(std::string& error) {
            uint64_t const at = m_parallel->StoppedAt();
            size_t const members = m_parallel->MembersRead();
            std::unique_ptr<ParallelGzip> pieces;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_parallel_counts = m_parallel->counts();
                pieces = std::move(m_parallel);
            }
            pieces.reset();
            if (::lseek(m_fd, static_cast<off_t>(at), SEEK_SET) < 0) {
                error = std::string("seek error: ") + std::strerror(errno);
                return false;
            }
            m_offset = at;
            m_read_pos = m_read_ahead = at;
            m_members = members;
            m_block_in.assign(size_t{4} << 20, 0);
            m_in_pos = m_in_end = 0;
            m_in_eof = false;
            m_gzip_members = true;
            unsigned const workers = StreamInflateThreads();
            m_switch_workers = workers;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_members_read = members;
                m_switched_to_members = true;
                m_switched_at = at;
                if (m_stop) return true;  // closing: the dispatcher ends at once
            }
            for (unsigned w = 0; w < workers; w++) m_worker_threads.emplace_back(&ThreadedGzStreambuf::Work, this);
            return true;
        }

        // With several inflating threads (a BGZF file): reads the blocks in order into batches (ReadBgzfJob), one per
        // free output block, and queues them for the workers (Work) and, in the same order, for the reader. If a gzip
        // member of another kind follows the BGZF blocks, it is inflated here, as Inflate does, into the next blocks.
        // What stops the reading (a cut file, a damaged header, a missing end-of-file block) is the reader's error once
        // it has read every block before it.
        void Dispatch() {
            bool end = false;
            while (!end) {
                size_t index;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return !m_free.empty() || m_stop; });
                    if (m_stop) break;
                    index = m_free.front();
                    m_free.pop_front();
                }
                Job& job = m_jobs[index];  // this thread's until it is queued
                job.complete = false;
                job.size = 0;
                job.error.clear();
                std::string error;
                bool const frames = m_bgzf || m_zstd_frames || m_gzip_members;  // whole blocks (frames, members) for the workers
                job.zstd = m_zstd_frames;
                job.gzip = m_gzip_members;
                if (m_bgzf) {
                    end = ReadBgzfJob(job, error);
                } else if (m_zstd_frames) {
                    end = ReadZstdJob(job, m_blocks[index], error);
                } else if (m_gzip_members) {
                    end = ReadGzipMembersJob(job, m_blocks[index], error);
                } else {
                    size_t size = 0;
                    end = FillStream(m_blocks[index].data(), size, error);
                    job.size = size;
                    job.complete = true;
                }
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    if (frames && !job.ends.empty()) {
                        m_ready.push_back(index);
                        m_todo.push_back(index);
                    } else if (!frames && job.size > 0) {
                        m_ready.push_back(index);
                    } else {
                        m_free.push_back(index);
                    }
                    if (!error.empty() && m_dispatch_error.empty()) m_dispatch_error = error;
                }
                m_cv.notify_all();
            }
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_done = true;
            }
            m_cv.notify_all();
        }

        // A worker: inflates queued batches block by block into their output blocks (bgzf::DecompressBlock, ISA-L, a
        // state per thread), until the dispatcher is done and nothing is queued, or the file is closed. A corrupt
        // block ends its batch there, with the error and the byte where the block starts.
        void Work() {
            while (true) {
                size_t index;
                {
                    std::unique_lock<std::mutex> lock(m_mutex);
                    m_cv.wait(lock, [this] { return !m_todo.empty() || m_done || m_stop; });
                    if (m_stop || m_todo.empty()) return;
                    index = m_todo.front();
                    m_todo.pop_front();
                }
                Job& job = m_jobs[index];
                char* out = m_blocks[index].data();
                size_t const room = m_blocks[index].size();  // a zstd job's block may have grown for a large frame
                size_t size = 0, begin = 0, members = 0;
                uint64_t offset = job.offset;
                std::string error;
                // A gzip job's members are inflated one after the other from the batch's start, whatever headers the
                // dispatcher cut at within it: only where the batch begins and ends must be a member's start and end.
                if (job.gzip && !InflateGzipMembers(job.in.data(), job.in.size(), out, room, size, members, job.last, error)) {
                    error += " in the gzip members from byte " + std::to_string(offset);
                }
                for (size_t i = 0; i < job.ends.size() && !job.gzip; i++) {
                    size_t const end = job.ends[i];
                    size_t n = 0;
                    if (job.zstd ? !DecompressZstdFrame(job.in.data() + begin, end - begin, out + size, room - size, n, error,
                                                        job.first + i)
                                 : !bgzf::DecompressBlock(job.in.data() + begin, end - begin, out + size, room - size, n,
                                                          error)) {
                        error += " at byte " + std::to_string(offset);
                        break;
                    }
                    size += n;
                    offset += end - begin;
                    begin = end;
                }
                {
                    std::lock_guard<std::mutex> lock(m_mutex);
                    job.size = size;
                    job.members = members;
                    job.error = std::move(error);
                    job.complete = true;
                }
                m_cv.notify_all();
            }
        }

        // FillBgzf without the inflating: reads whole BGZF blocks into job (job.in, each block's end in job.ends)
        // while their contents, as their footers give them, leave room for one more block of 64 KB in an output block;
        // true at the end of the file or on an error (then set). A gzip member of another kind after the blocks
        // switches the file to FillStream (StartGzip), and false is returned.
        bool ReadBgzfJob(Job& job, std::string& error) {
            job.in.clear();
            job.ends.clear();
            job.offset = m_offset;
            size_t out = 0;
            while (out + bgzf::kMaxBlock <= kBlockSize) {
                unsigned char* head = m_block_in.data();
                size_t const got = ReadUpTo(head, bgzf::kHeaderBytes, error);
                if (!error.empty()) return true;
                if (got == 0) {
                    if (!m_last_empty) error = "the BGZF end-of-file block is missing (truncated file?)";
                    return true;
                }
                if (got < bgzf::kHeaderBytes || !bgzf::IsBlockHeader(head)) {
                    if (got >= 2 && head[0] == 0x1f && head[1] == 0x8b) return !StartGzip(got, error);
                    error = got < bgzf::kHeaderBytes ? "the file ends inside a BGZF block header (truncated file?)"
                                                     : "no BGZF block at byte " + std::to_string(m_offset) +
                                                       " (a damaged block header?)";
                    return true;
                }
                size_t const block_size = bgzf::BlockSize(head);
                if (block_size < bgzf::kHeaderBytes + bgzf::kFooterBytes) {
                    error = "an invalid BGZF block size at byte " + std::to_string(m_offset);
                    return true;
                }
                size_t const start = job.in.size();
                job.in.resize(start + block_size);
                std::memcpy(job.in.data() + start, head, bgzf::kHeaderBytes);
                size_t const rest = block_size - bgzf::kHeaderBytes;
                if (ReadUpTo(job.in.data() + start + bgzf::kHeaderBytes, rest, error) != rest) {
                    job.in.resize(start);
                    if (error.empty()) error = "the file ends inside a BGZF block (truncated file?)";
                    return true;
                }
                uint32_t const isize = bgzf::LE32(job.in.data() + start + block_size - 4);
                job.ends.push_back(job.in.size());
                m_last_empty = isize == 0;
                m_members++;
                m_offset += block_size;
                out += isize;
            }
            return false;
        }

        // ReadBgzfJob for a zstd file (m_zstd_frames): whole frames into job (job.in, each frame's end in job.ends) while
        // their contents, as their headers give them, fit the output block `block` (one frame at least: the block grows
        // to hold a larger one, up to kMaxZstdFrame); true at the end of the file or on an error (then set). A frame whose
        // header gives no size (a stream's, written without knowing its length) or a larger one switches the rest of the
        // file to FillStream on this thread (the frames read before it first go to the workers): false is returned.
        // Skippable frames hold nothing.
        bool ReadZstdJob(Job& job, std::vector<char>& block, std::string& error) {
            job.in.clear();
            job.ends.clear();
            job.offset = m_offset;
            job.first = m_members + 1;
            size_t planned = 0;
            while (true) {
                if (!Available(kZstdFrameHeaderMax, error)) return true;
                size_t const left = m_in_end - m_in_pos;
                if (left == 0) break;  // the end of the file, after a whole frame
                unsigned char const* in = m_block_in.data() + m_in_pos;
                if (left < 4 || !IsZstdFrame(in)) {
                    error = left < 4 ? "the file ends inside zstd frame " + std::to_string(m_members + 1) + " (truncated file?)"
                                     : "the data at byte " + std::to_string(m_offset) + ", after zstd frame " +
                                           std::to_string(m_members) + ", is no zstd frame (corrupt file?)";
                    return true;
                }
                unsigned long long const content = ZSTD_getFrameContentSize(in, left);
                if (content == ZSTD_CONTENTSIZE_ERROR) {
                    error = "a corrupt header of zstd frame " + std::to_string(m_members + 1) + " at byte " +
                            std::to_string(m_offset) + " (corrupt or truncated file?)";
                    return true;
                }
                if (content == ZSTD_CONTENTSIZE_UNKNOWN || content > kMaxZstdFrame) {
                    if (!job.ends.empty()) break;  // the frames before it go to the workers first
                    if (!StartZstd(error)) return true;
                    m_zstd_frames = false;  // the rest on this thread, streamed (FillStream, Mode::Zstd)
                    return false;
                }
                if (!job.ends.empty() && planned + content > kBlockSize) break;  // the next job's
                // The whole frame in m_block_in: its compressed size is known only once all of it is there.
                size_t frame;
                while (true) {
                    frame = ZSTD_findFrameCompressedSize(m_block_in.data() + m_in_pos, m_in_end - m_in_pos);
                    if (!ZSTD_isError(frame)) break;
                    size_t const had = m_in_end - m_in_pos;
                    if (ZSTD_getErrorCode(frame) != ZSTD_error_srcSize_wrong || m_in_eof) {
                        error = m_in_eof && ZSTD_getErrorCode(frame) == ZSTD_error_srcSize_wrong
                                    ? "the file ends inside zstd frame " + std::to_string(m_members + 1) + " (truncated file?)"
                                    : std::string(ZSTD_getErrorName(frame)) + " in zstd frame " +
                                          std::to_string(m_members + 1) + " (corrupt file?) at byte " +
                                          std::to_string(m_offset);
                        return true;
                    }
                    if (m_block_in.size() < 2 * had) m_block_in.resize(2 * had);
                    if (!Available(2 * had, error)) return true;
                }
                job.in.insert(job.in.end(), m_block_in.data() + m_in_pos, m_block_in.data() + m_in_pos + frame);
                job.ends.push_back(job.in.size());
                m_in_pos += frame;
                m_offset += frame;
                m_members++;
                planned += static_cast<size_t>(content);
            }
            if (block.size() < planned) block.resize(planned);
            return job.ends.empty();
        }

        // A whole zstd frame (number `number` in the file) of `size` bytes decompressed into out (at most `room` bytes) on
        // this thread's context; false on an error (then set), as a frame that is not whole or holds more than its header
        // says.
        static bool DecompressZstdFrame(unsigned char const* frame, size_t size, char* out, size_t room, size_t& made,
                                        std::string& error, size_t number) {
            struct Free {
                void operator()(ZSTD_DCtx* d) const { ZSTD_freeDCtx(d); }
            };
            thread_local std::unique_ptr<ZSTD_DCtx, Free> dctx;
            if (!dctx) {
                dctx.reset(ZSTD_createDCtx());
                if (!dctx) {
                    error = "cannot allocate a zstd decompression context";
                    return false;
                }
                ZSTD_DCtx_setParameter(dctx.get(), ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
            }
            size_t const n = ZSTD_decompressDCtx(dctx.get(), out, room, frame, size);
            if (ZSTD_isError(n)) {
                error = std::string(ZSTD_getErrorName(n)) + " in zstd frame " + std::to_string(number) + " (corrupt file?)";
                return false;
            }
            made = n;
            return true;
        }

        static constexpr size_t kMaxGzipBatch = size_t{4} << 20;  // a batch's most compressed bytes (of more than one member)

        // Whether the regular file of `size` bytes being opened has a second gzip member within m_max_member bytes of its
        // start (gzip_member::Next): then it is read by members.
        bool SecondMemberWithin(size_t size) {
            size_t const n = std::min(size, m_max_member + gzip_member::kMaxHeader + 8);
            std::unique_ptr<unsigned char[]> head(new (std::nothrow) unsigned char[n]);
            if (!head) return false;
            size_t got = 0;
            while (got < n) {
                ssize_t const r = ::pread(m_fd, head.get() + got, n - got, static_cast<off_t>(got));
                if (r < 0 && errno == EINTR) continue;
                if (r <= 0) break;
                got += static_cast<size_t>(r);
            }
            return gzip_member::Next(head.get(), got, 0, m_max_member) > 0;  // a header cut off at `got` counts as none
        }

        // The bytes of the gzip member at the unread input (m_in_pos), up to the next member's header
        // (gzip_member::StartsAt) or the file's end, all in m_block_in then; kNone if no header follows within m_max_member
        // bytes (the member is longer). False on a read error (then set).
        bool MemberAt(size_t& member, std::string& error) {
            size_t from = gzip_member::kMinBytes;
            while (true) {
                size_t const left = m_in_end - m_in_pos;
                unsigned char const* const base = m_block_in.data() + m_in_pos;
                size_t const end = std::min(left, m_max_member + 1);
                bool told = true;  // every 0x1f before `from` is no header
                while (from < end) {
                    auto const* hit = static_cast<unsigned char const*>(std::memchr(base + from, 0x1f, end - from));
                    if (!hit) {
                        from = end;
                        break;
                    }
                    size_t const at = static_cast<size_t>(hit - base);
                    int const verdict = gzip_member::StartsAt(base, left, at, m_in_eof);
                    if (verdict > 0) {
                        member = at;
                        return true;
                    }
                    from = verdict < 0 ? at : at + 1;
                    if (verdict < 0) {
                        told = false;
                        break;
                    }
                }
                if (m_in_eof) {
                    member = left <= m_max_member ? left : kNone;
                    return true;
                }
                if (told && left > m_max_member) {
                    member = kNone;
                    return true;
                }
                // More of the file: m_block_in grows to hold the member and the header after it.
                size_t const want = std::min(left + std::max(left, size_t{1} << 20), m_max_member + gzip_member::kMaxHeader + 8);
                if (m_block_in.size() < want) m_block_in.resize(want);
                if (!Available(want, error)) return false;
            }
        }

        // ReadZstdJob for a gzip file of several members (m_gzip_members): whole members into job (job.in, each one's end
        // in job.ends) while their contents, as their trailers give them, fit the output block `block` (one member at least:
        // the block grows to hold a larger one, up to kMaxZstdFrame). A member ends where the next member's header begins
        // (MemberAt), so a header found by chance cuts a member wrongly: its batch then fails, and the reader reads the file
        // from that batch on on one thread (UnderflowFile). True at the end of the file or on a read error (then set). A
        // member longer than m_max_member, or holding more than kMaxZstdFrame, switches the rest of the file to FillStream
        // on this thread (the members before it first go to the workers): false is returned.
        bool ReadGzipMembersJob(Job& job, std::vector<char>& block, std::string& error) {
            job.in.clear();
            job.ends.clear();
            job.offset = m_offset;
            job.first = m_members + 1;
            job.last = false;
            job.members = 0;
            size_t planned = 0;
            while (true) {
                if (!Available(1, error)) return true;
                if (m_in_end == m_in_pos) break;  // the end of the file
                size_t member = 0;
                if (!MemberAt(member, error)) return true;
                unsigned char const* const start = m_block_in.data() + m_in_pos;
                bool const last = member == m_in_end - m_in_pos && m_in_eof;  // to the file's end
                uint32_t const isize = member != kNone && member >= 8 ? bgzf::LE32(start + member - 4) : 0;
                if (member == kNone || isize > kMaxZstdFrame) {
                    if (!job.ends.empty()) break;  // the members before it go to the workers first
                    return !StreamTheRest(error);
                }
                if (!job.ends.empty() && (planned + isize > kBlockSize || job.in.size() + member > kMaxGzipBatch)) break;
                job.in.insert(job.in.end(), start, start + member);
                job.ends.push_back(job.in.size());
                m_in_pos += member;
                m_offset += member;
                m_members++;
                planned += isize;
                if (last) {
                    job.last = true;
                    break;
                }
            }
            if (block.size() < planned) block.resize(planned);
            return job.ends.empty();
        }

        // The rest of the file, from the unread input on, inflated as a stream on this thread (FillStream).
        bool StreamTheRest(std::string& error) {
            try {
                m_inflate = std::make_unique<inflate_state>();  // ~85 KB
            } catch (std::bad_alloc const&) {
                error = "out of memory for the gzip inflate state";
                return false;
            }
            isal_inflate_init(m_inflate.get());
            m_mode = Mode::Look;
            m_gzip_members = false;
            return true;
        }

        // A batch of whole gzip members, one after the other, inflated into out (at most `room` bytes) on this thread's
        // state, each member's CRC-32 and length checked; false (error set) unless they end exactly at the batch's end: a
        // member that is damaged, that goes on past the batch (cut at a header found by chance), or that holds more than
        // there is room for. In the file's last batch zero bytes may follow the members.
        static bool InflateGzipMembers(unsigned char const* in, size_t size, char* out, size_t room, size_t& made,
                                       size_t& members, bool last, std::string& error) {
            thread_local std::unique_ptr<inflate_state> const state = std::make_unique<inflate_state>();  // ~85 KB
            inflate_state& st = *state;
            made = 0;
            members = 0;
            size_t at = 0;
            while (at < size) {
                if (last && std::all_of(in + at, in + size, [](unsigned char c) { return c == 0; })) break;
                isal_inflate_init(&st);
                st.crc_flag = ISAL_GZIP_NO_HDR_VER;
                isal_gzip_header header;
                isal_gzip_header_init(&header);  // no buffers: the name, comment and extra field are skipped
                st.next_in = const_cast<uint8_t*>(in + at);
                st.avail_in = static_cast<uint32_t>(size - at);
                if (isal_read_gzip_header(&st, &header) != ISAL_DECOMP_OK) {
                    error = "no gzip member header";
                    return false;
                }
                st.next_out = reinterpret_cast<uint8_t*>(out + made);
                st.avail_out = static_cast<uint32_t>(room - made);
                int const ret = isal_inflate(&st);
                made = static_cast<size_t>(reinterpret_cast<char*>(st.next_out) - out);
                if (ret != ISAL_DECOMP_OK) {
                    error = InflateError(ret);
                    return false;
                }
                if (st.block_state != ISAL_BLOCK_FINISH) {
                    error = made == room ? "a member larger than its trailer says" : "a member that goes on past its batch";
                    return false;
                }
                at = size - st.avail_in;  // ISA-L gives back what it read past the trailer
                members++;
            }
            return true;
        }

        // After a batch of gzip members that did not inflate as the file was cut: the dispatcher and the workers stopped,
        // the batches queued dropped, and the file read on one thread from `offset` on (the batch's start, a member's),
        // as Inflate reads a gzip file. False if that cannot be had (error then set).
        bool ReadOnOneThreadFrom(uint64_t offset, std::string& error) {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = true;
            }
            m_cv.notify_all();
            if (m_thread.joinable()) m_thread.join();
            for (auto& worker : m_worker_threads) worker.join();
            m_worker_threads.clear();
            if (::lseek(m_fd, static_cast<off_t>(offset), SEEK_SET) < 0) {
                error = std::string("seek error: ") + std::strerror(errno);
                return false;
            }
            m_gzip_members = false;
            m_workers = 1;
            m_offset = offset;
            m_read_pos = m_read_ahead = offset;
            m_members = m_members_read;
            m_restart_at = offset;
            if (!StartGzip(0, error)) return false;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stop = false;
                m_done = false;
                m_dispatch_error.clear();
                m_ready.clear();
                m_todo.clear();
                m_free.clear();
                for (size_t i = 0; i < m_jobs.size(); i++) {
                    m_jobs[i] = Job{};
                    m_free.push_back(i);
                }
                m_current = kNone;
            }
            m_member_restarts++;
            m_thread = std::thread(&ThreadedGzStreambuf::Inflate, this);
            return true;
        }

        // Reads up to size bytes at the file's current offset; the bytes read (fewer only at its end). Of a regular file,
        // the kReadAhead bytes after them are asked for (POSIX_FADV_WILLNEED) kReadAheadStep at a time: a file not in the
        // page cache (a sample's first read, from a network file system) is then read ahead of the reads.
        size_t ReadUpTo(unsigned char* dst, size_t size, std::string& error) {
            size_t done = 0;
            while (done < size) {
                ssize_t const n = ::read(m_fd, dst + done, size - done);
                if (n < 0 && errno == EINTR) continue;
                if (n < 0) {
                    error = std::string("read error: ") + std::strerror(errno);
                    break;
                }
                if (n == 0) break;
                done += static_cast<size_t>(n);
            }
            m_read_pos += done;
            if (m_regular && m_read_pos + kReadAhead >= m_read_ahead + kReadAheadStep) {
                uint64_t const from = std::max(m_read_ahead, m_read_pos);
#ifdef POSIX_FADV_WILLNEED
                ::posix_fadvise(m_fd, static_cast<off_t>(from), static_cast<off_t>(m_read_pos + kReadAhead - from), POSIX_FADV_WILLNEED);
#endif
                m_read_ahead = m_read_pos + kReadAhead;
            }
            return done;
        }

        // Inflates BGZF blocks into data (from size) while a whole block still fits; true at the end
        // of the file or on an error (then set).
        bool FillBgzf(char* data, size_t& size, std::string& error) {
            while (size + bgzf::kMaxBlock <= kBlockSize) {
                unsigned char* block = m_block_in.data();
                size_t const got = ReadUpTo(block, bgzf::kHeaderBytes, error);
                if (!error.empty()) return true;
                if (got == 0) {
                    if (!m_last_empty) error = "the BGZF end-of-file block is missing (truncated file?)";
                    return true;
                }
                if (got < bgzf::kHeaderBytes || !bgzf::IsBlockHeader(block)) {
                    // A gzip member of another kind after the BGZF blocks (cat a.bgzf.gz b.gz): on as gzip.
                    if (got >= 2 && block[0] == 0x1f && block[1] == 0x8b) {
                        return !StartGzip(got, error) || FillStream(data, size, error);
                    }
                    error = got < bgzf::kHeaderBytes ? "the file ends inside a BGZF block header (truncated file?)"
                                                     : "no BGZF block at byte " + std::to_string(m_offset) +
                                                       " (a damaged block header?)";
                    return true;
                }
                size_t const block_size = bgzf::BlockSize(block);
                if (block_size < bgzf::kHeaderBytes + bgzf::kFooterBytes) {
                    error = "an invalid BGZF block size at byte " + std::to_string(m_offset);
                    return true;
                }
                size_t const rest = block_size - bgzf::kHeaderBytes;
                if (ReadUpTo(block + bgzf::kHeaderBytes, rest, error) != rest) {
                    if (error.empty()) error = "the file ends inside a BGZF block (truncated file?)";
                    return true;
                }
                size_t n = 0;
                if (!bgzf::DecompressBlock(block, block_size, data + size, kBlockSize - size, n, error)) {
                    error += " at byte " + std::to_string(m_offset);
                    return true;
                }
                m_last_empty = n == 0;
                m_members++;
                size += n;
                m_offset += block_size;
            }
            return false;
        }

        // Copies the file's next inflated bytes from its pieces (m_parallel) into data (from size) until the block is
        // full; true at the end of the file or on an error (then set).
        bool FillParallel(char* data, size_t& size, std::string& error) {
            while (size < kBlockSize) {
                size_t const n = m_parallel->Read(data + size, kBlockSize - size);
                if (n == 0) {
                    if (!m_parallel->Error().empty()) error = "corrupt gzip data (corrupt or truncated file?): " + m_parallel->Error();
                    return true;
                }
                size += n;
            }
            return false;
        }

        // Reads the rest of the file with FillStream, from its first `got` bytes in m_block_in (read
        // there by FillBgzf, else none); false if the inflate state cannot be had (error then set).
        bool StartGzip(size_t got, std::string& error) {
            m_bgzf = false;
            m_block_in.resize(std::max<size_t>(m_block_in.size(), kInputBuffer));
            m_in_pos = 0;
            m_in_end = got;
            m_in_eof = got > 0 && got < bgzf::kHeaderBytes;  // FillBgzf read fewer only at the end
            m_mode = Mode::Look;
            try {
                m_inflate = std::make_unique<inflate_state>();  // ~85 KB
            } catch (std::bad_alloc const&) {
                error = "out of memory for the gzip inflate state";
                return false;
            }
            isal_inflate_init(m_inflate.get());
            return true;
        }

        // A gzip member begins at the unread input: its header is read next, then its deflate data and trailer,
        // which ISA-L checks (CRC-32 and length). The header is read with a state of its own (m_gz_header):
        // isal_inflate's own header reading (crc_flag ISAL_GZIP) fails on a header with a header CRC, or with two
        // of name, comment and extra field, that is split between two reads of the file (ISA-L 2.32.1).
        void StartMember() {
            isal_inflate_reset(m_inflate.get());
            m_inflate->crc_flag = ISAL_GZIP_NO_HDR_VER;
            isal_gzip_header_init(&m_gz_header);  // no buffers: the name, comment and extra field are skipped
            m_header_read = false;
        }

        // What an ISA-L error code says, as zlib words it.
        static char const* InflateError(int ret) {
            switch (ret) {
                case ISAL_INVALID_BLOCK: return "invalid block type";
                case ISAL_INVALID_SYMBOL: return "invalid code";
                case ISAL_INVALID_LOOKBACK: return "invalid distance too far back";
                case ISAL_INVALID_WRAPPER: return "incorrect header check";
                case ISAL_UNSUPPORTED_METHOD: return "unknown compression method";
                case ISAL_INCORRECT_CHECKSUM: return "incorrect data or length check";
                default: return "corrupt gzip data";
            }
        }

        // A zstd frame's magic number, or a skippable frame's (pzstd and the seekable format write them).
        static bool IsZstdFrame(unsigned char const* in) {
            return (in[0] == 0x28 && in[1] == 0xb5 && in[2] == 0x2f && in[3] == 0xfd) ||
                   ((in[0] & 0xf0) == 0x50 && in[1] == 0x2a && in[2] == 0x4d && in[3] == 0x18);
        }

        // FillZstd's decompression context; false if it cannot be made (error then set).
        bool StartZstd(std::string& error) {
            m_zstd = ZSTD_createDCtx();
            if (!m_zstd) {
                error = "cannot allocate a zstd decompression context";
                return false;
            }
            // The largest windows (zstd --long up to 2 GB); the memory follows the frame's own.
            ZSTD_DCtx_setParameter(m_zstd, ZSTD_d_windowLogMax, ZSTD_dParam_getBounds(ZSTD_d_windowLogMax).upperBound);
            m_zstd_in_frame = false;
            return true;
        }

        // Makes at least `want` unread bytes available in m_block_in unless the file ends first; false
        // on a read error (then set). Unread bytes move to the front first.
        bool Available(size_t want, std::string& error) {
            if (m_in_end - m_in_pos >= want || m_in_eof) return true;
            std::memmove(m_block_in.data(), m_block_in.data() + m_in_pos, m_in_end - m_in_pos);
            m_in_end -= m_in_pos;
            m_in_pos = 0;
            while (m_in_end < want && !m_in_eof) {
                size_t const got = ReadUpTo(m_block_in.data() + m_in_end, m_block_in.size() - m_in_end, error);
                if (!error.empty()) return false;
                if (got == 0) m_in_eof = true;
                m_in_end += got;
            }
            return true;
        }

        // The rest of the file is zero bytes (padding, which gzip ignores too).
        bool RestIsZero(std::string& error) {
            while (true) {
                for (size_t i = m_in_pos; i < m_in_end; i++) {
                    if (m_block_in[i] != 0) return false;
                }
                m_offset += m_in_end - m_in_pos;
                m_in_pos = m_in_end;
                if (!Available(1, error) || m_in_end == m_in_pos) return error.empty();
            }
        }

        // Inflates gzip members (decompresses a zstd file, or copies a file that is neither) into data
        // (from size) until the block is full; true at the end of the file or on an error (then set).
        // After a gzip member, the file must end, hold zero padding, or go on with another member:
        // anything else is taken for a damaged member, not ignored as zlib's gzread does.
        bool FillStream(char* data, size_t& size, std::string& error) {
            while (size < kBlockSize) {
                if (m_mode == Mode::Look) {
                    if (!Available(4, error)) return true;
                    size_t const left = m_in_end - m_in_pos;
                    if (left == 0) return true;  // the end of the file
                    unsigned char const* in = m_block_in.data() + m_in_pos;
                    if (left >= 2 && in[0] == 0x1f && in[1] == 0x8b) {
                        StartMember();
                        m_members++;
                        m_mode = Mode::Gzip;
                        if (m_format == 0) m_format = 'g';
                    } else if (m_members == 0 && left >= 4 && IsZstdFrame(in)) {
                        if (!StartZstd(error)) return true;
                        m_mode = Mode::Zstd;
                        if (m_format == 0) m_format = 'z';
                    } else if (m_members == 0) {
                        m_mode = Mode::Copy;  // neither gzip nor zstd: read as it is
                        if (m_format == 0) m_format = 'p';
                    } else {
                        uint64_t const at = m_offset;
                        if (!RestIsZero(error) && error.empty()) {
                            error = "the data at byte " + std::to_string(at) + ", after gzip member " + std::to_string(m_members) +
                                    ", is no gzip member (a damaged member header?)";
                        }
                        return true;
                    }
                }
                if (m_mode == Mode::Zstd) return FillZstd(data, size, error);
                if (m_mode == Mode::Copy) {
                    if (!Available(1, error)) return true;
                    size_t const n = std::min(m_in_end - m_in_pos, kBlockSize - size);
                    if (n == 0) return true;
                    std::memcpy(data + size, m_block_in.data() + m_in_pos, n);
                    m_in_pos += n;
                    m_offset += n;
                    size += n;
                    continue;
                }
                // Mode::Gzip: the member's header, then its deflate data and trailer. At the member's end ISA-L
                // gives back the bytes it read past the trailer, so the unread input starts where the member ends.
                if (!Available(1, error)) return true;
                if (m_in_end == m_in_pos) {
                    error = "the file ends inside gzip member " + std::to_string(m_members) + " (truncated file?)";
                    return true;
                }
                inflate_state& state = *m_inflate;
                size_t const had = m_in_end - m_in_pos;
                size_t const size_before = size;
                state.next_in = m_block_in.data() + m_in_pos;
                state.avail_in = static_cast<uint32_t>(had);
                int ret;
                if (!m_header_read) {
                    ret = isal_read_gzip_header(&state, &m_gz_header);
                    m_header_read = ret == ISAL_DECOMP_OK;
                    if (ret > 0) ret = ISAL_DECOMP_OK;  // ISAL_END_INPUT: more of the header to come
                } else {
                    state.next_out = reinterpret_cast<uint8_t*>(data + size);
                    state.avail_out = static_cast<uint32_t>(kBlockSize - size);
                    ret = isal_inflate(&state);
                    size = kBlockSize - state.avail_out;
                }
                size_t const used = had - state.avail_in;
                m_in_pos += used;
                m_offset += used;
                if (ret != ISAL_DECOMP_OK) {
                    error = std::string(InflateError(ret)) + " in gzip member " + std::to_string(m_members) + " (corrupt file?)";
                    return true;
                }
                if (m_header_read && state.block_state == ISAL_BLOCK_FINISH) {
                    m_mode = Mode::Look;
                } else if (used == 0 && size == size_before && size < kBlockSize) {
                    // Nothing consumed or made although there is room and input: more input is needed.
                    std::string more;
                    if (!Available(had + 1, more)) { error = more; return true; }
                    if (m_in_end - m_in_pos == had) {
                        error = "the file ends inside gzip member " + std::to_string(m_members) + " (truncated file?)";
                        return true;
                    }
                }
            }
            return false;
        }

        // Decompresses zstd frames into data (from size) until the block is full; true at the end of the
        // file or on an error (then set). ZSTD_decompressStream reads frame after frame, skippable ones
        // too, and returns 0 where one ends: the file must end there, not inside a frame.
        bool FillZstd(char* data, size_t& size, std::string& error) {
            while (size < kBlockSize) {
                if (!Available(1, error)) return true;
                size_t const had = m_in_end - m_in_pos;
                if (had == 0 && !m_zstd_in_frame) return true;  // the end of the file, after a whole frame
                ZSTD_inBuffer in {m_block_in.data() + m_in_pos, had, 0};
                ZSTD_outBuffer out {data + size, kBlockSize - size, 0};
                size_t const ret = ZSTD_decompressStream(m_zstd, &out, &in);
                m_in_pos += in.pos;
                m_offset += in.pos;
                size += out.pos;
                if (ZSTD_isError(ret)) {
                    error = std::string(ZSTD_getErrorName(ret)) + " in zstd frame " + std::to_string(m_members + 1) +
                            " (corrupt file?)";
                    return true;
                }
                if (in.pos == 0 && out.pos == 0) {
                    // No progress: the frame needs more input than there is unread (or, with none, it is cut).
                    std::string more;
                    if (!Available(had + 1, more)) { error = more; return true; }
                    if (m_in_end - m_in_pos == had) {
                        error = "the file ends inside zstd frame " + std::to_string(m_members + 1) + " (truncated file?)";
                        return true;
                    }
                    continue;
                }
                m_zstd_in_frame = ret != 0;
                if (ret == 0) m_members++;  // a frame ended
            }
            return false;
        }

        enum class Mode { Look, Gzip, Zstd, Copy };
        std::vector<std::string> m_paths;                // the stream's files, read one after the other
        size_t m_path = 0;                               // the one open (or that could not be opened)
        size_t m_files_begun = 0;                        // files opened so far (TakeLines ends a batch at a new one)
        char m_last_byte = '\n';                         // the last byte the reader took of the open file
        char m_newline[1] = { '\n' };                    // the get area of the '\n' a file's open last line gets
        OnOpen m_on_open;
        bool m_bgzf = false;                             // read with FillBgzf, else FillStream
        int m_fd = -1;
        std::string m_open_error;
        std::vector<unsigned char> m_block_in;           // one compressed BGZF block, or FillStream's input
        uint64_t m_offset = 0;                           // of the next unread byte of the file
        size_t m_in_pos = 0, m_in_end = 0;               // FillStream: unread input in m_block_in
        bool m_in_eof = false;                           // FillStream: the file has no more bytes
        size_t m_members = 0;                            // gzip members begun (BGZF blocks read; zstd frames ended)
        Mode m_mode = Mode::Look;
        std::unique_ptr<inflate_state> m_inflate;        // a gzip member's inflate state (ISA-L, ~85 KB)
        isal_gzip_header m_gz_header {};                 // its header's, while it is read
        bool m_header_read = false;                      // the member's header has been read
        ZSTD_DCtx* m_zstd = nullptr;                     // a zstd file's decompression
        bool m_zstd_in_frame = false;                    // a zstd frame has begun and not ended
        bool m_last_empty = false;                       // the last block read was empty (the EOF marker)
        static std::atomic<unsigned>& InflateThreadsSetting() {
            static std::atomic<unsigned> threads{ 1 };
            return threads;
        }
        static std::atomic<unsigned>& GzipThreadsSetting() {
            static std::atomic<unsigned> threads{ 1 };
            return threads;
        }
        static std::atomic<size_t>& MaxGzipMemberSetting() {
            static std::atomic<size_t> bytes{ kMaxGzipMember };
            return bytes;
        }
        // SDM_INFLATE_THREADS's value (at most kMaxInflateThreads), 0 if it is not set.
        static unsigned InflateThreadsOfTheEnvironment() {
            if (char const* env = std::getenv("SDM_INFLATE_THREADS")) {
                long const n = std::strtol(env, nullptr, 10);
                if (n >= 1) return static_cast<unsigned>(std::min<long>(n, kMaxInflateThreads));
            }
            return 0;
        }
        // The threads for the files this stream opens: set_threads' unless 0, else the settings' (the environment's first).
        unsigned StreamInflateThreads() const {
            if (unsigned const n = InflateThreadsOfTheEnvironment()) return n;
            return m_inflate_threads > 0 ? m_inflate_threads : InflateThreads();
        }
        unsigned StreamGzipThreads() const {
            if (unsigned const n = InflateThreadsOfTheEnvironment()) return n;
            return m_gzip_threads > 0 ? m_gzip_threads : GzipThreads();
        }
        unsigned m_inflate_threads = 0, m_gzip_threads = 0;  // set_threads' (0: the settings')

        std::thread m_thread;                            // Inflate, or with several inflating threads Dispatch
        unsigned m_workers = 1;                          // threads inflating BGZF blocks (1: m_thread does)
        std::unique_ptr<ParallelGzip> m_parallel;        // a gzip file inflated in pieces on several threads (Inflate copies)
        unsigned m_parallel_threads = 0;                 // its threads (0: none)
        std::string m_parallel_error;                    // why a gzip file is not inflated on several threads after all
        bool m_opened_bgzf = false;                      // the file began as BGZF
        bool m_zstd_frames = false;                      // a zstd file's frames read whole for the workers (ReadZstdJob)
        bool m_opened_zstd_frames = false;               // the file began so
        bool m_gzip_members = false;                     // a gzip file's members read whole for the workers (ReadGzipMembersJob)
        bool m_opened_gzip_members = false;              // the file began so
        size_t m_max_member = kMaxGzipMember;            // the open file's (MaxGzipMemberSetting when it was opened)
        size_t m_members_read = 0;                       // gzip members of the batches the reader has taken
        size_t m_member_restarts = 0;                    // batches of members after which the file was read on one thread
        bool m_switched_to_members = false;              // read in pieces, the file went on by members (SwitchToMembers)
        ParallelGzip::Counts m_parallel_counts;          // the pieces' counts once their reader is gone
        static constexpr uint64_t kNowhere = UINT64_MAX;
        uint64_t m_switched_at = kNowhere;               // the byte from which it went on by members
        unsigned m_switch_workers = 0;                   // on how many workers
        uint64_t m_restart_at = kNowhere;                // the byte from which it was read on one thread (ReadOnOneThreadFrom)
        char m_format = 0;                               // what its first bytes are: 'g' gzip, 'z' zstd, 'p' plain (0: not yet known)
        static constexpr uint64_t kReadAhead = uint64_t{64} << 20;      // a regular file's bytes asked to be read ahead (ReadUpTo)
        static constexpr uint64_t kReadAheadStep = uint64_t{16} << 20;  // at a time
        bool m_regular = false;                          // the open file is a regular one (no pipe), read ahead (ReadAheadWanted)
        uint64_t m_read_pos = 0;                         // where its next read() begins
        uint64_t m_read_ahead = 0;                       // the bytes before this asked to be read ahead
        unsigned m_opened_threads = 1;                   // the threads it was opened to inflate on (inflate_threads())
        std::vector<FileRead> m_files_read;              // the stream's files read so far (close)
        std::vector<std::thread> m_worker_threads;       // Work, with several
        std::vector<std::vector<char>> m_blocks;
        std::vector<Job> m_jobs;                         // per output block
        std::deque<size_t> m_free;                       // blocks the inflating thread (dispatcher) may fill
        std::deque<size_t> m_ready;                      // blocks in file order; the reader takes the first once complete
        std::deque<size_t> m_todo;                       // batches for the workers
        size_t m_current = kNone;                        // the block the get area points into
        bool m_stop = false;                             // close() asks the inflating threads to end
        bool m_done = false;                             // no more blocks will be queued
        bool m_failed = false;                           // the reader reached a corrupt block: nothing more is read
        std::string m_dispatch_error;                    // what stopped the dispatcher, the reader's once it gets there
        std::string m_error;
        mutable std::mutex m_mutex;
        std::condition_variable m_cv;
    };

    // The istream for ThreadedGzStreambuf, a drop-in for igzstream when reading.
    class ThreadedGzIstream : public std::istream {
    public:
        explicit ThreadedGzIstream(char const* path) : ThreadedGzIstream(std::vector<std::string>{ path }) {}

        // Several files read one after the other as one stream (ThreadedGzStreambuf::open); on_open, if given, is
        // called as each is opened.
        explicit ThreadedGzIstream(std::vector<std::string> paths, ThreadedGzStreambuf::OnOpen on_open = {}) : std::istream(nullptr) {
            std::istream::rdbuf(&m_buf);  // clears the badbit of the null buffer
            m_buf.set_on_open(std::move(on_open));
            if (!m_buf.open(std::move(paths))) setstate(std::ios::badbit);
        }

        ThreadedGzStreambuf* rdbuf() { return &m_buf; }

        void close() { m_buf.close(); }

    private:
        ThreadedGzStreambuf m_buf;
    };
}
