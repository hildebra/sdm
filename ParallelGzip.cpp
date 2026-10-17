// ParallelGzip.cpp - ParallelGzip: pieces of a gzip member decoded on several threads with ISA-L's own inflate
// (the stock library's API), each with the 32 KB before it labelled until it is known
// (protal docs/claude/2026-10-12-isal-parallel-inflate). Ported from protal (github.com/hildebra/protalFastHit,
// branch audit-fixes, 956e18b, src/IO/ParallelGzip.cpp).

#include "ParallelGzip.h"

#include <isa-l/crc.h>
#include <isa-l/igzip_lib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <new>
#include <thread>
#include <utility>
#include <vector>

#include "GzipMember.h"

// SDM_CLONE_V3 compiles a function twice, for plain x86-64 and for x86-64-v3 (AVX2, BMI1/2, ...), and the
// loader runs the copy the CPU supports (GCC's target_clones, through an ifunc; static glibc binaries resolve
// it too), so one sdm binary runs on any x86-64 CPU. Needs GCC 12 or later on x86-64 Linux with glibc;
// elsewhere, with -DSDM_NO_CLONES, or under ThreadSanitizer (the ifunc resolver runs before its runtime
// is set up, and the program crashes at load), the function is compiled once, for the build's target
// (protal's src/Utilities/TargetClones.h).
#if defined(__x86_64__) && defined(__linux__) && defined(__GLIBC__) && defined(__GNUC__) && !defined(__clang__) && \
    __GNUC__ >= 12 && !defined(SDM_NO_CLONES) && !defined(__SANITIZE_THREAD__)
#define SDM_CLONE_V3 __attribute__((target_clones("default", "arch=x86-64-v3")))
#else
#define SDM_CLONE_V3
#endif

namespace sdm {

    namespace {

        constexpr uint32_t kWindow = 32768;                           // deflate's window
        constexpr uint64_t kNone = UINT64_MAX;                        // no bit position (no end: the member's final block)
        constexpr uint64_t kUnknown = UINT64_MAX - 1;                 // a slot's block start not looked for yet
        constexpr size_t kInputChunk = size_t{1} << 30;               // ISA-L's counts are 32-bit
        constexpr size_t kMaxPieceOutput = size_t{256} << 20;         // a piece decoded whole holds at most this
        constexpr size_t kStreamSegment = size_t{4} << 20;            // a range decoded to its end gives this at a time
        constexpr size_t kReleaseStep = size_t{32} << 20;             // mapped bytes behind the reader released at a time
        constexpr size_t kReleaseMargin = size_t{1} << 20;            // and kept behind where it decodes next
        constexpr size_t kReadAhead = size_t{64} << 20;               // asked to be read ahead of the workers' furthest piece

        // Bytes that grow without being zeroed (a piece's output, decoded into place).
        struct Bytes {
            uint8_t* data = nullptr;
            size_t size = 0;
            size_t capacity = 0;
            Bytes() = default;
            Bytes(Bytes const&) = delete;
            Bytes& operator=(Bytes const&) = delete;
            ~Bytes() { std::free(data); }
            bool Reserve(size_t n) {
                if (n <= capacity) return true;
                size_t const c = std::max(n, capacity + capacity / 2);
                void* p = std::realloc(data, c);
                if (!p) return false;
                data = static_cast<uint8_t*>(p);
                capacity = c;
                return true;
            }
        };

        // The labels a piece's unknown window is decoded with, twice: a byte copied from window position p comes out
        // as A(p) in one pass and B(p) in the other, which differ for every p and give p back; where the passes agree,
        // the byte is the piece's own. A(p) = p & 0xff, B(p) = A(p) + 1 + (p >> 8) (mod 256), so
        // p = A | ((B - A - 1) & 0xff) << 8.
        struct Labels {
            std::array<uint8_t, kWindow> a{}, b{};
            Labels() {
                for (uint32_t p = 0; p < kWindow; p++) {
                    a[p] = static_cast<uint8_t>(p & 0xff);
                    b[p] = static_cast<uint8_t>((a[p] + 1 + (p >> 8)) & 0xff);
                }
            }
        };
        Labels const& TheLabels() {
            static Labels const labels;
            return labels;
        }


        uint32_t LoadLe32(uint8_t const* p) {
            return uint32_t{p[0]} | uint32_t{p[1]} << 8 | uint32_t{p[2]} << 16 | uint32_t{p[3]} << 24;
        }

        // ---- finding deflate blocks ---------------------------------------------------------------------------------
        // A dynamic Huffman block's header is checked as a decoder checks it: not the final block, at most 286
        // literal/length and 30 distance codes, a complete code-length code, the code lengths it gives fitting their
        // count (a repeat not first), an end-of-block code, and complete literal/length and distance codes (or one of
        // length 1, or no distance code at all, as zlib allows). Random bits pass this about never; a start found is
        // still checked by decoding from it (a piece's own end, and the piece before ending there).

        // The 64 bits of the file from byte `byte` on (zeros past its end).
        uint64_t Load64(uint8_t const* data, size_t size, size_t byte) {
            if (byte + 8 <= size) {
                uint64_t w;
                std::memcpy(&w, data + byte, 8);
                return w;
            }
            uint64_t w = 0;
            for (size_t i = 0; byte + i < size && i < 8; i++) w |= uint64_t{data[byte + i]} << (8 * i);
            return w;
        }

        class BitReader {
        public:
            BitReader(uint8_t const* data, size_t size, uint64_t bit) : m_data(data), m_size(size), m_bit(bit) {}
            uint32_t Peek(unsigned n) const {
                uint64_t const w = Load64(m_data, m_size, static_cast<size_t>(m_bit >> 3)) >> (m_bit & 7);
                return static_cast<uint32_t>(w & ((uint64_t{1} << n) - 1));  // n <= 32 < 57
            }
            uint32_t Get(unsigned n) {
                uint32_t const v = Peek(n);
                m_bit += n;
                return v;
            }
            void Skip(unsigned n) { m_bit += n; }
            bool Past() const { return (m_bit >> 3) > m_size; }

        private:
            uint8_t const* m_data;
            size_t m_size;
            uint64_t m_bit;
        };

        // Kraft's sum of code lengths (max_bits at most 15): 1 << max_bits for a complete code, more if
        // over-subscribed; the number of codes and the longest in `count`, `longest`.
        uint32_t Kraft(uint8_t const* lengths, size_t n, unsigned max_bits, size_t& count, unsigned& longest) {
            uint32_t sum = 0;
            count = 0;
            longest = 0;
            for (size_t i = 0; i < n; i++) {
                if (lengths[i] == 0) continue;
                sum += uint32_t{1} << (max_bits - lengths[i]);
                count++;
                longest = std::max<unsigned>(longest, lengths[i]);
            }
            return sum;
        }

        // Whether a code is one a decoder takes: complete, or a single code of length 1 (or none, if may_be_empty).
        bool CodeOk(uint8_t const* lengths, size_t n, unsigned max_bits, bool may_be_empty) {
            size_t count;
            unsigned longest;
            uint32_t const sum = Kraft(lengths, n, max_bits, count, longest);
            if (count == 0) return may_be_empty;
            if (sum == uint32_t{1} << max_bits) return true;
            return sum < (uint32_t{1} << max_bits) && count == 1 && longest == 1;
        }

        // Kraft's sum (7-bit code lengths) of four 3-bit lengths.
        constexpr std::array<uint16_t, 4096> kKraft4 = [] {
            std::array<uint16_t, 4096> sums{};
            for (uint32_t x = 0; x < 4096; x++) {
                for (uint32_t j = 0; j < 4; j++) {
                    uint32_t const l = (x >> (3 * j)) & 7;
                    if (l) sums[x] = static_cast<uint16_t>(sums[x] + (1u << (7 - l)));
                }
            }
            return sums;
        }();

        bool IsDynamicBlockHeader(uint8_t const* data, size_t size, uint64_t bit) {
            static constexpr uint8_t kOrder[19] = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };
            // The header's first 74 bits at most (17 then 3 per code length) from two loads: most positions are
            // refused within them, before the code lengths are decoded.
            size_t const byte = static_cast<size_t>(bit >> 3);
            unsigned const shift = bit & 7;
            uint64_t const w = Load64(data, size, byte) >> shift;  // 57 bits
            if ((w & 7) != 4) return false;                        // BFINAL 0, BTYPE 2 (dynamic)
            uint32_t const hlit = (w >> 3) & 31, hdist = (w >> 8) & 31;
            if (hlit > 29 || hdist > 29) return false;
            uint32_t const hclen = static_cast<uint32_t>((w >> 13) & 15) + 4;
            // The code-length code's lengths (3 bits each, hclen of them), and its Kraft sum by four at a time
            // (kKraft4), which must be 128 (complete): without a branch per length, as nearly all positions fail here.
            uint64_t const cl_bits = ((w >> 17) & ((uint64_t{1} << 39) - 1)) |                  // the first 13
                                     (((Load64(data, size, byte + 7) >> shift) & ((uint64_t{1} << 18) - 1)) << 39);  // the rest
            uint64_t const used = cl_bits & ((uint64_t{1} << (3 * hclen)) - 1);
            uint32_t const sum = kKraft4[used & 4095] + kKraft4[(used >> 12) & 4095] + kKraft4[(used >> 24) & 4095] +
                                 kKraft4[(used >> 36) & 4095] + kKraft4[(used >> 48) & 4095];
            if (sum != 128) return false;  // the code-length code must be complete
            uint8_t code_lengths[19] = {};
            for (uint32_t i = 0; i < hclen; i++) code_lengths[kOrder[i]] = static_cast<uint8_t>((used >> (3 * i)) & 7);
            BitReader r(data, size, bit + 17 + 3 * hclen);
            // Its decoding table, by the next 7 bits as they are read (codes are reversed in the stream).
            uint8_t table_symbol[128], table_length[128];
            {
                uint32_t count[8] = {}, next[8] = {};
                for (uint8_t l : code_lengths) count[l]++;
                count[0] = 0;
                uint32_t code = 0;
                for (unsigned l = 1; l <= 7; l++) {
                    code = (code + count[l - 1]) << 1;
                    next[l] = code;
                }
                for (unsigned s = 0; s < 19; s++) {
                    unsigned const l = code_lengths[s];
                    if (l == 0) continue;
                    uint32_t c = next[l]++, reversed = 0;
                    for (unsigned i = 0; i < l; i++) {
                        reversed = (reversed << 1) | (c & 1);
                        c >>= 1;
                    }
                    for (uint32_t k = reversed; k < 128; k += uint32_t{1} << l) {
                        table_symbol[k] = static_cast<uint8_t>(s);
                        table_length[k] = static_cast<uint8_t>(l);
                    }
                }
            }
            size_t const total = hlit + 257 + hdist + 1;
            uint8_t lengths[286 + 30];
            size_t n = 0;
            while (n < total) {
                uint32_t const k = r.Peek(7);
                unsigned const symbol = table_symbol[k];
                r.Skip(table_length[k]);
                if (symbol < 16) {
                    lengths[n++] = static_cast<uint8_t>(symbol);
                    continue;
                }
                size_t repeat;
                uint8_t value = 0;
                if (symbol == 16) {
                    if (n == 0) return false;
                    value = lengths[n - 1];
                    repeat = 3 + r.Get(2);
                } else if (symbol == 17) {
                    repeat = 3 + r.Get(3);
                } else {
                    repeat = 11 + r.Get(7);
                }
                if (n + repeat > total) return false;
                std::memset(lengths + n, value, repeat);
                n += repeat;
            }
            if (r.Past()) return false;
            if (lengths[256] == 0) return false;  // no end-of-block code
            return CodeOk(lengths, hlit + 257, 15, false) && CodeOk(lengths + hlit + 257, hdist + 1, 15, true);
        }

        // The first dynamic block header at a bit in [from, to) of data, or kNone. Most positions are refused by their
        // first three bits (BFINAL 0, BTYPE 2: bits 0, 0, 1), which are tested for 56 positions at once.
        uint64_t FindBlock(uint8_t const* data, size_t size, uint64_t from, uint64_t to) {
            to = std::min<uint64_t>(to, uint64_t{size} * 8);
            for (uint64_t at = from; at < to; at += 56) {
                uint64_t const w = Load64(data, size, static_cast<size_t>(at >> 3)) >> (at & 7);  // >= 57 valid bits
                uint64_t candidates = ~w & ~(w >> 1) & (w >> 2) & ((uint64_t{1} << 56) - 1);
                while (candidates) {
                    uint64_t const bit = at + static_cast<unsigned>(__builtin_ctzll(candidates));
                    candidates &= candidates - 1;
                    if (bit >= to) return kNone;
                    if (IsDynamicBlockHeader(data, size, bit)) return bit;
                }
            }
            return kNone;
        }

        // ---- decoding a range of a member ---------------------------------------------------------------------------

        enum class RangeEnd : uint8_t {
            Running,     // stopped at the output limit (Run again)
            AtBoundary,  // ended at the range's end, a block boundary
            Final,       // the member's final block ended (final_bit) before the range's end
            Failed,      // a decoding error, or the range's end is no block boundary (or the input ends: error says)
            TooBig,      // more output than the limit allows (a piece decoded whole)
        };

        // Decodes raw deflate from a block start, with the 32 KB before it (the true one, or labels) as ISA-L's
        // dictionary, up to a later block start (`end`), or through the member's final block (end kNone). The input
        // given ends at `end`: its last byte with the bits from `end` on zeroed, so that the decoder cannot read past
        // it (zero bits begin a stored block, whose header then waits for more input). The range ended at `end` if
        // ISA-L then waits for a block header that begins at `end`: when it runs out of input in a header, it keeps
        // the header's bytes and the bits before them (tmp_in_size, read_in_length) as they were where the header
        // began. Uses the stock ISA-L's API only.
        class RangeDecoder {
        public:
            RangeDecoder() : m_state(std::make_unique<inflate_state>()) {}

            void Start(uint8_t const* data, size_t size, uint64_t start, uint64_t end, uint8_t const* dict,
                       uint32_t dict_len) {
                m_data = data;
                m_size = size;
                m_end = end;
                m_error.clear();
                m_final_bit = kNone;
                inflate_state& st = *m_state;
                isal_inflate_init(&st);
                st.crc_flag = ISAL_DEFLATE;
                m_bad_start = false;
                if (dict_len > 0 && isal_inflate_set_dict(&st, const_cast<uint8_t*>(dict), dict_len) != COMP_OK) {
                    m_bad_start = true;
                }
                size_t byte = static_cast<size_t>(start >> 3);
                unsigned const shift = start & 7;
                if (shift) {
                    st.read_in = uint64_t{data[byte]} >> shift;
                    st.read_in_length = static_cast<int32_t>(8 - shift);
                    byte++;
                }
                m_next = byte;
                m_main_end = end == kNone ? size : static_cast<size_t>(end >> 3);
                if (end != kNone && (end <= start || m_main_end < byte || m_main_end > size)) m_bad_start = true;
                m_tail_bits = end == kNone ? 0 : static_cast<unsigned>(end & 7);
                m_tail_given = false;
                m_in_tail = false;
                if (m_tail_bits && !m_bad_start) {
                    if (m_main_end >= size) m_bad_start = true;
                    else m_tail = static_cast<uint8_t>(data[m_main_end] & ((1u << m_tail_bits) - 1));
                }
                st.next_in = const_cast<uint8_t*>(data + m_next);
                st.avail_in = 0;
            }

            // Decodes into out (appended) until out holds `limit` bytes, or the range ends. With `whole`, reaching
            // the limit is TooBig (a piece that must be decoded whole); else Running.
            RangeEnd Run(Bytes& out, size_t limit, bool whole) {
                if (m_bad_start) {
                    m_error = "a range that cannot be decoded";
                    return RangeEnd::Failed;
                }
                inflate_state& st = *m_state;
                for (;;) {
                    if (st.avail_in == 0) Feed();
                    if (out.size >= limit) return whole ? RangeEnd::TooBig : RangeEnd::Running;
                    size_t const room = std::min({ limit - out.size, std::max(out.size, size_t{4} << 20), kInputChunk });
                    if (!out.Reserve(out.size + room)) {
                        m_error = "out of memory";
                        return RangeEnd::Failed;
                    }
                    st.next_out = out.data + out.size;
                    st.avail_out = static_cast<uint32_t>(room);
                    uint32_t const avail_before = st.avail_in;
                    int32_t const bits_before = st.read_in_length;
                    auto const state_before = st.block_state;
                    int const ret = isal_inflate(&st);
                    size_t const produced = static_cast<size_t>(st.next_out - (out.data + out.size));
                    out.size += produced;
                    if (ret < 0) {
                        m_error = ret == ISAL_INVALID_LOOKBACK ? "invalid distance too far back"
                                  : ret == ISAL_INVALID_SYMBOL ? "invalid code"
                                  : ret == ISAL_INVALID_BLOCK  ? "invalid block type"
                                                               : "corrupt deflate data";
                        return RangeEnd::Failed;
                    }
                    if (st.block_state == ISAL_BLOCK_FINISH) {
                        m_final_bit = InputBit() * 8 - static_cast<uint64_t>(st.read_in_length);
                        if (m_end != kNone && m_final_bit > m_end) {
                            m_error = "a final block past the range's end";
                            return RangeEnd::Failed;
                        }
                        return RangeEnd::Final;
                    }
                    bool const progressed = produced > 0 || st.avail_in != avail_before || st.read_in_length != bits_before ||
                                            st.block_state != state_before;
                    if (progressed) continue;
                    if (st.avail_in == 0 && MoreInput()) continue;
                    break;  // all input given, nothing more decoded
                }
                if (m_end == kNone) {
                    m_error = "the file ends inside a gzip member (truncated file?)";
                    return RangeEnd::Failed;
                }
                if ((st.block_state == ISAL_BLOCK_HDR || st.block_state == ISAL_BLOCK_NEW_HDR) &&
                    st.tmp_out_valid == st.tmp_out_processed &&
                    InputBit() * 8 - uint64_t(st.tmp_in_size) * 8 - static_cast<uint64_t>(st.read_in_length) == m_end) {
                    return RangeEnd::AtBoundary;
                }
                m_error = "no block boundary at the range's end";
                return RangeEnd::Failed;
            }

            uint64_t final_bit() const { return m_final_bit; }
            std::string const& error() const { return m_error; }
            // The byte of the file the decoder reads next.
            uint64_t input_byte() const { return InputBit(); }

        private:
            bool MoreInput() const { return m_next < m_main_end || (m_tail_bits && !m_tail_given); }

            void Feed() {
                inflate_state& st = *m_state;
                if (m_next < m_main_end) {
                    size_t const n = std::min(m_main_end - m_next, kInputChunk);
                    st.next_in = const_cast<uint8_t*>(m_data + m_next);
                    st.avail_in = static_cast<uint32_t>(n);
                    m_next += n;
                    m_in_tail = false;
                } else if (m_tail_bits && !m_tail_given) {
                    st.next_in = &m_tail;
                    st.avail_in = 1;
                    m_tail_given = true;
                    m_in_tail = true;
                }
            }

            // The byte of the file next_in stands for (the masked last byte stands for its own).
            uint64_t InputBit() const {
                inflate_state const& st = *m_state;
                if (m_in_tail) return uint64_t{m_main_end} + static_cast<uint64_t>(st.next_in - &m_tail);
                return static_cast<uint64_t>(st.next_in - m_data);
            }

            std::unique_ptr<inflate_state> m_state;
            uint8_t const* m_data = nullptr;
            size_t m_size = 0;
            uint64_t m_end = kNone;
            size_t m_next = 0, m_main_end = 0;
            unsigned m_tail_bits = 0;
            uint8_t m_tail = 0;
            bool m_tail_given = false, m_in_tail = false, m_bad_start = false;
            uint64_t m_final_bit = kNone;
            std::string m_error;
        };

    }  // namespace

    namespace parallel_gzip {

        // A speculative piece's output (pass A, resolved in place; pass B) with the true window before it. Indices
        // into window ++ 0..255: a marker's window position, any other byte its own value + 32768. Cloned for
        // x86-64-v3 (SDM_CLONE_V3), and so not in the anonymous namespace.
        SDM_CLONE_V3 void Resolve(uint8_t* a, uint8_t const* b, size_t n, uint8_t const* window) {
            using U8x16 = uint8_t __attribute__((vector_size(16)));
            using U16x16 = uint16_t __attribute__((vector_size(32)));
            alignas(64) uint8_t table[kWindow + 256];
            std::memcpy(table, window, kWindow);
            for (unsigned v = 0; v < 256; v++) table[kWindow + v] = static_cast<uint8_t>(v);
            size_t i = 0;
            for (; i + 16 <= n; i += 16) {
                U8x16 va, vb;
                std::memcpy(&va, a + i, 16);
                std::memcpy(&vb, b + i, 16);
                auto const same8 = va == vb;
                uint64_t lo, hi;
                std::memcpy(&lo, &same8, 8);
                std::memcpy(&hi, reinterpret_cast<char const*>(&same8) + 8, 8);
                if ((lo & hi) == ~uint64_t{0}) continue;  // no marker among the 16
                U16x16 const wa = __builtin_convertvector(va, U16x16), wb = __builtin_convertvector(vb, U16x16);
                U16x16 const position = (wa | (((wb - wa - 1) & 0xff) << 8)) & 0x7fff;
                U16x16 const same = __builtin_convertvector(same8, U16x16) & 1;  // 1 where no marker
                U16x16 const index = same ? wa + static_cast<uint16_t>(kWindow) : position;
                alignas(32) uint16_t idx[16];
                std::memcpy(idx, &index, sizeof(idx));
                for (int j = 0; j < 16; j++) a[i + j] = table[idx[j]];
            }
            for (; i < n; i++) {
                if (a[i] != b[i]) a[i] = window[(a[i] | (((b[i] - a[i] - 1) & 0xffu) << 8)) & 0x7fffu];
            }
        }

    }  // namespace parallel_gzip

    namespace {

        // The last 32 KB of output after `n` more bytes (data) on top of window.
        void PushWindow(uint8_t* window, uint8_t const* data, size_t n) {
            if (n >= kWindow) {
                std::memcpy(window, data + n - kWindow, kWindow);
            } else {
                std::memmove(window, window + n, kWindow - n);
                std::memcpy(window + kWindow - n, data, n);
            }
        }

    }  // namespace

    // ---- the reader ---------------------------------------------------------------------------------------------------
    //
    // The file is cut into slots of kPieceBytes of compressed data; a slot's piece begins at the first dynamic block
    // header in it (FindBlock; none: the slot belongs to the piece before), the first slot's at the first member's
    // deflate data, and ends where the next piece begins (or at the member's final block). Workers decode pieces
    // ahead of the reader with labelled windows (two passes) and resolve them once their window is known; the
    // reader's own thread (the lead) decodes in order, with the true window, every piece no worker has taken, chains
    // the windows (from a speculative piece's last 32 KB alone), checks each member's CRC-32 and length, and reads
    // the members' trailers and headers. A speculative piece is used only if the piece before ends exactly where it
    // begins and it ends where it should (its own end, or a final block); otherwise the lead decodes it, and a range
    // whose end is no block boundary, or too large to hold, is decoded on to the member's end and given as it comes.
    struct ParallelGzip::Impl {
        enum class State : uint8_t { Free, Decoding, Decoded, Taken };

        struct Piece {
            State state = State::Free;
            uint64_t start = kNone, end = kNone;
            RangeEnd result = RangeEnd::Failed;
            uint64_t final_bit = kNone;
            std::unique_ptr<Bytes> a, b;
        };

        // A run of output in file order: a range the lead decoded (bytes), a speculative piece being resolved (slot),
        // or the lead's decoding to a member's end (stream), given as it comes; after it, maybe the member's end
        // (its trailer), the file's end, or an error.
        struct Item {
            std::unique_ptr<Bytes> bytes;
            size_t slot = SIZE_MAX;
            bool resolved = true;           // a speculative piece: its bytes resolved (under the mutex)
            std::unique_ptr<uint8_t[]> window;  // the true window before it, for its resolution
            bool resolving = false;         // taken by a thread to resolve
            bool streamed = false;          // the lead's stream (stream: not at its end yet)
            bool stream = false;
            bool member_end = false;
            uint32_t trailer_crc = 0, trailer_isize = 0;
            bool eof = false;
            std::string error;
        };

        uint8_t const* data = nullptr;
        size_t size = 0;
        size_t max_member = 0;         // small members after a member stop the reader (0: never)
        uint64_t stopped_at = ParallelGzip::kNotStopped;  // the member the reader stopped before (the lead's thread)
        size_t members_read = 0;       // and the members before it
        size_t piece_bytes = ParallelGzip::kChunkBytes;
        size_t slots = 0;
        unsigned threads = 1;
        size_t gap = 0, ahead = 0;

        std::mutex mutex;
        std::condition_variable changed;
        std::vector<uint64_t> starts;  // per slot: its piece's start (kUnknown: not looked for; kNone: no piece)
        std::vector<Piece> pieces;
        std::deque<Item> items;        // made by the chain, given by Read, in file order
        size_t chain_slot = 0;         // the slot the chain is in
        size_t emit_slot = 0;          // the slot of the item being given
        size_t next_find = 1;          // the next slot whose start the workers look for
        size_t next_spec = 1;          // the next slot to speculate on
        bool stop = false;
        std::vector<std::thread> workers;

        // The chain (the lead's thread only).
        uint64_t pos = 0;              // the next bit to decode: a block start, or a member's deflate start
        std::array<uint8_t, kWindow> window{};  // the last 32 KB of output
        uint64_t member_bytes = 0;     // of the current member so far (the true window's length, up to 32 KB)
        size_t member = 1;
        bool chain_done = false;
        RangeDecoder lead;             // the lead's decoder
        bool streaming = false;        // lead decodes to the member's end (the last item is its stream)
        size_t dropped = 0;            // false block starts dropped in a row
        ParallelGzip::Counts counts;   // how the ranges were decoded (the lead's thread)

        // Giving the output (Read).
        Bytes stream_bytes;            // a stream item's current segment
        size_t given = 0;              // of the front item's (or segment's) bytes
        uint32_t crc = 0;
        uint64_t member_out = 0;
        std::string error;
        size_t released = 0;           // the mapped file's bytes before this are released (Release)
        size_t read_ahead = 0;         // and those before this asked to be read ahead (Release)

        // The mapped file's pages before what the lead decodes next (less kReleaseMargin), released kReleaseStep at a time
        // (MADV_DONTNEED): every page read would otherwise stay in the resident set until the file is closed (a
        // paired-end run at r226 held 7.3 GB more, protal docs/claude/2026-10-13-gzip-members-regression). Nothing needs them
        // again; a worker that still looks at one (a block start it scans late) reads it from the page cache again, as
        // the mapping is private and never written. And the pages up to kReadAhead past the workers' furthest piece asked
        // for (MADV_WILLNEED), kReleaseStep at a time: a file not in the page cache (a sample's first read, from a network
        // file system) is then read ahead of the faults rather than at them. The lead's thread only.
        void Release() {
            size_t const page = static_cast<size_t>(::sysconf(_SC_PAGESIZE));
            uint64_t const next = streaming ? lead.input_byte() : pos >> 3;
            size_t const wanted = static_cast<size_t>(std::min<uint64_t>(size, next + (ahead + 2) * piece_bytes + kReadAhead));
            if (ReadAheadWanted() && (wanted >= read_ahead + kReleaseStep || (wanted == size && read_ahead < size))) {
                size_t const from = read_ahead / page * page;
                ::madvise(const_cast<uint8_t*>(data) + from, wanted - from, MADV_WILLNEED);
                read_ahead = wanted;
            }
            if (next < released + kReleaseStep + kReleaseMargin) return;
            size_t const upto = static_cast<size_t>(next - kReleaseMargin) / page * page;
            if (upto <= released) return;
            ::madvise(const_cast<uint8_t*>(data) + released, upto - released, MADV_DONTNEED);
            released = upto;
        }

        ~Impl() {
            {
                std::lock_guard<std::mutex> lock(mutex);
                stop = true;
            }
            changed.notify_all();
            for (auto& t : workers) t.join();
            if (data) ::munmap(const_cast<uint8_t*>(data), size);
        }

        // Output buffers, kept for the next pieces once given (a new one's pages would be faulted in again).
        std::mutex pool_mutex;
        std::vector<std::unique_ptr<Bytes>> pool;

        std::unique_ptr<Bytes> TakeBuffer() {
            {
                std::lock_guard<std::mutex> lock(pool_mutex);
                if (!pool.empty()) {
                    auto bytes = std::move(pool.back());
                    pool.pop_back();
                    bytes->size = 0;
                    return bytes;
                }
            }
            return std::make_unique<Bytes>();
        }
        void GiveBack(std::unique_ptr<Bytes>& bytes) {
            if (!bytes) return;
            std::lock_guard<std::mutex> lock(pool_mutex);
            if (pool.size() < 4 * ahead + 8 && bytes->capacity <= 4 * kMaxPieceOutput / 64) pool.push_back(std::move(bytes));
            bytes.reset();
        }

        size_t SlotOf(uint64_t bit) const { return static_cast<size_t>((bit >> 3) / piece_bytes); }

        // Looks for slot s's piece start (outside the mutex) and records it.
        uint64_t StartOf(size_t s) {
            {
                std::lock_guard<std::mutex> lock(mutex);
                if (starts[s] != kUnknown) return starts[s];
            }
            uint64_t const from = uint64_t{s} * piece_bytes * 8;
            uint64_t const found = FindBlock(data, size, std::max(from, starts[0] + 1), from + piece_bytes * 8);
            std::lock_guard<std::mutex> lock(mutex);
            if (starts[s] == kUnknown) starts[s] = found;
            return starts[s];
        }

        // The first piece start after `bit` within the look-ahead (kNone: none, the range goes to the member's end): a
        // file without dynamic blocks (zlib-ng's level 1 writes fixed codes only, which have no header to find) is not
        // searched to its end before its first bytes are given.
        uint64_t NextStartAfter(uint64_t bit) {
            size_t const first = SlotOf(bit);
            for (size_t s = first; s < slots && s <= first + ahead; s++) {
                uint64_t const start = StartOf(s);
                if (start != kNone && start > bit) return start;
            }
            return kNone;
        }

        void Work() {
            RangeDecoder decoder;
            Labels const& labels = TheLabels();
            std::unique_lock<std::mutex> lock(mutex);
            while (!stop) {
                // 1. A piece the reader waits for, or soon will: its resolution.
                Item* resolve = nullptr;
                for (auto& item : items) {
                    if (!item.resolved && !item.resolving) {
                        resolve = &item;
                        break;
                    }
                }
                if (resolve) {
                    resolve->resolving = true;
                    Piece& piece = pieces[resolve->slot];
                    lock.unlock();
                    parallel_gzip::Resolve(piece.a->data, piece.b->data, piece.a->size, resolve->window.get());
                    lock.lock();
                    GiveBack(piece.b);
                    resolve->resolved = true;
                    resolve->window.reset();
                    changed.notify_all();
                    continue;
                }
                // 2. Piece starts ahead of the chain (cheap; the lead needs them for its ends).
                size_t const limit = std::min(slots, emit_slot + ahead + 2);
                if (next_find < limit) {
                    size_t const s = next_find++;
                    lock.unlock();
                    StartOf(s);
                    lock.lock();
                    continue;
                }
                // 3. A piece far enough ahead that the lead will not reach it first.
                size_t const s = std::max(next_spec, chain_slot + gap);
                if (s < slots && s < emit_slot + ahead) {
                    next_spec = s + 1;
                    if (pieces[s].state != State::Free || starts[s] == kNone) continue;
                    Piece& piece = pieces[s];
                    piece.state = State::Decoding;
                    lock.unlock();
                    Speculate(piece, s, decoder, labels);
                    lock.lock();
                    piece.state = State::Decoded;
                    changed.notify_all();
                    continue;
                }
                changed.wait(lock);
            }
        }

        // Decodes piece s from its start to the next piece's (or the member's end) twice, with labelled windows.
        void Speculate(Piece& piece, size_t s, RangeDecoder& decoder, Labels const& labels) {
            piece.start = StartOf(s);
            if (piece.start == kNone) {
                piece.result = RangeEnd::Failed;
                return;
            }
            piece.end = NextStartAfter(piece.start);
            piece.result = RangeEnd::Failed;
            try {
                piece.a = TakeBuffer();
                piece.b = TakeBuffer();
                decoder.Start(data, size, piece.start, piece.end, labels.a.data(), kWindow);
                RangeEnd const ra = decoder.Run(*piece.a, kMaxPieceOutput / 4, true);
                uint64_t const final_a = decoder.final_bit();
                if (ra != RangeEnd::AtBoundary && ra != RangeEnd::Final) {
                    GiveBack(piece.a);
                    GiveBack(piece.b);
                    piece.result = ra;
                    return;
                }
                piece.b->Reserve(piece.a->size);
                decoder.Start(data, size, piece.start, piece.end, labels.b.data(), kWindow);
                RangeEnd const rb = decoder.Run(*piece.b, kMaxPieceOutput / 4, true);
                if (rb != ra || piece.b->size != piece.a->size || decoder.final_bit() != final_a) {
                    GiveBack(piece.a);
                    GiveBack(piece.b);
                    return;
                }
                piece.result = ra;
                piece.final_bit = final_a;
            } catch (std::bad_alloc const&) {
                GiveBack(piece.a);
                GiveBack(piece.b);
                piece.result = RangeEnd::Failed;
            }
        }

        // ---- the chain (the lead) -------------------------------------------------------------------------------------

        void Push(Item item) {
            std::lock_guard<std::mutex> lock(mutex);
            items.push_back(std::move(item));
            changed.notify_all();
        }

        // After the member's final block ended at `final_bit`: its trailer into `item`, and what follows: another
        // member's deflate data (pos), or the file's end (or zero bytes of padding), else an error.
        void MemberEnd(uint64_t final_bit, Item& item) {
            size_t const trailer = static_cast<size_t>((final_bit + 7) >> 3);
            if (trailer + 8 > size) {
                item.error = "the file ends inside gzip member " + std::to_string(member) + " (truncated file?)";
                chain_done = true;
                return;
            }
            item.member_end = true;
            item.trailer_crc = LoadLe32(data + trailer);
            item.trailer_isize = LoadLe32(data + trailer + 4);
            size_t const next = trailer + 8;
            if (next == size || std::all_of(data + next, data + size, [](uint8_t c) { return c == 0; })) {
                item.eof = true;
                chain_done = true;
                return;
            }
            uint64_t start = 0;
            if (!ReadHeader(next, start)) {
                item.error = "the data at byte " + std::to_string(next) + ", after gzip member " + std::to_string(member) +
                             ", is no gzip member";
                chain_done = true;
                return;
            }
            // Small members from here on, this one and the next each within max_member bytes: their blocks are mostly
            // final ones, which pieces cannot start at (FindBlock would scan whole slots in vain, the lead streaming a
            // member at a time). Read stops before this member, and the caller reads the rest by members.
            if (max_member > 0) {
                size_t const first = gzip_member::Next(data, size, next, max_member);
                if (first > 0 && gzip_member::Next(data, size, next + first, max_member) > 0) {
                    stopped_at = next;
                    members_read = member;
                    chain_done = true;
                    return;
                }
            }
            pos = start;
            member++;
            member_bytes = 0;
        }

        // A gzip header at byte `at`: its deflate data's first bit into start.
        bool ReadHeader(size_t at, uint64_t& start) const {
            if (at + 2 > size || data[at] != 0x1f || data[at + 1] != 0x8b) return false;
            auto st = std::make_unique<inflate_state>();
            isal_inflate_init(st.get());
            st->crc_flag = ISAL_GZIP_NO_HDR_VER;
            isal_gzip_header header;
            isal_gzip_header_init(&header);  // no buffers: name, comment and extra field skipped
            st->next_in = const_cast<uint8_t*>(data + at);
            st->avail_in = static_cast<uint32_t>(std::min(size - at, kInputChunk));
            if (isal_read_gzip_header(st.get(), &header) != ISAL_DECOMP_OK) return false;
            start = static_cast<uint64_t>(st->next_in - data) * 8 - static_cast<uint64_t>(st->read_in_length);
            return true;
        }

        void Chained(uint8_t const* bytes, size_t n) {
            PushWindow(window.data(), bytes, n);
            member_bytes += n;
        }

        // One step of the chain: the range from pos, as an item. False when there is nothing more to chain (the file's
        // end, an error, or the lead's stream, which Read drives), or, without may_wait, when a worker is decoding the
        // piece at pos.
        bool Advance(bool may_wait = true) {
            if (chain_done || streaming) return false;
            uint64_t const end = NextStartAfter(pos);
            size_t const s = SlotOf(pos);
            Piece* piece = nullptr;
            {
                std::unique_lock<std::mutex> lock(mutex);
                chain_slot = s;
                changed.notify_all();
                if (s < slots && starts[s] == pos) {
                    piece = &pieces[s];
                    if (piece->state == State::Decoding && !may_wait) return false;
                    changed.wait(lock, [&] { return piece->state != State::Decoding; });
                    if (piece->state == State::Free) {
                        piece->state = State::Taken;
                        piece = nullptr;
                    }
                }
            }
            if (piece && piece->state == State::Decoded && piece->end == end &&
                (piece->result == RangeEnd::AtBoundary || piece->result == RangeEnd::Final)) {
                // The window after it: its last 32 KB resolved with the window before it.
                Item item;
                item.slot = s;
                item.resolved = false;
                item.window.reset(new uint8_t[kWindow]);
                std::memcpy(item.window.get(), window.data(), kWindow);
                size_t const n = piece->a->size, tail = std::min<size_t>(n, kWindow);
                std::unique_ptr<uint8_t[]> last(new uint8_t[tail]);
                std::memcpy(last.get(), piece->a->data + n - tail, tail);
                parallel_gzip::Resolve(last.get(), piece->b->data + n - tail, tail, window.data());
                Chained(last.get(), tail);
                member_bytes += n - tail;
                RangeEnd const result = piece->result;
                uint64_t const final_bit = piece->final_bit;
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    piece->state = State::Taken;
                }
                counts.ahead++;
                if (result == RangeEnd::Final) MemberEnd(final_bit, item);
                else pos = end;
                Push(std::move(item));
                return true;
            }
            if (piece) {
                if (piece->state == State::Decoded) counts.unused++;
                std::lock_guard<std::mutex> lock(mutex);
                piece->state = State::Taken;
                GiveBack(piece->a);
                GiveBack(piece->b);
            }
            // The lead decodes the range itself, with the true window.
            uint32_t const dict_len = static_cast<uint32_t>(std::min<uint64_t>(member_bytes, kWindow));
            if (end != kNone) {
                auto bytes = TakeBuffer();
                lead.Start(data, size, pos, end, window.data() + kWindow - dict_len, dict_len);
                RangeEnd const result = lead.Run(*bytes, kMaxPieceOutput, true);
                if (result == RangeEnd::AtBoundary || result == RangeEnd::Final) {
                    Item item;
                    Chained(bytes->data, bytes->size);
                    item.bytes = std::move(bytes);
                    if (result == RangeEnd::Final) MemberEnd(lead.final_bit(), item);
                    else pos = end;
                    Push(std::move(item));
                    dropped = 0;
                    counts.in_order++;
                    return true;
                }
                if (result == RangeEnd::Failed && dropped < 2) {
                    // A block start found by FindBlock that is none (decoded up to it, or the zero bits after it taken
                    // for codes): the piece from it is dropped, and the range decoded again to the next start. Twice in
                    // a row: a corrupt member, decoded to its end below, where the error shows.
                    std::lock_guard<std::mutex> lock(mutex);
                    starts[SlotOf(end)] = kNone;
                    dropped++;
                    counts.dropped_starts++;
                    return true;
                }
                // A corrupt member, or too large to hold: decoded again to the member's end, given as it comes.
            }
            lead.Start(data, size, pos, kNone, window.data() + kWindow - dict_len, dict_len);
            streaming = true;
            counts.streamed++;
            Item item;
            item.streamed = item.stream = true;
            Push(std::move(item));
            return true;
        }

        // The next segment of the lead's stream (the front item) into stream_bytes; at the member's end, its trailer and
        // what follows into the item, and the chain goes on.
        void StreamSegment(Item& item) {
            stream_bytes.size = 0;
            RangeEnd const result = lead.Run(stream_bytes, kStreamSegment, false);
            if (result == RangeEnd::Running) {
                Chained(stream_bytes.data, stream_bytes.size);
                return;
            }
            Chained(stream_bytes.data, stream_bytes.size);
            item.stream = false;
            streaming = false;
            if (result == RangeEnd::Final) {
                MemberEnd(lead.final_bit(), item);
            } else {
                item.error = lead.error() + " in gzip member " + std::to_string(member);
                chain_done = true;
            }
        }

        // ---- giving the output --------------------------------------------------------------------------------------

        // Copies up to `size` bytes into out; 0 at the end or on an error (error then set).
        size_t Read(char* out, size_t want) {
            size_t copied = 0;
            while (copied < want && error.empty()) {
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (items.empty() && chain_done) break;
                }
                // The chain one item ahead of the reader (the lead decodes in order what the reader needs next), and
                // further, up to `ahead` slots, over pieces the workers have decoded (cheap: their windows chained,
                // their resolution queued) and, while the next item waits for a worker to resolve it, over pieces no
                // worker has (decoded in order meanwhile).
                for (;;) {
                    int step;  // 0: none, 1: needed (waits for a worker decoding the piece), 2: while waiting
                    {
                        std::lock_guard<std::mutex> lock(mutex);
                        size_t const s = SlotOf(pos);
                        bool const room = !chain_done && !streaming && chain_slot < emit_slot + ahead;
                        bool const decoded = room && s < slots && starts[s] == pos && pieces[s].state == State::Decoded;
                        bool const front_waits = room && !items.empty() && !items.front().resolved;
                        step = items.size() < 2 ? 1 : decoded || front_waits ? 2 : 0;
                    }
                    if (step == 0 || !Advance(step == 1)) break;
                }
                Item* front;
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    if (items.empty()) {
                        if (chain_done) break;
                        continue;
                    }
                    front = &items.front();
                    emit_slot = front->slot != SIZE_MAX ? front->slot : chain_slot;
                    if (!front->resolved && !front->resolving) {
                        front->resolving = true;  // the reader resolves it itself
                        Piece& piece = pieces[front->slot];
                        lock.unlock();
                        parallel_gzip::Resolve(piece.a->data, piece.b->data, piece.a->size, front->window.get());
                        lock.lock();
                        GiveBack(piece.b);
                        front->resolved = true;
                        front->window.reset();
                    }
                    changed.wait(lock, [&] { return front->resolved; });
                }
                uint8_t const* bytes;
                size_t n;
                if (front->streamed) {
                    if (front->stream && given == stream_bytes.size) {
                        StreamSegment(*front);
                        given = 0;
                        Release();
                    }
                    bytes = stream_bytes.data;
                    n = stream_bytes.size;
                } else if (front->bytes) {
                    bytes = front->bytes->data;
                    n = front->bytes->size;
                } else {
                    Bytes const& a = *pieces[front->slot].a;
                    bytes = a.data;
                    n = a.size;
                }
                size_t const take = std::min(want - copied, n - given);
                if (take > 0) {
                    std::memcpy(out + copied, bytes + given, take);
                    crc = crc32_gzip_refl(crc, bytes + given, take);
                    member_out += take;
                    copied += take;
                    given += take;
                }
                if (given < n || front->stream) continue;
                // The item is given whole: what follows it.
                if (front->member_end) {
                    if (crc != front->trailer_crc) {
                        error = "gzip member's CRC-32 does not match its trailer (corrupt data)";
                    } else if (static_cast<uint32_t>(member_out) != front->trailer_isize) {
                        error = "gzip member holds " + std::to_string(member_out) + " bytes, its trailer says " +
                                std::to_string(front->trailer_isize) + " (modulo 2^32)";
                    }
                    crc = 0;
                    member_out = 0;
                }
                if (error.empty() && !front->error.empty()) error = front->error;
                bool const eof = front->eof;
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (front->slot != SIZE_MAX) GiveBack(pieces[front->slot].a);
                    GiveBack(front->bytes);
                    if (front->streamed) stream_bytes.size = 0;
                    items.pop_front();
                    given = 0;
                    changed.notify_all();
                }
                Release();
                if (eof) break;
            }
            return copied;
        }
    };

    std::atomic<size_t>& ParallelGzip::PieceBytesSetting() {
        static std::atomic<size_t> bytes{ 0 };
        return bytes;
    }

    std::atomic<size_t>& ParallelGzip::GapSetting() {
        static std::atomic<size_t> gap{ 0 };
        return gap;
    }

    ParallelGzip::Counts ParallelGzip::counts() const { return m_impl ? m_impl->counts : Counts{}; }
    uint64_t ParallelGzip::StoppedAt() const { return m_impl ? m_impl->stopped_at : kNotStopped; }
    size_t ParallelGzip::MembersRead() const { return m_impl ? m_impl->members_read : 0; }
    size_t ParallelGzip::MembersBegun() const { return m_impl ? m_impl->member : 0; }

    ParallelGzip::ParallelGzip() = default;
    ParallelGzip::~ParallelGzip() = default;

    std::unique_ptr<ParallelGzip> ParallelGzip::Open(int fd, unsigned threads, std::string& error, size_t max_member) {
        struct stat st {};
        if (::fstat(fd, &st) != 0) {
            error = std::strerror(errno);
            return nullptr;
        }
        if (!S_ISREG(st.st_mode) || st.st_size < 18) {
            error = "not a regular gzip file";
            return nullptr;
        }
        size_t const size = static_cast<size_t>(st.st_size);
        void* map = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
        if (map == MAP_FAILED) {
            error = std::string("it cannot be mapped (") + std::strerror(errno) + ")";
            return nullptr;
        }
        std::unique_ptr<ParallelGzip> gz;
        try {
            gz.reset(new ParallelGzip());
            gz->m_impl = std::make_unique<Impl>();
        } catch (std::exception const& e) {
            ::munmap(map, size);
            error = e.what();
            return nullptr;
        }
        try {
            Impl& impl = *gz->m_impl;  // unmaps the file when it goes
            impl.data = static_cast<uint8_t const*>(map);
            impl.size = size;
            uint64_t start = 0;
            if (!impl.ReadHeader(0, start)) {
                error = "its first gzip header cannot be read";
                return nullptr;
            }
            impl.max_member = max_member;
            impl.threads = std::max(threads, 2u);
            impl.gap = GapSetting().load(std::memory_order_relaxed);
            // The workers from threads + 2 pieces ahead of the reader's thread, up to threads further: each piece in
            // that look-ahead holds up to two output buffers of ~4-16 MB (on FASTQ), so it sets the memory. Against
            // 2 * threads + 4 and 2 * threads (until 2026-10-13): -29 to -38% of the memory on 2-4 threads, the fastest
            // runs within 2% (protal docs/claude/2026-10-13-reading-input).
            if (impl.gap == 0) impl.gap = size_t{impl.threads} + 2;
            impl.ahead = impl.gap + size_t{impl.threads};
            impl.piece_bytes = PieceBytesSetting().load(std::memory_order_relaxed);
            if (impl.piece_bytes == 0) impl.piece_bytes = kChunkBytes;
            impl.slots = (size + impl.piece_bytes - 1) / impl.piece_bytes;
            impl.starts.assign(impl.slots, kUnknown);
            impl.starts[0] = start;
            impl.pieces.resize(impl.slots);
            impl.pieces[0].state = Impl::State::Taken;
            impl.pos = start;
            // Read in order, mostly: aggressive read-ahead where a page is faulted in, and the first window asked for.
            if (ReadAheadWanted()) ::madvise(map, size, MADV_SEQUENTIAL);
            impl.Release();
            TheLabels();
            for (unsigned w = 1; w < impl.threads; w++) impl.workers.emplace_back(&Impl::Work, &impl);
        } catch (std::exception const& e) {
            error = e.what();
            return nullptr;
        }
        return gz;
    }

    size_t ParallelGzip::Read(char* out, size_t size) {
        if (!m_error.empty() || !m_impl) return 0;
        size_t const n = m_impl->Read(out, size);
        if (!m_impl->error.empty() && n == 0) m_error = m_impl->error;
        return n;
    }

}  // namespace sdm
