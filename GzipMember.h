// GzipMember.h - where a gzip member's header begins, told from the bytes alone, for reading a gzip file of several
// members on several threads: ThreadedGzStreambuf cuts such a file at its members' headers and inflates whole members
// on its workers; ParallelGzip, reading a file in pieces, stops where small members follow so that ThreadedGzStreambuf
// reads the rest by members (protal docs/claude/2026-10-13-gzip-members-regression).
#pragma once

#include <isa-l/crc.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace sdm::gzip_member {

    inline constexpr size_t kMinBytes = 20;               // a member's least bytes: header 10, deflate 2, trailer 8
    inline constexpr size_t kMaxHeader = size_t{1} << 17;  // a header's most bytes looked at (its optional fields)

    // Whether a gzip member's header begins at base[at], `at` bytes after the start of the member before (at least
    // kMinBytes): its magic and method, flags without reserved bits, extra flags and system as writers set them, its
    // optional fields within kMaxHeader bytes (and their CRC right, if the header has one), its first deflate block's
    // header possible, and the length in the trailer before it no more than `at` bytes of deflate data hold (1032 to 1).
    // Deflate data rarely passes this by chance, and a reader must still check that the member before ends there; a
    // member it misses is only read with the one before. 1: yes; 0: no; -1: more than base's `size` bytes are needed
    // to tell (never with `all`, the bytes up to the file's end).
    inline int StartsAt(unsigned char const* base, size_t size, size_t at, bool all) {
        unsigned char const* const p = base + at;
        size_t const n = size - at;
        int const more = all ? 0 : -1;
        if (n < 3) return p[0] == 0x1f && (n < 2 || p[1] == 0x8b) ? more : 0;
        if (p[0] != 0x1f || p[1] != 0x8b || p[2] != 8) return 0;
        if (n < 10) return more;
        unsigned const flags = p[3];
        if ((flags & 0xe0) != 0 || (p[8] != 0 && p[8] != 2 && p[8] != 4) || (p[9] > 13 && p[9] != 255)) return 0;
        size_t h = 10;
        if (flags & 4) {  // an extra field
            if (n < h + 2) return more;
            h += 2 + (size_t{ p[h] } | size_t{ p[h + 1] } << 8);
        }
        for (unsigned const field : { 8u, 16u }) {  // a name, a comment: zero-terminated
            if (!(flags & field)) continue;
            size_t const limit = std::min(n, kMaxHeader);
            auto const* zero = h < limit ? static_cast<unsigned char const*>(std::memchr(p + h, 0, limit - h)) : nullptr;
            if (!zero) return n < kMaxHeader ? more : 0;
            h = static_cast<size_t>(zero - p) + 1;
        }
        if (h + 7 > kMaxHeader) return 0;
        if (flags & 2) {  // the header's CRC: the low 16 bits of its CRC-32
            if (n < h + 2) return more;
            if ((crc32_gzip_refl(0, p, h) & 0xffff) != (uint32_t{ p[h] } | uint32_t{ p[h + 1] } << 8)) return 0;
            h += 2;
        }
        if (n < h + 5) return more;
        unsigned const type = (p[h] >> 1) & 3;
        if (type == 3) return 0;
        if (type == 0 && ((p[h + 1] ^ p[h + 3]) != 0xff || (p[h + 2] ^ p[h + 4]) != 0xff)) return 0;  // LEN, NLEN
        if (type == 2 && (((p[h] >> 3) & 31) > 29 || (p[h + 1] & 31) > 29)) return 0;              // HLIT, HDIST
        uint32_t const isize = uint32_t{ p[-4] } | uint32_t{ p[-3] } << 8 | uint32_t{ p[-2] } << 16 | uint32_t{ p[-1] } << 24;
        return isize <= 1032 * uint64_t{ at } ? 1 : 0;
    }

    // Where the next member's header begins after the member that begins at data[from], at most `limit` bytes on
    // (StartsAt; data's `size` bytes reach the file's end): its distance from `from`, 0 if none does.
    inline size_t Next(unsigned char const* data, size_t size, size_t from, size_t limit) {
        if (from >= size) return 0;
        unsigned char const* const base = data + from;
        size_t const n = size - from;
        size_t const end = std::min(n, limit + 1);
        for (size_t at = kMinBytes; at < end;) {
            auto const* hit = static_cast<unsigned char const*>(std::memchr(base + at, 0x1f, end - at));
            if (!hit) break;
            at = static_cast<size_t>(hit - base);
            if (StartsAt(base, n, at, true) > 0) return at;
            at++;
        }
        return 0;
    }

}  // namespace sdm::gzip_member
