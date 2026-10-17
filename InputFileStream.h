#pragma once

#include <cstddef>
#include <ios>
#include <istream>
#include <memory>
#include <string>

// Opens an input file through one selection point. The format is told from the
// file's first bytes, not its name: gzip (BGZF or not) and zstd are decompressed
// with ISA-L and libzstd in threads of their own (ThreadedGzStream.h), anything
// else is read as it is. With parallelGzip and enough gzip worker capacity, BGZF
// blocks, zstd frames and the members of multi-member gzip are decompressed on
// several threads (single-member gzip on one; docs/compression.md); the bytes are
// the same either way. A compressed file that ends early (truncated, a BGZF file
// without its end-of-file block, data after a member that is no member) is read up to
// the damage with a warning on stderr; corrupt compressed data (a wrong CRC-32 or
// length, invalid deflate or zstd data) throws std::runtime_error, naming the file,
// when reading reaches it.
std::unique_ptr<std::istream> openInputFile(
    const std::string& filename,
    std::ios_base::openmode mode = std::ios_base::in,
    size_t bufferSize = 0,
    bool parallelGzip = true);

// Exposed for focused verification and diagnostics: whether the stream is a BGZF
// file whose blocks are inflated on several threads.
bool inputFileUsesParallelBgzf(const std::istream& stream);

// The compression libraries and their versions, for sdm's build information
// (e.g. "gzip/BGZF: ISA-L 2.32.1 (.gz output: BGZF); zstd: libzstd 1.5.7").
std::string compressionLibraries();
