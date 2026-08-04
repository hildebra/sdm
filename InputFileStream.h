#pragma once

#include <cstddef>
#include <ios>
#include <istream>
#include <memory>
#include <string>

// Opens plain, gzip, or indexed-BGZF input through one selection point.
// Indexed BGZF uses ordered parallel block inflation when parallelGzip is true
// and sufficient gzip worker capacity is available.
std::unique_ptr<std::istream> openInputFile(
    const std::string& filename,
    std::ios_base::openmode mode = std::ios_base::in,
    size_t bufferSize = 0,
    bool parallelGzip = true);

// Exposed for focused verification and diagnostics.
bool inputFileUsesParallelBgzf(const std::istream& stream);
