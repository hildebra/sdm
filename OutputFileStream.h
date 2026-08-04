#pragma once

#include <cstddef>
#include <ios>
#include <memory>
#include <ostream>
#include <string>

// Open and finalize output files through one backend-selection point.  Files
// ending in .gz or .gzip are written as BGZF and receive a standard .gzi
// sidecar index. Parallel block compression is used when parallelGzip is true
// and sufficient gzip worker capacity is available. Other files use std::ofstream.
std::unique_ptr<std::ostream> openOutputFile(
    const std::string& filename,
    std::ios_base::openmode mode,
    size_t bufferSize = 0,
    bool parallelGzip = true);

// Flushes, closes, and validates an output stream.  This is deliberately
// explicit because destructors cannot safely report delayed write failures.
void finalizeOutputFile(
    std::unique_ptr<std::ostream>& stream,
    const std::string& filename);

