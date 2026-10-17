#include "Benchmark.h"

// The production definition is kept with the executable entry point. Unit-test
// binaries deliberately omit Demultipl.cpp (and its main), so provide the inert
// implementation needed to link code paths that reference the global tracker.
void Instr::SimpleMemoryTracker::sample() {
}
