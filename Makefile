# Default: one fully static sdm executable with HTSlib enabled.
# make WITH_HTS=0 builds sdm without alignment input support.
# make STATIC=0 builds sdm with shared libraries.
CXX ?= g++
CC ?= gcc
PKG_CONFIG ?= pkg-config
WITH_HTS ?= 1
STATIC ?= 1
CXXFLAGS ?= -O3
CPPFLAGS += -I. -D__STDC_CONSTANT_MACROS
SDM_CXXFLAGS = -std=c++20 -pthread $(isa_FLAGS)
SDM_LDFLAGS = $(if $(filter 1,$(STATIC)),-static)
SDM_LIBS = -lz -pthread $(isa_FLAGS2)
ifeq ($(shell uname),Linux)
SDM_LIBS += -lrt
ifeq ($(STATIC),1)
# Before glibc 2.34 the thread functions live in libpthread.a and libstdc++ refers to them only
# weakly, so a static link can leave them out and std::thread/std::mutex break at run time
# (e.g. with conda's glibc 2.17/2.28 sysroots). From glibc 2.34 on libpthread.a is empty.
SDM_LIBS += -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif
endif

ifeq ($(WITH_HTS),1)
HTS_CFLAGS ?= $(shell $(PKG_CONFIG) --cflags htslib 2>/dev/null)
HTS_LIBS ?= $(shell $(PKG_CONFIG) --libs $(if $(filter 1,$(STATIC)),--static) htslib 2>/dev/null)
ifeq ($(strip $(HTS_LIBS)),)
HTS_LIBS = -lhts
endif
SDM_CPPFLAGS = $(HTS_CFLAGS) -DSDM_WITH_HTSLIB
SDM_HTS_LIBS = $(HTS_LIBS)
HTS_CHECK = check-htslib
endif

# libdeflate: faster BGZF compression, decompression and CRC-32 than zlib (same
# decompressed content). WITH_LIBDEFLATE=auto uses it when a test program links
# (for STATIC=1 that needs libdeflate.a), 1 requires it, 0 keeps zlib.
WITH_LIBDEFLATE ?= auto
LIBDEFLATE_LIBS ?= -ldeflate
ifeq ($(WITH_LIBDEFLATE),auto)
# A '#' inside $(shell) needs escaping before GNU make 4.3 and must not be escaped after.
HASH := \#
USE_LIBDEFLATE := $(shell printf '$(HASH)include <libdeflate.h>\nint main() { return libdeflate_crc32(0, 0, 0) != 0; }\n' | \
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(SDM_LDFLAGS) -x c++ - $(LIBDEFLATE_LIBS) -o /dev/null > /dev/null 2>&1 && echo 1 || echo 0)
else
USE_LIBDEFLATE := $(WITH_LIBDEFLATE)
endif
ifeq ($(USE_LIBDEFLATE),1)
SDM_CPPFLAGS += -DSDM_WITH_LIBDEFLATE
SDM_LIBS := $(LIBDEFLATE_LIBS) $(SDM_LIBS)
endif

SOURCES = $(wildcard *.cpp)
C_SOURCES = $(wildcard *.c)
BUILD_DIR = .build/$(if $(filter 1,$(WITH_HTS)),hts,portable)-$(STATIC)
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o) $(C_SOURCES:%.c=$(BUILD_DIR)/%.o)
TEST_OBJECTS = $(filter-out $(BUILD_DIR)/Demultipl.o,$(OBJECTS))
TEST_PROGRAMS = $(addprefix $(BUILD_DIR)/,TestAuditFixes TestInputHelpers TestIfBufStream TestGzipInput \
	TestBgzfOutput TestRegression TestSearch TestSeedOutputRouting TestCoarseDerep TestPartitionWorkPool)
ifeq ($(WITH_HTS),1)
TEST_PROGRAMS += $(BUILD_DIR)/TestHtsInput
endif
# Extra runner options, e.g. make test TEST_ARGS=--quick or TEST_ARGS="test_regressions_cli TestAuditFixes"
TEST_ARGS ?=

.DEFAULT_GOAL := all
.PHONY: all portable hts sdm clean distclean check-htslib test test-quick test-portable test-hts
all: sdm

# Convenience aliases select the configuration of the same executable.
portable:
	$(MAKE) WITH_HTS=0 sdm
hts:
	$(MAKE) WITH_HTS=1 sdm

# Separate object directories prevent macro/ABI mixing between configurations.
# Always relink sdm so switching configurations also replaces the executable.
sdm: $(OBJECTS) | $(HTS_CHECK)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SDM_LDFLAGS) $^ $(SDM_HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp Makefile | $(HTS_CHECK)
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(SDM_CPPFLAGS) $(CXXFLAGS) $(SDM_CXXFLAGS) -MMD -MP -c $< -o $@
$(BUILD_DIR)/%.o: %.c Makefile | $(HTS_CHECK)
	@mkdir -p $(@D)
	$(CC) $(CPPFLAGS) $(SDM_CPPFLAGS) $(CFLAGS) -MMD -MP -c $< -o $@

check-htslib:
	@mkdir -p $(BUILD_DIR)
	@printf '#include <htslib/hts.h>\nint main() { return hts_version() == nullptr; }\n' | \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(HTS_CFLAGS) $(LDFLAGS) $(SDM_LDFLAGS) -x c++ - \
	  -x none $(HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $(BUILD_DIR)/check-htslib > $(BUILD_DIR)/check-htslib.log 2>&1 || \
	  { cat $(BUILD_DIR)/check-htslib.log; \
	    echo 'HTSlib headers and libraries are required for the default sdm build.'; \
	    echo 'Install HTSlib development files and pkg-config, or set HTS_CFLAGS/HTS_LIBS.'; \
	    echo 'Static linking also requires static archives for all HTSlib dependencies.'; \
	    echo 'Use make STATIC=0 for shared linking, or make WITH_HTS=0 to omit HTSlib.'; exit 1; }

$(BUILD_DIR)/Test%: tests/Test%.cpp tests/TestSupport.cpp $(TEST_OBJECTS)
	$(CXX) $(CPPFLAGS) $(SDM_CPPFLAGS) $(CXXFLAGS) $(SDM_CXXFLAGS) $^ $(LDFLAGS) $(SDM_LDFLAGS) $(SDM_HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $@
# Benchmarks/experiments cited in docs/ live outside tests/ and are not part of make test;
# build one on request, e.g. make WITH_HTS=0 STATIC=0 .build/portable-0/TestCoarseBenchmark
BENCH_DIR := .claude/archive/benchmarks
$(BUILD_DIR)/Test%: $(BENCH_DIR)/Test%.cpp tests/TestSupport.cpp $(TEST_OBJECTS)
	$(CXX) $(CPPFLAGS) -Itests $(SDM_CPPFLAGS) $(CXXFLAGS) $(SDM_CXXFLAGS) $^ $(LDFLAGS) $(SDM_LDFLAGS) $(SDM_HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $@

# Runs every C++ and pipeline suite, continues past failures, prints a summary and
# exits non-zero if any suite failed. Logs: $(BUILD_DIR)/test-logs/<suite>.log
test: sdm $(TEST_PROGRAMS)
	python3 tests/run_tests.py --sdm ./sdm --bin-dir $(BUILD_DIR) $(if $(filter 1,$(WITH_HTS)),--hts) $(TEST_ARGS)
# Skips the slow parity and seed-extension suites (about a minute).
test-quick: sdm $(TEST_PROGRAMS)
	python3 tests/run_tests.py --sdm ./sdm --bin-dir $(BUILD_DIR) $(if $(filter 1,$(WITH_HTS)),--hts) --quick $(TEST_ARGS)
test-portable:
	$(MAKE) WITH_HTS=0 test
test-hts:
	$(MAKE) WITH_HTS=1 test

clean:
	$(RM) -r .build
	$(RM) sdm sdm-hts $(SOURCES:.cpp=.o) $(C_SOURCES:.c=.o)
distclean: clean

# The BGZF stream objects depend on the libdeflate setting: rebuild them when it changes.
DEFLATE_STAMP = $(BUILD_DIR)/libdeflate-$(USE_LIBDEFLATE).stamp
$(DEFLATE_STAMP):
	@mkdir -p $(@D)
	@rm -f $(BUILD_DIR)/libdeflate-*.stamp
	@touch $@
$(BUILD_DIR)/OutputFileStream.o $(BUILD_DIR)/InputFileStream.o: $(DEFLATE_STAMP)

-include $(OBJECTS:.o=.d)
