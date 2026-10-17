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
SDM_LIBS = -pthread $(isa_FLAGS2)
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

# ISA-L (igzip: every gzip input, BGZF or not, and the BGZF of .gz output) and zstd
# (.zst input and output) are required. Their flags come from ISAL_CFLAGS/ISAL_LIBS and
# ZSTD_CFLAGS/ZSTD_LIBS if given, else from the static libraries ./build_deps.sh built
# into DEPS_PREFIX (make deps), else from pkg-config (libisal, libzstd), else -lisal and
# -lzstd. STATIC=1 needs libisal.a and libzstd.a; Debian/Ubuntu's libisal-dev has no
# libisal.a, so build it with make deps.
DEPS_PREFIX ?= .deps
ifneq ($(wildcard $(DEPS_PREFIX)/lib/libisal.a),)
ISAL_CFLAGS ?= -I$(DEPS_PREFIX)/include
ISAL_LIBS ?= $(DEPS_PREFIX)/lib/libisal.a
endif
ifneq ($(wildcard $(DEPS_PREFIX)/lib/libzstd.a),)
ZSTD_CFLAGS ?= -I$(DEPS_PREFIX)/include
ZSTD_LIBS ?= $(DEPS_PREFIX)/lib/libzstd.a -pthread
endif
ifeq ($(origin ISAL_CFLAGS),undefined)
ISAL_CFLAGS := $(shell $(PKG_CONFIG) --cflags libisal 2>/dev/null)
endif
ifeq ($(origin ISAL_LIBS),undefined)
ISAL_LIBS := $(or $(strip $(shell $(PKG_CONFIG) --libs $(if $(filter 1,$(STATIC)),--static) libisal 2>/dev/null)),-lisal)
endif
ifeq ($(origin ZSTD_CFLAGS),undefined)
ZSTD_CFLAGS := $(shell $(PKG_CONFIG) --cflags libzstd 2>/dev/null)
endif
ifeq ($(origin ZSTD_LIBS),undefined)
ZSTD_LIBS := $(or $(strip $(shell $(PKG_CONFIG) --libs $(if $(filter 1,$(STATIC)),--static) libzstd 2>/dev/null)),-lzstd)
endif
SDM_CPPFLAGS += $(ISAL_CFLAGS) $(ZSTD_CFLAGS)
SDM_LIBS := $(ISAL_LIBS) $(ZSTD_LIBS) $(SDM_LIBS)

SOURCES = $(wildcard *.cpp)
C_SOURCES = $(wildcard *.c)
BUILD_DIR = .build/$(if $(filter 1,$(WITH_HTS)),hts,portable)-$(STATIC)
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o) $(C_SOURCES:%.c=$(BUILD_DIR)/%.o)
TEST_OBJECTS = $(filter-out $(BUILD_DIR)/Demultipl.o,$(OBJECTS))
TEST_PROGRAMS = $(addprefix $(BUILD_DIR)/,TestAuditFixes TestInputHelpers TestIfBufStream TestGzipInput \
	TestBgzfOutput TestThreadedGzStream TestRegression TestSearch TestSeedOutputRouting TestCoarseDerep TestPartitionWorkPool)
ifeq ($(WITH_HTS),1)
TEST_PROGRAMS += $(BUILD_DIR)/TestHtsInput
endif
# Extra runner options, e.g. make test TEST_ARGS=--quick or TEST_ARGS="test_regressions_cli TestAuditFixes"
TEST_ARGS ?=

.DEFAULT_GOAL := all
.PHONY: all portable hts sdm deps clean distclean check-htslib check-compression test test-quick test-portable test-hts
all: sdm

# Convenience aliases select the configuration of the same executable.
portable:
	$(MAKE) WITH_HTS=0 sdm
hts:
	$(MAKE) WITH_HTS=1 sdm

# Separate object directories prevent macro/ABI mixing between configurations.
# Always relink sdm so switching configurations also replaces the executable.
sdm: $(OBJECTS) | $(HTS_CHECK) check-compression
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SDM_LDFLAGS) $^ $(SDM_HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp Makefile | $(HTS_CHECK) check-compression
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(SDM_CPPFLAGS) $(CXXFLAGS) $(SDM_CXXFLAGS) -MMD -MP -c $< -o $@
$(BUILD_DIR)/%.o: %.c Makefile | $(HTS_CHECK) check-compression
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

# ISA-L and zstd, linked as sdm links them.
check-compression:
	@mkdir -p $(BUILD_DIR)
	@printf '#include <isa-l/igzip_lib.h>\n#include <zstd.h>\nint main() { inflate_state s; isal_inflate_init(&s); return ZSTD_versionNumber() == 0; }\n' | \
	  $(CXX) $(CPPFLAGS) $(SDM_CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(SDM_LDFLAGS) -x c++ - \
	  -x none $(LDLIBS) $(SDM_LIBS) -o $(BUILD_DIR)/check-compression > $(BUILD_DIR)/check-compression.log 2>&1 || \
	  { cat $(BUILD_DIR)/check-compression.log; \
	    echo 'ISA-L and zstd headers and libraries are required (ISA-L: all gzip input and the BGZF of .gz output;'; \
	    echo 'zstd: .zst input and output). Debian/Ubuntu: apt install libisal-dev libzstd-dev; conda: isa-l zstd.'; \
	    echo 'A static build (STATIC=1, the default) needs libisal.a and libzstd.a, which Debian/Ubuntu do not'; \
	    echo 'all ship: run make deps (./build_deps.sh: builds both into $(DEPS_PREFIX), which make then uses),'; \
	    echo 'or set ISAL_CFLAGS/ISAL_LIBS and ZSTD_CFLAGS/ZSTD_LIBS, or use make STATIC=0.'; exit 1; }

# Static ISA-L and zstd (and nasm if needed) from their release tarballs into DEPS_PREFIX.
deps:
	bash ./build_deps.sh --prefix $(DEPS_PREFIX)

# The test programs are rebuilt when a header of tests/ (e.g. GzipTestUtil.h) changes too.
TEST_HEADERS := $(wildcard tests/*.h)
$(BUILD_DIR)/Test%: tests/Test%.cpp tests/TestSupport.cpp $(TEST_OBJECTS) $(TEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(SDM_CPPFLAGS) $(CXXFLAGS) $(SDM_CXXFLAGS) $(filter-out %.h,$^) $(LDFLAGS) $(SDM_LDFLAGS) $(SDM_HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $@
# Benchmarks/experiments cited in docs/ live outside tests/ and are not part of make test;
# build one on request, e.g. make WITH_HTS=0 STATIC=0 .build/portable-0/TestCoarseBenchmark
BENCH_DIR := .claude/archive/benchmarks
$(BUILD_DIR)/Test%: $(BENCH_DIR)/Test%.cpp tests/TestSupport.cpp $(TEST_OBJECTS)
	$(CXX) $(CPPFLAGS) -Itests $(SDM_CPPFLAGS) $(CXXFLAGS) $(SDM_CXXFLAGS) $^ $(LDFLAGS) $(SDM_LDFLAGS) $(SDM_HTS_LIBS) $(LDLIBS) $(SDM_LIBS) -o $@

# Runs every C++ and pipeline suite, continues past failures, prints a summary and
# exits non-zero if any suite failed. Logs: $(BUILD_DIR)/test-logs/<suite>.log
test: sdm $(TEST_PROGRAMS)
	perl tests/run_tests.pl --sdm ./sdm --bin-dir $(BUILD_DIR) $(if $(filter 1,$(WITH_HTS)),--hts) $(TEST_ARGS)
# Skips the slow parity and seed-extension suites (about a minute).
test-quick: sdm $(TEST_PROGRAMS)
	perl tests/run_tests.pl --sdm ./sdm --bin-dir $(BUILD_DIR) $(if $(filter 1,$(WITH_HTS)),--hts) --quick $(TEST_ARGS)
test-portable:
	$(MAKE) WITH_HTS=0 test
test-hts:
	$(MAKE) WITH_HTS=1 test

clean:
	$(RM) -r .build
	$(RM) sdm sdm-hts $(SOURCES:.cpp=.o) $(C_SOURCES:.c=.o)
distclean: clean
	$(RM) -r $(DEPS_PREFIX)

-include $(OBJECTS:.o=.d)
