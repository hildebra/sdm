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

SOURCES = $(wildcard *.cpp)
C_SOURCES = $(wildcard *.c)
BUILD_DIR = .build/$(if $(filter 1,$(WITH_HTS)),hts,portable)-$(STATIC)
OBJECTS = $(SOURCES:%.cpp=$(BUILD_DIR)/%.o) $(C_SOURCES:%.c=$(BUILD_DIR)/%.o)
TEST_OBJECTS = $(filter-out $(BUILD_DIR)/Demultipl.o,$(OBJECTS))
TEST_PROGRAMS = $(BUILD_DIR)/TestInputHelpers $(BUILD_DIR)/TestRegression $(BUILD_DIR)/TestSearch $(BUILD_DIR)/TestCoarseDerep
ifeq ($(WITH_HTS),1)
TEST_PROGRAMS += $(BUILD_DIR)/TestHtsInput
endif

.DEFAULT_GOAL := all
.PHONY: all portable hts sdm clean distclean check-htslib test test-portable test-hts
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

test: sdm $(TEST_PROGRAMS)
	$(BUILD_DIR)/TestInputHelpers
	$(BUILD_DIR)/TestRegression
	$(BUILD_DIR)/TestSearch
	$(BUILD_DIR)/TestCoarseDerep
	python3 tests/test_ont_pipeline.py ./sdm
	python3 tests/test_diff_quality_pipeline.py ./sdm
	python3 tests/test_paired_diff_quality_pipeline.py ./sdm
ifeq ($(WITH_HTS),1)
	$(BUILD_DIR)/TestHtsInput
	python3 tests/test_hts_pipeline.py ./sdm $(BUILD_DIR)/TestHtsInput
endif
test-portable:
	$(MAKE) WITH_HTS=0 test
test-hts:
	$(MAKE) WITH_HTS=1 test

clean:
	$(RM) -r .build
	$(RM) sdm sdm-hts $(SOURCES:.cpp=.o) $(C_SOURCES:.c=.o)
distclean: clean

-include $(OBJECTS:.o=.d)
