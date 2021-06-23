# copyright: Michael Safyan
# modified by Falk Hildebrand


#CXX := $(if $(GCC),$(GCC),$(CXX))


program_NAME := sdm
program_C_SRCS := $(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o} 
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := ./
program_LIBRARY_DIRS := ${CPATH}
program_LIBRARIES :=

UNAME := $(shell uname)
CPPFLAGS += -O3 -std=c++17 -pthread
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
ifeq ($(UNAME), Darwin)
# 'shared_mutex' is available from macOS 10.12
CPPFLAGS += -mmacosx-version-min=10.12
endif
CXXFLAGS=-D__STDC_CONSTANT_MACROS
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))
LDLIBS += -lz
ifeq ($(UNAME), Linux)
LDLIBS += -lrt -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
else
LDLIBS += -lpthread
endif

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) $(LDLIBS) -static

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	@- $(RM) sdm

distclean: clean

#say_link: echo $(LINK.cc)
#	echo "hallo"
