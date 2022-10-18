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


CPPFLAGS += -O3 -std=c++17 -lrt
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CXXFLAGS=-D__STDC_CONSTANT_MACROS
LDFLAGS +=   $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir)) -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))
LDLIBS += -lz

.PHONY: all clean distclean

all: $(program_NAME)
#-static
$(program_NAME): $(program_OBJS)
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) $(LDLIBS) -static
clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	@- $(RM) sdm

distclean: clean

#say_link: echo $(LINK.cc)
#	echo "hallo"
