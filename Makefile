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


CPPFLAGS +=-Wall -O3 -lz -D__USE_XOPEN2K8 -std=c++17 -pthread
CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

.PHONY: all clean distclean

all: $(program_NAME)
	echo ${GCC}
	echo ${CPP}
	echo ${CXX}
	echo ${CC}


$(program_NAME): $(program_OBJS)
#	$(LINK.cc) $(program_OBJS) -o $(program_NAME) -lz --verbose
	$(LINK.cc) $(program_OBJS) -o $(program_NAME) -lz -static --verbose

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)
	rm sdm

distclean: clean

#say_link: echo $(LINK.cc)
#	echo "hallo"
