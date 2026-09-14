# modified by Falk Hildebrand
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
CPPFLAGS +=  -std=c++20 -pthread $(isa_FLAGS)
#DEBUG - don't use -static in linker!
#CPPFLAGS += -g -O1 -fno-omit-frame-pointer -fsanitize=address,undefined
#CPPFLAGS += -O1 -g

CPPFLAGS += -O3

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
CXXFLAGS=-D__STDC_CONSTANT_MACROS
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDLIBS += $(foreach library,$(program_LIBRARIES),-l$(library))
LDLIBS += -lz $(isa_FLAGS2)
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
