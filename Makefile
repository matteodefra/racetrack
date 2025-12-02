DIR := ${CURDIR}

CXX			= g++ -std=c++17 
CXXFLAGS  	= -g # -DBLOCKING_MODE -DFF_BOUNDED_BUFFER
 
LDFLAGS 	= -pthread
OPTFLAGS	= -fopenmp -finline-functions -w -DNDEBUG -O3

# PRELDFLAGS 	= LD_PRELOAD=${DIR}/jemalloc/lib/libjemalloc.so.2

TARGETS 	=	output/fulltest \
			output/scenario
						
.PHONY: all clean

.SUFFIXES: .cpp 

all: $(TARGETS)

output/fulltest: fulltest.cpp
	$(CXX) $^ $(CXXFLAGS) $(INCLUDES) $(OPTFLAGS) -o $@ $(LDFLAGS)

output/scenario: scenario.cpp
	$(CXX) $^ $(CXXFLAGS) $(INCLUDES) $(OPTFLAGS) -o $@ $(LDFLAGS)


clean:
	rm -rf output/*
