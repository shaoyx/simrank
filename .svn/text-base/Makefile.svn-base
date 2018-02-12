# Build targets
SimRank_O= core/SRBenchmark.cpp \
           core/srgs.cpp \
		   core/accuratesimrank.cpp \
		   core/FingerPrint.cpp \
		   core/SimRankMOD14.cpp \
		   core/SimRankMOD14Index.cpp \
		   core/SimRankWWW05.cpp \
		   core/topsim.cpp \
		   core/logmanager.cpp \
		   core/onewaygraph.cpp \
		   core/dynamicowgmanager.cpp \
           core/mytime.cpp \

USER_FLAGS=-Wno-unused-result -Wno-unused-but-set-variable -Wno-sign-compare
USER_LIBS=

# Compiler flags
CC  = gcc
CXX = g++
CFLAGS=$(USER_FLAGS) -O3 -fPIC -Wall -g -I. -I./include
CXXFLAGS=$(CFLAGS)
LDFLAGS=-lpthread

# The name of the excution that will be built
SimRank=SRbenchmark

all:
	$(CXX) $(CXXFLAGS) -pthread $(LDFLAGS) -o $(SimRank) $(SimRank_O) $(USER_LIBS)
	
clean:
	rm -f $(SimRank)
	find . -name '*.o' -print | xargs rm -f
