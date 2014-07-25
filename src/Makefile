CXXINCLUDE=-I ../boost-numeric-bindings -I/usr/include/atlas
CXXOPT=-O0 -g3
CCXFLAGS=-std=c++1y -DBOOST_LOG_DYN_LINK -pthread -Wall -Wextra -Wl,--no-as-needed -pedantic $(CXXINCLUDE) $(CXXOPT)
CXX=g++-4.8

all: myriaworld.h.gch

myriaworld.h.gch: myriaworld.hpp
	$(CXX) -x c++-header $< $(CXXFLAGS) -c -o $@

clean:
	rm -f myriaworld.h.gch
