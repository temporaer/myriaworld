CXXINCLUDE=-I ../boost-numeric-bindings -I/usr/include/atlas
CXXOPT=-O0 -g3
CXXFLAGS=-std=c++1y -DBOOST_LOG_DYN_LINK -pthread -Wall -Wextra -Wl,--no-as-needed -pedantic $(CXXINCLUDE) $(CXXOPT)
CXX=g++-4.8

all: myriaworld.h.gch geo.o read_shapefile_data.o SPHEmesh.o

myriaworld.h.gch: myriaworld.h
	$(CXX) -x c++-header $< $(CXXFLAGS) -c -o $@

geo.o: geo.cpp geo.hpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c geo.cpp -o $@

SPHEmesh.o: SPHEmesh.cpp SPHEmesh.hpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c geo.cpp -o $@

read_shapefile_data.o: read_shapefile_data.cpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c geo.cpp -o $@

clean:
	rm -f myriaworld.h.gch
	rm -f myriaworld.hpp.gch
