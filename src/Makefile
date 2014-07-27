CXXINCLUDE=-I ../boost-numeric-bindings -I/usr/include/atlas
CXXOPT=-O0 -g3
CXXFLAGS=-std=c++1y -DBOOST_LOG_DYN_LINK -pthread -Wall -Wextra -Wl,--no-as-needed -pedantic $(CXXINCLUDE) $(CXXOPT)
CXX=g++-4.8
LDFLAGS=-pthread -Wl,--no-as-needed -lboost_log -lboost_system -lboost_thread -llapack -llapack_atlas

all: myriaworld
pch: myriaworld.h.gch
objs: myriaworld.o geo.o read_shapefile_data.o SPHEmesh.o country2tria.o cutting.o

myriaworld: pch objs
	g++ $(LDFLAGS) -o myriaworld geo.o read_shapefile_data.o SPHEmesh.o country2tria.o cutting.o myriaworld.o

myriaworld.h.gch: myriaworld.h
	$(CXX) -x c++-header $< $(CXXFLAGS) -c -o $@

geo.o: geo.cpp geo.hpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c geo.cpp -o $@

SPHEmesh.o: SPHEmesh.cpp SPHEmesh.hpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c SPHEmesh.cpp -o $@

read_shapefile_data.o: read_shapefile_data.cpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c read_shapefile_data.cpp -o $@

country2tria.o: country2tria.cpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c country2tria.cpp -o $@

cutting.o: cutting.cpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c cutting.cpp -o $@

myriaworld.o: myriaworld.cpp myriaworld.h.gch
	$(CXX) $(CXXFLAGS) -c myriaworld.cpp -o $@

clean:
	rm -f *.o
	rm -f myriaworld.h.gch
	rm -f myriaworld.hpp.gch
