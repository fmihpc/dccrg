FLAGS = -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: load

load: load.cpp ../../dccrg.hpp ../../dccrg_index.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) load.cpp $(INCLUDES) -o load

c: clean
clean:
	rm -f load load.visit load_*vtk

