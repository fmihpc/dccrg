FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: load

load: load.cpp ../../dccrg.hpp ../../dccrg_index.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI load.cpp $(INCLUDES) -o load

d: data
data:
	rm -rf load.visit load_*vtk

c: clean
clean: data
	rm -f load
