# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: init

init: init.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI init.cpp $(INCLUDES) -o init

c: clean
clean:
	rm -f init
