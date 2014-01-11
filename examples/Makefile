# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
MPICXX = mpic++ $(FLAGS)

HEADERS = ../dccrg.hpp ../dccrg_cartesian_geometry.hpp

all: simple_game_of_life game_of_life game_of_life_with_output dc2vtk

game_of_life: game_of_life.cpp $(HEADERS) Makefile
	$(MPICXX) -DDCCRG_TRANSFER_USING_BOOST_MPI game_of_life.cpp $(INCLUDES) -o game_of_life

simple_game_of_life: simple_game_of_life.cpp $(HEADERS) Makefile
	$(MPICXX) -DDCCRG_TRANSFER_USING_BOOST_MPI simple_game_of_life.cpp $(INCLUDES) -o simple_game_of_life

game_of_life_with_output: game_of_life_with_output.cpp $(HEADERS) Makefile
	$(MPICXX) -DDCCRG_TRANSFER_USING_BOOST_MPI game_of_life_with_output.cpp $(INCLUDES) -o game_of_life_with_output

dc2vtk: dc2vtk.cpp ../dccrg_cartesian_geometry.hpp Makefile
	$(MPICXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

d: data
data:
	rm -f *.dc *.vtk

c: clean
clean: data
	rm -f simple_game_of_life game_of_life game_of_life_with_output dc2vtk
