FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lzoltan
MPICXX = mpic++ $(FLAGS)

HEADERS = ../dccrg.hpp ../dccrg_cartesian_geometry.hpp

all: simple_game_of_life game_of_life game_of_life_with_output dc2vtk

game_of_life: game_of_life.cpp $(HEADERS) Makefile
	$(MPICXX) game_of_life.cpp $(INCLUDES) -o game_of_life

simple_game_of_life: simple_game_of_life.cpp $(HEADERS) Makefile
	$(MPICXX) simple_game_of_life.cpp $(INCLUDES) -o simple_game_of_life

game_of_life_with_output: game_of_life_with_output.cpp $(HEADERS) Makefile
	$(MPICXX) game_of_life_with_output.cpp $(INCLUDES) -o game_of_life_with_output

dc2vtk: dc2vtk.cpp ../dccrg_cartesian_geometry.hpp Makefile
	$(MPICXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

d: data
data:
	rm -f *.dc *.vtk

c: clean
clean: data
	rm -f simple_game_of_life game_of_life game_of_life_with_output dc2vtk
