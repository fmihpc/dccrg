EXTRA_WARNINGS = -ansi -Winit-self -Wswitch-enum -Wunused-parameter -Wformat=2 -Wshadow -Wconversion -Wcast-align -Wsign-conversion -Wlogical-op -Wmissing-declarations -Wpacked -Wpadded -Wredundant-decls -Wunreachable-code -Winline -Winvalid-pch -Wdisabled-optimization
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x #$(EXTRA_WARNINGS)
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)

HEADERS = ../dccrg.hpp ../dccrg_constant_geometry.hpp

all: simple_game_of_life game_of_life game_of_life_with_output dc2vtk

game_of_life: game_of_life.cpp $(HEADERS) Makefile
	$(MPICXX) game_of_life.cpp $(INCLUDES) -o game_of_life

simple_game_of_life: simple_game_of_life.cpp $(HEADERS) Makefile
	$(MPICXX) simple_game_of_life.cpp $(INCLUDES) -o simple_game_of_life

game_of_life_with_output: game_of_life_with_output.cpp $(HEADERS) Makefile
	$(MPICXX) game_of_life_with_output.cpp $(INCLUDES) -o game_of_life_with_output

dc2vtk: dc2vtk.cpp ../dccrg_constant_geometry.hpp Makefile
	$(CXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

d: data
data:
	rm -f *.dc *.vtk

c: clean
clean: data
	rm -f simple_game_of_life game_of_life game_of_life_with_output dc2vtk
