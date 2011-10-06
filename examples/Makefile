EXTRA_WARNINGS = -ansi -Winit-self -Wswitch-enum -Wunused-parameter  -Wformat=2  -Wshadow -Wconversion -Wcast-align -Wsign-conversion -Wlogical-op -Wmissing-declarations -Wpacked -Wpadded -Wredundant-decls -Wunreachable-code -Winline -Winvalid-pch -Wdisabled-optimization -Wfatal-errors 
# Use -Wno-missing-braces until GCC stops warning about valid C++ (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=25137)
FLAGS = -O3 -W -Wall -Wextra -pedantic -Wno-missing-braces -std=c++0x #$(EXTRA_WARNINGS)
#FLAGS = -O0 -g -W -Wall -Wextra -pedantic -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)

all: simple_game_of_life game_of_life game_of_life_with_output dc2vtk

game_of_life: game_of_life.cpp ../dccrg.hpp Makefile
	$(MPICXX) game_of_life.cpp $(INCLUDES) -o game_of_life

simple_game_of_life: simple_game_of_life.cpp ../dccrg.hpp Makefile
	$(MPICXX) simple_game_of_life.cpp $(INCLUDES) -o simple_game_of_life

game_of_life_with_output: game_of_life_with_output.cpp ../dccrg.hpp Makefile
	$(MPICXX) game_of_life_with_output.cpp $(INCLUDES) -o game_of_life_with_output

dc2vtk: dc2vtk.cpp ../dccrg_geometry.hpp Makefile
	$(CXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

c: clean
clean:
	rm -f simple_game_of_life game_of_life game_of_life_with_output dc2vtk *dc *.vtk
