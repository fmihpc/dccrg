# use -fno-builtin-sprintf to avoid warning about gnu specific formatting in sprintf by boost lexical_cast
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -fno-builtin-sprintf
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = neighbor_list_length game_of_life advection

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_constant_geometry.hpp

all: $(PROGRAMS)

neighbor_list_length: neighbor_list_length.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG neighbor_list_length.cpp $(INCLUDES) -o neighbor_list_length

game_of_life: game_of_life.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG game_of_life.cpp $(INCLUDES) -o game_of_life

advection: advection.cpp $(HEADERS) ../advection/cell.hpp ../advection/save.hpp Makefile
	$(CXX) -DDEBUG advection.cpp $(INCLUDES) -o advection

d: data
data:
	rm -f *.dc *.vtk *.visit

c: clean
clean: data
	rm -f $(PROGRAMS)

