# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = simple

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_cartesian_geometry.hpp \
	../../dccrg_length.hpp \
	../../dccrg_mapping.hpp \
	../../dccrg_topology.hpp \
	../../dccrg_types.hpp \
	cell.hpp \

all: $(PROGRAMS)

simple: simple.cpp $(HEADERS) Makefile
	$(CXX) simple.cpp $(INCLUDES) -o simple

d: data
data:
	rm -f *.dc *.vtk *.visit

c: clean
clean: data
	rm -f $(PROGRAMS)

