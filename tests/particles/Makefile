FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lzoltan
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
	$(CXX) -DDEBUG simple.cpp $(INCLUDES) -o simple

d: data
data:
	rm -f *.dc *.vtk *.visit

c: clean
clean: data
	rm -f $(PROGRAMS)

