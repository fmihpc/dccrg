FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = simple

HEADERS = \
	cell.hpp \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_cartesian_geometry.hpp

all: $(PROGRAMS)

simple: simple.cpp $(HEADERS) Makefile
	$(CXX) simple.cpp $(INCLUDES) -o simple

d: data
data:
	rm -f *.dc *.vtk *.visit

c: clean
clean: data
	rm -f $(PROGRAMS)

