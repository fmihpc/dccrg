# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_cartesian_geometry.hpp \
	../../dccrg_length.hpp \
	../../dccrg_mapping.hpp \
	../../dccrg_stretched_cartesian_geometry.hpp \
	../../dccrg_types.hpp \
	../../dccrg_topology.hpp

all: scalability

scalability: scalability.cpp $(HEADERS) Makefile
	$(CXX) scalability.cpp $(INCLUDES) -o scalability

c: clean
clean:
	rm -f scalability
