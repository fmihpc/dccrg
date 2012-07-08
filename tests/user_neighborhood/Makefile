FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = neighbor_list_length

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_constant_geometry.hpp

all: $(PROGRAMS)

neighbor_list_length: neighbor_list_length.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG neighbor_list_length.cpp $(INCLUDES) -o neighbor_list_length

c: clean
clean:
	rm -f $(PROGRAMS)

