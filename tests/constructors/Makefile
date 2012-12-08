FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = simple copy

all: $(PROGRAMS)

simple: simple.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -o $@

copy: copy.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -o $@

c: clean
clean:
	rm -f $(PROGRAMS)

