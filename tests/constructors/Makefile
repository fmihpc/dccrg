FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: simple

simple: simple.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) simple.cpp $(INCLUDES) -o simple

c: clean
clean:
	rm -f simple

