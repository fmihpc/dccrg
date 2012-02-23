FLAGS = -O3 -W -Wall -Wextra -pedantic
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: load_balancing_test

load_balancing_test: load_balancing_test.cpp ../../dccrg.hpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) load_balancing_test.cpp $(INCLUDES) -o load_balancing_test

c: clean
clean:
	rm -f load_balancing_test *vtk *visit
