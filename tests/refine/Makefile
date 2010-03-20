#FLAGS = -O3 -W -Wall -pedantic -DNDEBUG
FLAGS = -O0 -g -W -Wall -pedantic
INCLUDES = -I. -I$$HOME/include -L. -L$$HOME/lib -lboost_mpi-mt -lboost_serialization-mt -lzoltan
CXX = mpic++ $(FLAGS)

all: refine_simple

refine_simple: refine_simple.cpp ../../dccrg.hpp Makefile
	$(CXX) refine_simple.cpp $(INCLUDES) -o refine_simple

c: clean
clean:
	rm -f refine_simple
