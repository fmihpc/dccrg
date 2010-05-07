#FLAGS = -O3 -W -Wall -pedantic -DNDEBUG
FLAGS = -O0 -g -W -Wall -pedantic
INCLUDES = -I. -I$$HOME/include -L. -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: init

init: init.cpp ../../dccrg.hpp Makefile
	$(CXX) init.cpp $(INCLUDES) -o init

c: clean
clean:
	rm -f init
