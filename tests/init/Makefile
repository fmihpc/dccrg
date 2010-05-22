#FLAGS = -O3 -W -Wall -Wextra -pedantic -DNDEBUG
FLAGS = -O0 -g -W -Wall -Wextra -pedantic
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: init

init: init.cpp ../../dccrg.hpp Makefile
	$(CXX) init.cpp $(INCLUDES) -o init

c: clean
clean:
	rm -f init
