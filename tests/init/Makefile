FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: init

init: init.cpp ../../dccrg.hpp ../../dccrg_arbitrary_geometry.hpp Makefile
	$(CXX) init.cpp $(INCLUDES) -o init

c: clean
clean:
	rm -f init
