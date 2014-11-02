FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)

all: init

init: init.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) init.cpp $(INCLUDES) -o init

c: clean
clean:
	rm -f init
