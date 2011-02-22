FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
#FLAGS = -O0 -g -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: load_balancing_test hierarchial

load_balancing_test: load_balancing_test.cpp ../../dccrg.hpp Makefile
	$(CXX) load_balancing_test.cpp $(INCLUDES) -o load_balancing_test

hierarchial: hierarchial.cpp ../../dccrg.hpp Makefile
	$(CXX) hierarchial.cpp $(INCLUDES) -o hierarchial

c: clean
clean:
	rm -f load_balancing_test hierarchial *vtk *visit
