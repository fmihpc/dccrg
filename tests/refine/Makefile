FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I../.. -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)

all: refine_simple unrefine_simple scalability

refine_simple: refine_simple.cpp ../../dccrg.hpp Makefile
	$(CXX) refine_simple.cpp $(INCLUDES) -o refine_simple

unrefine_simple: unrefine_simple.cpp ../../dccrg.hpp Makefile
	$(CXX) unrefine_simple.cpp $(INCLUDES) -o unrefine_simple

scalability: scalability.cpp ../../dccrg.hpp Makefile
	$(CXX) scalability.cpp $(INCLUDES) -o scalability

c: clean
clean:
	rm -f refine_simple unrefine_simple scalability *.visit *.vtk
