FLAGS = -O3 -W -Wall -Wextra -Wno-missing-braces -pedantic -std=c++0x
#FLAGS = -O0 -g -W -Wall -Wextra -pedantic -std=c++0x -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
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
	rm -f refine_simple unrefine_simple scalability
