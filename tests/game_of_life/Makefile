#FLAGS = -O0 -g -W -Wall -pedantic -fprofile-arcs -ftest-coverage
#FLAGS = -O3 -pg -W -Wall -pedantic
FLAGS = -O3 -W -Wall -pedantic -DNDEBUG
#FLAGS = -O0 -g -W -Wall -pedantic
INCLUDES = -I. -I$$HOME/include -L. -L$$HOME/lib -lboost_mpi-mt -lboost_serialization-mt -lzoltan
CXX = mpic++ $(FLAGS)

all: game_of_life_test scalability scalability3d

game_of_life_test: game_of_life_test.cpp ../../dccrg.hpp Makefile
	$(CXX) game_of_life_test.cpp $(INCLUDES) -o game_of_life_test

scalability: scalability.cpp ../../dccrg.hpp Makefile
	$(CXX) scalability.cpp $(INCLUDES) -o scalability

scalability3d: scalability3d.cpp ../../dccrg.hpp Makefile
	$(CXX) scalability3d.cpp $(INCLUDES) -o scalability3d

c: clean
clean:
	rm -f game_of_life_test scalability scalability3d
