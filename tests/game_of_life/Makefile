#FLAGS = -O0 -g -W -Wall -Wextra -pedantic -fprofile-arcs -ftest-coverage
#FLAGS = -O3 -pg -W -Wall -Wextra -pedantic
FLAGS = -O3 -W -Wall -Wextra -pedantic -DNDEBUG
#FLAGS = -O0 -g -W -Wall -Wextra -pedantic
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: game_of_life_test scalability scalability3d refined refined2d refined_scalability3d scalability1d

game_of_life_test: game_of_life_test.cpp ../../dccrg.hpp Makefile
	$(CXX) game_of_life_test.cpp $(INCLUDES) -o game_of_life_test

scalability: scalability.cpp ../../dccrg.hpp Makefile
	$(CXX) scalability.cpp $(INCLUDES) -o scalability

scalability3d: scalability3d.cpp ../../dccrg.hpp Makefile
	$(CXX) scalability3d.cpp $(INCLUDES) -o scalability3d

scalability1d: scalability1d.cpp ../../dccrg.hpp Makefile
	$(CXX) scalability1d.cpp $(INCLUDES) -o scalability1d

refined: refined.cpp ../../dccrg.hpp Makefile
	$(CXX) refined.cpp $(INCLUDES) -o refined

refined2d: refined2d.cpp ../../dccrg.hpp Makefile
	$(CXX) refined2d.cpp $(INCLUDES) -o refined2d

refined_scalability3d: refined_scalability3d.cpp ../../dccrg.hpp Makefile
	$(CXX) refined_scalability3d.cpp $(INCLUDES) -o refined_scalability3d

c: clean
clean:
	rm -f game_of_life_test scalability scalability3d refined refined2d refined_scalability3d scalability1d
