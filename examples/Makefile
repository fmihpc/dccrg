FLAGS = -O3 -W -Wall -Wextra -pedantic -DNDEBUG
INCLUDES = -I. -I$$HOME/include -L. -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: game_of_life

game_of_life: game_of_life.cpp ../dccrg.hpp Makefile
	$(CXX) game_of_life.cpp $(INCLUDES) -o game_of_life

c: clean
clean:
	rm -f game_of_life
