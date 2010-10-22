EXTRA_WARNINGS = -ansi -Winit-self -Wswitch-enum -Wunused-parameter  -Wformat=2  -Wshadow -Wconversion -Wcast-align -Wconversion -Wsign-conversion -Wlogical-op -Wmissing-declarations -Wpacked -Wpadded -Wredundant-decls -Wunreachable-code -Winline -Winvalid-pch -Wdisabled-optimization -Wfatal-errors 
FLAGS = -O3 -W -Wall -Wextra -pedantic -DNDEBUG -std=c++0x #$(EXTRA_WARNINGS)
#FLAGS = -O0 -W -Wall -Wextra -pedantic
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: simple_game_of_life game_of_life

game_of_life: game_of_life.cpp ../dccrg.hpp Makefile
	$(CXX) game_of_life.cpp $(INCLUDES) -o game_of_life

simple_game_of_life: simple_game_of_life.cpp ../dccrg.hpp Makefile
	$(CXX) simple_game_of_life.cpp $(INCLUDES) -o simple_game_of_life

c: clean
clean:
	rm -f simple_game_of_life game_of_life
