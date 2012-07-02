FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = g++ $(FLAGS)

all: constant_grid_speed arbitrary_stretch_speed

constant_grid_speed: constant_grid_speed.cpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) constant_grid_speed.cpp $(INCLUDES) -o constant_grid_speed

arbitrary_stretch_speed: arbitrary_stretch_speed.cpp ../../dccrg_arbitrary_geometry.hpp Makefile
	$(CXX) arbitrary_stretch_speed.cpp $(INCLUDES) -o arbitrary_stretch_speed

c: clean
clean:
	rm -f constant_grid_speed arbitrary_stretch_speed
