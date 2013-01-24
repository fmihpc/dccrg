FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = g++ $(FLAGS)

all: cartesian_grid_speed stretched_cartesian_grid_speed

cartesian_grid_speed: cartesian_grid_speed.cpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) cartesian_grid_speed.cpp $(INCLUDES) -o cartesian_grid_speed

stretched_cartesian_grid_speed: stretched_cartesian_grid_speed.cpp ../../dccrg_stretched_cartesian_geometry.hpp Makefile
	$(CXX) stretched_cartesian_grid_speed.cpp $(INCLUDES) -o stretched_cartesian_grid_speed

c: clean
clean:
	rm -f cartesian_grid_speed stretched_cartesian_grid_speed
