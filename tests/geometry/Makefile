FLAGS = -O3 -W -Wall -Wextra -pedantic -DNDEBUG
#FLAGS = -O0 -pg -g -W -Wall -Wextra -pedantic
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: constant_grid_speed_1d arbitrary_stretch_speed_1d

constant_grid_speed_1d: constant_grid_speed_1d.cpp ../../dccrg_cell_geometry.hpp Makefile
	$(CXX) constant_grid_speed_1d.cpp $(INCLUDES) -o constant_grid_speed_1d

arbitrary_stretch_speed_1d: arbitrary_stretch_speed_1d.cpp ../../dccrg_cell_geometry.hpp Makefile
	$(CXX) arbitrary_stretch_speed_1d.cpp $(INCLUDES) -o arbitrary_stretch_speed_1d

c: clean
clean:
	rm -f constant_grid_speed_1d arbitrary_stretch_speed_1d
