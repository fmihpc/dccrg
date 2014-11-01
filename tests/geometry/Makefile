FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lzoltan
MPICXX = mpic++ $(FLAGS)

all: cartesian_grid_speed stretched_cartesian_grid_speed

cartesian_grid_speed: cartesian_grid_speed.cpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(MPICXX) cartesian_grid_speed.cpp $(INCLUDES) -o cartesian_grid_speed

stretched_cartesian_grid_speed: stretched_cartesian_grid_speed.cpp ../../dccrg_stretched_cartesian_geometry.hpp Makefile
	$(MPICXX) stretched_cartesian_grid_speed.cpp $(INCLUDES) -o stretched_cartesian_grid_speed

c: clean
clean:
	rm -f cartesian_grid_speed stretched_cartesian_grid_speed
