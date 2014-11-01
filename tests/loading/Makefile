FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)

all: load

load: load.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) -DDEBUG load.cpp $(INCLUDES) -o load

d: data
data:
	rm -rf load.visit load_*vtk

c: clean
clean: data
	rm -f load
