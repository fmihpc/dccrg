FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lzoltan
CXX = mpic++ $(FLAGS)

all: variable_data_size variable_neighbour_data

variable_data_size: variable_data_size.cpp ../../dccrg.hpp ../../dccrg_stretched_cartesian_geometry.hpp Makefile
	$(CXX) -DDEBUG variable_data_size.cpp $(INCLUDES) -o variable_data_size

variable_neighbour_data: variable_neighbour_data.cpp ../../dccrg.hpp ../../dccrg_stretched_cartesian_geometry.hpp Makefile
	$(CXX) -DDEBUG variable_neighbour_data.cpp $(INCLUDES) -o variable_neighbour_data

c: clean
clean:
	rm -f variable_data_size variable_neighbour_data
