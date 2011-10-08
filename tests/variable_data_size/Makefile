# Use -Wno-missing-braces until GCC stops warning about valid C++ (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=25137)
FLAGS = -O3 -W -Wall -Wextra -Wno-missing-braces -pedantic -std=c++0x
#FLAGS = -O0 -g -W -Wall -Wextra -pedantic -std=c++0x -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: variable_data_size variable_neighbour_data

variable_data_size: variable_data_size.cpp ../../dccrg.hpp ../../dccrg_arbitrary_geometry.hpp Makefile
	$(CXX) variable_data_size.cpp $(INCLUDES) -o variable_data_size

variable_neighbour_data: variable_neighbour_data.cpp ../../dccrg.hpp ../../dccrg_arbitrary_geometry.hpp Makefile
	$(CXX) variable_neighbour_data.cpp $(INCLUDES) -o variable_neighbour_data

c: clean
clean:
	rm -f variable_data_size variable_neighbour_data
