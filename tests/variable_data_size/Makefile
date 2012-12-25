FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: variable_data_size variable_neighbour_data

variable_data_size: variable_data_size.cpp ../../dccrg.hpp ../../dccrg_arbitrary_geometry.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI variable_data_size.cpp $(INCLUDES) -o variable_data_size

variable_neighbour_data: variable_neighbour_data.cpp ../../dccrg.hpp ../../dccrg_arbitrary_geometry.hpp Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI variable_neighbour_data.cpp $(INCLUDES) -o variable_neighbour_data

c: clean
clean:
	rm -f variable_data_size variable_neighbour_data
