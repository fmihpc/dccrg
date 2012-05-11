FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I../game_of_life -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = restart_test

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_arbitrary_geometry.hpp \
	../game_of_life/cell.hpp \
	../game_of_life/initialize.hpp \
	../game_of_life/solve.hpp

all: $(PROGRAMS)

restart_test: restart_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER -DDCCRG_USER_MPI_DATA_TYPE restart_test.cpp $(INCLUDES) -o restart_test

d: data
data:
	rm -f *.dc

c: clean
clean: data
	rm -f $(PROGRAMS)

