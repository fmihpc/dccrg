FLAGS = -g -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I../game_of_life -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)

PROGRAMS = restart_test dc2vtk

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_constant_geometry.hpp \
	../game_of_life/initialize.hpp \
	../game_of_life/solve.hpp \
	cell.hpp \
	IO.hpp

all: $(PROGRAMS)

restart_test: restart_test.cpp $(HEADERS) Makefile
	$(MPICXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER -DDCCRG_USER_MPI_DATA_TYPE restart_test.cpp $(INCLUDES) -o restart_test

dc2vtk: dc2vtk.cpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

d: data
data:
	rm -f *.dc *.vtk

c: clean
clean: data
	rm -f $(PROGRAMS)

