# use -fno-builtin-sprintf to avoid warning about gnu specific formatting in sprintf by boost lexical_cast
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -fno-builtin-sprintf
BOOST_MPI_LIBS = -lboost_serialization -lboost_mpi
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = 2d_boost_mpi 2d_mpi 2d_boost_mpi_debug 2d_mpi_debug dc2vtk

HEADERS = \
	cell.hpp \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_constant_geometry.hpp \
	../../dccrg_mpi_support.hpp

all: $(PROGRAMS)

2d_boost_mpi: 2d.cpp $(HEADERS) Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI 2d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o 2d_boost_mpi

2d_mpi: 2d.cpp $(HEADERS) Makefile
	$(CXX) 2d.cpp $(INCLUDES) -o 2d_mpi

2d_boost_mpi_debug: 2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI 2d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o 2d_boost_mpi_debug

2d_mpi_debug: 2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG 2d.cpp $(INCLUDES) -o 2d_mpi_debug

dc2vtk: dc2vtk.cpp ../../dccrg_constant_geometry.hpp Makefile
	$(CXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

d: data
data:
	rm -f *.dc *.vtk

c: clean
clean: data
	rm -f $(PROGRAMS)

