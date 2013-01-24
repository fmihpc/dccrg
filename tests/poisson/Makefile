FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I.. -I../advection -I$$HOME/include -L$$HOME/lib
LIBS = -lboost_program_options -lzoltan
BOOST_MPI_LIBS = -lboost_serialization -lboost_mpi
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)
MPIRUN = mpirun -np 2
RUN =

PROGRAMS = \
	reference_poisson_example \
	poisson1d

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_stretched_cartesian_geometry.hpp \
	../../dccrg_cartesian_geometry.hpp \
	../advection/grid_support.hpp \
	poisson_solve.hpp \
	reference_poisson_solve.hpp

all: $(PROGRAMS)

poisson1d: poisson1d.cpp $(SOURCES) $(HEADERS) Makefile
	@echo "MPICXX  $@" && $(MPICXX) $(FLAGS) $(INCLUDES) $< $(LIBS) $(BOOST_MPI_LIBS) -o $@

reference_poisson_example: reference_poisson_example.cpp reference_poisson_solve.hpp Makefile
	@echo "CXX $@" && $(CXX) $(FLAGS) $(INCLUDES) $< $(LIBS) -o $@

c: clean
clean:
	rm -f $(PROGRAMS)

