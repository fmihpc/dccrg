FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
BOOST_MPI_LIBS = -lboost_serialization -lboost_mpi
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)
MPIRUN = mpirun -np 2
RUN =

PROGRAMS = reference_poisson_example

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_arbitrary_geometry.hpp \
	../../dccrg_constant_geometry.hpp \
	reference_poisson_solve.hpp

all: $(PROGRAMS)

reference_poisson_example: reference_poisson_example.cpp $(HEADERS) Makefile
	@echo "CXX $@" && $(CXX) $(FLAGS) $(INCLUDES) $< -o $@

c: clean
clean:
	rm -f $(PROGRAMS)

