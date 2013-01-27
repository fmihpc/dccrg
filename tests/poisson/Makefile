# use this many mpi processes in parallel tests
PROCS := 6

FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I../.. -I../advection -I$$HOME/include -L$$HOME/lib
LIBS = -lboost_program_options -lzoltan
BOOST_MPI_LIBS = -lboost_serialization -lboost_mpi
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)
MPIRUN = mpirun -np 2
RUN =

PROGRAMS = \
	reference_poisson_example \
	reference_poisson_test \
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

reference_poisson_test: reference_poisson_test.cpp reference_poisson_solve.hpp Makefile
	@echo "CXX $@" && $(CXX) $(FLAGS) $(INCLUDES) $< $(LIBS) -o $@

t: test
test: reference_poisson_test poisson1d
	@echo -n "Testing reference poisson solver... " && ./reference_poisson_test
	@echo -n "Testing poisson solver using 1 process... " && ./poisson1d

pt: parallel-test
parallel-test: poisson1d
	@echo -n "Testing poisson solver using $(PROCS) processes... " && mpirun -n $(PROCS) ./poisson1d

c: clean
clean:
	rm -f $(PROGRAMS)

