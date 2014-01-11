# use this many mpi processes in parallel tests
PROCS := 6

# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs
INCLUDES = -I../.. -I../advection -I$$HOME/include -L$$HOME/lib
LIBS = -lboost_program_options -lzoltan
BOOST_MPI_LIBS = -lboost_serialization -lboost_mpi
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)
MPIRUN = mpirun -n $(PROCS)
RUN =

PROGRAMS = \
	reference_poisson_example \
	reference_poisson_test \
	poisson1d \
	poisson1d_amr \
	poisson1d_skip_cells \
	poisson1d_stretched \
	poisson2d

TEST_PROGRAMS = \
	reference_poisson_test \
	poisson1d \
	poisson1d_amr \
	poisson1d_skip_cells \
	poisson1d_stretched \
	poisson2d

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_cartesian_geometry.hpp \
	../../dccrg_length.hpp \
	../../dccrg_mapping.hpp \
	../../dccrg_stretched_cartesian_geometry.hpp \
	../../dccrg_topology.hpp \
	../../dccrg_types.hpp \
	poisson_solve.hpp \
	reference_poisson_solve.hpp

all: $(PROGRAMS)

poisson1d_amr: poisson1d_amr.cpp $(HEADERS) Makefile
	@echo "MPICXX $@" && $(MPICXX) $(FLAGS) $(INCLUDES) $< $(LIBS) $(BOOST_MPI_LIBS) -o $@

poisson1d_skip_cells: poisson1d_skip_cells.cpp $(HEADERS) Makefile
	@echo "MPICXX $@" && $(MPICXX) $(FLAGS) $(INCLUDES) $< $(LIBS) $(BOOST_MPI_LIBS) -o $@

poisson1d_stretched: poisson1d_stretched.cpp $(HEADERS) Makefile
	@echo "MPICXX $@" && $(MPICXX) $(FLAGS) $(INCLUDES) $< $(LIBS) $(BOOST_MPI_LIBS) -o $@

poisson1d: poisson1d.cpp $(HEADERS) Makefile
	@echo "MPICXX $@" && $(MPICXX) $(FLAGS) $(INCLUDES) $< $(LIBS) $(BOOST_MPI_LIBS) -o $@

poisson2d: poisson2d.cpp $(HEADERS) Makefile
	@echo "MPICXX $@" && $(MPICXX) $(FLAGS) $(INCLUDES) $< $(LIBS) $(BOOST_MPI_LIBS) -o $@

reference_poisson_example: reference_poisson_example.cpp reference_poisson_solve.hpp Makefile
	@echo "CXX $@" && $(CXX) $(FLAGS) $(INCLUDES) $< $(LIBS) -o $@

reference_poisson_test: reference_poisson_test.cpp reference_poisson_solve.hpp Makefile
	@echo "CXX $@" && $(CXX) $(FLAGS) $(INCLUDES) $< $(LIBS) -o $@

plots:
	@echo "Plotting 1d stretched poisson test results with gnuplot... " && for f in *.gnuplot; do gnuplot $$f; done

t: test
test: $(TEST_PROGRAMS)
	@echo -n "Testing reference poisson solver... " && ./reference_poisson_test
	@echo -n "Testing 1d poisson solver in serial... " && ./poisson1d
	@echo -n "Testing 1d poisson solver on a refined grid in serial... " && ./poisson1d_amr
	@echo -n "Testing 1d poisson solver with skipped cells in serial... " && ./poisson1d_skip_cells
	@echo -n "Testing 1d poisson solver on a stretched grid in serial... " && ./poisson1d_stretched
	@echo -n "Testing 2d poisson solver in serial... " && ./poisson2d

pt: parallel-test
parallel-test: $(TEST_PROGRAMS)
	@echo -n "Testing 1d poisson solver using $(PROCS) processes... " && $(MPIRUN) ./poisson1d
	@echo -n "Testing 1d amr poisson solver using $(PROCS) processes... " && $(MPIRUN) ./poisson1d_amr
	@echo -n "Testing 1d poisson solver with skipped cells using $(PROCS) processes... " && $(MPIRUN) ./poisson1d_skip_cells
	@echo -n "Testing 1d stretched poisson solver using $(PROCS) processes... " && $(MPIRUN) ./poisson1d_stretched
	@echo -n "Testing 2d poisson solver using $(PROCS) processes... " && $(MPIRUN) ./poisson2d

d: data
data:
	rm -f poisson1d_stretched*gnuplot poisson1d_stretched*svg

c: clean
clean: data
	rm -f $(PROGRAMS)

