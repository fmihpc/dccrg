FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib
MPICXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np

all: some_reduce

some_reduce: some_reduce.cpp ../../dccrg_mpi_support.hpp Makefile
	$(MPICXX) some_reduce.cpp $(INCLUDES) -o some_reduce

t: test
test: some_reduce
	@for i in 1 2 3 4 5 6 7 8 9 10; do echo -n "MPIRUN some_reduce with $$i processes... " && $(MPIRUN) $$i ./some_reduce; done

c: clean
clean:
	rm -f some_reduce

