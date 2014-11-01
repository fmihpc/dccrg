FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib
MPICXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np

all: all_gather some_reduce

all_gather: all_gather.cpp ../../dccrg_mpi_support.hpp Makefile
	$(MPICXX) all_gather.cpp $(INCLUDES) -o all_gather

some_reduce: some_reduce.cpp ../../dccrg_mpi_support.hpp Makefile
	$(MPICXX) some_reduce.cpp $(INCLUDES) -o some_reduce

t: test
test: all_gather some_reduce
	@for i in 1 2 3 4 5 6 7 8 9 10; do echo -n "MPIRUN all_gather with $$i processes... " && $(MPIRUN) $$i ./all_gather; done
	@for i in 1 2 3 4 5 6 7 8 9 10; do echo -n "MPIRUN some_reduce with $$i processes... " && $(MPIRUN) $$i ./some_reduce; done

c: clean
clean:
	rm -f all_gather some_reduce

