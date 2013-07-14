FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lboost_serialization -lboost_mpi -lzoltan
MPICXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np 5
RUN = mpirun -np 1

PROGRAMS = test1

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_cartesian_geometry.hpp \
	../../dccrg_length.hpp \
	../../dccrg_mapping.hpp \
	../../dccrg_topology.hpp \
	../../dccrg_types.hpp

all: $(PROGRAMS)

test1: test1.cpp $(HEADERS) Makefile
	@echo "MPICXX $@" && $(MPICXX) -DDEBUG test1.cpp $(INCLUDES) -o test1

c: clean
clean:
	rm -f $(PROGRAMS)

t: test
test: $(PROGRAMS)
	@for i in $(PROGRAMS); do echo -n "RUN $$i... " && $(RUN) ./$$i; done
	@for i in $(PROGRAMS); do echo -n "MPIRUN $$i... " && $(MPIRUN) ./$$i; done

