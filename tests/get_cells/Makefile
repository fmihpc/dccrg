FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lzoltan
MPICXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np 3

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
	@echo "MPICXX $@" && $(MPICXX) -DDEBUG  $< $(INCLUDES) -o $@

c: clean
clean:
	rm -f $(PROGRAMS)

t: test
test: $(PROGRAMS)
	@for i in $(PROGRAMS); do echo -n "MPIRUN $$i... " && $(MPIRUN) $$i; done

