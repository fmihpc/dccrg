# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np 2
RUN =

PROGRAMS = simple copy
TEST_PROGRAMS = copy

all: $(PROGRAMS)

simple: simple.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -o $@

copy: copy.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -o $@

c: clean
clean:
	rm -f $(PROGRAMS)

t: test
test: $(TEST_PROGRAMS)
	@for i in $(TEST_PROGRAMS); do echo -n "RUN $$i... " && $(RUN) ./$$i; done
	@for i in $(TEST_PROGRAMS); do echo -n "MPIRUN $$i... " && $(MPIRUN) $$i; done

