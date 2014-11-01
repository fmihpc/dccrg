FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np 2
RUN =

PROGRAMS = simple copy
TEST_PROGRAMS = copy

all: $(PROGRAMS)

simple: simple.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -DDEBUG -o $@

copy: copy.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -DDEBUG -o $@

c: clean
clean:
	rm -f $(PROGRAMS)

t: test
test: $(TEST_PROGRAMS)
	@for i in $(TEST_PROGRAMS); do echo -n "RUN $$i... " && $(RUN) ./$$i; done
	@for i in $(TEST_PROGRAMS); do echo -n "MPIRUN $$i... " && $(MPIRUN) $$i; done

