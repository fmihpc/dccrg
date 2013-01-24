FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -DDEBUG
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = load_balancing_test multi_stage_load_balancing

all: $(PROGRAMS)

load_balancing_test: load_balancing_test.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< -DDCCRG_TRANSFER_USING_BOOST_MPI $(INCLUDES) -o $@

multi_stage_load_balancing: multi_stage_load_balancing.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -o $@

d: data
data:
	rm -f *vtk *visit

c: clean
clean: data
	rm -f $(PROGRAMS)

