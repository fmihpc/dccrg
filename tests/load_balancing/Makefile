FLAGS = -W -Wall -Wextra -pedantic -std=c++11 -O3
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS = load_balancing_test multi_stage_load_balancing

all: $(PROGRAMS)

load_balancing_test: load_balancing_test.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -DDEBUG -o $@

multi_stage_load_balancing: multi_stage_load_balancing.cpp ../../dccrg.hpp ../../dccrg_cartesian_geometry.hpp Makefile
	$(CXX) $< $(INCLUDES) -DDEBUG -o $@

d: data
data:
	rm -f *vtk *visit

c: clean
clean: data
	rm -f $(PROGRAMS)

