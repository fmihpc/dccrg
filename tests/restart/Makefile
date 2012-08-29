# use this many mpi processes in parallel tests
PROCS := 2

FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_program_options -lboost_serialization -lzoltan
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)

PROGRAMS = restart_test dc2vtk

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_constant_geometry.hpp \
	initialize.hpp \
	refine.hpp \
	solve.hpp \
	cell.hpp \
	IO.hpp

all: $(PROGRAMS)

restart_test: restart_test.cpp $(HEADERS) Makefile
	@echo "MPICXX restart_test " && $(MPICXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER -DDCCRG_USER_MPI_DATA_TYPE restart_test.cpp $(INCLUDES) -o restart_test

dc2vtk: dc2vtk.cpp ../../dccrg_constant_geometry.hpp Makefile
	@echo "CXX dc2vtk " && $(CXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

# test restart from there directories
DIRS := 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23

# remove generated data
d: data
data:
	@rm -f *.dc *.vtk
	@$(foreach dir, $(DIRS), rm -rf $(dir);)

# remove executables
c: clean
clean: data
	@rm -f $(PROGRAMS)

t: test
test: restart_test
	@echo "Running serial restart tests..."
	@echo -n "GoL without restart: " && ./restart_test
	@$(foreach dir, $(DIRS), mkdir -p $(dir);)
	@$(foreach dir, $(DIRS), echo -n "Restarting from step $(dir): " && cd $(dir) && ../restart_test --restart ../gol_$(dir).dc && cd ..;)

pt: parallel-test
parallel-test: restart_test
	@echo "Running parallel restart tests using $(PROCS) processes..."
	@echo -n "GoL without restart: " && mpirun -np 2 ./restart_test
	@$(foreach dir, $(DIRS), mkdir -p $(dir);)
	@$(foreach dir, $(DIRS), echo -n "Restarting from step $(dir): " && cd $(dir) && mpirun -np $(PROCS) ../restart_test --restart ../gol_$(dir).dc && cd ..;)

# convert results to vtk files
v: vtk
vtk: dc2vtk
	@./dc2vtk gol*.dc
	@$(foreach dir, $(DIRS), cd $(dir) && ../dc2vtk gol*.dc && cd ..;)
