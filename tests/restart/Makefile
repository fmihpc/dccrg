# use this many mpi processes in parallel tests
PROCS := 6

# boost causes too many warnings at least on gcc-4.8.2 without -Wno-...
FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x -Wno-unused-local-typedefs
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
MPICXX = mpic++ $(FLAGS)
CXX = g++ $(FLAGS)

PROGRAMS = restart_test restart_test2 variable_cell_data dc2vtk

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_cartesian_geometry.hpp \
	../../dccrg_mapping.hpp \
	../../dccrg_mpi_support.hpp \
	../../dccrg_topology.hpp \
	../../dccrg_types.hpp \
	initialize.hpp \
	refine.hpp \
	solve.hpp \
	cell.hpp \
	IO.hpp

all: $(PROGRAMS)

restart_test: restart_test.cpp $(HEADERS) Makefile
	@echo "MPICXX restart_test " && $(MPICXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS restart_test.cpp $(INCLUDES) -o restart_test

restart_test2: restart_test2.cpp $(HEADERS) Makefile
	@echo "MPICXX restart_test2 " && $(MPICXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS restart_test2.cpp $(INCLUDES) -o restart_test2

variable_cell_data: variable_cell_data.cpp $(HEADERS) Makefile
	@echo "MPICXX variable_cell_data " && $(MPICXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS variable_cell_data.cpp $(INCLUDES) -o variable_cell_data

dc2vtk: dc2vtk.cpp ../../dccrg_cartesian_geometry.hpp Makefile
	@echo "MPICXX dc2vtk " && $(MPICXX) dc2vtk.cpp $(INCLUDES) -o dc2vtk

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
test: restart_test restart_test2 variable_cell_data
	@echo "Running serial restart tests..."
	@echo -n "GoL1 without restart: " && ./restart_test
	@$(foreach dir, $(DIRS), mkdir -p $(dir);)
	@$(foreach dir, $(DIRS), echo -n "Restarting from step $(dir): " && cd $(dir) && ../restart_test --restart ../gol_$(dir).dc && cd ..;)
	@$(foreach dir, $(DIRS), rm -rf $(dir);)
	@echo -n "GoL2 without restart: " && ./restart_test2
	@$(foreach dir, $(DIRS), mkdir -p $(dir);)
	@$(foreach dir, $(DIRS), echo -n "Restarting from step $(dir): " && cd $(dir) && ../restart_test2 --restart ../gol_$(dir).dc && cd ..;)
	@echo -n "Variable cell data without restart: " && ./variable_cell_data
	@echo -n "Variable cell data with restart: " && ./variable_cell_data --restart variable_cell_data.dc

pt: parallel-test
parallel-test: restart_test restart_test2 variable_cell_data
	@echo "Running parallel restart tests using $(PROCS) processes..."
	@echo -n "GoL1 without restart: " && mpirun -np 2 ./restart_test
	@$(foreach dir, $(DIRS), mkdir -p $(dir);)
	@$(foreach dir, $(DIRS), echo -n "Restarting from step $(dir): " && cd $(dir) && mpirun -np $(PROCS) ../restart_test --restart ../gol_$(dir).dc && cd ..;)
	@$(foreach dir, $(DIRS), rm -rf $(dir);)
	@echo -n "GoL2 without restart: " && mpirun -np 2 ./restart_test2
	@$(foreach dir, $(DIRS), mkdir -p $(dir);)
	@$(foreach dir, $(DIRS), echo -n "Restarting from step $(dir): " && cd $(dir) && mpirun -np $(PROCS) ../restart_test2 --restart ../gol_$(dir).dc && cd ..;)
	@echo -n "Variable cell data without restart: " && mpirun -n $(PROCS) ./variable_cell_data
	@echo -n "Variable cell data with restart: " && mpirun -n $(PROCS) ./variable_cell_data --restart variable_cell_data.dc

# convert results to vtk files
v: vtk
vtk: dc2vtk
	@./dc2vtk gol*.dc
	@$(foreach dir, $(DIRS), cd $(dir) && ../dc2vtk gol*.dc && cd ..;)

