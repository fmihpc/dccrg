FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
BOOST_MPI_LIBS = -lboost_serialization -lboost_mpi
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_program_options -lzoltan
CXX = mpic++ $(FLAGS)
MPIRUN = mpirun -np 2
RUN =

PROGRAMS = \
	game_of_life_test_all_cells_boost_mpi \
	game_of_life_test_all_cells_mpi \
	game_of_life_test_boost_mpi \
	game_of_life_test_mpi \
	unrefined2d_all_cells_boost_mpi \
	unrefined2d_all_cells_mpi \
	unrefined2d_boost_mpi \
	unrefined2d_mpi \
	scalability \
	scalability3d \
	refined \
	refined2d \
	refined_scalability3d \
	scalability1d \
	hierarchical_test \
	pinned_cells

TEST_PROGRAMS = \
	game_of_life_test_all_cells_boost_mpi \
	game_of_life_test_all_cells_mpi \
	game_of_life_test_boost_mpi \
	game_of_life_test_mpi \
	unrefined2d_all_cells_boost_mpi \
	unrefined2d_all_cells_mpi \
	unrefined2d_boost_mpi \
	unrefined2d_mpi

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_arbitrary_geometry.hpp \
	../../dccrg_constant_geometry.hpp \
	cell.hpp initialize.hpp save.hpp solve.hpp

all: $(PROGRAMS)

game_of_life_test_all_cells_boost_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI game_of_life_test.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o game_of_life_test_all_cells_boost_mpi

game_of_life_test_all_cells_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_all_cells_mpi

game_of_life_test_boost_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI -DDCCRG_SEND_SINGLE_CELLS game_of_life_test.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o game_of_life_test_boost_mpi

game_of_life_test_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_mpi

scalability: scalability.cpp $(HEADERS) Makefile
	$(CXX) scalability.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o scalability

scalability3d: scalability3d.cpp $(HEADERS) Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI scalability3d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o scalability3d

scalability1d: scalability1d.cpp $(HEADERS) Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI scalability1d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o scalability1d

refined: refined.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI refined.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o refined

refined2d: refined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI refined2d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o refined2d

unrefined2d_all_cells_boost_mpi: unrefined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI unrefined2d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o unrefined2d_all_cells_boost_mpi

unrefined2d_all_cells_mpi: unrefined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG unrefined2d.cpp $(INCLUDES) -o unrefined2d_all_cells_mpi

unrefined2d_boost_mpi: unrefined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI -DDCCRG_SEND_SINGLE_CELLS unrefined2d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o unrefined2d_boost_mpi

unrefined2d_mpi: unrefined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS unrefined2d.cpp $(INCLUDES) -o unrefined2d_mpi

refined_scalability3d: refined_scalability3d.cpp $(HEADERS) Makefile
	$(CXX) -DDCCRG_TRANSFER_USING_BOOST_MPI refined_scalability3d.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o refined_scalability3d

hierarchical_test: hierarchical_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI hierarchical_test.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o hierarchical_test

pinned_cells: pinned_cells.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_TRANSFER_USING_BOOST_MPI pinned_cells.cpp $(INCLUDES) $(BOOST_MPI_LIBS) -o pinned_cells

d: data
data:
	rm -f *.vtk *.visit

c: clean
clean: data
	rm -f $(PROGRAMS)

t: test
test: $(TEST_PROGRAMS)
	@for i in $(TEST_PROGRAMS); do echo -n "RUN $$i... " && $(RUN) ./$$i; done
	@for i in $(TEST_PROGRAMS); do echo -n "MPIRUN $$i... " && $(MPIRUN) $$i; done
