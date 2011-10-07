# Use -Wno-missing-braces until GCC stops warning about valid C++ (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=25137)
FLAGS = -O3 -W -Wall -Wextra -pedantic -Wno-missing-braces -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

PROGRAMS =						\
	game_of_life_test_all_cells_boost_mpi		\
	game_of_life_test_all_cells_mpi			\
	game_of_life_test_all_cells_mpi_datatype	\
	game_of_life_test_boost_mpi			\
	game_of_life_test_mpi				\
	game_of_life_test_mpi_datatype			\
	scalability					\
	scalability3d					\
	refined						\
	refined2d					\
	unrefined2d					\
	refined_scalability3d				\
	scalability1d					\
	hierarchial_test				\
	pinned_cells

HEADERS = \
	../../dccrg.hpp \
	../../dccrg_index.hpp \
	../../dccrg_geometry.hpp

all: $(PROGRAMS)

game_of_life_test_all_cells_boost_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_all_cells_boost_mpi

game_of_life_test_all_cells_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_CELL_DATA_SIZE_FROM_USER game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_all_cells_mpi

game_of_life_test_all_cells_mpi_datatype: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_CELL_DATA_SIZE_FROM_USER -DDCCRG_USER_MPI_DATA_TYPE game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_all_cells_mpi_datatype

game_of_life_test_boost_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_boost_mpi

game_of_life_test_mpi: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_mpi

game_of_life_test_mpi_datatype: game_of_life_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER -DDCCRG_USER_MPI_DATA_TYPE game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_mpi_datatype

scalability: scalability.cpp $(HEADERS) Makefile
	$(CXX) scalability.cpp $(INCLUDES) -o scalability

scalability3d: scalability3d.cpp $(HEADERS) Makefile
	$(CXX) scalability3d.cpp $(INCLUDES) -o scalability3d

scalability1d: scalability1d.cpp $(HEADERS) Makefile
	$(CXX) scalability1d.cpp $(INCLUDES) -o scalability1d

refined: refined.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG refined.cpp $(INCLUDES) -o refined

refined2d: refined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG refined2d.cpp $(INCLUDES) -o refined2d

unrefined2d: unrefined2d.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG unrefined2d.cpp $(INCLUDES) -o unrefined2d

refined_scalability3d: refined_scalability3d.cpp $(HEADERS) Makefile
	$(CXX) refined_scalability3d.cpp $(INCLUDES) -o refined_scalability3d

hierarchial_test: hierarchial_test.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG hierarchial_test.cpp $(INCLUDES) -o hierarchial_test

pinned_cells: pinned_cells.cpp $(HEADERS) Makefile
	$(CXX) -DDEBUG pinned_cells.cpp $(INCLUDES) -o pinned_cells

c: clean
clean:
	rm -f $(PROGRAMS) *.vtk *visit

