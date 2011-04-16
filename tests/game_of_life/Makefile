FLAGS = -O3 -W -Wall -Wextra -pedantic -std=c++0x
INCLUDES = -I$$HOME/include -L$$HOME/lib -lboost_mpi -lboost_serialization -lzoltan
CXX = mpic++ $(FLAGS)

all: game_of_life_test_all_cells game_of_life_test_single_cells game_of_life_user_data_size scalability scalability3d refined refined2d unrefined2d refined_scalability3d scalability1d hierarchial_test pinned_cells

game_of_life_test_all_cells: game_of_life_test.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_all_cells

game_of_life_test_single_cells: game_of_life_test.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS game_of_life_test.cpp $(INCLUDES) -o game_of_life_test_single_cells

game_of_life_user_data_size: game_of_life_user_data_size.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG -DDCCRG_SEND_SINGLE_CELLS -DDCCRG_CELL_DATA_SIZE_FROM_USER game_of_life_user_data_size.cpp $(INCLUDES) -o game_of_life_user_data_size

scalability: scalability.cpp ../../dccrg.hpp Makefile
	$(CXX) -DNDEBUG scalability.cpp $(INCLUDES) -o scalability

scalability3d: scalability3d.cpp ../../dccrg.hpp Makefile
	$(CXX) -DNDEBUG scalability3d.cpp $(INCLUDES) -o scalability3d

scalability1d: scalability1d.cpp ../../dccrg.hpp Makefile
	$(CXX) -DNDEBUG scalability1d.cpp $(INCLUDES) -o scalability1d

refined: refined.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG refined.cpp $(INCLUDES) -o refined

refined2d: refined2d.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG refined2d.cpp $(INCLUDES) -o refined2d

unrefined2d: unrefined2d.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG unrefined2d.cpp $(INCLUDES) -o unrefined2d

refined_scalability3d: refined_scalability3d.cpp ../../dccrg.hpp Makefile
	$(CXX) -DNDEBUG refined_scalability3d.cpp $(INCLUDES) -o refined_scalability3d

hierarchial_test: hierarchial_test.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG hierarchial_test.cpp $(INCLUDES) -o hierarchial_test

pinned_cells: pinned_cells.cpp ../../dccrg.hpp Makefile
	$(CXX) -DDEBUG pinned_cells.cpp $(INCLUDES) -o pinned_cells

c: clean
clean:
	rm -f game_of_life_test_all_cells game_of_life_test_single_cells game_of_life_user_data_size scalability scalability3d refined refined2d unrefined2d refined_scalability3d scalability1d hierarchial_test pinned_cells *vtk *visit
