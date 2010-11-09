/*
The simplest 2 D game of life program that demonstrates the basic usage of dccrg
The program doesn't scale, takes no input and produces no output.
For good scaling see the file game_of_life.cpp, for output see the files game_of_life_with_output.cpp and dc2vtk.cpp
*/

#include "boost/mpi.hpp"
#include "cstdlib"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;
using namespace boost::mpi;


// store in every cell of the grid whether the cell is alive and the number of live neighbours it has
struct game_of_life_cell {

	// boost requires this from user data
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
		/* live_neighbour_count from neighbouring cells is not used
		ar & live_neighbour_count;*/
	}

	bool is_alive;
	unsigned int live_neighbour_count;
};


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	// create the grid
	#define GRID_X_SIZE 10	// in unrefined cells
	#define GRID_Y_SIZE 10
	#define GRID_Z_SIZE 1
	#define CELL_SIZE 1.0
	#define STENCIL_SIZE 1	// the cells that share a vertex are considered neighbours
	#define MAX_REFINEMENT_LEVEL 0
	dccrg<game_of_life_cell> game_grid(comm, "RCB", 0, 0, 0, CELL_SIZE, CELL_SIZE, CELL_SIZE, GRID_X_SIZE, GRID_Y_SIZE, GRID_Z_SIZE, STENCIL_SIZE, MAX_REFINEMENT_LEVEL);	// use the recursive coordinate bisection method for load balancing (http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg_rcb.html)

	// since the grid doesn't change (isn't refined / unrefined) during the game, workload can be balanced just once in the beginning
	game_grid.balance_load();

	// get the cells on this process just once, since the grid doesn't change during the game
	vector<uint64_t> cells = game_grid.get_cells();


	// initialize the game
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbour_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	#define TURNS 100
	for (int turn = 0; turn < TURNS; turn++) {

		game_grid.update_remote_neighbour_data();

		// get the neighbour counts for every local cell
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			cell_data->live_neighbour_count = 0;
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);

			for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++) {
				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data->is_alive) {
					cell_data->live_neighbour_count++;
				}
			}
		}

		// calculate the next turn
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbour_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbour_count != 2) {
				cell_data->is_alive = false;
			}
		}

	}

	return EXIT_SUCCESS;
}
