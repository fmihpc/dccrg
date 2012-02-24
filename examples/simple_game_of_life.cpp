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
using namespace dccrg;

// store in every cell of the grid whether the cell is alive and the number of live neighbors it has
struct game_of_life_cell {

	// boost requires this from user data
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
		/* live_neighbor_count from neighboring cells is not used
		ar & live_neighbor_count;*/
	}

	bool is_alive;
	unsigned int live_neighbor_count;
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

	Dccrg<game_of_life_cell> game_grid;

	#define X_LENGTH 10	// in unrefined cells
	#define Y_LENGTH 10
	#define Z_LENGTH 1
	#define CELL_SIZE 1.0
	game_grid.set_geometry(X_LENGTH, Y_LENGTH, Z_LENGTH, 0, 0, 0, CELL_SIZE, CELL_SIZE, CELL_SIZE);

	// the cells that share a vertex are considered neighbors
	#define NEIGHBORHOOD_SIZE 1
	#define MAX_REFINEMENT_LEVEL 0
	// use the recursive coordinate bisection method for load balancing (http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg_rcb.html)
	game_grid.initialize(comm, "RCB", NEIGHBORHOOD_SIZE, MAX_REFINEMENT_LEVEL);

	// since the grid doesn't change (isn't refined / unrefined) during the game, workload can be balanced just once in the beginning
	game_grid.balance_load();

	// get the cells on this process just once, since the grid doesn't change during the game
	vector<uint64_t> cells = game_grid.get_cells();


	// initialize the game
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	#define TURNS 100
	for (int turn = 0; turn < TURNS; turn++) {

		game_grid.update_remote_neighbor_data();

		// get the neighbor counts for every local cell
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			cell_data->live_neighbor_count = 0;
			const vector<uint64_t>* neighbors = game_grid.get_neighbors(*cell);

			for (vector<uint64_t>::const_iterator neighbor = neighbors->begin(); neighbor != neighbors->end(); neighbor++) {

				/*
				neighbors that would be outside of the grid
				are recorded as 0, skip them
				*/
				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// calculate the next turn
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbor_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbor_count != 2) {
				cell_data->is_alive = false;
			}
		}

	}

	return EXIT_SUCCESS;
}
