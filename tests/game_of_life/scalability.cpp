/*
Tests the scalability of the grid in 2 D
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "../../dccrg.hpp"
#include "fstream"
#include "iostream"
#include "zoltan.h"


struct game_of_life_cell {

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
	}

	bool is_alive;
	unsigned int live_neighbour_count;
};


using namespace std;
using namespace boost::mpi;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	clock_t before, after, total = 0;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}


	#define STARTING_CORNER 0.0
	#define GRID_SIZE 1000	// in unrefined cells
	#define CELL_SIZE (1.0 / GRID_SIZE)
	#define STENCIL_SIZE 1
	#define MAX_REFINEMENT_LEVEL 0
	dccrg<game_of_life_cell> game_grid(comm, "RCB", STARTING_CORNER, STARTING_CORNER, STARTING_CORNER, CELL_SIZE, GRID_SIZE, GRID_SIZE, 1, STENCIL_SIZE, MAX_REFINEMENT_LEVEL);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << game_grid.get_max_refinement_level() << endl;
		cout << "Number of cells: " << GRID_SIZE * GRID_SIZE * 1 << endl << endl;
	}

	game_grid.balance_load();
	vector<uint64_t> cells = game_grid.get_cells();
	// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
	sort(cells.begin(), cells.end());

	// initialize the game with a line of living cells in the x direction in the middle
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbour_count = 0;

		double y = game_grid.get_cell_y(*cell);
		if (fabs(0.5 + 0.1 * game_grid.get_cell_size(*cell) - y) < 0.5 * game_grid.get_cell_size(*cell)) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	if (comm.rank() == 0) {
		cout << "step: ";
	}

	#define TIME_STEPS 100
	before = clock();
	for (int step = 0; step < TIME_STEPS; step++) {
		comm.barrier();

		if (comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		// get the neighbour counts of every cell
		game_grid.update_remote_neighbour_data();
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			cell_data->live_neighbour_count = 0;
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			if (neighbours == NULL) {
				cout << "Process " << comm.rank() << ": neighbour list for cell " << *cell << " not available" << endl;
				exit(EXIT_FAILURE);
			}

			for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++) {

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data == NULL) {
					cout << "Process " << comm.rank() << ": neighbour " << *neighbour << " data for cell " << *cell << " not available" << endl;
					exit(EXIT_FAILURE);
				}

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
	after = clock();
	total += after - before;
	if (comm.rank() == 0) {
		cout << endl;
	}
	comm.barrier();

	cout << "Process " << comm.rank() << ": " << cells.size() * TIME_STEPS << " cells processed at the speed of " << double(cells.size() * TIME_STEPS) * CLOCKS_PER_SEC / total << " cells / second"<< endl;

	return EXIT_SUCCESS;
}
