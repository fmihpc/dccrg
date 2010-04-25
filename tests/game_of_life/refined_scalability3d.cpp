/*
Tests the scalability of the grid in 3 D with refined grid
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "../../dccrg.hpp"
#include "fstream"
#include "functional"
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

	#define STARTING_CORNER 0.0
	#define GRID_SIZE 50
	#define CELL_SIZE (1.0 / GRID_SIZE)
	#define STENCIL_SIZE 1
	dccrg<game_of_life_cell> game_grid(comm, "RCB", STARTING_CORNER, STARTING_CORNER, STARTING_CORNER, CELL_SIZE, GRID_SIZE, GRID_SIZE, GRID_SIZE, STENCIL_SIZE, 2);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << game_grid.get_max_refinement_level() << endl;
	}
	game_grid.balance_load();

	vector<uint64_t> cells = game_grid.get_cells();

	// refine random cells until every process has enough cells
	#define MAX_CELLS 1000000.0
	vector<uint64_t> new_cells;
	do {
		random_shuffle(cells.begin(), cells.end());

		cout << "Proc " << comm.rank() << " refining" << endl;

		// refine 1/10 of the cells needed to fill the quota every round
		for (int i = 0; i < int(cells.size()) && i < int((MAX_CELLS / comm.size() - cells.size()) / 10); i++) {
			game_grid.refine_completely(cells[i]);
		}

		new_cells = game_grid.stop_refining();
		cells = game_grid.get_cells();
	} while (all_reduce(comm, int(new_cells.size()), plus<int>()) > 0);

	game_grid.balance_load();
	cells = game_grid.get_cells();

	// initialize the game with random cells alive
	srand(time(NULL));
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbour_count = 0;

		if (double(rand()) / RAND_MAX <= 0.2) {
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
	after = clock();
	total += double(after) - before;
	if (comm.rank() == 0) {
		cout << endl;
	}
	comm.barrier();

	cout << "Process " << comm.rank() << ": " << cells.size() * TIME_STEPS << " cells processed at the speed of " << double(cells.size() * TIME_STEPS) * CLOCKS_PER_SEC / total << " cells / second"<< endl;

	return EXIT_SUCCESS;
}
