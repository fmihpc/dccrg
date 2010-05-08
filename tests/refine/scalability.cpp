/*
Tests the speed of refining the grid in 3-d by refining random cells until enough cell exist
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "../../dccrg.hpp"
#include "fstream"
#include "functional"
#include "iostream"
#include "zoltan.h"


using namespace std;
using namespace boost::mpi;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}


	#define STARTING_CORNER 0.0
	#define GRID_SIZE 21
	#define CELL_SIZE (1.0 / GRID_SIZE)
	#define STENCIL_SIZE 1
	dccrg<int> grid(comm, "RCB", STARTING_CORNER, STARTING_CORNER, STARTING_CORNER, CELL_SIZE, GRID_SIZE, GRID_SIZE, GRID_SIZE, STENCIL_SIZE);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << grid.get_max_refinement_level() << endl;
	}
	grid.balance_load();

	vector<uint64_t> cells = grid.get_cells();

	clock_t before, after, total = 0;

	// refine random cells until every process has enough cells
	#define MAX_CELLS (100 * GRID_SIZE * GRID_SIZE * GRID_SIZE)
	vector<uint64_t> new_cells;
	do {
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		before = clock();
		for (int i = 0; i < int(cells.size() / 15); i++) {
			grid.refine_completely(cells[i]);
		}
		new_cells = grid.stop_refining();
		after = clock();
		cout << "Proc " << comm.rank() << ": " << new_cells.size() << " new cells" << endl;

		for (vector<uint64_t>::const_iterator new_cell = new_cells.begin(); new_cell != new_cells.end(); new_cell++) {
			cells.push_back(*new_cell);
		}
		new_cells.clear();

		total += after - before;
	} while (all_reduce(comm, int(cells.size()), plus<int>()) < MAX_CELLS);

	cout << "Process " << comm.rank() << ": " << (cells.size() - GRID_SIZE * GRID_SIZE * GRID_SIZE) << " new cells created in " << double(total) / CLOCKS_PER_SEC << " s, (" << (cells.size() - GRID_SIZE * GRID_SIZE * GRID_SIZE) / (double(total) / CLOCKS_PER_SEC) << " new cells / s)" << endl;

	return EXIT_SUCCESS;
}

