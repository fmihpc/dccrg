/*
Tests the speed of refining the grid in 3-d by refining random cells until enough cell exist
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "functional"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg_arbitrary_geometry.hpp"
#include "../../dccrg.hpp"


using namespace std;
using namespace boost::mpi;
using namespace dccrg;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	Dccrg<int, ArbitraryGeometry> grid;

	#define GRID_SIZE 21
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
		y_coordinates.push_back(i * CELL_SIZE);
		z_coordinates.push_back(i * CELL_SIZE);
	}
	grid.set_geometry(x_coordinates, y_coordinates, z_coordinates);

	#define NEIGHBORHOOD_SIZE 1
	grid.initialize(comm, "RCB", NEIGHBORHOOD_SIZE);
	grid.balance_load();
	vector<uint64_t> cells = grid.get_cells();

	clock_t before, after, total = 0;

	// refine random cells until there are enough cells in total
	#define MAX_CELLS (100 * GRID_SIZE * GRID_SIZE * GRID_SIZE)
	do {
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		before = clock();
		for (int i = 0; i < int(cells.size() / 15) && i < 10000; i++) {
			grid.refine_completely(cells[i]);
		}
		vector<uint64_t> new_cells = grid.stop_refining();
		after = clock();
		total += after - before;
		cout << "Proc " << comm.rank() << ": " << new_cells.size() << " new cells" << endl;

		cells = grid.get_cells();
	} while (all_reduce(comm, int(cells.size()), plus<int>()) < MAX_CELLS);

	cout << "Process " << comm.rank() << ": " << (cells.size() - GRID_SIZE * GRID_SIZE * GRID_SIZE) << " new cells created in " << double(total) / CLOCKS_PER_SEC << " s, (" << (cells.size() - GRID_SIZE * GRID_SIZE * GRID_SIZE) / (double(total) / CLOCKS_PER_SEC) << " new cells / s)" << endl;

	return EXIT_SUCCESS;
}

