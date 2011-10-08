/*
Tests the grid with variable amount of data in cells
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "iostream"
#include "unistd.h"
#include "vector"
#include "zoltan.h"

#include "../../dccrg_arbitrary_geometry.hpp"
#include "../../dccrg.hpp"

using namespace std;
using namespace boost::mpi;
using namespace dccrg;

struct CellData {

	vector<double> variables;

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & variables;
	}
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

	Dccrg<CellData, ArbitraryGeometry> grid;

	#define GRID_SIZE 3
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
		y_coordinates.push_back(i * CELL_SIZE);
	}
	z_coordinates.push_back(0);
	z_coordinates.push_back(1);
	grid.set_geometry(x_coordinates, y_coordinates, z_coordinates);

	grid.initialize(comm, "RANDOM", 1, 0);

	// populate the grid, number of variables in a cell is equal to its id
	vector<uint64_t> cells = grid.get_cells();
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		CellData* cell_data = grid[*cell];
		for (uint64_t i = 0; i < *cell; i++) {
			cell_data->variables.push_back(*cell + i);
		}
	}

	// print cell data
	for (int proc = 0; proc < comm.size(); proc++) {
		comm.barrier();
		if (proc != comm.rank()) {
			continue;
		}

		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			cout << "Cell " << *cell << " data (on process " << comm.rank() << "): ";

			CellData* cell_data = grid[*cell];
			for (vector<double>::const_iterator variable = cell_data->variables.begin(); variable != cell_data->variables.end(); variable++) {
				cout << *variable << " ";
			}
			cout << endl;
		}
		cout.flush();
		sleep(3);
	}

	grid.balance_load();

	if (comm.rank() == 0) {
		cout << endl;
	}

	cells = grid.get_cells();

	// print cell data again
	for (int proc = 0; proc < comm.size(); proc++) {
		comm.barrier();
		if (proc != comm.rank()) {
			continue;
		}

		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			cout << "Cell " << *cell << " data (on process " << comm.rank() << "): ";

			CellData* cell_data = grid[*cell];
			for (vector<double>::const_iterator variable = cell_data->variables.begin(); variable != cell_data->variables.end(); variable++) {
				cout << *variable << " ";
			}
			cout << endl;
		}
		cout.flush();
		sleep(3);
	}

	return EXIT_SUCCESS;
}
