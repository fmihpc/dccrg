/*
Tests the grid with variable amount of data in cells
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "../../dccrg.hpp"
#include "iostream"
#include "unistd.h"
#include "vector"
#include "zoltan.h"

using namespace std;

struct CellData {

	vector<double> variables;

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & variables;
	}
};


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

	#define GRID_SIZE 3
	dccrg<CellData> grid(comm, "RANDOM", 0.0, 0.0, 0.0, 1, GRID_SIZE, GRID_SIZE, 1, 1, 0);

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
