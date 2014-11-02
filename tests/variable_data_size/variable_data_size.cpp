/*
Tests the grid with variable amount of data in cells using serialization
*/

#include "cstdlib"
#include "ctime"
#include "iostream"
#include "unistd.h"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg_stretched_cartesian_geometry.hpp"
#include "../../dccrg.hpp"

using namespace std;
using namespace dccrg;

struct CellData {

	std::vector<double> variables;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(this->variables.data(), int(this->variables.size()), MPI_DOUBLE);
	}
};


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	Dccrg<CellData> grid;

	const std::array<uint64_t, 3> grid_length = {{3, 1, 1}};
	grid.initialize(grid_length, comm, "RANDOM", 1, 0);

	// populate the grid, number of variables in a cell is equal to its id
	vector<uint64_t> cells = grid.get_cells();
	for (auto cell: cells) {
		auto* const cell_data = grid[cell];
		for (uint64_t i = 0; i < cell; i++) {
			cell_data->variables.push_back(cell + i);
		}
	}

	// print cell data
	for (int proc = 0; proc < comm_size; proc++) {
		MPI_Barrier(comm);
		if (proc != rank) {
			continue;
		}

		for (auto cell: cells) {
			cout << "Cell " << cell << " data (on process " << rank << "): ";

			const auto* const cell_data = grid[cell];
			for (auto variable: cell_data->variables) {
				cout << variable << " ";
			}
			cout << endl;
		}
		cout.flush();
		sleep(3);
	}

	grid.initialize_balance_load(true);
	// make room for incoming cell data
	const auto& cells_to_receive = grid.get_cells_to_receive();
	for (const auto& sender: cells_to_receive) {
		for (const auto item: sender.second) {
			const auto cell = item.first;
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) { abort(); }
			cell_data->variables.resize(cell);
		}
	}
	grid.continue_balance_load();
	grid.finish_balance_load();

	if (rank == 0) {
		cout << endl;
	}

	cells = grid.get_cells();

	// print cell data again
	for (int proc = 0; proc < comm_size; proc++) {
		MPI_Barrier(comm);
		if (proc != rank) {
			continue;
		}

		for (auto cell: cells) {

			cout << "Cell " << cell << " data (on process " << rank << "): ";

			const auto* const cell_data = grid[cell];
			for (auto variable: cell_data->variables) {
				cout << variable << " ";
			}
			cout << endl;
		}
		cout.flush();
		sleep(3);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

