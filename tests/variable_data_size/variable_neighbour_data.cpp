/*
Tests the grid with variable amount of data in cells and
variable amount of data sent during neighbor data updates
using serialization.
*/

#include "array"
#include "cstdlib"
#include "ctime"
#include "iostream"
#include "unistd.h"
#include "vector"
#include "zoltan.h"

#include "../../dccrg_stretched_cartesian_geometry.hpp"
#include "../../dccrg.hpp"

using namespace std;
using namespace dccrg;


static bool send_variables2;

class CellData {
public:
	vector<int> variables1, variables2;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		if (not send_variables2) {
			return std::make_tuple(this->variables1.data(), this->variables1.size(), MPI_INT);
		} else {
			array<int, 2> counts{int(this->variables1.size()), int(this->variables2.size())};
			array<MPI_Aint, 2> displacements{
				0,
				(uint8_t*) this->variables2.data() - (uint8_t*) this->variables1.data()
			};
			array<MPI_Datatype, 2> datatypes{MPI_INT, MPI_INT};

			MPI_Datatype final_datatype = MPI_DATATYPE_NULL;
			if (
				MPI_Type_create_struct(
					2,
					counts.data(),
					displacements.data(),
					datatypes.data(),
					&final_datatype
				) != MPI_SUCCESS
			) {
				return std::make_tuple((void*) this->variables1.data(), 0, MPI_BYTE);
			}
			return std::make_tuple((void*) this->variables1.data(), 1, final_datatype);
		}
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
	grid.allocate_copies_of_remote_neighbors();

	// populate the grid, number of variables in a cell is equal to its id
	auto cells = grid.get_cells();
	for (auto cell: cells) {
		auto* const cell_data = grid[cell];
		if (cell_data == NULL) { abort(); }

		for (uint64_t i = 0; i < cell; i++) {
			cell_data->variables1.push_back(cell + i);
			cell_data->variables2.push_back(-(cell + i));
		}
	}

	// make room for incoming cell data
	auto remote_cells = grid.get_remote_cells_on_process_boundary();
	for (auto cell: remote_cells) {
		auto* const cell_data = grid[cell];
		if (cell_data == NULL) { abort(); }

		cell_data->variables1.resize(cell);
		cell_data->variables2.resize(cell);
	}

	send_variables2 = true;
	grid.update_copies_of_remote_neighbors();


	// print cell variables1 and neighbor variables2
	for (int proc = 0; proc < comm_size; proc++) {
		MPI_Barrier(comm);
		if (proc != rank) {
			continue;
		}

		for (auto cell: cells) {
			cout << "Cell " << cell << " data (on process " << rank << "): ";

			const auto* const cell_data = grid[cell];
			if (cell_data == NULL) { abort(); }

			for (auto variable: cell_data->variables1) {
				cout << variable << " ";
			}

			const auto* neighbors = grid.get_neighbors_of(cell);
			for (auto neighbor: *neighbors) {
				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = grid[neighbor];
				if (neighbor_data == NULL) { abort(); }

				for (auto variable: neighbor_data->variables2) {
					cout << variable << " ";
				}
			}
			cout << endl;
		}
		cout.flush();
		sleep(2);
	}

	grid.initialize_balance_load(true);
	// make room for incoming cell data
	const auto& cells_to_receive = grid.get_cells_to_receive();
	for (const auto& sender: cells_to_receive) {
		for (const auto item: sender.second) {
			const auto cell = item.first;
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) { abort(); }

			cell_data->variables1.resize(cell);
			cell_data->variables2.resize(cell);
		}
	}
	grid.continue_balance_load();
	grid.finish_balance_load();

	grid.allocate_copies_of_remote_neighbors();

	remote_cells = grid.get_remote_cells_on_process_boundary();
	for (auto cell: remote_cells) {
		auto* const cell_data = grid[cell];
		if (cell_data == NULL) { abort(); }

		cell_data->variables1.resize(cell);
		cell_data->variables2.resize(cell);
	}
	grid.update_copies_of_remote_neighbors();

	cells = grid.get_cells();

	if (rank == 0) {
		cout << endl;
	}

	// print cell data again
	for (int proc = 0; proc < comm_size; proc++) {
		MPI_Barrier(comm);
		if (proc != rank) {
			continue;
		}

		for (auto cell: cells) {
			cout << "Cell " << cell << " data (on process " << rank << "): ";

			const auto* const cell_data = grid[cell];
			for (auto variable: cell_data->variables1) {
				cout << variable << " ";
			}

			const auto* neighbors = grid.get_neighbors_of(cell);
			for (auto neighbor: *neighbors) {
				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = grid[neighbor];

				for (auto variable: neighbor_data->variables2) {
					cout << variable << " ";
				}
			}
			cout << endl;
		}
		cout.flush();
		sleep(2);
	}

	grid.initialize_balance_load(true);
	for (const auto& sender: cells_to_receive) {
		for (const auto item: sender.second) {
			const auto cell = item.first;
			auto* const cell_data = grid[cell];
			if (cell_data == NULL) { abort(); }

			cell_data->variables1.resize(cell);
			cell_data->variables2.resize(cell);
		}
	}
	grid.continue_balance_load();
	grid.finish_balance_load();

	grid.allocate_copies_of_remote_neighbors();
	remote_cells = grid.get_remote_cells_on_process_boundary();
	for (auto cell: remote_cells) {
		auto* const cell_data = grid[cell];
		if (cell_data == NULL) { abort(); }

		cell_data->variables1.resize(cell);
		cell_data->variables2.resize(cell);
	}

	cells = grid.get_cells();

	if (rank == 0) {
		cout << endl;
	}

	send_variables2 = false;
	grid.update_copies_of_remote_neighbors();

	// print cell data
	for (int proc = 0; proc < comm_size; proc++) {
		MPI_Barrier(comm);
		if (proc != rank) {
			continue;
		}

		for (auto cell: cells) {
			cout << "Cell " << cell << " data (on process " << rank << "): ";

			const auto* const cell_data = grid[cell];
			for (auto variable: cell_data->variables1) {
				cout << variable << " ";
			}

			const auto* const neighbors = grid.get_neighbors_of(cell);
			for (auto neighbor: *neighbors) {
				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = grid[neighbor];

				for (auto variable: neighbor_data->variables2) {
					cout << variable << " ";
				}
			}
			cout << endl;
		}
		cout.flush();
		sleep(2);
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

