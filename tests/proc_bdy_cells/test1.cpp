/*
Tests dccrg logic for cells on process boundaries.

Copyright 2018 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "cstdlib"
#include "iostream"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;
using namespace dccrg;

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

	Dccrg<int> grid;
	grid
		.set_initial_length({3, 1, 1})
		.set_neighborhood_length(2)
		.set_maximum_refinement_level(1)
		.set_load_balancing_method("RCB")
		.initialize(comm)
		.balance_load();

	grid.refine_completely(3);
	const auto new_cells = grid.stop_refining(true);
	grid.balance_load();

	grid.update_copies_of_remote_neighbors();

	bool fail = false;
	for (const auto cell_id: grid.get_remote_cells_on_process_boundary()) {
		const auto* const cell_data = grid[cell_id];
		if (cell_data == nullptr) {
			for (const auto& cell: grid.local_cells()) {
				if (cell.id == cell_id) {
					std::cerr << __FILE__  "(" << __LINE__ << ")" << std::endl;
					fail = true;
				}
				for (const auto& neighbor: cell.neighbors_of) {
					if (neighbor.id == cell_id) {
						std::cerr << __FILE__  "(" << __LINE__ << ")" << std::endl;
						fail = true;
					}
				}
			}
			// otherwise nullptr is not an error
		}
	}

	MPI_Finalize();

	if (fail) {
		return EXIT_FAILURE;
	} else {
		return EXIT_SUCCESS;
	}
}
