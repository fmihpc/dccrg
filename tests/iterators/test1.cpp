/*
Program for testing dccrg iterators.

Copyright 2013, 2014, 2015, 2016, 2018 Finnish Meteorological Institute
Copyright 2018 Ilja Honkonen

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
#include "set"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;
using namespace dccrg;

struct Cell {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) this, 0, MPI_BYTE);
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

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}

	// initialize grid
	Dccrg<Cell> grid;
	const std::array<uint64_t, 3> grid_length = {{1000, 1, 1}};
	grid.initialize(grid_length, comm, "RANDOM", 3);

	// do a few iterations with random load balancing
	for (int i = 0; i < 5; i++) {
		set<uint64_t> inner_ref, outer_ref;
		for (const auto& cell: grid.cells) {
			if (not grid.is_local(cell.id)) {
				continue;
			}

			const auto* const neighbors_of = grid.get_neighbors_of(cell.id);
			if (neighbors_of == nullptr) {
				cout << "No neighbors_of for cell " << cell.id << endl;
				abort();
			}

			bool local = true;
			for (const auto& neighbor: *neighbors_of) {
				if (neighbor.first == dccrg::error_cell) {
					continue;
				}
				if (not grid.is_local(neighbor.first)) {
					local = false;
					break;
				}
			}

			const auto* const neighbors_to = grid.get_neighbors_to(cell.id);
			if (neighbors_to == nullptr) {
				cout << "No neighbors_to for cell " << cell.id << endl;
				abort();
			}

			for (const auto& neighbor: *neighbors_to) {
				if (neighbor.first == dccrg::error_cell) {
					continue;
				}
				if (not grid.is_local(neighbor.first)) {
					local = false;
					break;
				}
			}

			if (local) {
				inner_ref.insert(cell.id);
			} else {
				outer_ref.insert(cell.id);
			}
		}

		// check that inner cell iterator works
		for (const auto& cell: grid.inner_cells) {
			if (inner_ref.count(cell.id) == 0) {
				cout << "FAILED on turn " << i << endl;
				cerr << "Cell " << cell.id
					<< " is not in inner cells"
					<< endl;
				abort();
			}

			if (outer_ref.count(cell.id) > 0) {
				cout << "FAILED on turn " << i << endl;
				cerr << "Cell " << cell.id
					<< " is in outer cells"
					<< endl;
				abort();
			}
		}

		// check that outer cell iterator works
		for (const auto& cell: grid.outer_cells) {
			if (inner_ref.count(cell.id) > 0) {
				cout << "FAILED on turn " << i << endl;
				cerr << "Cell " << cell.id
					<< " is in inner cells"
					<< endl;
				abort();
			}

			if (outer_ref.count(cell.id) == 0) {
				cout << "FAILED on turn " << i << endl;
				cerr << "Cell " << cell.id
					<< " is not in outer cells"
					<< endl;
				abort();
			}
		}

		grid.balance_load();
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

