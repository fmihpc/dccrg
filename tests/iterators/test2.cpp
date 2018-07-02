/*
Program for testing dccrg neighbor iterators.

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
#include "vector"

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
	grid
		.set_initial_length({1000, 1, 1})
		.set_neighborhood_length(3)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	// do a few iterations with random load balancing
	for (int i = 0; i < 5; i++) {
		// check that local cell neighbor iterators work
		for (const auto& cell: grid.local_cells) {
			vector<uint64_t> ref_neighbors_of;
			const auto
				start_cell = [&]()->uint64_t{
					if (cell.id < 4) {
						return 1;
					} else {
						return cell.id - 3;
					}
				}(),
				end_cell = [&]()->uint64_t{
					if (cell.id > 997) {
						return 1000;
					} else {
						return cell.id + 3;
					}
				}();
			for (uint64_t c = start_cell; c <= end_cell; c++) {
				if (c == cell.id) {
					continue;
				}
				ref_neighbors_of.push_back(c);
			}

			vector<uint64_t> neighbors_of;
			for (const auto& neighbor: cell.neighbors_of) {
				neighbors_of.push_back(neighbor.id);
			}
			if (ref_neighbors_of != neighbors_of) {
				cerr << "FAILED" << endl;
				cout << "Wrong neighbors_of for cell " << cell.id
					<< "\nResult: ";
				for (const auto& n: neighbors_of) cout << n << ", ";
				cout << "\nReference: ";
				for (const auto& n: ref_neighbors_of) cout << n << ", ";
				cout << endl;
				return EXIT_FAILURE;
			}

			// check offset
			for (const auto& neighbor: cell.neighbors_of) {
				const int ref_diff = neighbor.id - cell.id;
				if (ref_diff != neighbor.x) {
					cerr << "FAILED" << endl;
					cout << "Wrong offset for neighbor_of " << neighbor.id
						<< " of cell " << cell.id << ": " << neighbor.x
						<< ", should be " << ref_diff << endl;
					return EXIT_FAILURE;
				}
			}
		}

		grid.balance_load();
	}

	/* TODO: test user-defined neighborhood
	if (not grid.add_neighborhood(0, {{1, 0, 0}, {2, 0, 0}, {3, 0, 0}})) {
		cout << "FAILED" << endl;
		cerr << "Failed to add neighborhood with id 0" << endl;
	}*/

	MPI_Finalize();

	return EXIT_SUCCESS;
}
