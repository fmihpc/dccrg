/*
Program for testing dccrg neighbor offsets with AMR.

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
	const array<uint64_t, 3> grid_length = {{2, 1, 1}};
	grid.initialize(grid_length, comm, "RANDOM", 1);

	if (grid.is_local(1)) {
		grid.refine_completely(1);
	}
	grid.stop_refining();
	grid.clear_refined_unrefined_data();

	if (grid.is_local(3)) {
		grid.refine_completely(3);
	}
	grid.stop_refining();
	grid.clear_refined_unrefined_data();

	// do a few iterations with random load balancing
	for (int i = 0; i < 5; i++) {
		// check that local cell neighbor iterators work
		for (const auto& cell: grid.local_cells) {
			const auto cell_i = grid.mapping.get_indices(cell.id);
			const auto cell_length = grid.mapping.get_cell_length_in_indices(cell.id);
			for (const auto& neighbor: cell.neighbors_of) {
				const auto neigh_i = grid.mapping.get_indices(neighbor.id);
				const auto neigh_len = grid.mapping.get_cell_length_in_indices(neighbor.id);

				std::array<int, 3> ref_offsets{
					(int)neigh_i[0] - (int)cell_i[0],
					(int)neigh_i[1] - (int)cell_i[1],
					(int)neigh_i[2] - (int)cell_i[2]
				};
				for (size_t i = 0; i < ref_offsets.size(); i++) {
					if (
						neigh_len < cell_length
						and grid.indices_overlap(cell_i[i], cell_length, neigh_i[i], neigh_len)
					) {
						ref_offsets[i] = 0;
						continue;
					}
					if (neigh_len < cell_length) {
						ref_offsets[i] *= 2;
					}
					ref_offsets[i] /= (int)cell_length;
				}

				if (
					neighbor.x != ref_offsets[0]
					or neighbor.y != ref_offsets[1]
					or neighbor.z != ref_offsets[2]
				) {
					cerr << "FAILED" << endl;
					cout << "Wrong neighbor offset: " << neighbor.x << ", "
						<< neighbor.y << ", " << neighbor.z << ", should be "
						<< ref_offsets[0] << ", " << ref_offsets[1] << ", "
						<< ref_offsets[2] << " (excluding denominator) for neighbor "
						<< neighbor.id << " of cell " << cell.id
						<< " with respective indices " << neigh_i[0] << ", "
						<< neigh_i[1] << ", " << neigh_i[2] << " and "
						<< cell_i[0] << ", " << cell_i[1] << ", " << cell_i[2] << endl;
					return EXIT_FAILURE;
				}
			}
		}

		grid.balance_load();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
