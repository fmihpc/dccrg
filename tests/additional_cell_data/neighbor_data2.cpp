/*
Tests additional neighbor data functionality of dccrg.

Copyright 2018 Ilja Honkonen
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
#include "tuple"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"
#include "../../dccrg_no_geometry.hpp"

using namespace std;
using namespace dccrg;

struct Cell {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) this, 0, MPI_BYTE);
	}
};

struct Additional_Neighbor_Data_Item1 {
	uint64_t
		additional_neighbor_data1 = dccrg::error_cell,
		additional_neighbor_data2 = dccrg::error_cell;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid&, const Cell_Item& cell, const Neighbor_Item& neighbor, const int&, const Additional_Neighbor_Data_Item1&
	) {
		additional_neighbor_data1 = neighbor.id;
		additional_neighbor_data2 = cell.id;
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
	Dccrg<
		Cell,
		No_Geometry,
		std::tuple<>,
		std::tuple<Additional_Neighbor_Data_Item1>
	> grid; grid
		.set_initial_length({2, 2, 2})
		.set_neighborhood_length(3)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	for (const auto& cell: grid.local_cells()) {
		for (const auto& neighbor: cell.neighbors_of) {
			if (neighbor.additional_neighbor_data1 != neighbor.id) {
				std::cerr << __LINE__ << ": " << neighbor.id << " != " << neighbor.additional_neighbor_data1 << std::endl;
				abort();
			}
			if (neighbor.additional_neighbor_data2 != cell.id) {
				std::cerr << __LINE__ << ": " << cell.id << " != " << neighbor.additional_neighbor_data2 << std::endl;
				abort();
			}
		}
	}

	grid.balance_load();

	for (const auto& cell: grid.local_cells()) {
		for (const auto& neighbor: cell.neighbors_of) {
			if (neighbor.additional_neighbor_data1 != neighbor.id) {
				std::cerr << __LINE__ << ": " << neighbor.id << " != " << neighbor.additional_neighbor_data1 << std::endl;
				abort();
			}
			if (neighbor.additional_neighbor_data2 != cell.id) {
				std::cerr << __LINE__ << ": " << cell.id << " != " << neighbor.additional_neighbor_data2 << std::endl;
				abort();
			}
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
