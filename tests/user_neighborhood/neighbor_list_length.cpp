/*
Program for testing user defined neighborhoods of dccrg.

Copyright 2012, 2013, 2014, 2015, 2016 Finnish Meteorological Institute

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

#include "array"
#include "cstdlib"
#include "iostream"
#include "tuple"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;

struct Cell
{
	int data;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) &(this->data), 1, MPI_INT);
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
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}

	dccrg::Dccrg<Cell> grid;

	const std::array<uint64_t, 3> grid_length = {{10, 10, 10}};
	grid
		.set_initial_length(grid_length)
		.set_periodic(true, true, true)
		.set_maximum_refinement_level(0)
		.set_neighborhood_length(2)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	// default neighbor lists should have 5^3 - 1 = 124 neighbors
	for (const auto& cell: grid.local_cells()) {
		const auto nofs = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		if (nofs != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< ": " << nofs
				<< std::endl;
			abort();
		}
		const auto ntos = std::distance(cell.neighbors_to.begin(), cell.neighbors_to.end());
		if (ntos != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< ": " << ntos
				<< std::endl;
			abort();
		}
	}

	typedef dccrg::Types<3>::neighborhood_item_t neigh_t;

	// create a neighborhood of one cell at (-2, -2, -2)
	{
	const int hood_id = 1;
	const std::vector<neigh_t> neighborhood{{-2, -2, -2}};
	if (!grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood failed"
				<< std::endl;
			abort();
	}
	for (const auto& cell: grid.local_cells(hood_id)) {
		const auto nofs = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		if (nofs != 1) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << nofs
				<< std::endl;
			abort();
		}
		const auto ntos = std::distance(cell.neighbors_to.begin(), cell.neighbors_to.end());
		if (ntos != 1) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << ntos
				<< std::endl;
			abort();
		}
	}
	grid.remove_neighborhood(hood_id);
	}

	// create a neighborhood of two cells at (-1, -1, -1) and (2, 2, 2)
	{
	const int hood_id = 2;
	const std::vector<neigh_t> neighborhood{
		{-1, -1, -1},
		{ 2,  2,  2}
	};
	if (!grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood failed"
				<< std::endl;
			abort();
	}
	for (const auto& cell: grid.local_cells(hood_id)) {
		const auto nofs = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		if (nofs != 2) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << nofs
				<< std::endl;
			abort();
		}
		const auto ntos = std::distance(cell.neighbors_to.begin(), cell.neighbors_to.end());
		if (ntos != 2) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << ntos
				<< std::endl;
			abort();
		}
	}
	grid.remove_neighborhood(hood_id);
	}

	// create a neighborhood of 24 cells at z offset == 0
	{
	const int hood_id = -3;
	std::vector<neigh_t> neighborhood;
	for (int i = -2; i <= 2; i++)
	for (int j = -2; j <= 2; j++) {
		if (i == 0 && j == 0) {
			continue;
		}
		neigh_t temp = {{i, j, 0}};
		neighborhood.push_back(temp);
	}
	if (!grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood failed"
				<< std::endl;
			abort();
	}
	for (const auto& cell: grid.local_cells(hood_id)) {
		const auto nofs = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		if (nofs != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << nofs
				<< std::endl;
			abort();
		}
		const auto ntos = std::distance(cell.neighbors_to.begin(), cell.neighbors_to.end());
		if (ntos != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << ntos
				<< std::endl;
			abort();
		}
	}
	grid.remove_neighborhood(hood_id);
	}

	// create a neighborhood of 24 cells at x offset == 0
	{
	const int hood_id = -4;
	std::vector<neigh_t> neighborhood;
	for (int j = -2; j <= 2; j++)
	for (int k = -2; k <= 2; k++) {
		if (j == 0 && k == 0) {
			continue;
		}
		neigh_t temp = {{0, j, k}};
		neighborhood.push_back(temp);
	}
	if (!grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood failed"
				<< std::endl;
			abort();
	}
	for (const auto& cell: grid.local_cells(hood_id)) {
		const auto nofs = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		if (nofs != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << nofs
				<< std::endl;
			abort();
		}
		const auto ntos = std::distance(cell.neighbors_to.begin(), cell.neighbors_to.end());
		if (ntos != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << ntos
				<< std::endl;
			abort();
		}
	}
	grid.remove_neighborhood(hood_id);
	}

	// create a full neighborhood identical to the default one
	{
	const int hood_id = 0;
	std::vector<neigh_t> neighborhood;
	for (int k = -2; k <= 2; k++)
	for (int j = -2; j <= 2; j++)
	for (int i = -2; i <= 2; i++) {
		if (i == 0 && j == 0 && k == 0) {
			continue;
		}
		neigh_t temp = {{i, j, k}};
		neighborhood.push_back(temp);
	}
	if (!grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood failed"
				<< std::endl;
			abort();
	}
	for (const auto& cell: grid.local_cells(hood_id)) {
		// check number of neighbors
		const auto nofs = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		if (nofs != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << nofs
				<< std::endl;
			abort();
		}
		const auto ntos = std::distance(cell.neighbors_to.begin(), cell.neighbors_to.end());
		if (ntos != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell.id
				<< " and neighborhood id " << hood_id
				<< ": " << ntos
				<< std::endl;
			abort();
		}
	}
	}

	// try to create a neighborhood too far away
	{
	const int hood_id = 5;
	const std::vector<neigh_t> neighborhood{{-1, 2, -3}};
	if (grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood succeeded"
				<< std::endl;
			abort();
	}
	}

	// try to use an existing neighborhood id
	{
	const int hood_id = 0;
	const std::vector<neigh_t> neighborhood{{1, 2, -2}};
	if (grid.add_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_neighborhood succeeded"
				<< std::endl;
			abort();
	}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
