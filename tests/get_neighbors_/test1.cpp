/*
Program for testing get_neighbors_() function of dccrg.

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

#include "dccrg.hpp"

int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "Coudln't initialize MPI." << std::endl;
		abort();
	}
	MPI_Comm comm = MPI_COMM_WORLD;

	{dccrg::Dccrg<int> grid; grid
		.set_initial_length({1, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(0)
		.initialize(comm);
	const auto& neighs = grid.get_neighbors_();
	for (const auto& neighs_i: neighs) {
		for (const auto& neigh: neighs_i.second) {
			if (neigh != dccrg::error_cell) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh << std::endl;
				abort();
			}
		}
	}}

	{dccrg::Dccrg<int> grid; grid
		.set_periodic(true, true, true)
		.set_initial_length({1, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(0)
		.initialize(comm);
	const auto& neighs = grid.get_neighbors_();
	for (const auto& neighs_i: neighs) {
		for (const auto& neigh: neighs_i.second) {
			if (neigh != 1) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh << std::endl;
				abort();
			}
		}
	}}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
