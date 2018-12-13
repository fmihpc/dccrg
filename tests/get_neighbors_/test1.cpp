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

	for (int neigh_len: {0, 1, 2, 3}) {
		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != dccrg::error_cell
				or neighs_i.second[1] != dccrg::error_cell
				or neighs_i.second[2] != dccrg::error_cell
				or neighs_i.second[3] != dccrg::error_cell
				or neighs_i.second[4] != dccrg::error_cell
				or neighs_i.second[5] != dccrg::error_cell
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_periodic(true, true, true)
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != 1
				or neighs_i.second[1] != 1
				or neighs_i.second[2] != 1
				or neighs_i.second[3] != 1
				or neighs_i.second[4] != 1
				or neighs_i.second[5] != 1
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_periodic(false, true, true)
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != dccrg::error_cell
				or neighs_i.second[1] != dccrg::error_cell
				or neighs_i.second[2] != 1
				or neighs_i.second[3] != 1
				or neighs_i.second[4] != 1
				or neighs_i.second[5] != 1
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_periodic(true, false, true)
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != 1
				or neighs_i.second[1] != 1
				or neighs_i.second[2] != dccrg::error_cell
				or neighs_i.second[3] != dccrg::error_cell
				or neighs_i.second[4] != 1
				or neighs_i.second[5] != 1
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_periodic(true, true, false)
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != 1
				or neighs_i.second[1] != 1
				or neighs_i.second[2] != 1
				or neighs_i.second[3] != 1
				or neighs_i.second[4] != dccrg::error_cell
				or neighs_i.second[5] != dccrg::error_cell
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}


		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({3, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[2] != dccrg::error_cell
				or neighs_i.second[3] != dccrg::error_cell
				or neighs_i.second[4] != dccrg::error_cell
				or neighs_i.second[5] != dccrg::error_cell
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch (neighs_i.first) {
			case 1:
				if (
					neighs_i.second[0] != dccrg::error_cell
					or neighs_i.second[1] != 2
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			case 2:
				if (
					neighs_i.second[0] != 1
					or neighs_i.second[1] != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			case 3:
				if (
					neighs_i.second[0] != 2
					or neighs_i.second[1] != dccrg::error_cell
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 3, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != dccrg::error_cell
				or neighs_i.second[1] != dccrg::error_cell
				or neighs_i.second[4] != dccrg::error_cell
				or neighs_i.second[5] != dccrg::error_cell
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch (neighs_i.first) {
			case 1:
				if (
					neighs_i.second[2] != dccrg::error_cell
					or neighs_i.second[3] != 2
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			case 2:
				if (
					neighs_i.second[2] != 1
					or neighs_i.second[3] != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			case 3:
				if (
					neighs_i.second[2] != 2
					or neighs_i.second[3] != dccrg::error_cell
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 3})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != dccrg::error_cell
				or neighs_i.second[1] != dccrg::error_cell
				or neighs_i.second[2] != dccrg::error_cell
				or neighs_i.second[3] != dccrg::error_cell
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch (neighs_i.first) {
			case 1:
				if (
					neighs_i.second[4] != dccrg::error_cell
					or neighs_i.second[5] != 2
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			case 2:
				if (
					neighs_i.second[4] != 1
					or neighs_i.second[5] != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			case 3:
				if (
					neighs_i.second[4] != 2
					or neighs_i.second[5] != dccrg::error_cell
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
					abort();
				}
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
				abort();
			}
		}}


		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.initialize(comm);
		grid.refine_completely(1);
		grid.stop_refining();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			constexpr auto err = dccrg::error_cell;
			std::array<uint64_t, 6> ref{err, err, err, err, err, err};

			if (neighs_i.first == 2) {
				ref[1] = 3;
				ref[3] = 4;
				ref[5] = 6;
			}
			if (neighs_i.first == 3) {
				ref[0] = 2;
				ref[3] = 5;
				ref[5] = 7;
			}
			if (neighs_i.first == 4) {
				ref[1] = 5;
				ref[2] = 2;
				ref[5] = 8;
			}
			if (neighs_i.first == 5) {
				ref[0] = 4;
				ref[2] = 3;
				ref[5] = 9;
			}
			if (neighs_i.first == 6) {
				ref[1] = 7;
				ref[3] = 8;
				ref[4] = 2;
			}
			if (neighs_i.first == 7) {
				ref[0] = 6;
				ref[3] = 9;
				ref[4] = 3;
			}
			if (neighs_i.first == 8) {
				ref[1] = 9;
				ref[2] = 6;
				ref[4] = 4;
			}
			if (neighs_i.first == 9) {
				ref[0] = 8;
				ref[2] = 7;
				ref[4] = 5;
			}

			if (neighs_i.second != ref) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
				abort();
			}
		}}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
