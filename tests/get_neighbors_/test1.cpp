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

	constexpr auto err = dccrg::error_cell;

	for (int neigh_len: {0, 1, 2, 3}) {
		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != err
				or neighs_i.second[1] != err
				or neighs_i.second[2] != err
				or neighs_i.second[3] != err
				or neighs_i.second[4] != err
				or neighs_i.second[5] != err
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
				neighs_i.second[0] != err
				or neighs_i.second[1] != err
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
				or neighs_i.second[2] != err
				or neighs_i.second[3] != err
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
				or neighs_i.second[4] != err
				or neighs_i.second[5] != err
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}


		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({3, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[2] != err
				or neighs_i.second[3] != err
				or neighs_i.second[4] != err
				or neighs_i.second[5] != err
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch (neighs_i.first) {
			case 1:
				if (
					neighs_i.second[0] != err
					or neighs_i.second[1] != 2
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.second[0] << ", " << neighs_i.second[1] << std::endl;
					abort();
				}
				break;
			case 2:
				if (
					neighs_i.second[0] != 1
					or neighs_i.second[1] != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.second[0] << ", " << neighs_i.second[1] << std::endl;
					abort();
				}
				break;
			case 3:
				if (
					neighs_i.second[0] != 2
					or neighs_i.second[1] != err
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.second[0] << ", " << neighs_i.second[1] << std::endl;
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
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != err
				or neighs_i.second[1] != err
				or neighs_i.second[4] != err
				or neighs_i.second[5] != err
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch (neighs_i.first) {
			case 1:
				if (
					neighs_i.second[2] != err
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
					or neighs_i.second[3] != err
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
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			if (
				neighs_i.second[0] != err
				or neighs_i.second[1] != err
				or neighs_i.second[2] != err
				or neighs_i.second[3] != err
			) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch (neighs_i.first) {
			case 1:
				if (
					neighs_i.second[4] != err
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
					or neighs_i.second[5] != err
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
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			std::array<uint64_t, 6> ref{err, err, err, err, err, err};

			const auto get_cell = [&](uint64_t x, uint64_t y, uint64_t z){
				return grid.mapping.get_cell_from_indices({x, y, z}, 1);
			};
			if (neighs_i.first == 2) { // index: 0, 0, 0
				ref[1] = get_cell(1, 0, 0);
				ref[3] = get_cell(0, 1, 0);
				ref[5] = get_cell(0, 0, 1);
			}
			if (neighs_i.first == 3) { // 1, 0, 0
				ref[0] = get_cell(0, 0, 0);
				ref[3] = get_cell(1, 1, 0);
				ref[5] = get_cell(1, 0, 1);
			}
			if (neighs_i.first == 4) { // 0, 1, 0
				ref[1] = get_cell(1, 1, 0);
				ref[2] = get_cell(0, 0, 0);
				ref[5] = get_cell(0, 1, 1);
			}
			if (neighs_i.first == 5) { // 1, 1, 0
				ref[0] = get_cell(0, 1, 0);
				ref[2] = get_cell(1, 0, 0);
				ref[5] = get_cell(1, 1, 1);
			}
			if (neighs_i.first == 6) { // 0, 0, 1
				ref[1] = get_cell(1, 0, 1);
				ref[3] = get_cell(0, 1, 1);
				ref[4] = get_cell(0, 0, 0);
			}
			if (neighs_i.first == 7) { // 1, 0, 1
				ref[0] = get_cell(0, 0, 1);
				ref[3] = get_cell(1, 1, 1);
				ref[4] = get_cell(1, 0, 0);
			}
			if (neighs_i.first == 8) { // 0, 1, 1
				ref[1] = get_cell(1, 1, 1);
				ref[2] = get_cell(0, 0, 1);
				ref[4] = get_cell(0, 1, 0);
			}
			if (neighs_i.first == 9) { // 1, 1, 1
				ref[0] = get_cell(0, 1, 1);
				ref[2] = get_cell(1, 0, 1);
				ref[4] = get_cell(1, 1, 0);
			}

			if (neighs_i.second != ref) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neighs_i.first << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 1})
			.set_periodic(true, true, true)
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			std::array<uint64_t, 6> ref{err, err, err, err, err, err};

			const auto get_cell = [&](uint64_t x, uint64_t y, uint64_t z){
				return grid.mapping.get_cell_from_indices({x, y, z}, 1);
			};
			if (neighs_i.first == 2) { // 0, 0, 0
				ref[0] = get_cell(1, 0, 0);
				ref[1] = get_cell(1, 0, 0);
				ref[2] = get_cell(0, 1, 0);
				ref[3] = get_cell(0, 1, 0);
				ref[4] = get_cell(0, 0, 1);
				ref[5] = get_cell(0, 0, 1);
			}
			if (neighs_i.first == 3) { // 1, 0, 0
				ref[0] = get_cell(0, 0, 0);
				ref[1] = get_cell(0, 0, 0);
				ref[2] = get_cell(1, 1, 0);
				ref[3] = get_cell(1, 1, 0);
				ref[4] = get_cell(1, 0, 1);
				ref[5] = get_cell(1, 0, 1);
			}
			if (neighs_i.first == 4) { // 0, 1, 0
				ref[0] = get_cell(1, 1, 0);
				ref[1] = get_cell(1, 1, 0);
				ref[2] = get_cell(0, 0, 0);
				ref[3] = get_cell(0, 0, 0);
				ref[4] = get_cell(0, 1, 1);
				ref[5] = get_cell(0, 1, 1);
			}
			if (neighs_i.first == 5) { // 1, 1, 0
				ref[0] = get_cell(0, 1, 0);
				ref[1] = get_cell(0, 1, 0);
				ref[2] = get_cell(1, 0, 0);
				ref[3] = get_cell(1, 0, 0);
				ref[4] = get_cell(1, 1, 1);
				ref[5] = get_cell(1, 1, 1);
			}
			if (neighs_i.first == 6) { // 0, 0, 1
				ref[0] = get_cell(1, 0, 1);
				ref[1] = get_cell(1, 0, 1);
				ref[2] = get_cell(0, 1, 1);
				ref[3] = get_cell(0, 1, 1);
				ref[4] = get_cell(0, 0, 0);
				ref[5] = get_cell(0, 0, 0);
			}
			if (neighs_i.first == 7) { // 1, 0, 1
				ref[0] = get_cell(0, 0, 1);
				ref[1] = get_cell(0, 0, 1);
				ref[2] = get_cell(1, 1, 1);
				ref[3] = get_cell(1, 1, 1);
				ref[4] = get_cell(1, 0, 0);
				ref[5] = get_cell(1, 0, 0);
			}
			if (neighs_i.first == 8) { // 0, 1, 1
				ref[0] = get_cell(1, 1, 1);
				ref[1] = get_cell(1, 1, 1);
				ref[2] = get_cell(0, 0, 1);
				ref[3] = get_cell(0, 0, 1);
				ref[4] = get_cell(0, 1, 0);
				ref[5] = get_cell(0, 1, 0);
			}
			if (neighs_i.first == 9) { // 1, 1, 1
				ref[0] = get_cell(0, 1, 1);
				ref[1] = get_cell(0, 1, 1);
				ref[2] = get_cell(1, 0, 1);
				ref[3] = get_cell(1, 0, 1);
				ref[4] = get_cell(1, 1, 0);
				ref[5] = get_cell(1, 1, 0);
			}

			if (neighs_i.first != 1 and neighs_i.second != ref) {
				std::cerr << __FILE__ "(" << __LINE__ << ") "
				<< neighs_i.first << ": "
				<< ref[0] << "," << ref[1] << "," << ref[2] << ","
				<< ref[3] << "," << ref[4] << "," << ref[5] << " != "
				<< neighs_i.second[0] << "," << neighs_i.second[1] << ","
				<< neighs_i.second[2] << "," << neighs_i.second[3] << ","
				<< neighs_i.second[4] << "," << neighs_i.second[5]
				<< std::endl;
				abort();
			}
		}}


		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({2, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			std::array<uint64_t, 6> ref{err, err, err, err, err, err};

			const auto get_cell = [&](uint64_t x, uint64_t y, uint64_t z){
				return grid.get_existing_cell({x, y, z}, 0, 1);
			};
			switch (neighs_i.first) {
			case 1:
				break;
			case 2: // 2, 0, 0
				ref[0] = get_cell(1, 0, 0);
				break;
			case 3: // 0, 0, 0
				ref[1] = get_cell(1, 0, 0);
				ref[3] = get_cell(0, 1, 0);
				ref[5] = get_cell(0, 0, 1);
				break;
			case 4: // 1, 0, 0
				ref[0] = get_cell(0, 0, 0);
				ref[1] = get_cell(2, 0, 0);
				ref[3] = get_cell(1, 1, 0);
				ref[5] = get_cell(1, 0, 1);
				break;
			case 7: // 0, 1, 0
				ref[1] = get_cell(1, 1, 0);
				ref[2] = get_cell(0, 0, 0);
				ref[5] = get_cell(0, 1, 1);
				break;
			case 8: // 1, 1, 0
				ref[0] = get_cell(0, 1, 0);
				ref[1] = get_cell(2, 1, 0);
				ref[2] = get_cell(1, 0, 0);
				ref[5] = get_cell(1, 1, 1);
				break;
			case 11: // 0, 0, 1
				ref[1] = get_cell(1, 0, 1);
				ref[3] = get_cell(0, 1, 1);
				ref[4] = get_cell(0, 0, 0);
				break;
			case 12: // 1, 0, 1
				ref[0] = get_cell(0, 0, 1);
				ref[1] = get_cell(2, 0, 1);
				ref[3] = get_cell(1, 1, 1);
				ref[4] = get_cell(1, 0, 0);
				break;
			case 15: // 0, 1, 1
				ref[1] = get_cell(1, 1, 1);
				ref[2] = get_cell(0, 0, 1);
				ref[4] = get_cell(0, 1, 0);
				break;
			case 16: // 1, 1, 1
				ref[0] = get_cell(0, 1, 1);
				ref[1] = get_cell(2, 1, 1);
				ref[2] = get_cell(1, 0, 1);
				ref[4] = get_cell(1, 1, 0);
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			if (neighs_i.first != 1 and neighs_i.second != ref) {
				std::cerr << __FILE__ "(" << __LINE__ << ") "
				<< neighs_i.first << ": "
				<< ref[0] << "," << ref[1] << "," << ref[2] << ","
				<< ref[3] << "," << ref[4] << "," << ref[5] << " != "
				<< neighs_i.second[0] << "," << neighs_i.second[1] << ","
				<< neighs_i.second[2] << "," << neighs_i.second[3] << ","
				<< neighs_i.second[4] << "," << neighs_i.second[5]
				<< std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 2, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		const auto& neighs = grid.get_neighbors_();
		for (const auto& neighs_i: neighs) {
			std::array<uint64_t, 6> ref{err, err, err, err, err, err};

			const auto get_cell = [&](uint64_t x, uint64_t y, uint64_t z){
				return grid.get_existing_cell({x, y, z}, 0, 1);
			};
			switch (neighs_i.first) {
			case 1:
				break;
			case 2: // 0, 2, 0
				ref[2] = get_cell(0, 1, 0);
				break;
			case 3: // 0, 0, 0
				ref[1] = get_cell(1, 0, 0);
				ref[3] = get_cell(0, 1, 0);
				ref[5] = get_cell(0, 0, 1);
				break;
			case 4: // 1, 0, 0
				ref[0] = get_cell(0, 0, 0);
				ref[3] = get_cell(1, 1, 0);
				ref[5] = get_cell(1, 0, 1);
				break;
			case 5: // 0, 1, 0
				ref[1] = get_cell(1, 1, 0);
				ref[2] = get_cell(0, 0, 0);
				ref[3] = get_cell(0, 2, 0);
				ref[5] = get_cell(0, 1, 1);
				break;
			case 6: // 1, 1, 0
				ref[0] = get_cell(0, 1, 0);
				ref[2] = get_cell(1, 0, 0);
				ref[3] = get_cell(1, 2, 0);
				ref[5] = get_cell(1, 1, 1);
				break;
			case 11: // 0, 0, 1
				ref[1] = get_cell(1, 0, 1);
				ref[3] = get_cell(0, 1, 1);
				ref[4] = get_cell(0, 0, 0);
				break;
			case 12: // 1, 0, 1
				ref[0] = get_cell(0, 0, 1);
				ref[3] = get_cell(1, 1, 1);
				ref[4] = get_cell(1, 0, 0);
				break;
			case 13: // 0, 1, 1
				ref[1] = get_cell(1, 1, 1);
				ref[2] = get_cell(0, 0, 1);
				ref[3] = get_cell(0, 2, 1);
				ref[4] = get_cell(0, 1, 0);
				break;
			case 14: // 1, 1, 1
				ref[0] = get_cell(0, 1, 1);
				ref[2] = get_cell(1, 0, 1);
				ref[3] = get_cell(1, 2, 1);
				ref[4] = get_cell(1, 1, 0);
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
				abort();
			}

			if (neighs_i.first != 1 and neighs_i.second != ref) {
				std::cerr << __FILE__ "(" << __LINE__ << ") "
				<< neighs_i.first << ": "
				<< ref[0] << "," << ref[1] << "," << ref[2] << ","
				<< ref[3] << "," << ref[4] << "," << ref[5] << " != "
				<< neighs_i.second[0] << "," << neighs_i.second[1] << ","
				<< neighs_i.second[2] << "," << neighs_i.second[3] << ","
				<< neighs_i.second[4] << "," << neighs_i.second[5]
				<< std::endl;
				abort();
			}
		}}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
