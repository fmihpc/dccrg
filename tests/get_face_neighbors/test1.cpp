/*
Program for testing get_face_neighbors() function of dccrg.

Copyright 2020 Finnish Meteorological Institute

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

//	constexpr auto err = dccrg::error_cell;

	for (int neigh_len: {0, 1, 2, 3}) {
		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.initialize(comm);
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() > 0) {
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
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 6) {
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
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 4) {
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
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 4) {
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
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 4) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}


		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({2, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.initialize(comm);
		if (grid.get_comm_size() > 1) {
			grid.pin(1, 0);
			grid.pin(2, 1);
		}
		grid.balance_load();

		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 1) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch(cell.id) {
			case 1:
				if (neighs[0].first != 2) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
					abort();
				}
				break;
			case 2:
				if (neighs[0].first != 1) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
					abort();
				}
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 2, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.initialize(comm);
		if (grid.get_comm_size() > 1) {
			grid.pin(1, 0);
			grid.pin(2, 1);
		}
		grid.balance_load();

		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 1) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch(cell.id) {
			case 1:
				if (neighs[0].first != 2) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
					abort();
				}
				break;
			case 2:
				if (neighs[0].first != 1) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
					abort();
				}
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 2})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(0)
			.set_load_balancing_method("RANDOM")
			.initialize(comm);
		if (grid.get_comm_size() > 1) {
			grid.pin(1, 0);
			grid.pin(2, 1);
		}
		grid.balance_load();

		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 1) {
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
				abort();
			}

			switch(cell.id) {
			case 1:
				if (neighs[0].first != 2) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
					abort();
				}
				break;
			case 2:
				if (neighs[0].first != 1) {
					std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
					abort();
				}
				break;
			default:
				std::cerr << __FILE__ "(" << __LINE__ << "): " << neigh_len << std::endl;
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
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 3) {
				std::cerr << __FILE__ "(" << __LINE__ << "): "
					<< neigh_len << ", " << cell.id << std::endl;
				abort();
			}
		}
		grid.balance_load();
		for (const auto& cell: grid.all_cells()) {
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (neighs.size() != 3) {
				std::cerr << __FILE__ "(" << __LINE__ << "): "
					<< neigh_len << ", " << cell.id << std::endl;
				abort();
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({3, 1, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm);
		if (grid.get_comm_size() > 1) {
			grid.pin(1, 0);
			grid.pin(2, 0);
			grid.pin(3, 1);
		}
		grid.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		if (grid.get_comm_size() > 1) {
			grid.pin(7, 1);
		}
		grid.balance_load();
		for (const auto& cell: grid.all_cells()) {
			if (cell.id != 7) {
				continue;
			}
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (grid.get_rank() == 0) {
				if (
					grid.get_comm_size() > 1
					and neigh_len <= 1
					and neighs.size() != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): "
						<< neigh_len << ", " << cell.id << std::endl;
					abort();
				} else if (neighs.size() != 4) {
					std::cerr << __FILE__ "(" << __LINE__ << "): "
						<< neigh_len << ", " << cell.id << std::endl;
					abort();
				}
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 3, 1})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm);
		if (grid.get_comm_size() > 1) {
			grid.pin(1, 0);
			grid.pin(2, 0);
			grid.pin(3, 1);
		}
		grid.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		if (grid.get_comm_size() > 1) {
			grid.pin(10, 1);
		}
		grid.balance_load();
		for (const auto& cell: grid.all_cells()) {
			if (cell.id != 10) {
				continue;
			}
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (grid.get_rank() == 0) {
				if (
					grid.get_comm_size() > 1
					and neigh_len <= 1
					and neighs.size() != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): "
						<< neigh_len << ", " << cell.id << std::endl;
					abort();
				} else if (neighs.size() != 4) {
					std::cerr << __FILE__ "(" << __LINE__ << "): "
						<< neigh_len << ", " << cell.id << std::endl;
					abort();
				}
			}
		}}

		{dccrg::Dccrg<int> grid; grid
			.set_initial_length({1, 1, 3})
			.set_neighborhood_length(neigh_len)
			.set_maximum_refinement_level(1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm);
		if (grid.get_comm_size() > 1) {
			grid.pin(1, 0);
			grid.pin(2, 0);
			grid.pin(3, 1);
		}
		grid.balance_load();
		grid.refine_completely(1);
		grid.stop_refining();
		if (grid.get_comm_size() > 1) {
			grid.pin(16, 1);
		}
		grid.balance_load();
		for (const auto& cell: grid.all_cells()) {
			if (cell.id != 16) {
				continue;
			}
			const auto& neighs = grid.get_face_neighbors_of(cell.id);
			if (grid.get_rank() == 0) {
				if (
					grid.get_comm_size() > 1
					and neigh_len <= 1
					and neighs.size() != 3
				) {
					std::cerr << __FILE__ "(" << __LINE__ << "): "
						<< neigh_len << ", " << cell.id << std::endl;
					abort();
				} else if (neighs.size() != 4) {
					std::cerr << __FILE__ "(" << __LINE__ << "): "
						<< neigh_len << ", " << cell.id << std::endl;
					abort();
				}
			}
		}}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
