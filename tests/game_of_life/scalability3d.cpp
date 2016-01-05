/*
Test for scalability of dccrg in 3 D

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute

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

#include "algorithm"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "iostream"
#include "unordered_set"

#include "mpi.h"
#include "zoltan.h"

#include "dccrg_stretched_cartesian_geometry.hpp"
#include "dccrg.hpp"


struct game_of_life_cell {
	unsigned int is_alive, live_neighbor_count;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->is_alive), 1, MPI_UNSIGNED);
	}
};


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

	time_t before, after, total = 0;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	Dccrg<game_of_life_cell, Stretched_Cartesian_Geometry> game_grid;

	const std::array<uint64_t, 3> grid_length = {{100, 100, 100}};
	const double cell_length = 1.0 / grid_length[0];

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	game_grid.set_geometry(geom_params);

	#define NEIGHBORHOOD_SIZE 1
	#define MAX_REFINEMENT_LEVEL 0
	game_grid.initialize(grid_length, comm, "RCB", NEIGHBORHOOD_SIZE, MAX_REFINEMENT_LEVEL);
	if (rank == 0) {
		cout << "Maximum refinement level of the grid: "
			<< game_grid.get_maximum_refinement_level()
			<< "\nNumber of cells: "
			<< (geom_params.coordinates[0].size() - 1)
				* (geom_params.coordinates[1].size() - 1)
				* (geom_params.coordinates[2].size() - 1)
			<< endl << endl;
	}

	game_grid.balance_load();
	MPI_Barrier(comm);

	auto inner_cells = game_grid.get_local_cells_not_on_process_boundary();
	auto outer_cells = game_grid.get_local_cells_on_process_boundary();

	// initialize the game with random cells alive
	srand(time(NULL));
	for (const auto& cell: inner_cells) {

		auto* const cell_data = game_grid[cell];
		cell_data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX <= 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}
	for (const auto& cell: outer_cells) {

		auto* const cell_data = game_grid[cell];
		cell_data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX <= 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	if (rank == 0) {
		cout << "step: ";
	}

	#define TIME_STEPS 100
	before = time(NULL);
	for (int step = 0; step < TIME_STEPS; step++) {
		MPI_Barrier(comm);

		if (step % 10 == 0 && rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		game_grid.start_remote_neighbor_copy_updates();

		/*
		Get the neighbor counts of every cell, starting with the cells whose neighbor
		data doesn't come from other processes while those transfers are carried out.
		*/
		for (const auto& cell: inner_cells) {

			auto* const cell_data = game_grid[cell];
			if (cell_data == NULL) { abort(); }

			cell_data->live_neighbor_count = 0;

			const auto* const neighbors = game_grid.get_neighbors_of(cell);
			for (const auto& neighbor: *neighbors) {

				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = game_grid[neighbor];
				if (neighbor_data == NULL) { abort(); }

				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// wait for neighbor data updates to finish and go through the rest of the cells
		game_grid.wait_remote_neighbor_copy_updates();
		for (const auto& cell: outer_cells) {

			auto* const cell_data = game_grid[cell];
			if (cell_data == NULL) { abort(); }

			cell_data->live_neighbor_count = 0;

			const auto* const neighbors = game_grid.get_neighbors_of(cell);
			for (const auto& neighbor: *neighbors) {

				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = game_grid[neighbor];
				if (neighbor_data == NULL) { abort(); }

				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// calculate the next turn
		for (const auto& cell: inner_cells) {

			auto* const cell_data = game_grid[cell];
			if (cell_data == NULL) { abort(); }

			if (cell_data->live_neighbor_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbor_count != 2) {
				cell_data->is_alive = false;
			}
		}
		for (const auto& cell: outer_cells) {

			auto* const cell_data = game_grid[cell];
			if (cell_data == NULL) { abort(); }

			if (cell_data->live_neighbor_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbor_count != 2) {
				cell_data->is_alive = false;
			}
		}
	}
	after = time(NULL);
	total += after - before;
	if (rank == 0) {
		cout << endl;
	}
	MPI_Barrier(comm);

	int number_of_cells = inner_cells.size() + outer_cells.size();
	cout << "Process " << rank << ": "
		<< number_of_cells * TIME_STEPS
		<< " cells processed at the speed of "
		<< double(number_of_cells * TIME_STEPS) / total << " cells / second"
		<< endl;

	MPI_Finalize();

	return EXIT_SUCCESS;
}
