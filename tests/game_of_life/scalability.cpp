/*
Test for scalability of dccrg in 2 D

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016, 2018 Finnish Meteorological Institute
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

#include "algorithm"
#include "array"
#include "chrono"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "iostream"
#include "tuple"
#include "unordered_set"

#include "zoltan.h"

#include "dccrg_stretched_cartesian_geometry.hpp"
#include "dccrg.hpp"

// TODO: move this to a separate file
struct game_of_life_cell {

	// data[0] == 1 if cell is alive, data[1] holds the number of live neighbors
	unsigned int data[2];

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->data), 1, MPI_INT);
	}
};


using namespace std;
using namespace std::chrono;
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

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	Dccrg<game_of_life_cell, Stretched_Cartesian_Geometry> grid;

	const std::array<uint64_t, 3> grid_length = {{1000, 1000, 1}};
	const double cell_length = 1.0 / grid_length[0];
	#define NEIGHBORHOOD_SIZE 1
	#define MAX_REFINEMENT_LEVEL 0
	grid.initialize(
		grid_length,
		comm,
		"RCB",
		NEIGHBORHOOD_SIZE,
		MAX_REFINEMENT_LEVEL,
		true, true, false
	);

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	grid.set_geometry(geom_params);

	if (rank == 0) {
		cout << "Maximum refinement level of the grid: " << grid.get_maximum_refinement_level() << endl;
		cout << "Number of cells: "
			<< (geom_params.coordinates[0].size() - 1)
				* (geom_params.coordinates[1].size() - 1)
				* (geom_params.coordinates[2].size() - 1)
			<< endl << endl;
	}

	grid.balance_load();

	cout << "Process " << rank
		<< ": number of cells with local neighbors: "
		<< std::distance(grid.inner_cells.begin(), grid.inner_cells.end())
		<< ", number of cells with a remote neighbor: "
		<< std::distance(grid.outer_cells.begin(), grid.outer_cells.end())
		<< endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (const auto& cell: grid.local_cells) {
		cell.data->data[1] = 0;

		const auto indices = grid.mapping.get_indices(cell.id);
		if (indices[1] == 500) {
			cell.data->data[0] = 1;
		} else {
			cell.data->data[0] = 0;
		}
	}

	if (rank == 0) {
		cout << "step: ";
	}

	MPI_Barrier(comm);

	constexpr size_t TIME_STEPS = 16;
	const auto before = high_resolution_clock::now();
	for (size_t step = 0; step < TIME_STEPS; step++) {

		if (rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		grid.start_remote_neighbor_copy_updates();
		/*
		Get the neighbor counts of every cell, starting with the cells whose neighbor data
		doesn't come from other processes
		*/
		for (const auto& cell: grid.inner_cells) {
			cell.data->data[1] = 0;

			for (const auto& neighbor: cell.neighbors_of) {
				if (neighbor.data->data[0] == 1) {
					cell.data->data[1]++;
				}
			}
		}

		// wait for neighbor data updates to this process to finish and go through the rest of the cells
		grid.wait_remote_neighbor_copy_update_receives();
		for (const auto& cell: grid.outer_cells) {
			cell.data->data[1] = 0;

			for (const auto& neighbor: cell.neighbors_of) {
				if (neighbor.data->data[0] == 1) {
					cell.data->data[1]++;
				}
			}
		}

		// calculate the next turn
		for (const auto& cell: grid.inner_cells) {
			if (cell.data->data[1] == 3) {
				cell.data->data[0] = 1;
			} else if (cell.data->data[1] != 2) {
				cell.data->data[0] = 0;
			}
		}
		/*
		Wait for neighbor data updates from this process to finish until
		updating live status of own cells with remote neighbors
		*/
		grid.wait_remote_neighbor_copy_update_sends();

		for (const auto& cell: grid.outer_cells) {
			if (cell.data->data[1] == 3) {
				cell.data->data[0] = 1;
			} else if (cell.data->data[1] != 2) {
				cell.data->data[0] = 0;
			}
		}
	}
	const auto after = high_resolution_clock::now();
	const auto total = duration_cast<duration<double>>(after - before).count();
	if (rank == 0) {
		cout << endl;
	}
	MPI_Barrier(comm);

	for (const auto& cell: grid.local_cells) {
		const auto indices = grid.mapping.get_indices(cell.id);
		if (indices[1] + TIME_STEPS == 500 or indices[1] - TIME_STEPS == 500) {
			if (cell.data->data[0] == 0) {
				std::cout << __FILE__ "(" << __LINE__ << "): "
					<< indices[0] << ", " << indices[1] << ", " << indices[2]
					<< std::endl;
				abort();
			}
		} else {
			if (cell.data->data[0] == 1) {
				std::cout << __FILE__ "(" << __LINE__ << "): "
					<< indices[0] << ", " << indices[1] << ", " << indices[2]
					<< std::endl;
				abort();
			}
		}
	}

	const auto number_of_cells = std::distance(grid.local_cells.begin(), grid.local_cells.end());
	cout << "Process " << rank
		<< ": " << number_of_cells * TIME_STEPS << " cells processed at the speed of "
		<< double(number_of_cells * TIME_STEPS) / total << " cells / second"
		<< endl;

	MPI_Finalize();

	return EXIT_SUCCESS;
}
