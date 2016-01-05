/*
Test for scalability of dccrg in 2 D

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
#include "array"
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

	const std::array<uint64_t, 3> grid_length = {{1000, 1000, 1}};
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
		cout << "Maximum refinement level of the grid: " << game_grid.get_maximum_refinement_level() << endl;
		cout << "Number of cells: "
			<< (geom_params.coordinates[0].size() - 1)
				* (geom_params.coordinates[1].size() - 1)
				* (geom_params.coordinates[2].size() - 1)
			<< endl << endl;
	}

	game_grid.balance_load();

	vector<uint64_t> inner_cells = game_grid.get_local_cells_not_on_process_boundary();
	vector<uint64_t> outer_cells = game_grid.get_local_cells_on_process_boundary();
	cout << "Process " << rank
		<< ": number of cells with local neighbors: " << inner_cells.size()
		<< ", number of cells with a remote neighbor: " << outer_cells.size()
		<< endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (vector<uint64_t>::const_iterator
		cell = inner_cells.begin();
		cell != inner_cells.end();
		cell++
	) {
		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->data[1] = 0;

		const std::array<double, 3>
			cell_center = game_grid.geometry.get_center(*cell),
			cell_length = game_grid.geometry.get_length(*cell);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell_data->data[0] = 1;
		} else {
			cell_data->data[0] = 0;
		}
	}
	for (vector<uint64_t>::const_iterator
		cell = outer_cells.begin();
		cell != outer_cells.end();
		cell++
	) {
		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->data[1] = 0;

		const std::array<double, 3>
			cell_center = game_grid.geometry.get_center(*cell),
			cell_length = game_grid.geometry.get_length(*cell);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell_data->data[0] = 1;
		} else {
			cell_data->data[0] = 0;
		}
	}

	if (rank == 0) {
		cout << "step: ";
	}

	MPI_Barrier(comm);

	#define TIME_STEPS 100
	before = time(NULL);
	for (int step = 0; step < TIME_STEPS; step++) {

		if (rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		game_grid.start_remote_neighbor_copy_updates();
		/*
		Get the neighbor counts of every cell, starting with the cells whose neighbor data
		doesn't come from other processes
		*/
		for (vector<uint64_t>::const_iterator
			cell = inner_cells.begin();
			cell != inner_cells.end();
			cell++
		) {
			game_of_life_cell* cell_data = game_grid[*cell];
			cell_data->data[1] = 0;

			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			for (vector<uint64_t>::const_iterator
				neighbor = neighbors->begin();
				neighbor != neighbors->end();
				neighbor++
			) {
				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data->data[0] == 1) {
					cell_data->data[1]++;
				}
			}
		}

		// wait for neighbor data updates to this process to finish and go through the rest of the cells
		game_grid.wait_remote_neighbor_copy_update_receives();
		for (vector<uint64_t>::const_iterator
			cell = outer_cells.begin();
			cell != outer_cells.end();
			cell++
		) {
			game_of_life_cell* cell_data = game_grid[*cell];
			cell_data->data[1] = 0;

			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			for (vector<uint64_t>::const_iterator
				neighbor = neighbors->begin();
				neighbor != neighbors->end();
				neighbor++
			) {
				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data->data[0] == 1) {
					cell_data->data[1]++;
				}
			}
		}
		/*
		Wait for neighbor data updates from this process to finish until
		updating live status of own cells
		*/
		game_grid.wait_remote_neighbor_copy_update_sends();

		// calculate the next turn
		for (vector<uint64_t>::const_iterator
			cell = inner_cells.begin();
			cell != inner_cells.end();
			cell++
		) {
			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->data[1] == 3) {
				cell_data->data[0] = 1;
			} else if (cell_data->data[1] != 2) {
				cell_data->data[0] = 0;
			}
		}
		for (vector<uint64_t>::const_iterator
			cell = outer_cells.begin();
			cell != outer_cells.end();
			cell++
		) {
			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->data[1] == 3) {
				cell_data->data[0] = 1;
			} else if (cell_data->data[1] != 2) {
				cell_data->data[0] = 0;
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
	cout << "Process " << rank
		<< ": " << number_of_cells * TIME_STEPS << " cells processed at the speed of "
		<< double(number_of_cells * TIME_STEPS) / total << " cells / second"
		<< endl;

	return EXIT_SUCCESS;
}

