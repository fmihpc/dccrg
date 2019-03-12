/*
Test for scalability of dccrg in 3 D with refined grid

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016, 2018 Finnish Meteorological Institute

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
#include "functional"
#include "iostream"
#include "unordered_map"
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

	Dccrg<game_of_life_cell, Stretched_Cartesian_Geometry> grid;

	const std::array<uint64_t, 3> grid_length = {{21, 21, 21}};
	const double cell_length = 1.0 / grid_length[0];

	#define NEIGHBORHOOD_SIZE 1
	grid
		.set_initial_length(grid_length)
		.set_neighborhood_length(NEIGHBORHOOD_SIZE)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	grid.set_geometry(geom_params);

	grid.balance_load();

	vector<uint64_t> cells;
	for (const auto& cell: grid.local_cells()) {
		cells.push_back(cell.id);
	}
	// refine random cells until every process has enough cells
	#define MAX_CELLS (100 * grid_length[0] * grid_length[1] * grid_length[2])
	unsigned int total_cells = 0;
	do {
		cout << "Process " << rank << ", number of cells: " << cells.size() << endl;
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		for (int i = 0; i < int(cells.size() / 15) && i < 10000; i++) {
			grid.refine_completely(cells[i]);
		}
		grid.stop_refining();
		cells = grid.get_cells();
		unsigned int local_cells = cells.size();
		MPI_Allreduce(&local_cells, &total_cells, 1, MPI_UNSIGNED, MPI_SUM, comm);
	} while (total_cells < MAX_CELLS);
	grid.balance_load();

	cout << "Process " << rank
		<< ": number of cells with local neighbors: "
		<< std::distance(grid.inner_cells().begin(), grid.inner_cells().end())
		<< ", number of cells with a remote neighbor: "
		<< std::distance(grid.outer_cells().begin(), grid.outer_cells().end())
		<< endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (const auto& cell: grid.inner_cells()) {
		cell.data->live_neighbor_count = 0;

		const std::array<double, 3>
			cell_center = grid.geometry.get_center(cell.id),
			cell_length = grid.geometry.get_length(cell.id);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell.data->is_alive = 1;
		} else {
			cell.data->is_alive = 0;
		}
	}
	for (const auto& cell: grid.outer_cells()) {
		cell.data->live_neighbor_count = 0;

		const std::array<double, 3>
			cell_center = grid.geometry.get_center(cell.id),
			cell_length = grid.geometry.get_length(cell.id);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell.data->is_alive = 1;
		} else {
			cell.data->is_alive = 0;
		}
	}

	// get some statistics
	unordered_map<int, int> neighbor_count_histogram;
	for (int i = 0; i < 8 * (9 + 8 + 9); i++) {
		neighbor_count_histogram[i] = 0;
	}
	double avg_neighbors = 0.0;
	int max_neighbors = 0, min_neighbors = 999;
	const auto number_of_cells = std::distance(grid.local_cells().begin(), grid.local_cells().end());

	for (const auto& cell: grid.local_cells()) {
		const double neighbors_size = std::distance(cell.neighbors_of.begin(), cell.neighbors_of.end());
		avg_neighbors += neighbors_size / number_of_cells;

		if (max_neighbors < int(neighbors_size)) {
			max_neighbors = neighbors_size;
		}

		if (min_neighbors > int(neighbors_size)) {
			min_neighbors = neighbors_size;
		}

		if (neighbors_size >= 7 && neighbors_size <= 8 * (9 + 8 + 9)) {
			neighbor_count_histogram[neighbors_size] += 1;
		} else {
			cout << "Impossible number of neighbors for cell " << cell.id
				<< ": " << neighbors_size
				<< endl;
		}
	}

	// print statistics
	int total_max_neighbors = 0;
	MPI_Allreduce(&max_neighbors, &total_max_neighbors, 1, MPI_INT, MPI_SUM, comm);
	int total_min_neighbors = 0;
	MPI_Allreduce(&min_neighbors, &total_min_neighbors, 1, MPI_INT, MPI_SUM, comm);
	double total_avg_neighbors = 0;
	MPI_Allreduce(&avg_neighbors, &total_avg_neighbors, 1, MPI_DOUBLE, MPI_SUM, comm);

	for (int i = 0; i <= 8 * (9 + 8 + 9); i++) {
		int total = 0;
		MPI_Allreduce(&(neighbor_count_histogram[i]), &total, 1, MPI_INT, MPI_SUM, comm);
		neighbor_count_histogram[i] = total;
	}

	if (rank == 0) {
		cout << "Max neighbors: " << total_max_neighbors << endl;
		cout << "Min neighbors: " << total_min_neighbors << endl;
		cout << "Avg. neighbors: " << total_avg_neighbors << endl;
		cout << "Neighbour count histogram: (neighbors, count)" << endl;
		for (int i = 0; i <= 8 * (9 + 8 + 9); i++) {
			cout << i << " " << neighbor_count_histogram[i] << endl;
		}
	}

	if (rank == 0) {
		cout << "step: ";
	}

	#define TIME_STEPS 100
	const auto before = high_resolution_clock::now();
	uint64_t processed_neighbors = 0;
	for (int step = 0; step < TIME_STEPS; step++) {

		if (rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		grid.start_remote_neighbor_copy_updates();
		// get the neighbor counts of every cell, starting with the
		// cells whose neighbor data doesn't come from other processes
		for (const auto& cell: grid.inner_cells()) {
			cell.data->live_neighbor_count = 0;
			for (const auto& neighbor: cell.neighbors_of) {
				if (neighbor.data->is_alive) {
					cell.data->live_neighbor_count++;
				}
				processed_neighbors++;
			}
		}

		// wait for neighbor data updates to finish and go through the rest of the cells
		grid.wait_remote_neighbor_copy_updates();
		for (const auto& cell: grid.outer_cells()) {
			cell.data->live_neighbor_count = 0;
			for (const auto& neighbor: cell.neighbors_of) {
				if (neighbor.data->is_alive) {
					cell.data->live_neighbor_count++;
				}
				processed_neighbors++;
			}
		}

		// calculate the next turn
		for (const auto& cell: grid.inner_cells()) {
			if (cell.data->live_neighbor_count == 3) {
				cell.data->is_alive = 1;
			} else if (cell.data->live_neighbor_count != 2) {
				cell.data->is_alive = 0;
			}
		}
		for (const auto& cell: grid.outer_cells()) {
			if (cell.data->live_neighbor_count == 3) {
				cell.data->is_alive = 1;
			} else if (cell.data->live_neighbor_count != 2) {
				cell.data->is_alive = 0;
			}
		}
	}
	const auto after = high_resolution_clock::now();
	const auto total = duration_cast<duration<double>>(after - before).count();
	if (rank == 0) {
		cout << endl;
	}
	MPI_Barrier(comm);

	cout << "processed neighbors: " << processed_neighbors << endl;
	cout << "Process " << rank << ": " << number_of_cells * TIME_STEPS
		<< " cells processed at the speed of " << double(number_of_cells * TIME_STEPS) / total
		<< " cells / second"
		<< endl;

	MPI_Finalize();

	return EXIT_SUCCESS;
}
