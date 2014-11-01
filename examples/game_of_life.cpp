/*
A simple 2 D game of life program to demonstrate the efficient parallel usage of dccrg.

Serial performance has not been optimized, for that the get_live_neighbor_counts() and
apply_rules() functions would have to be modified to use cached pointer to cell data
instead of fetching the pointer each time cell data is accessed.
*/

#include "chrono"
#include "cstdlib"
#include "ctime"
#include "tuple"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;

// store in every cell of the grid whether the cell is alive and the number of live neighbors it has
struct game_of_life_cell {
	unsigned int is_alive = 0, live_neighbor_count = 0;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) &(this->is_alive), 1, MPI_UNSIGNED);
	}
};


/*!
Initializes the given cells, all of which must be local
*/
void initialize_game(
	const vector<uint64_t>& cells,
	dccrg::Dccrg<game_of_life_cell>& game_grid
) {
	for (const auto& cell: cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == nullptr) {
			abort();
		}

		cell_data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell_data->is_alive = 1;
		} else {
			cell_data->is_alive = 0;
		}
	}
}


/*!
Calculates the number of live neihgbours for every cell given, all of which must be local
*/
void get_live_neighbor_counts(
	const vector<uint64_t>& cells,
	dccrg::Dccrg<game_of_life_cell>& game_grid
) {
	for (const auto& cell: cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == nullptr) {
			abort();
		}

		cell_data->live_neighbor_count = 0;

		const auto* const neighbors = game_grid.get_neighbors_of(cell);

		for (const auto& neighbor: *neighbors) {
			if (neighbor == dccrg::error_cell) {
				continue;
			}

			const auto* const neighbor_data = game_grid[neighbor];
			if (neighbor_data == nullptr) {
				abort();
			}

			if (neighbor_data->is_alive > 0) {
				cell_data->live_neighbor_count++;
			}
		}
	}
}


/*!
Applies the game of life rules to every given cell, all of which must be local.
*/
void apply_rules(
	const vector<uint64_t>& cells,
	dccrg::Dccrg<game_of_life_cell>& game_grid
) {
	for (const auto& cell: cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == nullptr) {
			abort();
		}

		if (cell_data->live_neighbor_count == 3) {
			cell_data->is_alive = 1;
		} else if (cell_data->live_neighbor_count != 2) {
			cell_data->is_alive = 0;
		}
	}
}


/*!
See the comments in simple_game_of_life.cpp for an explanation of the basics.
*/
int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	if (rank < 0) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}
	MPI_Comm_size(comm, &comm_size);
	if (comm_size < 0) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	dccrg::Dccrg<game_of_life_cell> game_grid;

	const int neighborhood_size = 1, maximum_refinement_level = 0;
	const std::array<uint64_t, 3> grid_length{{500, 500, 1}};
	game_grid.initialize(
		grid_length,
		comm,
		"RCB",
		neighborhood_size,
		maximum_refinement_level
	);

	game_grid.balance_load();

	/*
	Get the cells on this process just once, since
	the grid doesn't change during the game
	To make the game scale better, separate local cells into those
	without even one neighbor on another process and those that do.
	While updating cell data between processes, start calculating
	the next turn for cells which don't have neighbors on other processes
	*/
	const vector<uint64_t>
		inner_cells = game_grid.get_local_cells_not_on_process_boundary(),
		outer_cells = game_grid.get_local_cells_on_process_boundary();

	initialize_game(inner_cells, game_grid);
	initialize_game(outer_cells, game_grid);


	// time the game to examine its scalability
	const auto time_start = chrono::high_resolution_clock::now();

	const int turns = 100;
	for (int turn = 0; turn < turns; turn++) {

		// start updating cell data from other processes
		// and calculate the next turn for cells without
		// neighbors on other processes in the meantime
		game_grid.start_remote_neighbor_copy_updates();
		get_live_neighbor_counts(inner_cells, game_grid);

		// wait for neighbor data updates to finish and the
		// calculate the next turn for rest of the cells on this process
		game_grid.wait_remote_neighbor_copy_updates();
		get_live_neighbor_counts(outer_cells, game_grid);

		// update the state of life for all local cells
		apply_rules(inner_cells, game_grid);
		apply_rules(outer_cells, game_grid);
	}

	const auto time_end
		= chrono::high_resolution_clock::now();
	const auto total_time
		= chrono::duration_cast<
			chrono::duration<double>
		>(time_end - time_start).count();

	// calculate some timing statistics
	double
		total_cells = double(turns) * (inner_cells.size() + outer_cells.size()),
		min_speed_local = total_cells / total_time, min_speed = 0,
		max_speed_local = total_cells / total_time, max_speed = 0,
		avg_speed_local = total_cells / total_time, avg_speed = 0,
		total_global_cells = 0;

	MPI_Reduce(&min_speed_local, &min_speed, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&max_speed_local, &max_speed, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&avg_speed_local, &avg_speed, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	avg_speed /= comm_size;
	MPI_Reduce(&total_cells, &total_global_cells, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

	// print the statistics
	if (rank == 0) {
		cout << "Game played at " << avg_speed
			<< " cells / process / s (minimum: " << min_speed
			<< ", maximum: " << max_speed << ")"
			<< endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

