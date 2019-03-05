/*
A 2 D game of life program to demonstrate the efficient parallel usage of dccrg.

Serial performance has not been optimized, for that see game_of_life_optimized.cpp.
*/

#include "algorithm"
#include "chrono"
#include "cstdlib"
#include "ctime"
#include "tuple"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;

/*
Stores in every cell of the grid whether the cell is
alive as well as the number of live neighbors it has.
*/
struct game_of_life_cell {
	unsigned int is_alive = 0, live_neighbor_count = 0;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) &(this->is_alive), 1, MPI_UNSIGNED);
	}
};


/*!
Initializes game with random live cells.
*/
template<class Cell_Iterator> void initialize_game(const Cell_Iterator& cells) {
	for (const auto& cell: cells) {
		cell.data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell.data->is_alive = 1;
		} else {
			cell.data->is_alive = 0;
		}
	}
}


/*!
Calculates the number of live neihgbours for given local cells.
*/
template<class Cell_Iterator> void get_live_neighbor_counts(
	const Cell_Iterator& cells
) {
	for (const auto& cell: cells) {
		cell.data->live_neighbor_count = 0;
		for (const auto& neighbor: cell.neighbors_of) {
			if (neighbor.data->is_alive > 0) {
				cell.data->live_neighbor_count++;
			}
		}
	}
}


/*!
Applies the game of life rules to given local cells.
*/
template<class Cell_Iterator> void apply_rules(const Cell_Iterator& cells) {
	for (const auto& cell: cells) {
		if (cell.data->live_neighbor_count == 3) {
			cell.data->is_alive = 1;
		} else if (cell.data->live_neighbor_count != 2) {
			cell.data->is_alive = 0;
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

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	dccrg::Dccrg<game_of_life_cell> grid;

	grid
		.set_initial_length({500, 500, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(0)
		.initialize(comm)
		.balance_load();

	initialize_game(grid.local_cells);

	// time the game to examine its scalability
	const auto time_start = chrono::high_resolution_clock::now();

	const int turns = 100;
	for (int turn = 0; turn < turns; turn++) {
		// start updating cell data from other processes
		// and calculate the next turn for cells without
		// neighbors on other processes in the meantime
		grid.start_remote_neighbor_copy_updates();
		get_live_neighbor_counts(grid.inner_cells);

		// wait for neighbor data updates to finish and
		// calculate the next turn for rest of the cells
		// on this process and afterwards apply next turn
		// to cells without remove neighbors
		grid.wait_remote_neighbor_copy_update_receives();
		get_live_neighbor_counts(grid.outer_cells);
		apply_rules(grid.inner_cells);

		// apply next turn for cells with remote neighbors
		// after their state was sent to other processes
		grid.wait_remote_neighbor_copy_update_sends();
		apply_rules(grid.outer_cells);
	}

	const auto time_end
		= chrono::high_resolution_clock::now();
	const auto total_time
		= chrono::duration_cast<
			chrono::duration<double>
		>(time_end - time_start).count();

	// calculate some timing statistics
	const auto
		inner_size = std::distance(grid.inner_cells.begin(), grid.inner_cells.end()),
		outer_size = std::distance(grid.outer_cells.begin(), grid.outer_cells.end());
	double
		total_cells = double(turns) * (inner_size + outer_size),
		min_speed_local = total_cells / total_time, min_speed = 0,
		max_speed_local = total_cells / total_time, max_speed = 0,
		avg_speed_local = total_cells / total_time, avg_speed = 0,
		total_global_cells = 0;

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
