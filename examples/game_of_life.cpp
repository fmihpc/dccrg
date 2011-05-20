/*
A simple 2 D game of life program to demonstrate the efficient usage of dccrg
*/

#include "boost/mpi.hpp"
#include "cstdlib"
#include "ctime"
#include "vector"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;
using namespace boost::mpi;


// store in every cell of the grid whether the cell is alive and the number of live neighbours it has
struct game_of_life_cell {

	// boost requires this from user data
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
		/* live_neighbour_count from neighbouring cells is not used
		ar & live_neighbour_count;*/
	}

	bool is_alive;
	unsigned int live_neighbour_count;
};


/*!
Initializes the given cells, all of which must be local
*/
void initialize_game(const vector<uint64_t>* cells, dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];
		cell_data->live_neighbour_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}
}


/*!
Calculates the number of live neihgbours for every cell given, all of which must be local
*/
void get_live_neighbour_counts(const vector<uint64_t>* cells, dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];

		cell_data->live_neighbour_count = 0;
		const vector<uint64_t>* neighbours = game_grid->get_neighbours(*cell);

		for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++) {
			if (*neighbour == 0) {
				continue;
			}

			game_of_life_cell* neighbour_data = (*game_grid)[*neighbour];
			if (neighbour_data->is_alive) {
				cell_data->live_neighbour_count++;
			}
		}
	}
}


/*!
Applies the game of life rules to every given cell, all of which must be local
*/
void apply_rules(const vector<uint64_t>* cells, dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];

		if (cell_data->live_neighbour_count == 3) {
			cell_data->is_alive = true;
		} else if (cell_data->live_neighbour_count != 2) {
			cell_data->is_alive = false;
		}
	}
}


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	// create the grid
	#define GRID_X_SIZE 1000	// in unrefined cells
	#define GRID_Y_SIZE 1000
	#define GRID_Z_SIZE 1
	#define CELL_SIZE 1.0
	#define STENCIL_SIZE 1	// the cells that share a vertex are considered neighbours
	#define MAX_REFINEMENT_LEVEL 0
	dccrg<game_of_life_cell> game_grid(comm, "RCB", 0, 0, 0, CELL_SIZE, CELL_SIZE, CELL_SIZE, GRID_X_SIZE, GRID_Y_SIZE, GRID_Z_SIZE, STENCIL_SIZE, MAX_REFINEMENT_LEVEL);	// use the recursive coordinate bisection method for load balancing (http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg_rcb.html)

	// since the grid doesn't change (isn't refined / unrefined) during the game, workload can be balanced just once in the beginning
	game_grid.balance_load();

	/*
	Get the cells on this process just once, since the grid doesn't change during the game
	To make the game scale better, separate local cells into those without even one neighbour on another process and those that do.
	While updating cell data between processes, start calculating the next turn for cells which don't have neighbours on other processes
	*/
	vector<uint64_t> cells_with_local_neighbours = game_grid.get_cells_with_local_neighbours();
	vector<uint64_t> cells_with_remote_neighbour = game_grid.get_cells_with_remote_neighbour();

	initialize_game(&cells_with_local_neighbours, &game_grid);
	initialize_game(&cells_with_remote_neighbour, &game_grid);


	// time the game to examine its scalability
	time_t before = time(NULL);
	#define TURNS 100
	for (int turn = 0; turn < TURNS; turn++) {

		// start updating cell data from other processes and calculate the next turn for cells without neighbours on other processes in the meantime
		game_grid.start_remote_neighbour_data_update();
		get_live_neighbour_counts(&cells_with_local_neighbours, &game_grid);

		// wait for neighbour data updates to finish and the calculate the next turn for rest of the cells on this process
		game_grid.wait_neighbour_data_update();
		get_live_neighbour_counts(&cells_with_remote_neighbour, &game_grid);

		// update the state of life for all local cells
		apply_rules(&cells_with_local_neighbours, &game_grid);
		apply_rules(&cells_with_remote_neighbour, &game_grid);
	}
	time_t after = time(NULL);


	// calculate some timing statistics
	time_t total_time = after - before;
	uint64_t total_cells = TURNS * (cells_with_local_neighbours.size() + cells_with_remote_neighbour.size());

	double min_speed = all_reduce(comm, total_cells / total_time, minimum<double>());
	double max_speed = all_reduce(comm, total_cells / total_time, maximum<double>());
	double avg_speed = all_reduce(comm, total_cells / total_time, plus<double>()) / comm.size();

	uint64_t total_global_cells = all_reduce(comm, total_cells, plus<uint64_t>());
	double avg_global_speed = all_reduce(comm, total_global_cells / (all_reduce(comm, total_time, plus<double>()) / comm.size()), plus<double>()) / comm.size();

	// print the statistics
	if (comm.rank() == 0) {
		cout << "Game played at " << avg_speed << " cells / process / s (average speed, minimum: " << min_speed << ", maximum: " << max_speed << ")" << endl;
		cout << "Average total playing speed " << avg_global_speed << " cells / s" << endl;
	}

	return game_grid[cells_with_local_neighbours[0]]->is_alive;
}
