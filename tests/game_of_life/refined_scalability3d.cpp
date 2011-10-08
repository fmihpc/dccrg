/*
Tests the scalability of the grid in 3 D with refined grid
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "functional"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg_arbitrary_geometry.hpp"
#include "../../dccrg.hpp"


struct game_of_life_cell {

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
	}

	bool is_alive;
	unsigned int live_neighbour_count;
};


using namespace std;
using namespace boost;
using namespace boost::mpi;
using namespace dccrg;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	time_t before, after, total = 0;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	Dccrg<game_of_life_cell, ArbitraryGeometry> game_grid;

	#define GRID_SIZE 21
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
		y_coordinates.push_back(i * CELL_SIZE);
		z_coordinates.push_back(i * CELL_SIZE);
	}
	game_grid.set_geometry(x_coordinates, y_coordinates, z_coordinates);

	#define NEIGHBORHOOD_SIZE 1
	game_grid.initialize(comm, "RCB", NEIGHBORHOOD_SIZE);
	game_grid.balance_load();

	vector<uint64_t> cells = game_grid.get_cells();
	// refine random cells until every process has enough cells
	#define MAX_CELLS (100 * GRID_SIZE * GRID_SIZE * GRID_SIZE)
	do {
		cout << "Process " << comm.rank() << ", number of cells: " << cells.size() << endl;
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		for (int i = 0; i < int(cells.size() / 15) && i < 10000; i++) {
			game_grid.refine_completely(cells[i]);
		}
		game_grid.stop_refining();
		cells = game_grid.get_cells();
	} while (all_reduce(comm, int(cells.size()), plus<int>()) < MAX_CELLS);
	game_grid.balance_load();

	vector<uint64_t> cells_with_local_neighbours = game_grid.get_cells_with_local_neighbours();
	vector<uint64_t> cells_with_remote_neighbour = game_grid.get_cells_with_remote_neighbour();
	cout << "Process " << comm.rank() << ": number of cells with local neighbours: " << cells_with_local_neighbours.size() << ", number of cells with a remote neighbour: " << cells_with_remote_neighbour.size() << endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (vector<uint64_t>::const_iterator cell = cells_with_local_neighbours.begin(); cell != cells_with_local_neighbours.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbour_count = 0;

		double y = game_grid.get_cell_y(*cell);
		if (fabs(0.5 + 0.1 * game_grid.get_cell_y_size(*cell) - y) < 0.5 * game_grid.get_cell_y_size(*cell)) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}
	for (vector<uint64_t>::const_iterator cell = cells_with_remote_neighbour.begin(); cell != cells_with_remote_neighbour.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbour_count = 0;

		double y = game_grid.get_cell_y(*cell);
		if (fabs(0.5 + 0.1 * game_grid.get_cell_y_size(*cell) - y) < 0.5 * game_grid.get_cell_y_size(*cell)) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	// get some statistics
	unordered_map<int, int> neighbour_count_histogram;
	for (int i = 0; i < 8 * (9 + 8 + 9); i++) {
		neighbour_count_histogram[i] = 0;
	}
	double avg_neighbours = 0.0;
	int max_neighbours = 0, min_neighbours = 999;
	uint64_t number_of_cells = cells_with_local_neighbours.size() + cells_with_remote_neighbour.size();

	for (vector<uint64_t>::const_iterator cell = cells_with_local_neighbours.begin(); cell != cells_with_local_neighbours.end(); cell++) {

		const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
		avg_neighbours += double(neighbours->size()) / number_of_cells;

		if (max_neighbours < int(neighbours->size())) {
			max_neighbours = neighbours->size();
		}

		if (min_neighbours > int(neighbours->size())) {
			min_neighbours = neighbours->size();
		}

		if (neighbours->size() >= 7 && neighbours->size() <= 8 * (9 + 8 + 9)) {
			neighbour_count_histogram[neighbours->size()] += 1;
		} else {
			cout << "Impossible number of neighbours for cell " << *cell << ": " << neighbours->size() << endl;
		}
	}
	for (vector<uint64_t>::const_iterator cell = cells_with_remote_neighbour.begin(); cell != cells_with_remote_neighbour.end(); cell++) {

		const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
		avg_neighbours += double(neighbours->size()) / number_of_cells;

		if (max_neighbours < int(neighbours->size())) {
			max_neighbours = neighbours->size();
		}

		if (min_neighbours > int(neighbours->size())) {
			min_neighbours = neighbours->size();
		}

		if (neighbours->size() >= 7 && neighbours->size() <= 8 * (9 + 8 + 9)) {
			neighbour_count_histogram[neighbours->size()] += 1;
		} else {
			cout << "Impossible number of neighbours for cell " << *cell << ": " << neighbours->size() << endl;
		}
	}

	// print statistics
	int total_max_neighbours = all_reduce(comm, max_neighbours, plus<int>());
	int total_min_neighbours = all_reduce(comm, min_neighbours, plus<int>());
	double total_avg_neighbours = all_reduce(comm, avg_neighbours, plus<double>());

	for (int i = 0; i <= 8 * (9 + 8 + 9); i++) {
		int total = all_reduce(comm, neighbour_count_histogram[i], plus<int>());
		neighbour_count_histogram[i] = total;
	}

	if (comm.rank() == 0) {
		cout << "Max neighbours: " << total_max_neighbours << endl;
		cout << "Min neighbours: " << total_min_neighbours << endl;
		cout << "Avg. neighbours: " << total_avg_neighbours << endl;
		cout << "Neighbour count histogram: (neighbours, count)" << endl;
		for (int i = 0; i <= 8 * (9 + 8 + 9); i++) {
			cout << i << " " << neighbour_count_histogram[i] << endl;
		}
	}

	if (comm.rank() == 0) {
		cout << "step: ";
	}

	#define TIME_STEPS 100
	before = time(NULL);
	uint64_t processed_neighbours = 0;
	for (int step = 0; step < TIME_STEPS; step++) {

		if (comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		game_grid.start_remote_neighbour_data_update();
		// get the neighbour counts of every cell, starting with the cells whose neighbour data doesn't come from other processes
		for (vector<uint64_t>::const_iterator cell = cells_with_local_neighbours.begin(); cell != cells_with_local_neighbours.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];
			cell_data->live_neighbour_count = 0;

			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++, processed_neighbours++) {

				if (*neighbour == 0) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data->is_alive) {
					cell_data->live_neighbour_count++;
				}
			}
		}

		// wait for neighbour data updates to finish and go through the rest of the cells
		game_grid.wait_neighbour_data_update();
		for (vector<uint64_t>::const_iterator cell = cells_with_remote_neighbour.begin(); cell != cells_with_remote_neighbour.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];
			cell_data->live_neighbour_count = 0;

			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++, processed_neighbours++) {

				if (*neighbour == 0) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data->is_alive) {
					cell_data->live_neighbour_count++;
				}
			}
		}

		// calculate the next turn
		for (vector<uint64_t>::const_iterator cell = cells_with_local_neighbours.begin(); cell != cells_with_local_neighbours.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbour_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbour_count != 2) {
				cell_data->is_alive = false;
			}
		}
		for (vector<uint64_t>::const_iterator cell = cells_with_remote_neighbour.begin(); cell != cells_with_remote_neighbour.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbour_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbour_count != 2) {
				cell_data->is_alive = false;
			}
		}
	}
	after = time(NULL);
	total += after - before;
	if (comm.rank() == 0) {
		cout << endl;
	}
	comm.barrier();
	cout << "processed neighbours: " << processed_neighbours << endl;
	cout << "Process " << comm.rank() << ": " << number_of_cells * TIME_STEPS << " cells processed at the speed of " << double(number_of_cells * TIME_STEPS) / total << " cells / second"<< endl;

	return EXIT_SUCCESS;
}
