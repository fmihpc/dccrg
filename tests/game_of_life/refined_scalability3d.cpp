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

#include "../../dccrg_stretched_cartesian_geometry.hpp"
#include "../../dccrg.hpp"


struct game_of_life_cell {

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
	}

	bool is_alive;
	unsigned int live_neighbor_count;
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

	Dccrg<game_of_life_cell, Stretched_Cartesian_Geometry> game_grid;

	const boost::array<uint64_t, 3> grid_length = {{21, 21, 21}};
	const double cell_length = 1.0 / grid_length[0];

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	game_grid.set_geometry(geom_params);

	#define NEIGHBORHOOD_SIZE 1
	game_grid.initialize(grid_length, comm, "RCB", NEIGHBORHOOD_SIZE);

	game_grid.balance_load();

	vector<uint64_t> cells = game_grid.get_cells();
	// refine random cells until every process has enough cells
	#define MAX_CELLS (100 * grid_length[0] * grid_length[1] * grid_length[2])
	do {
		cout << "Process " << comm.rank() << ", number of cells: " << cells.size() << endl;
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		for (int i = 0; i < int(cells.size() / 15) && i < 10000; i++) {
			game_grid.refine_completely(cells[i]);
		}
		game_grid.stop_refining();
		cells = game_grid.get_cells();
	} while (all_reduce(comm, int(cells.size()), plus<int>()) < int(MAX_CELLS));
	game_grid.balance_load();

	vector<uint64_t> inner_cells = game_grid.get_local_cells_not_on_process_boundary();
	vector<uint64_t> outer_cells = game_grid.get_local_cells_on_process_boundary();
	cout << "Process " << comm.rank() << ": number of cells with local neighbors: " << inner_cells.size() << ", number of cells with a remote neighbor: " << outer_cells.size() << endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (vector<uint64_t>::const_iterator cell = inner_cells.begin(); cell != inner_cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbor_count = 0;

		const boost::array<double, 3>
			cell_center = game_grid.geometry.get_center(*cell),
			cell_length = game_grid.geometry.get_length(*cell);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}
	for (vector<uint64_t>::const_iterator cell = outer_cells.begin(); cell != outer_cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbor_count = 0;

		const boost::array<double, 3>
			cell_center = game_grid.geometry.get_center(*cell),
			cell_length = game_grid.geometry.get_length(*cell);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	// get some statistics
	unordered_map<int, int> neighbor_count_histogram;
	for (int i = 0; i < 8 * (9 + 8 + 9); i++) {
		neighbor_count_histogram[i] = 0;
	}
	double avg_neighbors = 0.0;
	int max_neighbors = 0, min_neighbors = 999;
	uint64_t number_of_cells = inner_cells.size() + outer_cells.size();

	for (vector<uint64_t>::const_iterator cell = inner_cells.begin(); cell != inner_cells.end(); cell++) {

		const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
		avg_neighbors += double(neighbors->size()) / number_of_cells;

		if (max_neighbors < int(neighbors->size())) {
			max_neighbors = neighbors->size();
		}

		if (min_neighbors > int(neighbors->size())) {
			min_neighbors = neighbors->size();
		}

		if (neighbors->size() >= 7 && neighbors->size() <= 8 * (9 + 8 + 9)) {
			neighbor_count_histogram[neighbors->size()] += 1;
		} else {
			cout << "Impossible number of neighbors for cell " << *cell << ": " << neighbors->size() << endl;
		}
	}
	for (vector<uint64_t>::const_iterator cell = outer_cells.begin(); cell != outer_cells.end(); cell++) {

		const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
		avg_neighbors += double(neighbors->size()) / number_of_cells;

		if (max_neighbors < int(neighbors->size())) {
			max_neighbors = neighbors->size();
		}

		if (min_neighbors > int(neighbors->size())) {
			min_neighbors = neighbors->size();
		}

		if (neighbors->size() >= 7 && neighbors->size() <= 8 * (9 + 8 + 9)) {
			neighbor_count_histogram[neighbors->size()] += 1;
		} else {
			cout << "Impossible number of neighbors for cell " << *cell << ": " << neighbors->size() << endl;
		}
	}

	// print statistics
	int total_max_neighbors = all_reduce(comm, max_neighbors, plus<int>());
	int total_min_neighbors = all_reduce(comm, min_neighbors, plus<int>());
	double total_avg_neighbors = all_reduce(comm, avg_neighbors, plus<double>());

	for (int i = 0; i <= 8 * (9 + 8 + 9); i++) {
		int total = all_reduce(comm, neighbor_count_histogram[i], plus<int>());
		neighbor_count_histogram[i] = total;
	}

	if (comm.rank() == 0) {
		cout << "Max neighbors: " << total_max_neighbors << endl;
		cout << "Min neighbors: " << total_min_neighbors << endl;
		cout << "Avg. neighbors: " << total_avg_neighbors << endl;
		cout << "Neighbour count histogram: (neighbors, count)" << endl;
		for (int i = 0; i <= 8 * (9 + 8 + 9); i++) {
			cout << i << " " << neighbor_count_histogram[i] << endl;
		}
	}

	if (comm.rank() == 0) {
		cout << "step: ";
	}

	#define TIME_STEPS 100
	before = time(NULL);
	uint64_t processed_neighbors = 0;
	for (int step = 0; step < TIME_STEPS; step++) {

		if (comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		game_grid.start_remote_neighbor_copy_updates();
		// get the neighbor counts of every cell, starting with the cells whose neighbor data doesn't come from other processes
		for (vector<uint64_t>::const_iterator cell = inner_cells.begin(); cell != inner_cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];
			cell_data->live_neighbor_count = 0;

			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			for (vector<uint64_t>::const_iterator neighbor = neighbors->begin(); neighbor != neighbors->end(); neighbor++, processed_neighbors++) {

				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// wait for neighbor data updates to finish and go through the rest of the cells
		game_grid.wait_remote_neighbor_copy_updates();
		for (vector<uint64_t>::const_iterator cell = outer_cells.begin(); cell != outer_cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];
			cell_data->live_neighbor_count = 0;

			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			for (vector<uint64_t>::const_iterator neighbor = neighbors->begin(); neighbor != neighbors->end(); neighbor++, processed_neighbors++) {

				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// calculate the next turn
		for (vector<uint64_t>::const_iterator cell = inner_cells.begin(); cell != inner_cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbor_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbor_count != 2) {
				cell_data->is_alive = false;
			}
		}
		for (vector<uint64_t>::const_iterator cell = outer_cells.begin(); cell != outer_cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbor_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbor_count != 2) {
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
	cout << "processed neighbors: " << processed_neighbors << endl;
	cout << "Process " << comm.rank() << ": " << number_of_cells * TIME_STEPS << " cells processed at the speed of " << double(number_of_cells * TIME_STEPS) / total << " cells / second"<< endl;

	return EXIT_SUCCESS;
}
