/*
Tests the grid with some simple game of life patters using a hierarchical load balancing method.
Returns EXIT_SUCCESS if everything went ok.
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg.hpp"


struct game_of_life_cell {

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
	}

	bool is_alive;
	unsigned int live_neighbor_count;
};


using namespace std;
using namespace boost::mpi;
using namespace dccrg;


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	Dccrg<game_of_life_cell> game_grid;

	const boost::array<uint64_t, 3> grid_length = {{34, 7, 1}};
	game_grid.initialize(grid_length, comm, "HIER", 1);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: "
			<< game_grid.get_maximum_refinement_level()
			<< endl;
	}

	game_grid.add_partitioning_level(12);
	game_grid.add_partitioning_option(0, "LB_METHOD", "HYPERGRAPH");
	game_grid.add_partitioning_option(0, "IMBALANCE_TOL", "1.05");

	game_grid.add_partitioning_level(1);
	game_grid.add_partitioning_option(1, "LB_METHOD", "HYPERGRAPH");

	game_grid.balance_load();
	vector<uint64_t> cells = game_grid.get_cells();
	sort(cells.begin(), cells.end());

	// initialize the game
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "hierarchical_test_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("hierarchical_test.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
		cout << "step: ";
	}

	#define TIME_STEPS 25
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.start_remote_neighbor_copy_updates();
		game_grid.wait_remote_neighbor_copy_updates();

		if (comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		// write the game state into a file named according to the current time step
		string current_output_name("");
		ostringstream step_string;
		step_string.fill('0');
		step_string.width(5);
		step_string << step;
		current_output_name += basename.str();
		current_output_name += step_string.str();
		current_output_name += suffix.str();

		// visualize the game with visit -o game_of_life_test.visit
		if (comm.rank() == 0) {
			for (int process = 0; process < comm.size(); process++) {
				visit_file << "hierarchical_test_" << process << "_" << step_string.str() << suffix.str() << endl;
			}
		}


		// write the grid into a file
		game_grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// go through the grids cells and write their state into the file
		outfile << "SCALARS is_alive float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->is_alive == true) {
				outfile << "1";
			} else {
				outfile << "0";
			}
			outfile << endl;

		}

		// write each cells live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			outfile << cell_data->live_neighbor_count << endl;

		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			outfile << neighbors->size() << endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << comm.rank() << endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << *cell << endl;
		}
		outfile.close();


		// get the neighbor counts of every cell, starting with the cells whose neighbor data doesn't come from other processes
		vector<uint64_t> inner_cells = game_grid.get_local_cells_not_on_process_boundary();
		for (vector<uint64_t>::const_iterator cell = inner_cells.begin(); cell != inner_cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			cell_data->live_neighbor_count = 0;
			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			if (neighbors == NULL) {
				cout << "Process " << comm.rank() << ": neighbor list for cell " << *cell << " not available" << endl;
				exit(EXIT_FAILURE);
			}

			for (vector<uint64_t>::const_iterator neighbor = neighbors->begin(); neighbor != neighbors->end(); neighbor++) {
				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data == NULL) {
					cout << "Process " << comm.rank() << ": neighbor " << *neighbor << " data for cell " << *cell << " not available" << endl;
					exit(EXIT_FAILURE);
				}
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		vector<uint64_t> outer_cells = game_grid.get_local_cells_on_process_boundary();
		for (vector<uint64_t>::const_iterator cell = outer_cells.begin(); cell != outer_cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			cell_data->live_neighbor_count = 0;
			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(*cell);
			if (neighbors == NULL) {
				cout << "Process " << comm.rank() << ": neighbor list for cell " << *cell << " not available" << endl;
				exit(EXIT_FAILURE);
			}

			for (vector<uint64_t>::const_iterator neighbor = neighbors->begin(); neighbor != neighbors->end(); neighbor++) {
				if (*neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[*neighbor];
				if (neighbor_data == NULL) {
					cout << "Process " << comm.rank() << ": neighbor " << *neighbor << " data for cell " << *cell << " not available" << endl;
					exit(EXIT_FAILURE);
				}
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

	if (comm.rank() == 0) {
		visit_file.close();
		cout << endl;
	}

	return EXIT_SUCCESS;
}
