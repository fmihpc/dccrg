/*
Tests the grid with a game of life on a refined grid in 3 D with neighbours only in the ? plane
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "../../dccrg.hpp"
#include "fstream"
#include "iostream"
#include "zoltan.h"


struct game_of_life_cell {

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
	}

	bool is_alive;
	unsigned int live_neighbour_count;
};


using namespace std;
using namespace boost::mpi;

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


	#define STARTING_CORNER 0.0
	#define GRID_SIZE 1
	#define CELL_SIZE (1.0 / GRID_SIZE)
	#define STENCIL_SIZE 1
	dccrg<game_of_life_cell> game_grid(comm, "RANDOM", STARTING_CORNER, STARTING_CORNER, STARTING_CORNER, CELL_SIZE, GRID_SIZE, GRID_SIZE, 3, STENCIL_SIZE);
	vector<uint64_t> cells = game_grid.get_cells();
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << game_grid.get_max_refinement_level() << endl;
	}
	cout << "Process " << comm.rank() << ": number of cells: " << cells.size() << endl;

	// refine the grid increasingly in the z direction a few times
	game_grid.balance_load();
	cells = game_grid.get_cells();
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		if (game_grid.get_cell_z(*cell) > 1 * CELL_SIZE) {
			game_grid.refine_completely(*cell);
		}
	}
	game_grid.stop_refining();
	game_grid.balance_load();
	cells = game_grid.get_cells();
	cout << "Process " << comm.rank() << ": number of cells after refining: " << cells.size() << endl;
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		if (game_grid.get_cell_z(*cell) > 1.5 * CELL_SIZE) {
			game_grid.refine_completely(*cell);
		}
	}
	game_grid.stop_refining();
	game_grid.balance_load();
	cells = game_grid.get_cells();
	cout << "Process " << comm.rank() << ": number of cells after refining: " << cells.size() << endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

		game_of_life_cell* cell_data = game_grid[*cell];
		cell_data->live_neighbour_count = 0;

		double y = game_grid.get_cell_y(*cell);
		if (fabs(0.5 + 0.1 * game_grid.get_cell_size(*cell) - y) < 0.5 * game_grid.get_cell_size(*cell)) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "refined_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("refined.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	#define TIME_STEPS 25
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.update_remote_neighbour_data();
		cells = game_grid.get_cells();
		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
		sort(cells.begin(), cells.end());

		if (comm.rank() == 0) {
			cout << "step: " << step << endl;
		}
		/*if (step == 20 || step == 21) {
			game_grid.print_cell_to_process_mappings();
			game_grid.print_remote_neighbours();
		}
		cout.flush();*/

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
				visit_file << "refined_" << process << "_" << step_string.str() << suffix.str() << endl;
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

		// write each cells live neighbour count
		outfile << "SCALARS live_neighbour_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			outfile << cell_data->live_neighbour_count << endl;

		}

		// write each cells neighbour count
		outfile << "SCALARS neighbours int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			outfile << neighbours->size() << endl;
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

		// get the neighbour counts of every cell
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			cell_data->live_neighbour_count = 0;
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			if (neighbours == NULL) {
				cout << "Process " << comm.rank() << ": neighbour list for cell " << *cell << " not available" << endl;
				exit(EXIT_FAILURE);
			}

			for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++) {

				// only consider neighbours in the same z plane
				if (game_grid.get_cell_z(*cell) != game_grid.get_cell_z(*neighbour)) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data == NULL) {
					cout << "Process " << comm.rank() << ": neighbour " << *neighbour << " data of cell " << *cell << " not available" << endl;
					exit(EXIT_FAILURE);
				}
				if (neighbour_data->is_alive) {
					cell_data->live_neighbour_count++;
				}
			}
		}

		// calculate the next turn
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			if (cell_data->live_neighbour_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbour_count != 2) {
				cell_data->is_alive = false;
			}
		}

	}

	if (comm.rank() == 0) {
		visit_file.close();
	}

	return EXIT_SUCCESS;
}
