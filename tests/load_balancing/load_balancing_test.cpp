/*
Tests a Zoltan load balancing algorithm on a refined grid in 3d.
Visualize the results using for example VisIt.
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg.hpp"

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


	#define X_LENGTH 10
	#define Y_LENGTH 10
	#define Z_LENGTH 10
	#define CELL_X_SIZE (1.0 / X_LENGTH)
	#define CELL_Y_SIZE (1.0 / Y_LENGTH)
	#define CELL_Z_SIZE (1.0 / Z_LENGTH)
	#define STENCIL_SIZE 1
	#define MAX_REFINEMENT 2
	dccrg<int> game_grid(comm, "HYPERGRAPH", -0.5, -0.5, -0.5, CELL_X_SIZE, CELL_Y_SIZE, CELL_Z_SIZE, X_LENGTH, Y_LENGTH, Z_LENGTH, STENCIL_SIZE, MAX_REFINEMENT);

	vector<uint64_t> cells;
	for (int i = 0; i < MAX_REFINEMENT; i++) {
		cells = game_grid.get_cells();
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			double x = game_grid.get_cell_x(*cell);
			double y = game_grid.get_cell_y(*cell);
			double z = game_grid.get_cell_z(*cell);

			if (sqrt(x * x + y * y + z * z) < 0.1) {
				game_grid.refine_completely(*cell);
			}
		}
		game_grid.stop_refining();
	}

	if (comm.rank() == 0) {
		cout << "Total cells after refining: " << all_reduce(comm, cells.size(), plus<uint64_t>()) << endl;
	} else {
		all_reduce(comm, cells.size(), plus<uint64_t>());
	}


	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "load_balancing_test_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("load_balancing_test.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	#define TIME_STEPS 1
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.start_remote_neighbour_data_update();
		game_grid.wait_neighbour_data_update();
		cells = game_grid.get_cells();

		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
		sort(cells.begin(), cells.end());

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
				visit_file << "load_balancing_test_" << process << "_" << step_string.str() << suffix.str() << endl;
			}
		}


		// write the grid into a file
		game_grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

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
	}

	return EXIT_SUCCESS;
}
