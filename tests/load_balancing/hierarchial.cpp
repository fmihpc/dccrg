/*
Tests hierarchial Zoltan load balancing on a 2d grid.
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


	#define X_LENGTH 100
	#define Y_LENGTH 100
	#define Z_LENGTH 1
	#define CELL_X_SIZE 1.0
	#define CELL_Y_SIZE 1.0
	#define CELL_Z_SIZE 1.0
	#define STENCIL_SIZE 1
	dccrg<int> game_grid(comm, "HIER", 0, 0, 0, CELL_X_SIZE, CELL_Y_SIZE, CELL_Z_SIZE, X_LENGTH, Y_LENGTH, Z_LENGTH, STENCIL_SIZE);

	// first divide the grid cells into two process blocks
	game_grid.add_partitioning_level(2);
	game_grid.add_partitioning_option(0, "LB_METHOD", "RCB");
	game_grid.add_partitioning_option(0, "IMBALANCE_TOL", "1.05");
	// check whether reserved options can be changed
	game_grid.add_partitioning_option(0, "NUM_GID_ENTRIES", "3");

	// then assign cells randomly to processes within each block
	game_grid.add_partitioning_level(1);
	game_grid.add_partitioning_option(1, "LB_METHOD", "RANDOM");

	game_grid.balance_load();

	vector<uint64_t> cells = game_grid.get_cells();
	sort(cells.begin(), cells.end());

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "hierarchial_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("hierarchial.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	// write the game state into a file named according to the current time step
	string current_output_name("");
	current_output_name += basename.str();
	current_output_name += suffix.str();

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		for (int process = 0; process < comm.size(); process++) {
			visit_file << "hierarchial_" << process << "_" << suffix.str() << endl;
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

	return EXIT_SUCCESS;
}
