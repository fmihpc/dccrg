/*
Tests the grid using simple refinement which should induce refinement also across processes
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "../../dccrg.hpp"
#include "fstream"
#include "iostream"
#include "unistd.h"
#include "zoltan.h"


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
	#define GRID_SIZE 2
	#define CELL_SIZE (1.0 / GRID_SIZE)
	#define STENCIL_SIZE 1
	dccrg<int> grid(comm, "RANDOM", STARTING_CORNER, STARTING_CORNER, STARTING_CORNER, CELL_SIZE, GRID_SIZE, 1, 1, STENCIL_SIZE);

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "refine_simple_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("refine_simple.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	#define TIME_STEPS 3
	for (int step = 0; step < TIME_STEPS; step++) {

		// refine the smallest cell that is closest to the starting corner
		grid.refine_completely_at(0.00000001 * CELL_SIZE, 0.00000001 * CELL_SIZE, 0.00000001 * CELL_SIZE);
		grid.stop_refining();
		grid.balance_load();
		vector<uint64_t> cells = grid.get_cells();
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
				visit_file << "refine_simple_" << process << "_" << step_string.str() << suffix.str() << endl;
			}
		}

		// write the grid into a file
		grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// write each cells neighbour count
		outfile << "SCALARS neighbours int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			const vector<uint64_t>* neighbours = grid.get_neighbours(*cell);
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

	vector<uint64_t> cells = grid.get_cells();
	for (int i = 0; i < comm.size(); i++) {
		comm.barrier();
		if (i != comm.rank()) {
			continue;
		}
		for (vector<uint64_t>::const_iterator c = cells.begin(); c != cells.end(); c++) {
			const vector<uint64_t>* neighbours = grid.get_neighbours(*c);
			vector<uint64_t> sorted_neighbours(neighbours->begin(), neighbours->end());
			sort(sorted_neighbours.begin(), sorted_neighbours.end());
			cout << "Cell " << *c << " neighbours (" << sorted_neighbours.size() << "): ";
			for (vector<uint64_t>::const_iterator n = sorted_neighbours.begin(); n != sorted_neighbours.end(); n++) {
				cout << *n << " ";
			}
			cout << endl;
		}
		cout.flush();
		sleep(3);
	}

	if (comm.rank() == 0) {
		visit_file.close();
	}

	return EXIT_SUCCESS;
}
