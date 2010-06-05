/*
Tests the grid using simple refining and unrefining which should induce unrefinement also across processes
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "../../dccrg.hpp"
#include "fstream"
#include "iostream"
#include "zoltan.h"


using namespace std;
using namespace boost::mpi;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	clock_t before, after;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}


	#define GRID_SIZE 2
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
	}
	y_coordinates.push_back(0);
	y_coordinates.push_back(1);
	z_coordinates.push_back(0);
	z_coordinates.push_back(1);
	#define STENCIL_SIZE 1
	dccrg<int> grid(comm, "RANDOM", x_coordinates, y_coordinates, z_coordinates, STENCIL_SIZE, 5);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << grid.get_max_refinement_level() << endl;
		cout << "Number of cells: " << (x_coordinates.size() - 1) * (y_coordinates.size() - 1) * (z_coordinates.size() - 1) << endl << endl;
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "unrefine_simple_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("unrefine_simple.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	#define TIME_STEPS 8
	for (int step = 0; step < TIME_STEPS; step++) {

		if (comm.rank() == 0) {
			cout << "step " << step << endl;
		}

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
				visit_file << "unrefine_simple_" << process << "_" << step_string.str() << suffix.str() << endl;
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

		before = clock();

		// refine / unrefine the smallest cell that is closest to the grid starting corner
		if (step < 4) {
			/*// unrefine all cells in the grid
			for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
				grid.unrefine_completely(*cell);
			}
			cout << "unrefined 1" << endl;*/

			grid.refine_completely_at(0.0001 * CELL_SIZE, 0.0001 * CELL_SIZE, 0.0001 * CELL_SIZE);

			/*for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
				grid.unrefine_completely(*cell);
			}
			cout << "unrefined 2" << endl;*/
		} else {
			grid.unrefine_completely_at(0.0001 * CELL_SIZE, 0.0001 * CELL_SIZE, 0.0001 * CELL_SIZE);
		}

		vector<uint64_t> new_cells = grid.stop_refining();

		after = clock();
		cout << "Process " << comm.rank() <<": Refining / unrefining took " << double(after - before) / CLOCKS_PER_SEC << " seconds, " << new_cells.size() << " new cells created" << endl;
	}

	if (comm.rank() == 0) {
		visit_file.close();
	}

	return EXIT_SUCCESS;
}
