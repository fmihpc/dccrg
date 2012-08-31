/*
As refined2d.cpp but refines / unrefines the grid constantly and randomly
*/

#include "algorithm"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "string"
#include "zoltan.h"

#include "../../dccrg_arbitrary_geometry.hpp"
#include "../../dccrg.hpp"

#include "cell.hpp"
#include "initialize.hpp"
#include "refine.hpp"
#include "save.hpp"
#include "solve.hpp"

using namespace std;
using namespace boost;
using namespace dccrg;

int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	/*
	Options
	*/
	char direction;
	bool save = false, verbose = false;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("direction",
			boost::program_options::value<char>(&direction)->default_value('z'),
			"Create a 2d grid with normal into direction arg (x, y or z)")
		("save", "Save the game to vtk files")
		("verbose", "Print information about the game");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Barrier(comm);
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	if (option_variables.count("save") > 0) {
		save = true;
	}

	// initialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose && rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	// initialize grids, reference grid doesn't refine/unrefine
	Dccrg<Cell, ArbitraryGeometry> game_grid, reference_grid;

	const int grid_size = 15;	// in unrefined cells
	const double cell_size = 1.0 / grid_size;
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	switch (direction) {
	case 'x':
		for (int i = 0; i <= grid_size; i++) {
			y_coordinates.push_back(i * cell_size);
			z_coordinates.push_back(i * cell_size);
		}
		x_coordinates.push_back(0);
		x_coordinates.push_back(1);
		break;

	case 'y':
		for (int i = 0; i <= grid_size; i++) {
			x_coordinates.push_back(i * cell_size);
			z_coordinates.push_back(i * cell_size);
		}
		y_coordinates.push_back(0);
		y_coordinates.push_back(1);
		break;

	case 'z':
		for (int i = 0; i <= grid_size; i++) {
			x_coordinates.push_back(i * cell_size);
			y_coordinates.push_back(i * cell_size);
		}
		z_coordinates.push_back(0);
		z_coordinates.push_back(1);
		break;

	default:
		cerr << "Unsupported direction given: " << direction << endl;
		break;
	}

	if (!game_grid.set_geometry(x_coordinates, y_coordinates, z_coordinates)) {
		cerr << "Couldn't set grid geometry" << endl;
		exit(EXIT_FAILURE);
	}

	if (!reference_grid.set_geometry(x_coordinates, y_coordinates, z_coordinates)) {
		cerr << "Couldn't set reference grid geometry" << endl;
		exit(EXIT_FAILURE);
	}

	const unsigned int neighborhood_size = 1;
	game_grid.initialize(comm, "RANDOM", neighborhood_size);
	// play complete reference game on each process
	reference_grid.initialize(MPI_COMM_SELF, "RANDOM", neighborhood_size);

	Initialize<ArbitraryGeometry>::initialize(game_grid, grid_size);
	Initialize<ArbitraryGeometry>::initialize(reference_grid, grid_size);

	// every process outputs the game state into its own file
	string basename("unrefined2d_");
	basename.append(1, direction).append("_").append(lexical_cast<string>(rank)).append("_");

	// visualize the game with visit -o game_of_life_test.visit
	ofstream visit_file;
	if (save && rank == 0) {
		string visit_file_name("unrefined2d_");
		visit_file_name += direction;
		visit_file_name += ".visit";
		visit_file.open(visit_file_name.c_str());
		visit_file << "!NBLOCKS " << comm_size << endl;
	}

	if (verbose && rank == 0) {
		cout << "step: ";
		cout.flush();
	}

	const int time_steps = 25;
	for (int step = 0; step < time_steps; step++) {

		Refine<ArbitraryGeometry>::refine(game_grid, grid_size, step, comm_size);

		game_grid.balance_load();
		game_grid.update_remote_neighbor_data();

		if (verbose && rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		if (save) {
			// write the game state into a file named according to the current time step
			string output_name(basename);
			output_name.append(lexical_cast<string>(step)).append(".vtk");
			Save<ArbitraryGeometry>::save(output_name, rank, game_grid);

			// visualize the game with visit -o game_of_life_test.visit
			if (rank == 0) {
				for (int process = 0; process < comm_size; process++) {
					visit_file << "unrefined2d_"
						<< direction << "_"
						<< process << "_"
						<< lexical_cast<string>(step)
						<< ".vtk"
						<< endl;
				}
			}
		}

		Solve<ArbitraryGeometry>::solve(game_grid);
		Solve<ArbitraryGeometry>::solve(reference_grid);

		// verify refined/unrefined game
		const vector<uint64_t> cells = game_grid.get_cells();
		BOOST_FOREACH(const uint64_t cell, cells) {

			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			uint64_t reference_cell;

			const int refinement_level = game_grid.get_refinement_level(cell);
			if (refinement_level > 0) {
				reference_cell = game_grid.get_parent_for_removed(cell);
			} else {
				reference_cell = cell;
			}

			Cell* reference_data = reference_grid[reference_cell];
			if (reference_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for reference cell " << reference_cell
					<< std::endl;
				abort();
			}

			if (cell_data->data[0] != reference_data->data[0]) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell's " << cell << " life doesn't agree with reference."
					<< std::endl;
				abort();
			}
		}
	}

	if (rank == 0) {
		if (save) {
			visit_file.close();
		}

		if (verbose) {
			cout << endl;
		}
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
