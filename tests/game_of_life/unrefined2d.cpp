/*
As refined2d.cpp but refines / unrefines the grid constantly and randomly
*/

#include "algorithm"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
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
#include "save.hpp"
#include "solve.hpp"

using namespace std;
using namespace boost;
using namespace boost::mpi;
using namespace dccrg;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

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
		if (comm.rank() == 0) {
			cout << options << endl;
		}
		comm.barrier();
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
	if (verbose && comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	// initialize grid
	Dccrg<Cell, ArbitraryGeometry> game_grid;

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

	const unsigned int neighborhood_size = 1;
	game_grid.initialize(comm, "RANDOM", neighborhood_size);

	Initialize<ArbitraryGeometry>::initialize(game_grid, grid_size);

	// every process outputs the game state into its own file
	string basename("unrefined2d_");
	basename.append(1, direction).append("_").append(lexical_cast<string>(comm.rank())).append("_");

	// visualize the game with visit -o game_of_life_test.visit
	ofstream visit_file;
	if (save && comm.rank() == 0) {
		string visit_file_name("unrefined2d_");
		visit_file_name += direction;
		visit_file_name += ".visit";
		visit_file.open(visit_file_name.c_str());
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	if (verbose && comm.rank() == 0) {
		cout << "step: ";
		cout.flush();
	}

	const int time_steps = 25;
	for (int step = 0; step < time_steps; step++) {

		// refine random unrefined cells and unrefine random refined cells
		vector<uint64_t> cells = game_grid.get_cells();
		random_shuffle(cells.begin(), cells.end());

		if (step % 2 == 0) {

			for (int i = 0, refined = 0;
				i < int(cells.size()) && refined <= grid_size * grid_size / (5 * comm.size());
				i++
			) {
				if (game_grid.get_refinement_level(cells[i]) == 0) {
					game_grid.refine_completely(cells[i]);
					refined++;
				}
			}

		} else {

			for (int i = 0, unrefined = 0;
				i < int(cells.size()) && unrefined <= grid_size * grid_size / (4 * comm.size());
				i++
			) {
				if (game_grid.get_refinement_level(cells[i]) > 0) {
					game_grid.unrefine_completely(cells[i]);
					unrefined++;
				}
			}
		}

		vector<uint64_t> new_cells = game_grid.stop_refining();

		// assign parents' state to children
		for (vector<uint64_t>::const_iterator new_cell = new_cells.begin(); new_cell != new_cells.end(); new_cell++) {
			Cell* new_cell_data = game_grid[*new_cell];
			if (new_cell_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " no data for created cell " << *new_cell
					<< std::endl;
				abort();
			}
			Cell* parent_data = game_grid[game_grid.get_parent(*new_cell)];
			if (parent_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " no data for parent cell " << game_grid.get_parent(*new_cell)
					<< std::endl;
				abort();
			}
			new_cell_data->data[0] = parent_data->data[0];
		}

		// "interpolate" parent cell's value from unrefined children
		vector<uint64_t> removed_cells = game_grid.get_removed_cells();
		for (vector<uint64_t>::const_iterator
			removed_cell = removed_cells.begin();
			removed_cell != removed_cells.end();
			removed_cell++
		) {
			Cell* removed_cell_data = game_grid[*removed_cell];
			if (removed_cell_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " no data for removed cell after unrefining: " << *removed_cell
					<< std::endl;
				abort();
			}
			Cell* parent_data = game_grid[game_grid.get_parent_for_removed(*removed_cell)];
			if (parent_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " no data for parent cell after unrefining: " << game_grid.get_parent_for_removed(*removed_cell)
					<< std::endl;
				abort();
			}
			parent_data->data[0] = removed_cell_data->data[0];
		}
		game_grid.clear_refined_unrefined_data();

		game_grid.balance_load();
		game_grid.update_remote_neighbor_data();

		if (verbose && comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		if (save) {
			// write the game state into a file named according to the current time step
			string output_name(basename);
			output_name.append(lexical_cast<string>(step)).append(".vtk");
			Save<ArbitraryGeometry>::save(output_name, comm.rank(), game_grid);

			// visualize the game with visit -o game_of_life_test.visit
			if (comm.rank() == 0) {
				for (int process = 0; process < comm.size(); process++) {
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
	}

	if (comm.rank() == 0) {
		if (save) {
			visit_file.close();
		}

		if (verbose) {
			cout << endl;
		}
	}

	return EXIT_SUCCESS;
}
