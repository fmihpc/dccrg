/*
Program for testing dccrg using Conway's game of life.

Copyright 2010, 2011, 2012 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "algorithm"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
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


/*!
Returns EXIT_SUCCESS if the state of the given game at given timestep is correct on this process, returns EXIT_FAILURE otherwise.
timestep == 0 means before any turns have been taken.
*/
int check_game_of_life_state(int timestep, const Dccrg<Cell, ArbitraryGeometry>& grid)
{
	vector<uint64_t> cells = grid.get_cells();
	for (vector<uint64_t>::const_iterator
		cell = cells.begin();
		cell != cells.end();
		cell++
	) {
		Cell* data = grid[*cell];
		if (data == NULL) {
			cerr << "No data for cell " << *cell << endl;
			return EXIT_FAILURE;
		}

		// check cells that are always supposed to be alive
		switch (*cell) {
		case 22:
		case 23:
		case 32:
		case 33:
		case 36:
		case 39:
		case 47:
		case 48:
		case 52:
		case 53:
		case 94:
		case 95:
		case 110:
		case 122:
		case 137:
		case 138:
		case 188:
		case 199:
		case 206:
			if (!data->data[0]) {
				cerr << "Cell " << *cell
					<< " isn't alive on timestep " << timestep
					<< endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		// these are supposed to be alive every other turn
		if (timestep % 2 == 0) {

		switch (*cell) {
		case 109:
		case 123:
		case 189:
		case 190:
		case 198:
		case 200:
		case 204:
		case 205:
			if (!data->data[0]) {
				cerr << "Cell " << *cell
					<< " isn't alive on timestep " << timestep
					<< endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		} else {

		switch (*cell) {
		case 174:
		case 184:
		case 214:
		case 220:
			if (!data->data[0]) {
				cerr << "Cell " << *cell
					<< " isn't alive on timestep " << timestep
					<< endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		}

		// check that the glider is moving correctly
		switch (timestep) {
		/* can't be bothered manually for cases 1-19, use an automatic method later */
		case 20:
			switch (*cell) {
			case 43:
			case 44:
			case 45:
			case 60:
			case 74:
				if (!data->data[0]) {
					cerr << "Cell " << *cell
						<< " isn't alive on timestep " << timestep
						<< endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 21:
			switch (*cell) {
			case 29:
			case 44:
			case 45:
			case 58:
			case 60:
				if (!data->data[0]) {
					cerr << "Cell " << *cell
						<< " isn't alive on timestep " << timestep
						<< endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 22:
			switch (*cell) {
			case 29:
			case 30:
			case 43:
			case 45:
			case 60:
				if (!data->data[0]) {
					cerr << "Cell " << *cell
						<< " isn't alive on timestep " << timestep
						<< endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 23:
			switch (*cell) {
			case 29:
			case 30:
			case 45:
			case 59:
				if (!data->data[0]) {
					cerr << "Cell " << *cell
						<< " isn't alive on timestep " << timestep
						<< endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 24:
			switch (*cell) {
			case 29:
			case 30:
			case 45:
				if (!data->data[0]) {
					cerr << "Cell " << *cell
						<< " isn't alive on timestep " << timestep
						<< endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;
		default:
			break;
		}
	}

	return EXIT_SUCCESS;
}


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

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
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

	if (verbose && comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << game_grid.get_maximum_refinement_level() << endl;
		cout << "Number of cells: " << (x_coordinates.size() - 1) * (y_coordinates.size() - 1) * (z_coordinates.size() - 1) << endl << endl;
	}

	Initialize<ArbitraryGeometry>::initialize(game_grid, grid_size);

	// every process outputs the game state into its own file
	string basename("game_of_life_test_");
	basename.append(1, direction).append("_").append(lexical_cast<string>(comm.rank())).append("_");

	// visualize the game with visit -o game_of_life_test.visit
	ofstream visit_file;
	if (save && comm.rank() == 0) {
		string visit_file_name("game_of_life_test_");
		visit_file_name += direction;
		visit_file_name += ".visit";
		visit_file.open(visit_file_name.c_str());
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	const int time_steps = 25;
	if (verbose && comm.rank() == 0) {
		cout << "step: ";
	}
	for (int step = 0; step < time_steps; step++) {

		game_grid.balance_load();
		game_grid.start_remote_neighbor_data_update();
		game_grid.wait_neighbor_data_update();
		vector<uint64_t> cells = game_grid.get_cells();

		int result = check_game_of_life_state(step, game_grid);
		if (grid_size != 15 || result != EXIT_SUCCESS) {
			cout << "Process " << comm.rank() << ": Game of Life test failed on timestep: " << step << endl;
			abort();
		}

		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
		sort(cells.begin(), cells.end());

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
					visit_file << "game_of_life_test_"
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
		if (verbose) {
			cout << endl;
		}

		cout << "Passed" << endl;

		if (save) {
			visit_file.close();
		}
	}

	return EXIT_SUCCESS;
}

