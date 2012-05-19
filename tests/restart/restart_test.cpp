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
#include "string"
#include "zoltan.h"

#include "../../dccrg.hpp"

// include restart cell before gol cell
#include "cell.hpp"
#include "../game_of_life/initialize.hpp"
#include "../game_of_life/solve.hpp"
#include "IO.hpp"

using namespace std;
using namespace boost;
using namespace boost::mpi;
using namespace dccrg;


/*!
Returns EXIT_SUCCESS if the state of the given game at given timestep is correct on this process, returns EXIT_FAILURE otherwise.
timestep == 0 means before any turns have been taken.
*/
int check_game_of_life_state(int timestep, const Dccrg<Cell, ConstantGeometry>& grid)
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
	string restart_name;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("restart",
			boost::program_options::value<std::string>(&restart_name)->default_value(""),
			"restart the game from file arg (don't restart if name is empty)");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (comm.rank() == 0) {
			cout << options << endl;
		}
		comm.barrier();
		return EXIT_SUCCESS;
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}

	// initialize grid
	Dccrg<Cell, ConstantGeometry> game_grid;

	const int grid_size = 15;	// in unrefined cells
	const double cell_size = 1.0 / grid_size;

	if (!game_grid.set_geometry(grid_size, grid_size, 1, 0, 0, 0, cell_size, cell_size, cell_size)) {
		cerr << "Couldn't set grid geometry" << endl;
		exit(EXIT_FAILURE);
	}

	const unsigned int neighborhood_size = 1;
	game_grid.initialize(comm, "RANDOM", neighborhood_size);
	game_grid.balance_load();

	const uint64_t time_steps = 25;
	uint64_t step = 0;

	// start a new game
	if (restart_name == "") {

		Initialize<ConstantGeometry>::initialize(game_grid, grid_size);

		// save initial state
		IO<ConstantGeometry>::save(
			"gol_0.dc",
			0,
			comm,
			game_grid
		);

	// restart from saved game
	} else {

		IO<ConstantGeometry>::load(
			restart_name,
			step,
			comm,
			game_grid
		);

		step++;
	}

	for ( ; step < time_steps; step++) {

		game_grid.balance_load();
		game_grid.update_remote_neighbor_data();
		vector<uint64_t> cells = game_grid.get_cells();

		int result = check_game_of_life_state(step, game_grid);
		if (grid_size != 15 || result != EXIT_SUCCESS) {
			cout << "Process " << comm.rank()
				<< ": Game of Life test failed on time step: " << step
				<< endl;
			abort();
		}

		Solve<ConstantGeometry>::solve(game_grid);

		IO<ConstantGeometry>::save(
			"gol_" + boost::lexical_cast<std::string>(step) + ".dc",
			step,
			comm,
			game_grid
		);
	}

	if (comm.rank() == 0) {
		cout << "Passed" << endl;
	}

	return EXIT_SUCCESS;
}

