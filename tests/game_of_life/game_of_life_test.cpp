/*
Program for testing dccrg using Conway's game of life.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016, 2018 Finnish Meteorological Institute
Copyright 2018 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "algorithm"
#include "array"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "unordered_set"

#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "mpi.h"
#include "zoltan.h"

#include "dccrg_stretched_cartesian_geometry.hpp"
#include "dccrg.hpp"

#include "cell.hpp"
#include "initialize.hpp"
#include "save.hpp"
#ifdef OPTIMIZED
#include "solve_optimized.hpp"
#else
#include "solve.hpp"
#endif

using namespace std;
using namespace boost;
using namespace dccrg;


/*!
Returns EXIT_SUCCESS if the state of the given game at given timestep is correct on this process, returns EXIT_FAILURE otherwise.
timestep == 0 means before any turns have been taken.
*/
int check_game_of_life_state(int timestep, const Dccrg<Cell, Stretched_Cartesian_Geometry>& grid)
{
	for (const auto& cell: grid.local_cells()) {
		// check cells that are always supposed to be alive
		switch (cell.id) {
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
			if (!cell.data->data[0]) {
				cerr << "Cell " << cell.id
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

		switch (cell.id) {
		case 109:
		case 123:
		case 189:
		case 190:
		case 198:
		case 200:
		case 204:
		case 205:
			if (!cell.data->data[0]) {
				cerr << "Cell " << cell.id
					<< " isn't alive on timestep " << timestep
					<< endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		} else {

		switch (cell.id) {
		case 174:
		case 184:
		case 214:
		case 220:
			if (!cell.data->data[0]) {
				cerr << "Cell " << cell.id
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
			switch (cell.id) {
			case 43:
			case 44:
			case 45:
			case 60:
			case 74:
				if (!cell.data->data[0]) {
					cerr << "Cell " << cell.id
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
			switch (cell.id) {
			case 29:
			case 44:
			case 45:
			case 58:
			case 60:
				if (!cell.data->data[0]) {
					cerr << "Cell " << cell.id
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
			switch (cell.id) {
			case 29:
			case 30:
			case 43:
			case 45:
			case 60:
				if (!cell.data->data[0]) {
					cerr << "Cell " << cell.id
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
			switch (cell.id) {
			case 29:
			case 30:
			case 45:
			case 59:
				if (!cell.data->data[0]) {
					cerr << "Cell " << cell.id
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
			switch (cell.id) {
			case 29:
			case 30:
			case 45:
				if (!cell.data->data[0]) {
					cerr << "Cell " << cell.id
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
	bool save_results = false, verbose = false;
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
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	if (option_variables.count("save") > 0) {
		save_results = true;
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}
	if (verbose && rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	// initialize grid
	Dccrg<Cell, Stretched_Cartesian_Geometry> grid;

	const uint64_t base_length = 15;
	const double cell_length = 1.0 / base_length;

	// set grid length in each dimension based on direction given by user
	std::array<uint64_t, 3> grid_length = {{0, 0, 0}};
	switch (direction) {
		case 'x':
			grid_length[0] = 1;
			grid_length[1] = base_length;
			grid_length[2] = base_length;
			break;

		case 'y':
			grid_length[0] = base_length;
			grid_length[1] = 1;
			grid_length[2] = base_length;
			break;

		case 'z':
			grid_length[0] = base_length;
			grid_length[1] = base_length;
			grid_length[2] = 1;
			break;

		default:
			cerr << "Unsupported direction given: " << direction << endl;
			break;
	}

	const unsigned int neighborhood_size = 1;
	grid
		.set_initial_length(grid_length)
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(0)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (uint64_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}

	grid.set_geometry(geom_params);

	#ifdef SEND_SINGLE_CELLS
	grid.set_send_single_cells(true);
	#endif

	if (verbose && rank == 0) {
		cout << "Maximum refinement level of the grid: " << grid.get_maximum_refinement_level()
			<< "\nNumber of cells: "
			<< (geom_params.coordinates[0].size() - 1)
				* (geom_params.coordinates[1].size() - 1)
				* (geom_params.coordinates[2].size() - 1)
			<< "\nSending single cells: " << boolalpha << grid.get_send_single_cells()
			<< endl << endl;
	}

	initialize(grid, base_length);

	// every process outputs the game state into its own file
	string basename("tests/game_of_life/game_of_life_test_");
	basename.append(1, direction).append("_").append(lexical_cast<string>(rank)).append("_");

	// visualize the game with visit -o game_of_life_test.visit
	ofstream visit_file;
	if (save_results && rank == 0) {
		string visit_file_name("tests/game_of_life/game_of_life_test_");
		visit_file_name += direction;
		visit_file_name += ".visit";
		visit_file.open(visit_file_name.c_str());
		visit_file << "!NBLOCKS " << comm_size << endl;
	}

	const int time_steps = 25;
	if (verbose && rank == 0) {
		cout << "step: ";
	}
	for (int step = 0; step < time_steps; step++) {

		grid.balance_load();
		grid.start_remote_neighbor_copy_updates();
		grid.wait_remote_neighbor_copy_updates();

		int result = check_game_of_life_state(step, grid);
		if ((grid_length[0] != 15 and grid_length[1] != 15) or result != EXIT_SUCCESS) {
			cout << "Process " << rank << ": Game of Life test failed on timestep: " << step << endl;
			abort();
		}

		if (verbose && rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		if (save_results) {
			// write the game state into a file named according to the current time step
			string output_name(basename);
			output_name.append(lexical_cast<string>(step)).append(".vtk");
			save(output_name, rank, grid);

			// visualize the game with visit -o game_of_life_test.visit
			if (rank == 0) {
				for (int process = 0; process < comm_size; process++) {
					visit_file << "game_of_life_test_"
						<< direction << "_"
						<< process << "_"
						<< lexical_cast<string>(step)
						<< ".vtk"
						<< endl;
				}
			}
		}

		get_live_neighbors(grid);
	}

	if (rank == 0 and save_results) {
		visit_file.close();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
