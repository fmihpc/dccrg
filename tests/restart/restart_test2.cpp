/*
Program for testing dccrg restart using Conway's game of life, stretched geometry version.

Copyright 2010, 2011, 2012, 2013, 2014,
2015 Finnish Meteorological Institute

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
#include "array"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "string"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"
#include "../../dccrg_stretched_cartesian_geometry.hpp"

// include restart cell before gol cell
#include "cell.hpp"
#include "IO.hpp"
#include "initialize.hpp"
#include "solve.hpp"
#include "refine.hpp"

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
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Barrier(comm);
		return EXIT_SUCCESS;
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}


	Dccrg<Cell, Stretched_Cartesian_Geometry> game_grid, reference_grid;


	/*
	Common grid parameters
	*/
	const std::array<uint64_t, 3> grid_length = {{15, 15, 1}};
	const unsigned int neighborhood_size = 1;

	Stretched_Cartesian_Geometry::Parameters geom_params;
	geom_params.coordinates[2].push_back(-2);
	geom_params.coordinates[2].push_back(-1);
	for (size_t i = 0; i <= grid_length[0]; i++) {
		geom_params.coordinates[0].push_back(i * 0.9);
		geom_params.coordinates[1].push_back(i * 1.1);
	}

	/*
	Setup reference grid and game
	*/

	reference_grid.initialize(grid_length, MPI_COMM_SELF, "RANDOM", neighborhood_size);

	if (!reference_grid.set_geometry(geom_params)) {
		cerr << "Couldn't set reference grid geometry" << endl;
		exit(EXIT_FAILURE);
	}

	Initialize<Stretched_Cartesian_Geometry>::initialize(reference_grid, grid_length[0]);


	/*
	Setup restart grid and game
	*/

	uint64_t step = 0;

	// either start a new game...
	if (restart_name == "") {

		game_grid.initialize(grid_length, comm, "RANDOM", neighborhood_size);

		if (!game_grid.set_geometry(geom_params)) {
			cerr << "Couldn't set grid geometry" << endl;
			exit(EXIT_FAILURE);
		}

		Initialize<Stretched_Cartesian_Geometry>::initialize(game_grid, grid_length[0]);

		// save initial state
		IO<Stretched_Cartesian_Geometry>::save(
			"gol_0.dc",
			0,
			game_grid
		);

	// ...or restart from saved game
	} else {
		step = IO<Stretched_Cartesian_Geometry>::load(
			comm,
			restart_name,
			game_grid
		);

		// play the reference game to the same step
		for (uint64_t i = 0; i < step; i++) {
			Solve<Stretched_Cartesian_Geometry>::solve(reference_grid);
		}
	}

	game_grid.balance_load();

	const uint64_t time_steps = 25;
	while (step < time_steps) {

		Refine<Stretched_Cartesian_Geometry>::refine(game_grid, grid_length[0], step, comm_size);

		game_grid.balance_load();
		game_grid.update_copies_of_remote_neighbors();
		const vector<uint64_t> cells = game_grid.get_cells();

		Solve<Stretched_Cartesian_Geometry>::solve(game_grid);
		Solve<Stretched_Cartesian_Geometry>::solve(reference_grid);

		step++;

		IO<Stretched_Cartesian_Geometry>::save(
			"gol_" + boost::lexical_cast<std::string>(step) + ".dc",
			step,
			game_grid
		);

		// verify refined/unrefined game
		for (const auto& cell: cells) {

			auto* const cell_data = game_grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			uint64_t reference_cell;

			const int refinement_level = game_grid.get_refinement_level(cell);
			if (refinement_level > 0) {
				reference_cell = game_grid.mapping.get_level_0_parent(cell);
			} else {
				reference_cell = cell;
			}

			auto* const reference_data = reference_grid[reference_cell];
			if (reference_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for reference cell " << reference_cell
					<< std::endl;
				abort();
			}

			if (cell_data->data[0] != reference_data->data[0]) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell's " << cell
					<< " life doesn't agree with reference at the beginning of step " << step
					<< std::endl;
				abort();
			}
		}
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

