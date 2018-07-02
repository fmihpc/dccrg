/*
Program for testing dccrg restart using Conway's game of life.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016, 2018 Finnish Meteorological Institute

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
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "string"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"
#include "../../dccrg_cartesian_geometry.hpp"

// include restart cell before gol cell
#include "cell.hpp"
#include "IO.hpp"
#include "initialize.hpp"
#include "solve.hpp"
#include "refine.hpp"

using namespace std;
using namespace dccrg;


/*
Migrates cells off process 0.
*/
template <class Grid_T> void migrate_cells(Grid_T& grid)
{
	if (grid.get_comm_size() == 1) {
		return;
	}

	if (grid.get_rank() == 0) {
		const std::vector<uint64_t> cells = grid.get_cells();
		for (const auto& cell: cells) {
			grid.pin(cell, 1);
		}
	}

	grid.balance_load();
	grid.unpin_local_cells();
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


	Dccrg<Cell, Cartesian_Geometry> grid, reference_grid;


	/*
	Common grid parameters
	*/
	const std::array<uint64_t, 3> grid_length = {{15, 15, 1}};
	const unsigned int neighborhood_size = 1;

	Cartesian_Geometry::Parameters geom_params;
	geom_params.start[0] =
	geom_params.start[1] =
	geom_params.start[2] = 0;
	geom_params.level_0_cell_length[0] =
	geom_params.level_0_cell_length[1] =
	geom_params.level_0_cell_length[2] = 1.0 / grid_length[0];


	/*
	Setup reference grid and game
	*/

	reference_grid
		.set_initial_length(grid_length)
		.set_neighborhood_length(neighborhood_size)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(MPI_COMM_SELF);

	reference_grid.set_geometry(geom_params);
	initialize(reference_grid, grid_length[0]);


	/*
	Setup restart grid and game
	*/

	uint64_t step = 0;

	// either start a new game...
	if (restart_name == "") {

		grid
			.set_initial_length(grid_length)
			.set_neighborhood_length(neighborhood_size)
			.set_maximum_refinement_level(-1)
			.set_load_balancing_method("RANDOM")
			.initialize(comm)
			.set_geometry(geom_params);

		migrate_cells(grid);

		initialize(grid, grid_length[0]);

		// save initial state
		save("gol_0.dc", 0, grid);

	// ...or restart from saved game
	} else {
		step = load(comm, restart_name, grid);

		migrate_cells(grid);

		// play the reference game to the same step
		for (uint64_t i = 0; i < step; i++) {
			solve(reference_grid);
		}
	}

	const uint64_t time_steps = 25;
	while (step < time_steps) {

		refine(grid, grid_length[0], step, comm_size);

		grid.balance_load();
		grid.update_copies_of_remote_neighbors();

		solve(grid);
		solve(reference_grid);

		step++;

		save(
			"gol_" + boost::lexical_cast<std::string>(step) + ".dc",
			step,
			grid
		);

		// verify refined/unrefined game
		for (const auto& cell: grid.local_cells) {

			const int refinement_level = grid.get_refinement_level(cell.id);
			uint64_t reference_cell;
			if (refinement_level > 0) {
				reference_cell = grid.mapping.get_level_0_parent(cell.id);
			} else {
				reference_cell = cell.id;
			}

			auto* const reference_data = reference_grid[reference_cell];
			if (reference_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for reference cell " << reference_cell
					<< std::endl;
				abort();
			}

			if (cell.data->data[0] != reference_data->data[0]) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell's " << cell.id
					<< " life doesn't agree with reference at the beginning of step " << step
					<< std::endl;
				abort();
			}
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

