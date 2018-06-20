/*
As refined2d.cpp but refines / unrefines the grid constantly and randomly

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute

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
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "string"
#include "unordered_set"

#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "mpi.h"
#include "zoltan.h"

#include "dccrg_stretched_cartesian_geometry.hpp"
#include "dccrg.hpp"
#include "cell.hpp"
#include "initialize.hpp"
#include "refine.hpp"
#include "save.hpp"
#ifdef OPTIMIZED
#include "solve_optimized.hpp"
#else
#include "solve.hpp"
#endif

using namespace std;
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

	// initialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		exit(EXIT_FAILURE);
	}

	// initialize grids, reference grid doesn't refine/unrefine
	Dccrg<Cell, Stretched_Cartesian_Geometry> grid, reference_grid;

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
	grid.initialize(grid_length, comm, "RANDOM", neighborhood_size, 1);
	// play complete reference game on each process
	reference_grid.initialize(grid_length, MPI_COMM_SELF, "RANDOM", neighborhood_size, 0);


	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (uint64_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	if (!grid.set_geometry(geom_params)) {
		cerr << "Couldn't set grid geometry" << endl;
		exit(EXIT_FAILURE);
	}
	if (!reference_grid.set_geometry(geom_params)) {
		cerr << "Couldn't set reference grid geometry" << endl;
		exit(EXIT_FAILURE);
	}


	initialize(grid, grid_length[0]);
	initialize(reference_grid, grid_length[0]);

	// every process outputs the game state into its own file
	string basename("tests/game_of_life/unrefined2d_");
	basename.append(1, direction).append("_").append(boost::lexical_cast<string>(rank)).append("_");

	// visualize the game with visit -o game_of_life_test.visit
	ofstream visit_file;
	if (save_results && rank == 0) {
		string visit_file_name("tests/game_of_life/unrefined2d_");
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

		Refine<Stretched_Cartesian_Geometry>::refine(grid, int(grid_length[0]), step, comm_size);

		grid.balance_load();
		grid.update_copies_of_remote_neighbors();

		if (verbose && rank == 0) {
			cout << step << " ";
			cout.flush();
		}

		if (save_results) {
			// write the game state into a file named according to the current time step
			string output_name(basename);
			output_name.append(boost::lexical_cast<string>(step)).append(".vtk");
			save(output_name, rank, grid);

			// visualize the game with visit -o game_of_life_test.visit
			if (rank == 0) {
				for (int process = 0; process < comm_size; process++) {
					visit_file << "unrefined2d_"
						<< direction << "_"
						<< process << "_"
						<< boost::lexical_cast<string>(step)
						<< ".vtk"
						<< endl;
				}
			}
		}

		get_live_neighbors(grid);
		get_live_neighbors(reference_grid);

		// verify refined/unrefined game
		for (const auto& cell: grid.cells) {
			uint64_t reference_cell_id = cell.id;
			const int refinement_level = grid.get_refinement_level(cell.id);
			if (refinement_level > 0) {
				reference_cell_id = grid.mapping.get_parent(cell.id);
			}

			const auto* const reference_data = reference_grid[reference_cell_id];
			if (reference_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for reference cell " << reference_cell_id
					<< std::endl;
				abort();
			}

			if (cell.data->data[0] != reference_data->data[0]) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell's " << cell.id << " life doesn't agree with reference."
					<< std::endl;
				abort();
			}
		}
	}

	if (rank == 0) {
		if (save_results) {
			visit_file.close();
		}

		if (verbose) {
			cout << endl;
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

