/*
User neighborhood version of advection test of dccrg.

Copyright 2012, 2013, 2014 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "algorithm"
#include "array"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "cstdlib"
#include "iomanip"
#include "iostream"
#include "string"
#include "unordered_set"
#include "utility"
#include "vector"

#include "mpi.h"
#include "zoltan.h"
#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

#include "adapter.hpp"
#include "cell.hpp"
#include "initialize.hpp"
#include "save.hpp"
#include "solve.hpp"


using namespace std;
using namespace dccrg;


#define GRID_START_X 0
#define GRID_START_Y 0
#define GRID_START_Z 0

#define GRID_END_X 1.0
#define GRID_END_Y 1.0
#define GRID_END_Z 1.0

#define NEIGHBORHOOD_ID 1


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
	//char direction;
	bool verbose = false;
	char direction = 'z';
	unsigned int cells;
	int max_ref_lvl, save_n, balance_n, adapt_n;
	string load_balancing_method;
	double tmax, relative_diff, unrefine_sensitivity, diff_threshold, cfl;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("cells",
			boost::program_options::value<unsigned int>(&cells)->default_value(400),
			"Total number of unrefined cells at the start of the simulation")
		("max-ref-lvl",
			boost::program_options::value<int>(&max_ref_lvl)->default_value(2),
			"Maximum refinement level of cells in the grid (0 means unrefined)")
		("relative-diff",
			boost::program_options::value<double>(&relative_diff)->default_value(0.025),
			"Maximum relative difference in variables for a cell which to keep at maximum refinement level")
		("diff-threshold",
			boost::program_options::value<double>(&diff_threshold)->default_value(0.25),
			"TODO")
		("unrefine-sensitivity",
			boost::program_options::value<double>(&unrefine_sensitivity)->default_value(0.5),
			"TODO")
		("save-n",
			boost::program_options::value<int>(&save_n)->default_value(0),
			"Save results every arg'th time step (0 = only save initial and final result,"
			" -1 = never save)")
		("tmax",
			boost::program_options::value<double>(&tmax)->default_value(25.5),
			"Duration of run in seconds")
		("load-balancing-method",
			boost::program_options::value<string>(&load_balancing_method)->default_value("RCB"),
			"Use arg as load balancing method")
		("balance-n",
			boost::program_options::value<int>(&balance_n)->default_value(25),
			"Balance computational load every argth time step (-1 == never balance load)")
		("adapt-n",
			boost::program_options::value<int>(&adapt_n)->default_value(1),
			"Check for grid adaptation every argth timestep")
		/*("direction",
			boost::program_options::value<char>(&direction)->default_value('z'),
			"Create a 2d grid with normal into direction arg (x, y or z)")*/
		("cfl",
			boost::program_options::value<double>(&cfl)->default_value(0.5),
			"Fraction of ... to use (0..1)")
		("verbose", "Print information during the simulation");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	// check simulation parameters
	if (save_n < -1) {
		cerr << "save_n must be >= -1" << endl;
		return EXIT_FAILURE;
	}

	if (balance_n < -1) {
		cerr << "balance_n must be >= -1" << endl;
		return EXIT_FAILURE;
	}

	if (cfl < 0 || cfl > 1) {
		cerr << "cfl must be >= 0 and <= 1" << endl;
		return EXIT_FAILURE;
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

	// transform user-given parameters to internal units
	cells = (unsigned int) round(sqrt(double(cells)));

	// initialize grid
	Dccrg<Cell, Cartesian_Geometry> grid;
	Cartesian_Geometry::Parameters geom_params;

	std::array<uint64_t, 3> grid_length = {{0, 0, 0}};

	switch (direction) {
	case 'x':
		grid_length[0] = 1;
		grid_length[1] = cells;
		grid_length[2] = cells;

		grid.initialize(
			grid_length,
			comm,
			load_balancing_method.c_str(),
			2,
			max_ref_lvl,
			false, true, true
		);

		geom_params.start[0] = GRID_START_X;
		geom_params.start[1] = GRID_START_Y;
		geom_params.start[2] = GRID_START_Z;
		geom_params.level_0_cell_length[0] = GRID_END_X / cells;
		geom_params.level_0_cell_length[1] = GRID_END_Y / cells;
		geom_params.level_0_cell_length[2] = GRID_END_Z / cells;
		break;

	case 'y':
		grid_length[0] = cells;
		grid_length[1] = 1;
		grid_length[2] = cells;

		grid.initialize(
			grid_length,
			comm,
			load_balancing_method.c_str(),
			2,
			max_ref_lvl,
			true, false, true
		);

		geom_params.start[0] = GRID_START_X;
		geom_params.start[1] = GRID_START_Y;
		geom_params.start[2] = GRID_START_Z;
		geom_params.level_0_cell_length[0] = GRID_END_X / cells;
		geom_params.level_0_cell_length[1] = GRID_END_Y / cells;
		geom_params.level_0_cell_length[2] = GRID_END_Z / cells;
		break;

	case 'z':
		grid_length[0] = cells;
		grid_length[1] = cells;
		grid_length[2] = 1;

		grid.initialize(
			grid_length,
			comm,
			load_balancing_method.c_str(),
			2,
			max_ref_lvl,
			true, true, false
		);

		geom_params.start[0] = GRID_START_X;
		geom_params.start[1] = GRID_START_Y;
		geom_params.start[2] = GRID_START_Z;
		geom_params.level_0_cell_length[0] = GRID_END_X / cells;
		geom_params.level_0_cell_length[1] = GRID_END_Y / cells;
		geom_params.level_0_cell_length[2] = GRID_END_Z / cells;
		break;

	default:
		cerr << "Unsupported direction given: " << direction << endl;
		break;
	}

	if (!grid.set_geometry(geom_params)) {
		cerr << __FILE__ << ":" << __LINE__ << ": Couldn't set grid geometry" << endl;
		abort();
	}

	// create the neighborhood
	typedef dccrg::Types<3>::neighborhood_item_t neigh_t;
	const std::vector<neigh_t> neighborhood{
		{ 0,  0, -1},
		{ 0, -1,  0},
		{-1,  0,  0},
		{ 1,  0,  0},
		{ 0,  1,  0},
		{ 0,  0,  1}
	};

	if (!grid.add_neighborhood(NEIGHBORHOOD_ID, neighborhood)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set neighborhood"
			<< std::endl;
		abort();
	}

	if (balance_n > -1) {
		grid.balance_load();
	}

	// apply initial condition 1st time for prerefining the grid
	Initialize()(grid);

	std::unordered_set<uint64_t> cells_to_refine, cells_not_to_unrefine, cells_to_unrefine;

	uint64_t created_cells = 0, removed_cells = 0;

	// prerefine up to maximum refinement level
	for (int ref_lvl = 0; ref_lvl < max_ref_lvl; ref_lvl++) {
		Adapter().check_for_adaptation(
			relative_diff / grid.get_maximum_refinement_level(),
			diff_threshold,
			unrefine_sensitivity,
			cells_to_refine,
			cells_not_to_unrefine,
			cells_to_unrefine,
			grid
		);

		const std::pair<uint64_t, uint64_t> adapted_cells
			= Adapter().adapt_grid(
				cells_to_refine,
				cells_not_to_unrefine,
				cells_to_unrefine,
				grid
			);
		created_cells += adapted_cells.first;
		removed_cells += adapted_cells.second;

		// apply initial condition on a finer grid
		Initialize()(grid);
	}

	double dt = Solver().max_time_step(comm, grid);
	if (verbose && rank == 0) {
		cout << "Initial timestep: " << dt << endl;
	}

	// save initial state
	#ifndef DEBUG
	const string base_output_name("2d_");
	#else
	const string base_output_name("2d_debug_");
	#endif
	unsigned int files_saved = 0;
	if (save_n > -1) {
		if (verbose && rank == 0) {
			cout << "Saving initial state of simulation" << endl;
		}
		Save()(File_Namer()(0, base_output_name), comm, grid);
		files_saved++;
	}

	if (verbose && rank == 0) {
		cout << "Starting simulation" << endl;
	}

	vector<uint64_t> inner_cells = grid.get_local_cells_not_on_process_boundary(NEIGHBORHOOD_ID);
	vector<uint64_t> outer_cells = grid.get_local_cells_on_process_boundary(NEIGHBORHOOD_ID);

	// record solution time for inner cells and amount of neighbor data received
	double inner_solve_time = 0, outer_solve_time = 0, neighbor_receive_size = 0;

	double time = 0;
	unsigned int step = 0;
	while (time < tmax) {

		if (verbose && rank == 0) {
			cout << "Simulation time: " << time << endl;
		}

		grid.start_remote_neighbor_copy_updates(NEIGHBORHOOD_ID);

		// solve inner cells
		const double inner_solve_start = MPI_Wtime();
		Solver().calculate_fluxes(cfl * dt, inner_cells, grid);
		inner_solve_time += MPI_Wtime() - inner_solve_start;

		// wait for remote neighbor data
		grid.wait_remote_neighbor_copy_update_receives(NEIGHBORHOOD_ID);

		// solve outer cells
		const double outer_solve_start = MPI_Wtime();
		Solver().calculate_fluxes(cfl * dt, outer_cells, grid);
		outer_solve_time += MPI_Wtime() - outer_solve_start;

		// wait until local data has been sent
		grid.wait_remote_neighbor_copy_update_sends();

		neighbor_receive_size +=
			4 * sizeof(double)
			* (grid.get_number_of_update_receive_cells() + grid.get_number_of_update_send_cells());

		/*
		Starting from this point local cells and copies of remote cells have
		data from the same timestep (flux and max_diff isn't transferred).
		*/

		// check where to adapt the grid
		if (adapt_n > 0 && step % adapt_n == 0) {

			if (verbose && rank == 0) {
				cout << "Checking which cells to adapt in the grid" << endl;
			}

			Adapter().check_for_adaptation(
				relative_diff / grid.get_maximum_refinement_level(),
				diff_threshold,
				unrefine_sensitivity,
				cells_to_refine,
				cells_not_to_unrefine,
				cells_to_unrefine,
				grid
			);
		}

		// save simulation state
		if (save_n > 0 && step % save_n == 0) {
			if (verbose && rank == 0) {
				cout << "Saving simulation at " << time << endl;
			}
			Save()(File_Namer()(time, base_output_name), comm, grid);
			files_saved++;
		}

		/*
		Up to this point local cells and copies of remote cells have
		data from the same timestep (variables not fluxes which aren't transferred anyway).
		*/

		// apply fluxes
		Solver().apply_fluxes(grid);

		// adapt the grid
		if (adapt_n > 0 && step % adapt_n == 0) {

			if (verbose && rank == 0) {
				cout << "Adapting grid" << endl;
			}

			const std::pair<uint64_t, uint64_t> adapted_cells
				= Adapter().adapt_grid(
					cells_to_refine,
					cells_not_to_unrefine,
					cells_to_unrefine,
					grid
				);
			created_cells += adapted_cells.first;
			removed_cells += adapted_cells.second;

			inner_cells = grid.get_local_cells_not_on_process_boundary(NEIGHBORHOOD_ID);
			outer_cells = grid.get_local_cells_on_process_boundary(NEIGHBORHOOD_ID);

			// update maximum allowed time step
			dt = Solver().max_time_step(comm, grid);
			if (verbose && rank == 0) {
				cout << "New timestep: " << dt << endl;
			}
		}

		// balance load
		if (balance_n > 0 && step % balance_n == 0) {

			if (verbose && rank == 0) {
				cout << "Balancing load" << endl;
			}

			grid.balance_load();

			inner_cells = grid.get_local_cells_not_on_process_boundary(NEIGHBORHOOD_ID);
			outer_cells = grid.get_local_cells_on_process_boundary(NEIGHBORHOOD_ID);
		}

		step++;
		time += dt;
	}

	if (save_n > -1) {
		if (verbose && rank == 0) {
			cout << "Saving final state of simulation" << endl;
		}

		Save()(File_Namer()(tmax, base_output_name), comm, grid);
		files_saved++;
	}

	// gather statistics about solving time and tranferred data
	double min_inner_solve_time = 0, max_inner_solve_time = 0, total_inner_solve_time = 0,
		min_outer_solve_time = 0, max_outer_solve_time = 0, total_outer_solve_time = 0,
		min_receive_size = 0, max_receive_size = 0, total_receive_size = 0,
		// fractions of the above
		min_fraction = 0, max_fraction = 0, total_fraction = 0;
	uint64_t total_created_cells = 0, total_removed_cells = 0;

	MPI_Reduce(&inner_solve_time, &min_inner_solve_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&inner_solve_time, &max_inner_solve_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&inner_solve_time, &total_inner_solve_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	MPI_Reduce(&outer_solve_time, &min_outer_solve_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&outer_solve_time, &max_outer_solve_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&outer_solve_time, &total_outer_solve_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	MPI_Reduce(&neighbor_receive_size, &min_receive_size, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&neighbor_receive_size, &max_receive_size, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&neighbor_receive_size, &total_receive_size, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	MPI_Reduce(&created_cells, &total_created_cells, 1, MPI_UINT64_T, MPI_SUM, 0, comm);
	MPI_Reduce(&removed_cells, &total_removed_cells, 1, MPI_UINT64_T, MPI_SUM, 0, comm);
	double fraction = neighbor_receive_size / inner_solve_time;
	MPI_Reduce(&fraction, &min_fraction, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&fraction, &max_fraction, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&fraction, &total_fraction, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

	if (rank == 0) {
		cout << endl;
		cout << "Processes used: " << comm_size << endl;
		cout << "Initial grid size: " << cells * cells << endl;
		cout << "Total timesteps calculated: " << tmax << endl;
		cout << "Total files saved: " << files_saved << endl;
		cout << "Total created and removed cells: "
			<< total_created_cells << ", "
			<< total_removed_cells << endl;
		cout << "Inner cell solution time / step (s, avg, max, min):          "
			<< total_inner_solve_time / comm_size / tmax << "\t"
			<< max_inner_solve_time / tmax << "\t"
			<< min_inner_solve_time / tmax << endl;
		cout << "Outer cell solution time / step (s, avg, max, min):          "
			<< total_outer_solve_time / comm_size / tmax << "\t"
			<< max_outer_solve_time / tmax << "\t"
			<< min_outer_solve_time / tmax << endl;
		cout << "Remote neighbor data receive size / step (B, avg, max, min): "
			<< total_receive_size / comm_size / tmax << "\t"
			<< max_receive_size / tmax << "\t"
			<< min_receive_size / tmax << endl;
		cout << "Per process fractions of the above (B / s, avg, max, min):   "
			<< total_fraction / comm_size << "\t"
			<< max_fraction << "\t"
			<< min_fraction << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

