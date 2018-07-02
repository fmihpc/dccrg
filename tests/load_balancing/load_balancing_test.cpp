/*
Tests a Zoltan load balancing algorithm on a refined grid in 3d.
Visualize the results using for example VisIt.
*/

#include "algorithm"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "vector"
#include "sstream"
#include "utility"

#include "boost/program_options.hpp"
#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"
#include "../../dccrg_cartesian_geometry.hpp"

using namespace std;
using namespace dccrg;

struct Cell {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) this, 0, MPI_BYTE);
	}
};

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

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	/*
	Options
	*/
	uint64_t x_length, y_length, z_length;
	unsigned int neighborhood_length, refine_n, iterations;
	string load_balancing_method;
	vector<string> hier_lb_methods;
	vector<int> hier_lb_procs_per_level;
	bool save;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("x_length",
			boost::program_options::value<uint64_t>(&x_length)->default_value(10),
			"Create a grid with arg number of unrefined cells in the x direction")
		("y_length",
			boost::program_options::value<uint64_t>(&y_length)->default_value(10),
			"Create a grid with arg number of unrefined cells in the y direction")
		("z_length",
			boost::program_options::value<uint64_t>(&z_length)->default_value(10),
			"Create a grid with arg number of unrefined cells in the z direction")
		("load_balancing_method",
			boost::program_options::value<string>(&load_balancing_method)->default_value("HYPERGRAPH"),
			"Use arg as the load balancing method (supported values: NONE, BLOCK, RANDOM, RCB, RIB, HSFC, GRAPH, HYPERGRAPH, HIER)")
		("iterations",
			boost::program_options::value<unsigned int>(&iterations)->default_value(10),
			"First randomize processes of cells then balance the load using given options arg times")
		("hier_lb_methods",
			boost::program_options::value<vector<string> >(&hier_lb_methods)->composing(),
			"Load balancing method used by HIER is specified here, once per number of hierarchies, default HYPERGRAPH, number of these must equal number of hier_lb_procs_per_level (for example --hier_... RCB --hier_... RIB --hier_... RANDOM for three hierarhies)")
		("hier_lb_procs_per_level",
			boost::program_options::value<vector<int> >(&hier_lb_procs_per_level)->composing(),
			"Number of processes per hierarchy with HIER load balancing method, default 1, number of these must equal number of hier_lb_methods (for example --hier_... 12 --hier_... 4 --hier_... 2 for three hierarchies)")
		("refine_n",
			boost::program_options::value<unsigned int>(&refine_n)->default_value(0),
			"Refine cells arg times")
		("save",
			boost::program_options::value<bool>(&save)->default_value(false),
			"Save the grid after every iteration")
		("neighborhood_length",
			boost::program_options::value<unsigned int>(&neighborhood_length)->default_value(1),
			"Size of a cell's neighborhood in cells of equal size (0 means only face neighbors are neighbors)");

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

	// check for invalid options
	if (load_balancing_method == "HIER") {
		if (hier_lb_methods.size() == 0) {
			if (rank == 0) {
				cerr << "At least one load balancing method for HIER partitioning must be given" << endl;
			}
			return EXIT_FAILURE;
		}

		if (hier_lb_methods.size() != hier_lb_procs_per_level.size()) {
			if (rank == 0) {
				cerr << "Number of load balancing methods for HIER must equal the number of processes per partition options" << endl;
			}
			return EXIT_FAILURE;
		}

		for (unsigned int i = 0; i < hier_lb_procs_per_level.size(); i++) {
			if (hier_lb_procs_per_level[i] <= 0) {
				if (rank == 0) {
					cerr << "Processes per partition must be a positive number" << endl;
				}
				return EXIT_FAILURE;
			}
		}
	}

	Dccrg<Cell, Cartesian_Geometry> grid;
	grid
		.set_initial_length({x_length, y_length, z_length})
		.set_neighborhood_length(neighborhood_length)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method(load_balancing_method)
		.initialize(comm)
		.set_geometry({
			{-0.5, -0.5, -0.5},
			{1.0 / x_length, 1.0 / y_length, 1.0 / z_length}
		});

	// set load balancing options
	if (load_balancing_method == "HIER") {
		for (unsigned int i = 0; i < hier_lb_methods.size(); i++) {
			grid.add_partitioning_level(hier_lb_procs_per_level[i]);
			grid.add_partitioning_option(i, "LB_METHOD", hier_lb_methods[i]);
		}
	}

	auto cells = grid.get_cells();
	for (unsigned int i = 0; i < refine_n; i++) {
		for (const auto& cell: cells) {
			const std::array<double, 3> cell_center = grid.geometry.get_center(cell);
			double
				x = cell_center[0],
				y = cell_center[1],
				z = cell_center[2];

			if (sqrt(x * x + y * y + z * z) < 0.1) {
				grid.refine_completely(cell);
			}
		}
		grid.stop_refining();
		cells = grid.get_cells();
	}

	uint64_t local_cells = cells.size(), total_cells = 0;
	MPI_Reduce(&local_cells, &total_cells, 1, MPI_UINT64_T, MPI_SUM, 0, comm);
	if (rank == 0) {
		cout << "Total cells after refining: " << total_cells << endl;
	}


	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "load_balancing_test_" << rank << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (rank == 0) {
		visit_file.open("load_balancing_test.visit");
		visit_file << "!NBLOCKS " << comm_size << endl;
	}

	for (unsigned int step = 0; step < iterations; step++) {

		// scatter cells to random processes
		cells = grid.get_cells();
		for (const auto& cell: cells) {
			grid.pin(cell, rand() % comm_size);
		}
		grid.balance_load(false);
		grid.unpin_all_cells();

		// balance the load using given options
		grid.balance_load();
		/*grid.start_remote_neighbor_data_update();
		grid.wait_neighbor_data_update();*/
		cells = grid.get_cells();

		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
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
		if (rank == 0) {
			for (int process = 0; process < comm_size; process++) {
				visit_file << "load_balancing_test_" << process << "_"
					<< step_string.str() << suffix.str()
					<< endl;
			}
		}


		// write the grid into a file
		grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			const auto* const neighbors = grid.get_neighbors_of(cell);
			outfile << neighbors->size() << endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (size_t i = 0; i < cells.size(); i++) {
			outfile << rank << endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			outfile << cell << endl;
		}
		outfile.close();
	}

	if (rank == 0) {
		cout << "Passed" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

