/*
Tests the speed of refining the grid in 3-d by refining random cells until enough cell exist
*/

#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "functional"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg.hpp"


using namespace std;
using namespace boost::mpi;
using namespace dccrg;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	/*
	Options
	*/
	uint64_t x_length, y_length, z_length, maximum_cells, refines_per_process;
	int maximum_refinement_level, neighborhood_size;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("x_length",
			boost::program_options::value<uint64_t>(&x_length)->default_value(22),
			"Create a grid with arg number of unrefined cells in the x direction")
		("y_length",
			boost::program_options::value<uint64_t>(&y_length)->default_value(22),
			"Create a grid with arg number of unrefined cells in the y direction")
		("z_length",
			boost::program_options::value<uint64_t>(&z_length)->default_value(22),
			"Create a grid with arg number of unrefined cells in the z direction")
		("maximum_refinement_level",
			boost::program_options::value<int>(&maximum_refinement_level)->default_value(-1),
			"Maximum refinement level of the grid (0 == not refined, -1 == maximum possible for given lengths)")
		("maximum_cells",
			boost::program_options::value<uint64_t>(&maximum_cells)->default_value(1061208),
			"Stop refining after grid has arg number of total cells")
		("refines_per_process",
			boost::program_options::value<uint64_t>(&refines_per_process)->default_value(1000),
			"Refine at most arg cells at a time per process (not including induced refines)")
		("neighborhood_size",
			boost::program_options::value<int>(&neighborhood_size)->default_value(1),
			"Size of a cell's neighborhood in cells of equal size (0 means only face neighbors are neighbors)");

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

	// initialize
	Dccrg<int> grid;

	if (!grid.set_geometry(
		x_length, y_length, z_length,
		0.0, 0.0, 0.0,
		1.0 / x_length, 1.0 / y_length, 1.0 / z_length
	)) {
		cerr << "Couldn't set grid geometry" << endl;
		return EXIT_FAILURE;
	}

	grid.initialize(comm, "RCB", neighborhood_size, maximum_refinement_level);
	grid.balance_load();

	vector<uint64_t> cells = grid.get_cells();

	clock_t before, after, total = 0;

	do {
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		before = clock();
		uint64_t refined = 0;
		for (int i = 0; i < int(cells.size() / 15) && refined < refines_per_process; i++) {
			grid.refine_completely(cells[i]);
			refined++;
		}
		vector<uint64_t> new_cells = grid.stop_refining();
		after = clock();
		total += after - before;
		cout << "Proc " << comm.rank() << ": " << new_cells.size() << " new cells" << endl;

		cells = grid.get_cells();
	} while (all_reduce(comm, cells.size(), plus<uint64_t>()) < maximum_cells);

	cout << "Process " << comm.rank()
		<< ": " << cells.size() - x_length * y_length * z_length
		<< " new cells created in " << double(total) / CLOCKS_PER_SEC
		<< " s, (" << (cells.size() - x_length * y_length * z_length) / (double(total) / CLOCKS_PER_SEC)
		<< " new cells / s)"
		<< endl;

	return EXIT_SUCCESS;
}

