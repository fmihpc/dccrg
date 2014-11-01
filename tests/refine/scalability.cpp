/*
Tests the speed of refining the grid in 3-d by refining random cells until enough cell exist
*/

#include "cstdlib"
#include "ctime"
#include "fstream"
#include "functional"
#include "iostream"
#include "sstream"

#include "boost/program_options.hpp"
#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"


using namespace std;
using namespace dccrg;

struct Cell {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(this, 0, MPI_BYTE);
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

	/*
	Options
	*/
	uint64_t x_length, y_length, z_length, maximum_cells, max_refines;
	int maximum_refinement_level, neighborhood_size;
	double max_refine_fraction;
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
		("max_refines",
			boost::program_options::value<uint64_t>(&max_refines)->default_value(1000),
			"Refine at most arg local cells at a time per process (not including induced refines)")
		("max_refine_fraction",
			boost::program_options::value<double>(&max_refine_fraction)->default_value(0.1),
			"Refine at most fraction arg of local cells at a time per process (not including induced refines)")
		("neighborhood_size",
			boost::program_options::value<int>(&neighborhood_size)->default_value(1),
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

	// initialize
	Dccrg<Cell> grid;

	const std::array<uint64_t, 3> grid_length = {{x_length, y_length, z_length}};
	grid.initialize(grid_length, comm, "RCB", neighborhood_size, maximum_refinement_level);
	grid.balance_load();

	vector<uint64_t> cells = grid.get_cells();
	const uint64_t initial_cells = cells.size();

	clock_t before, after, total = 0;

	size_t local_cells = cells.size(), total_cells = 0;
	do {
		random_shuffle(cells.begin(), cells.end());

		// refine a fraction of all cells each round
		before = clock();
		uint64_t refined = 0;
		for (int i = 0; i < int(cells.size() * max_refine_fraction) && refined < max_refines; i++) {
			grid.refine_completely(cells[i]);
			refined++;
		}
		vector<uint64_t> new_cells = grid.stop_refining();
		after = clock();
		total += after - before;
		cout << "Proc " << rank << ": " << new_cells.size() << " new cells" << endl;

		cells = grid.get_cells();

		local_cells = cells.size();
		total_cells = 0;
		MPI_Allreduce(&local_cells, &total_cells, 1, MPI_UINT64_T, MPI_SUM, comm);
	} while (total_cells < maximum_cells);

	cout << "Process " << rank
		<< ": " << cells.size() - initial_cells
		<< " new cells created in " << double(total) / CLOCKS_PER_SEC
		<< " s, (" << (cells.size() - initial_cells) / (double(total) / CLOCKS_PER_SEC)
		<< " new cells / s)"
		<< endl;

	MPI_Finalize();

	return EXIT_SUCCESS;
}

