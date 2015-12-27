/*
Tests the scalability of initializing the grid
*/

#include "cstdlib"
#include "ctime"

#include "boost/program_options.hpp"
#include "mpi.h"
#include "zoltan.h"

#include "dccrg_cartesian_geometry.hpp"
#include "dccrg.hpp"


struct CellData {

	double variables[3];

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->variables), 3, MPI_DOUBLE);
	}
};


using namespace std;
using namespace dccrg;


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = -1;
	MPI_Comm_rank(comm, &rank);

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	/*
	Options
	*/
	uint64_t x_length, y_length, z_length;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("x_length",
			boost::program_options::value<uint64_t>(&x_length)->default_value(100),
			"Create a grid with arg number of unrefined cells in the x direction")
		("y_length",
			boost::program_options::value<uint64_t>(&y_length)->default_value(100),
			"Create a grid with arg number of unrefined cells in the y direction")
		("z_length",
			boost::program_options::value<uint64_t>(&z_length)->default_value(100),
			"Create a grid with arg number of unrefined cells in the z direction");

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

	Dccrg<CellData> grid;
	const std::array<uint64_t, 3> grid_length = {{x_length, y_length, z_length}};

	clock_t before = clock();
	grid.initialize(grid_length, comm, "RCB", 1, 0);
	clock_t after = clock();

	cout << "Process " << rank
		<< ": grid initialization took " << double(after - before) / CLOCKS_PER_SEC
		<< " seconds (total grid size " << x_length * y_length * z_length << ")"
		<< endl;

	MPI_Finalize();

	return EXIT_SUCCESS;
}

