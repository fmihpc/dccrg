/*
Tests the scalability of initializing the grid
*/

#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "zoltan.h"

#include "../../dccrg_cartesian_geometry.hpp"
#include "../../dccrg.hpp"


struct CellData {

	double variables[3];

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & variables[0] & variables[1] & variables[2];
	}
};


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
		if (comm.rank() == 0) {
			cout << options << endl;
		}
		comm.barrier();
		return EXIT_SUCCESS;
	}

	Dccrg<CellData> grid;
	const boost::array<uint64_t, 3> grid_length = {{x_length, y_length, z_length}};

	clock_t before = clock();
	grid.initialize(grid_length, comm, "RCB", 1, 0);
	clock_t after = clock();

	cout << "Process " << comm.rank()
		<< ": grid initialization took " << double(after - before) / CLOCKS_PER_SEC
		<< " seconds (total grid size " << x_length * y_length * z_length << ")"
		<< endl;

	return EXIT_SUCCESS;
}

