/*
Tests loading of the grid
*/

#include "boost/assign/list_of.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "cstdlib"
#include "ctime"
#include "iostream"
#include "zoltan.h"

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

	// initialize grid
	Dccrg<CellData> grid;
	const boost::array<uint64_t, 3> grid_length = {{10, 10, 10}};
	grid.initialize(grid_length, comm, "RANDOM", 2, 5);
	grid.balance_load();

	const std::vector<uint64_t> cells_to_load = boost::assign::list_of
		// refine the corners twice
		(9001)
		(9040)
		(10561)
		(10600)
		(71401)
		(71440)
		(72961)
		(73000)
		// throw in some cells with refinement level > 5
		(37449000 + 100)
		(37449000 + 10000)
		(37449000 + 1000000)
		(37449000 + 100000000);

	grid.load_cells(cells_to_load);

	// save the grid
	const std::string
		base_output_name("load_"),
		output_name_suffix(".vtk");

	grid.write_vtk_file(
		(base_output_name
		+ boost::lexical_cast<std::string>(comm.rank())
		+ output_name_suffix).c_str()
	);

	// visualize with visit -o load.visit
	if (comm.rank() == 0) {
		ofstream visit_file("load.visit");

		visit_file << "!NBLOCKS " << comm.size() << endl;
		for (int i = 0; i < comm.size(); i++) {
			visit_file << base_output_name << i << output_name_suffix << "\n";
		}

		visit_file.close();
	}

	if (comm.rank() == 0) {
		cout << "PASSED" << endl;
	}

	return EXIT_SUCCESS;
}

