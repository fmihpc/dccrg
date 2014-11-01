/*
Tests loading of the grid
*/

#include "boost/lexical_cast.hpp"
#include "cstdlib"
#include "ctime"
#include "iostream"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"

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

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	// initialize grid
	Dccrg<CellData> grid;
	const std::array<uint64_t, 3> grid_length = {{10, 10, 10}};
	grid.initialize(grid_length, comm, "RANDOM", 2, 5);
	grid.balance_load();

	const std::vector<uint64_t> cells_to_load{
		// refine the corners twice
		9001,
		9040,
		10561,
		10600,
		71401,
		71440,
		72961,
		73000,
		// throw in some cells with refinement level > 5
		37449000 + 100,
		37449000 + 10000,
		37449000 + 1000000,
		37449000 + 100000000
	};

	grid.load_cells(cells_to_load);

	// save the grid
	const std::string
		base_output_name("load_"),
		output_name_suffix(".vtk");

	grid.write_vtk_file(
		base_output_name
		+ boost::lexical_cast<std::string>(rank)
		+ output_name_suffix
	);

	// visualize with visit -o load.visit
	if (rank == 0) {
		ofstream visit_file("load.visit");

		visit_file << "!NBLOCKS " << comm_size << endl;
		for (int i = 0; i < comm_size; i++) {
			visit_file << base_output_name << i << output_name_suffix << "\n";
		}

		visit_file.close();
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

