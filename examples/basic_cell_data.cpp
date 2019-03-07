/*
Example that uses basic types as cell data.
*/

#include "cstdlib"

#include "mpi.h"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;

int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	dccrg::Dccrg<int> grid;
	grid
		.set_initial_length({7, 13, 11})
		.set_maximum_refinement_level(1)
		.set_neighborhood_length(2)
		.set_load_balancing_method("RANDOM")
		.set_periodic(true, true, true)
		.initialize(MPI_COMM_WORLD)
		.balance_load();
	// set cell id as value for cell data
	for (const auto& cell: grid.local_cells()) {
		*cell.data = (int)cell.id;
	}
	// check that cell data is updated
	// correctly between processes
	grid.update_copies_of_remote_neighbors();
	for (const auto& cell: grid.local_cells()) {
		for (const auto& neighbor: cell.neighbors_of) {
			if (*neighbor.data != (int)neighbor.id) {
				std::cerr << __FILE__ "(" << __LINE__
					<< ") Wrong cell data for neighbor " << neighbor.id
					<< " of cell " << cell.id << ": " << *neighbor.data
					<< std::endl;
				abort();
			}
		}
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}
