/*
Tests dont_refine function.

Copyright 2010, 2011, 2012, 2013, 2014,
2018, 2019 Finnish Meteorological Institute
*/

#include "algorithm"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "sstream"
#include "unistd.h"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"
#include "../../dccrg_mpi_support.hpp"


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

	{
	Dccrg<int> grid; grid
		.set_initial_length({3, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	grid.dont_refine(1);
	grid.refine_completely(1);
	grid.stop_refining();

	uint64_t failed = 0;
	for (const auto& cell: grid.local_cells()) {
		if (cell.id == 4) {
			cerr << "Cell 1 was refined directly." << std::endl;
			failed = 1;
		}
	}
	if (All_Reduce()(failed, comm) > 0) {
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	}

	{
	Dccrg<int> grid; grid
		.set_initial_length({3, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	grid.pin(1, 0);
	grid.pin(2, 1);
	grid.balance_load(false);

	grid.refine_completely(2);
	grid.stop_refining();
	grid.dont_refine(1);
	grid.refine_completely(6);
	grid.stop_refining();

	uint64_t failed = 0;
	for (const auto& cell: grid.local_cells()) {
		if (cell.id == 4) {
			cerr << "Cell 1 was refined from face neighbor" << std::endl;
			failed = 1;
		}
	}
	if (All_Reduce()(failed, comm) > 0) {
		MPI_Finalize();
		return EXIT_FAILURE;
	}


	grid.dont_refine(1);
	grid.refine_completely(7);
	grid.stop_refining();

	failed = 0;
	for (const auto& cell: grid.local_cells()) {
		if (cell.id == 4) {
			cerr << "Cell 1 was refined from further neighbor" << std::endl;
			failed = 1;
		}
	}
	if (All_Reduce()(failed, comm) > 0) {
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
