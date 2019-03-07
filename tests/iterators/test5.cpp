/*
Program for testing dccrg neighbor iterators with AMR.

Copyright 2019 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "cstdlib"
#include "iostream"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"

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

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}

	{
	Dccrg<Cell> grid;
	grid
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm)
		.balance_load();

	if (grid.is_local(1)) {
		grid.refine_completely(1);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(3)) {
		grid.refine_completely(3);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(19)) {
		grid.refine_completely(19);
	}
	grid.stop_refining();
	}

	{
	Dccrg<Cell> grid;
	grid
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm)
		.balance_load();

	if (grid.is_local(2)) {
		grid.refine_completely(2);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(18)) {
		grid.refine_completely(18);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(146)) {
		grid.refine_completely(146);
	}
	grid.stop_refining();
	}


	{
	Dccrg<Cell> grid;
	grid
		.set_initial_length({1, 2, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm)
		.balance_load();

	if (grid.is_local(1)) {
		grid.refine_completely(1);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(3)) {
		grid.refine_completely(3);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(19)) {
		grid.refine_completely(19);
	}
	grid.stop_refining();
	}

	{
	Dccrg<Cell> grid;
	grid
		.set_initial_length({1, 2, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm)
		.balance_load();

	if (grid.is_local(2)) {
		grid.refine_completely(2);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(18)) {
		grid.refine_completely(18);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(146)) {
		grid.refine_completely(146);
	}
	grid.stop_refining();
	}


	{
	Dccrg<Cell> grid;
	grid
		.set_initial_length({1, 1, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm)
		.balance_load();

	if (grid.is_local(1)) {
		grid.refine_completely(1);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(3)) {
		grid.refine_completely(3);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(19)) {
		grid.refine_completely(19);
	}
	grid.stop_refining();
	}

	{
	Dccrg<Cell> grid;
	grid
		.set_initial_length({1, 1, 2})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm)
		.balance_load();

	if (grid.is_local(2)) {
		grid.refine_completely(2);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(18)) {
		grid.refine_completely(18);
	}
	grid.stop_refining();
	grid.balance_load();
	if (grid.is_local(146)) {
		grid.refine_completely(146);
	}
	grid.stop_refining();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
