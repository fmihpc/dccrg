/*
Tests copy construction of dccrg.

Copyright 2015, 2016, 2018, 2019 Finnish Meteorological Institute
Copyright 2015, 2016 Ilja Honkonen

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
#include "ctime"
#include "iostream"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_stretched_cartesian_geometry.hpp"

using namespace std;

/*!
Cell data in grid1.
*/
struct Cell1 {
	int data = -1;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) &(this->data), 1, MPI_INT);
	}
};

/*!
Cell data in grid2.
*/
struct Cell2 {
	double data = -2;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->data), 1, MPI_DOUBLE);
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
		cerr << "Zoltan_Initialize failed" << endl;
		return EXIT_FAILURE;
	}

	dccrg::Dccrg<Cell1, dccrg::Cartesian_Geometry> grid1;

	grid1
		.set_initial_length({10, 10, 10})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(0)
		.initialize(comm)
		.set_geometry(dccrg::Cartesian_Geometry::Parameters{{0,0,0}, {1,1,1}});

	// check that remote neighbor update works in original grid
	// data in grid1 == process rank
	const vector<uint64_t> cells1 = grid1.get_cells();
	for (const auto& cell: cells1) {
		auto* const cell_data = grid1[cell];
		cell_data->data = rank;
	}

	grid1.update_copies_of_remote_neighbors();

	const std::unordered_set<uint64_t>& remote_neighbors1
		= grid1.get_remote_cells_on_process_boundary_internal();

	const std::unordered_map<uint64_t, uint64_t>& cell_process1
		= grid1.get_cell_process();

	for (const auto& cell: remote_neighbors1) {
		auto* const cell_data = grid1[cell];
		if (cell_data == NULL) {
			cerr << "No data for cell " << cell << " in grid1" << endl;
			abort();
		}

		if (cell_data->data != (int) cell_process1.at(cell)) {
			cerr << "Data of cell " << cell
				<< " in grid1 incorrect: " << cell_data->data
				<< ", should be " << cell_process1.at(cell)
				<< endl;
			abort();
		}
	}

	// check copy constructor
	dccrg::Dccrg<Cell2, dccrg::Stretched_Cartesian_Geometry> grid2(grid1);
	const auto cells2 = grid2.get_cells();

	if (cells1.size() != cells2.size()) {
		cerr << "Rank " << rank << ": Number of cells doesn't match" << endl;
		abort();
	}

	for (const auto& cell: cells2) {
		auto* const cell_data = grid2[cell];
		if (cell_data->data != -2) {
			cerr << "Rank " << rank
				<< ": Wrong data in cell " << cell
				<< ": " << cell_data->data
				<< ", should be -2"
				<< endl;
			abort();
		}
	}

	// check that remote neighbor update works in original grid
	// data in grid2 == 2 * process rank
	for (const auto& cell: cells2) {
		auto* const cell_data = grid2[cell];
		cell_data->data = 2 * rank;
	}

	grid2.update_copies_of_remote_neighbors();

	const std::unordered_set<uint64_t>& remote_neighbors2
		= grid2.get_remote_cells_on_process_boundary_internal();

	const std::unordered_map<uint64_t, uint64_t>& cell_process2
		= grid2.get_cell_process();

	for (const auto& cell: remote_neighbors2) {
		auto* const cell_data = grid2[cell];
		if (cell_data == NULL) {
			cerr << "No data for cell " << cell << " in grid2" << endl;
			abort();
		}

		if (cell_data->data != 2 * cell_process2.at(cell)) {
			cerr << "Data of cell " << cell
				<< " in grid2 incorrect: " << cell_data->data
				<< ", should be " << 2 * cell_process2.at(cell)
				<< endl;
			abort();
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
