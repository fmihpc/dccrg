/*
Program for testing dccrg neighbor iterators with AMR.

Copyright 2013, 2014, 2015, 2016, 2018 Finnish Meteorological Institute
Copyright 2018 Ilja Honkonen

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

	// initialize grid
	Dccrg<Cell> grid;
	grid
		.set_initial_length({2, 1, 1})
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	if (grid.is_local(1)) {
		grid.refine_completely(1);
	}
	grid.stop_refining();
	grid.clear_refined_unrefined_data();
	set<uint64_t> all_cells{2, 3, 4, 7, 8, 11, 12, 15, 16};

	// do a few iterations with random load balancing
	for (int i = 0; i < 5; i++) {
		// check that local cell neighbor iterators work
		for (const auto& cell: grid.local_cells) {
			set<uint64_t> ref_neighbors_of{all_cells};
			ref_neighbors_of.erase(cell.id);
			if (cell.id == 3 or cell.id == 7 or cell.id == 11 or cell.id == 15) {
				ref_neighbors_of.erase(2);
			}

			set<uint64_t> neighbors_of;
			for (const auto& neighbor: cell.neighbors_of) {
				neighbors_of.insert(neighbor.id);
			}
			if (ref_neighbors_of != neighbors_of) {
				cerr << "FAILED" << endl;
				cout << "Wrong neighbors_of for cell " << cell.id
					<< "\nResult: ";
				for (const auto& n: neighbors_of) cout << n << ", ";
				cout << "\nReference: ";
				for (const auto& n: ref_neighbors_of) cout << n << ", ";
				cout << endl;
				return EXIT_FAILURE;
			}
		}

		grid.balance_load();
	}


	if (grid.is_local(3)) {
		grid.refine_completely(3);
	}
	grid.stop_refining();
	grid.clear_refined_unrefined_data();

	for (int i = 0; i < 5; i++) {
		// check that local cell neighbor iterators work
		for (const auto& cell: grid.local_cells) {
			set<uint64_t> ref_neighbors_of;

			if (cell.id > 18) {
				ref_neighbors_of = {19, 20, 27, 28, 51, 52, 59, 60};
				if (cell.id % 2 == 0) ref_neighbors_of.insert(4);
				if (cell.id > 50) ref_neighbors_of.insert(11);
				if (cell.id == 52 or cell.id == 60) ref_neighbors_of.insert(12);
				if (cell.id == 28 or cell.id == 60) ref_neighbors_of.insert(8);
				if (cell.id == 60) ref_neighbors_of.insert(16);
				if (cell.id == 59 or cell.id == 60) ref_neighbors_of.insert(15);
				if (
					cell.id == 27
					or cell.id == 28
					or cell.id == 59
					or cell.id == 60
				) ref_neighbors_of.insert(7);

			} else {

				ref_neighbors_of = {4, 5, 8, 9, 12, 13, 16, 17};
				if (cell.id % 4 == 1 or cell.id % 4 == 2) {
					for (const auto& c: {6, 10, 14, 18}) ref_neighbors_of.insert(c);
				}
				if (cell.id % 4 == 2) {
					for (const auto& c: {4, 8, 12, 16}) ref_neighbors_of.erase(c);
				}
				if (cell.id % 4 == 0 or cell.id % 4 == 3) {
					for (const auto& c: {7, 11, 15, 19, 20, 27, 28, 51, 52, 59, 60}) ref_neighbors_of.insert(c);
				}
				if (cell.id % 4 == 3) {
					for (const auto& c: {5, 9, 13, 17}) ref_neighbors_of.erase(c);
				}
			}
			ref_neighbors_of.erase(cell.id);

			set<uint64_t> neighbors_of;
			for (const auto& neighbor: cell.neighbors_of) {
				neighbors_of.insert(neighbor.id);
			}

			if (ref_neighbors_of != neighbors_of) {
				cerr << "FAILED" << endl;
				cout << "Wrong neighbors_of for cell " << cell.id
					<< "\nResult: ";
				for (const auto& n: neighbors_of) cout << n << ", ";
				cout << "\nReference: ";
				for (const auto& n: ref_neighbors_of) cout << n << ", ";
				cout << endl;
				return EXIT_FAILURE;
			}
		}

		grid.balance_load();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
