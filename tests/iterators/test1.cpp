/*
Program for testing dccrg iterators.

Copyright 2013, 2014 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "boost/tuple/tuple.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg.hpp"


using namespace std;
using namespace boost;
using namespace dccrg;


class Cell
{
public:
	int unused;

	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

	template<typename Archiver> void serialize(
		Archiver& ar,
		const unsigned int /*version*/
	) {
		ar & unused;
	}

	#else // ifdef DCCRG_TRANSFER_USING_BOOST_MPI

	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhoo_id*/
	) {
		return boost::make_tuple(&(this->unused), 1, MPI_INT);
	}

	#endif // ifdef DCCRG_TRANSFER_USING_BOOST_MPI

};


struct is_larger {
	is_larger(int given_x) : x(given_x) {}
	bool operator()(int given_x)
	{
		return given_x > this->x;
	}
	int x;
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
	const boost::array<uint64_t, 3> grid_length = {{1000, 1, 1}};
	grid.initialize(grid_length, comm, "RANDOM", 3);

	// do a few iterations with random load balancing
	for (int i = 0; i < 5; i++) {

		grid.balance_load();

		const vector<uint64_t>
			inner_temp = grid.get_local_cells_not_on_process_boundary(),
			outer_temp = grid.get_local_cells_on_process_boundary();

		const boost::unordered_set<uint64_t>
			inner_cells(inner_temp.begin(), inner_temp.end()),
			outer_cells(outer_temp.begin(), outer_temp.end());

		// check that inner cell iterator works
		for (auto item = grid.begin_inner(); item != grid.end_inner(); item++) {
			const uint64_t cell = item->first;

			if (inner_cells.count(cell) == 0) {
				cout << "FAILED" << endl;
				cerr << "Cell " << cell
					<< " is not in inner cells"
					<< endl;
				abort();
			}

			if (outer_cells.count(cell) > 0) {
				cout << "FAILED" << endl;
				cerr << "Cell " << cell
					<< " is in outer cells"
					<< endl;
				abort();
			}
		}

		// check that outer cell iterator works
		for (auto item = grid.begin_outer(); item != grid.end_outer(); item++) {
			const uint64_t cell = item->first;

			if (inner_cells.count(cell) > 0) {
				cout << "FAILED" << endl;
				cerr << "Cell " << cell
					<< " is in inner cells"
					<< endl;
				abort();
			}

			if (outer_cells.count(cell) == 0) {
				cout << "FAILED" << endl;
				cerr << "Cell " << cell
					<< " is not in outer cells"
					<< endl;
				abort();
			}
		}
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

