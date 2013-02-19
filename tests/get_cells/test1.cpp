/*
Program for testing the get_cells() function of dccrg.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute

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

#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "boost/mpi.hpp"
#include "cstdlib"
#include "iostream"
#include "zoltan.h"

#define DCCRG_TRANSFER_USING_BOOST_MPI
#include "../../dccrg.hpp"

using namespace std;
using namespace boost;
using namespace boost::assign;
using namespace boost::mpi;
using namespace dccrg;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	if (comm.size() < 2) {
		if (comm.rank() == 0) {
			cerr << "Must be run with at least 2 processes" << endl;
			cout << "FAILED" << endl;
		}
		return EXIT_FAILURE;
	}

	Dccrg<int> grid;

	if (!grid.set_geometry(7, 1, 1, 0, 0, 0, 1, 1, 1)) {
		cerr << "Couldn't set grid geometry" << endl;
		exit(EXIT_FAILURE);
	}

	grid.initialize(comm, "RANDOM", 2);

	// distribute cells to processes like this,
	// ______________  processes >= 2 don't get cells
	// |0|0|1|0|1|0|0|
	// --------------
	const std::vector<uint64_t> initial_cells = grid.get_cells();
	BOOST_FOREACH(const uint64_t cell, initial_cells) {
		if (cell == 3 || cell == 5) {
			grid.pin(cell, 1);
		} else {
			grid.pin(cell, 0);
		}
	}
	grid.balance_load(false);
	grid.unpin_all_cells();


	// use different neighborhoods for querying cells with different neighbor types
	typedef dccrg::Types<3>::neighborhood_item_t neigh_t;

	const std::vector<neigh_t>
		// neighborhood without neighbors
		neighborhood0,
		// one neighbor in positive x direction
		neighborhood1 = list_of<neigh_t>(list_of(2)(0)(0)),
		// two neighbors in +y and +z
		neighborhood2 = list_of<neigh_t>(list_of(0)(1)(0))(list_of(0)(0)(1));

	const int
		neighborhood_id0 = 0,
		neighborhood_id1 = 1,
		neighborhood_id2 = 2;

	if (!grid.add_remote_update_neighborhood(neighborhood_id0, neighborhood0)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set neighborhood " << neighborhood_id0
			<< std::endl;
		abort();
	}
	if (!grid.add_remote_update_neighborhood(neighborhood_id1, neighborhood1)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set neighborhood " << neighborhood_id1
			<< std::endl;
		abort();
	}
	if (!grid.add_remote_update_neighborhood(neighborhood_id2, neighborhood2)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set neighborhood " << neighborhood_id2
			<< std::endl;
		abort();
	}

	// cells returned for different neighbor types
	const std::vector<uint64_t>
		no_neigh_cells = grid.get_cells(list_of(0), true, neighborhood_id0),
		one_neigh_cells = grid.get_cells(std::vector<int>(), true, neighborhood_id1),
		two_neigh_cells = grid.get_cells(list_of(0), true, neighborhood_id2),
		local_neigh_of_cells = grid.get_cells(list_of(has_local_neighbor_of), true, neighborhood_id1),
		local_neigh_to_cells = grid.get_cells(list_of(has_local_neighbor_to), true, neighborhood_id1),
		remote_neigh_of_cells = grid.get_cells(list_of(has_remote_neighbor_of), true, neighborhood_id1),
		remote_neigh_to_cells = grid.get_cells(list_of(has_remote_neighbor_to), true, neighborhood_id1),
		local_of_remote_to_cells
			= grid.get_cells(
				list_of(has_local_neighbor_of | has_remote_neighbor_to),
				true,
				neighborhood_id1
			),
		local_to_remote_of_cells
			= grid.get_cells(
				list_of(has_local_neighbor_to | has_remote_neighbor_of),
				true,
				neighborhood_id1
			),
		local_neigh_both_cells = grid.get_cells(list_of(has_local_neighbor_both), true, neighborhood_id1),
		remote_neigh_both_cells = grid.get_cells(list_of(has_remote_neighbor_both), true, neighborhood_id1);

	// number of unfiltered cells
	const size_t nr_of_cells = one_neigh_cells.size();
	if (comm.rank() == 0) {
		if (nr_of_cells != 5) {
			cerr << "Proc 0: Incorrect number of unfiltered cells: " << nr_of_cells
				<< ", should be 5" << endl;
			cout << "FAILED" << endl;
			abort();
		}
	}
	if (comm.rank() == 1) {
		if (nr_of_cells != 2) {
			cerr << "Proc 1: Incorrect number of unfiltered cells: " << nr_of_cells
				<< ", should be 2" << endl;
			cout << "FAILED" << endl;
			abort();
		}
	}
	if (comm.rank() >= 2) {
		if (nr_of_cells != 0) {
			cerr << "Proc " << comm.rank() << ": Incorrect number of unfiltered cells: " << nr_of_cells
				<< ", should be 0" << endl;
			cout << "FAILED" << endl;
			abort();
		}
	}

	// no neighbors case
	if (no_neigh_cells.size() != nr_of_cells) {
		cerr << "Proc " << comm.rank()
			<< ": Incorrect number of cells without neighbors for neighborhood id 0: "
			<< no_neigh_cells.size() << ", should be " << nr_of_cells
			<< endl;
		cout << "FAILED" << endl;
		abort();
	}
	if (two_neigh_cells.size() != nr_of_cells) {
		cerr << "Proc " << comm.rank()
			<< ": Incorrect number of cells without neighbors for neighborhood id 2: "
			<< two_neigh_cells.size() << ", should be " << nr_of_cells
			<< endl;
		cout << "FAILED" << endl;
		abort();
	}

	// rest of tests for processes 0 and 1 only
	if (comm.rank() >= 2) {
		return EXIT_SUCCESS;
	}

	if (comm.rank() == 0) {
		// reminder from above, cell ids start from 1, left to right
		// 0010100
		if (local_neigh_of_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (local_neigh_of_cells[0] != 2) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_neigh_to_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (local_neigh_to_cells[0] != 6) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (remote_neigh_of_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (remote_neigh_of_cells[0] != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (remote_neigh_to_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (remote_neigh_to_cells[0] != 7) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_of_remote_to_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_to_remote_of_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_neigh_both_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (local_neigh_both_cells[0] != 4) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (remote_neigh_both_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}

	} else {

		// reminder from above, cell ids start from 1, left to right
		// 0010100
		if (local_neigh_of_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_neigh_to_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (remote_neigh_of_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (remote_neigh_to_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_of_remote_to_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (local_of_remote_to_cells[0] != 3) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_to_remote_of_cells.size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		} else if (local_to_remote_of_cells[0] != 5) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (local_neigh_both_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
		if (remote_neigh_both_cells.size() != 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << endl;
			abort();
		}
	}

	if (comm.rank() == 0) {
		cout << "PASSED" << endl;
	}

	return EXIT_SUCCESS;
}

