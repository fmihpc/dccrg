/*
Program for testing user defined neighborhoods of dccrg.

Copyright 2012 Finnish Meteorological Institute

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

#include "algorithm"
#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "vector"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;

struct Cell
{
	int data;
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	template<typename Archiver> void serialize(
		Archiver& ar,
		const unsigned int
	) {
		ar & data;
	}
	#else
	void* at()
	{
		return this;
	}
	const void* at() const
	{
		return this;
	}
	#ifdef DCCRG_USER_MPI_DATA_TYPE
	MPI_Datatype mpi_datatype() const
	{
		MPI_Datatype type;
		MPI_Type_contiguous(sizeof(this->data), MPI_BYTE, &type);
		return type;
	}
	#else
	static size_t size(void)
	{
		return sizeof(int);
	}
	#endif
	#endif	// ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
};


int main(int argc, char* argv[])
{
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}

	dccrg::Dccrg<Cell> grid;

	if (!grid.set_geometry(10, 10, 10, 0, 0, 0, 1, 1, 1)) {
		cerr << "Couldn't set grid geometry" << endl;
		exit(EXIT_FAILURE);
	}

	grid.initialize(comm, "RANDOM", 2, 0, true, true, true);

	const vector<uint64_t> cells = grid.get_cells();

	// default neighbor lists should have 5^3 - 1 = 124 neighbors
	BOOST_FOREACH(const uint64_t cell, cells) {
		const vector<uint64_t>* neighbors = grid.get_neighbors(cell);
		if (neighbors->size() != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< ": " << neighbors->size()
				<< std::endl;
			abort();
		}
		const vector<uint64_t>* neighbors_to = grid.get_neighbors2(cell);
		if (neighbors_to->size() != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< ": " << neighbors_to->size()
				<< std::endl;
			abort();
		}
	}

	typedef dccrg::Types<3>::neighborhood_item_t neigh_t;

	// create a neighborhood of one cell at (-2, -2, -2)
	{
	const int hood_id = 1;
	const std::vector<neigh_t> neighborhood
		= boost::assign::list_of<neigh_t>(boost::assign::list_of(-2)(-2)(-2));
	if (!grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood failed"
				<< std::endl;
			abort();
	}
	BOOST_FOREACH(const uint64_t cell, cells) {
		const vector<uint64_t>* neighbors = grid.get_neighbors(cell, hood_id);
		if (neighbors->size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors->size()
				<< std::endl;
			abort();
		}
		const vector<uint64_t>* neighbors_to = grid.get_neighbors2(cell, hood_id);
		if (neighbors_to->size() != 1) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors_to->size()
				<< std::endl;
			abort();
		}
	}
	grid.remove_remote_update_neighborhood(hood_id);
	}

	// create a neighborhood of two cells at (-1, -1, -1) and (2, 2, 2)
	{
	const int hood_id = 2;
	const std::vector<neigh_t> neighborhood
		= boost::assign::list_of<neigh_t>
			(boost::assign::list_of(-1)(-1)(-1))
			(boost::assign::list_of( 2)( 2)( 2));
	if (!grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood failed"
				<< std::endl;
			abort();
	}
	BOOST_FOREACH(const uint64_t cell, cells) {
		const vector<uint64_t>* neighbors = grid.get_neighbors(cell, hood_id);
		if (neighbors->size() != 2) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors->size()
				<< std::endl;
			abort();
		}
		const vector<uint64_t>* neighbors_to = grid.get_neighbors2(cell, hood_id);
		if (neighbors_to->size() != 2) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors_to->size()
				<< std::endl;
			abort();
		}
	}
	grid.remove_remote_update_neighborhood(hood_id);
	}

	// create a neighborhood of 24 cells at z offset == 0
	{
	const int hood_id = -3;
	std::vector<neigh_t> neighborhood;
	for (int i = -2; i <= 2; i++)
	for (int j = -2; j <= 2; j++) {
		if (i == 0 && j == 0) {
			continue;
		}
		neigh_t temp = {{i, j, 0}};
		neighborhood.push_back(temp);
	}
	if (!grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood failed"
				<< std::endl;
			abort();
	}
	BOOST_FOREACH(const uint64_t cell, cells) {
		const vector<uint64_t>* neighbors = grid.get_neighbors(cell, hood_id);
		if (neighbors->size() != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors->size()
				<< std::endl;
			abort();
		}
		const vector<uint64_t>* neighbors_to = grid.get_neighbors2(cell, hood_id);
		if (neighbors_to->size() != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors_to->size()
				<< std::endl;
			abort();
		}
	}
	grid.remove_remote_update_neighborhood(hood_id);
	}

	// create a neighborhood of 24 cells at x offset == 0
	{
	const int hood_id = -4;
	std::vector<neigh_t> neighborhood;
	for (int j = -2; j <= 2; j++)
	for (int k = -2; k <= 2; k++) {
		if (j == 0 && k == 0) {
			continue;
		}
		neigh_t temp = {{0, j, k}};
		neighborhood.push_back(temp);
	}
	if (!grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood failed"
				<< std::endl;
			abort();
	}
	BOOST_FOREACH(const uint64_t cell, cells) {
		const vector<uint64_t>* neighbors = grid.get_neighbors(cell, hood_id);
		if (neighbors->size() != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors->size()
				<< std::endl;
			abort();
		}
		const vector<uint64_t>* neighbors_to = grid.get_neighbors2(cell, hood_id);
		if (neighbors_to->size() != 24) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors_to->size()
				<< std::endl;
			abort();
		}
	}
	grid.remove_remote_update_neighborhood(hood_id);
	}

	// create a full neighborhood identical to the default one
	{
	const int hood_id = 0;
	std::vector<neigh_t> neighborhood;
	for (int k = -2; k <= 2; k++)
	for (int j = -2; j <= 2; j++)
	for (int i = -2; i <= 2; i++) {
		if (i == 0 && j == 0 && k == 0) {
			continue;
		}
		neigh_t temp = {{i, j, k}};
		neighborhood.push_back(temp);
	}
	if (!grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood failed"
				<< std::endl;
			abort();
	}
	BOOST_FOREACH(const uint64_t cell, cells) {
		// check number of neighbors
		const vector<uint64_t>* neighbors = grid.get_neighbors(cell, hood_id);
		if (neighbors->size() != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors->size()
				<< std::endl;
			abort();
		}
		const vector<uint64_t>* neighbors_to = grid.get_neighbors2(cell, hood_id);
		if (neighbors_to->size() != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors_to for cell " << cell
				<< " and neighborhood id " << hood_id
				<< ": " << neighbors_to->size()
				<< std::endl;
			abort();
		}
		// check ids and ordering of neighbors_of
		const vector<uint64_t>* default_neighbors = grid.get_neighbors(cell);
		if (default_neighbors->size() != 124) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of neighbors for cell " << cell
				<< " in default neighborhood"
				<< std::endl;
			abort();
		}
		if (
			!std::equal(
				neighbors->begin(),
				neighbors->end(),
				default_neighbors->begin()
			)
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " User neighbor list of cell " << cell << ":\n";
			BOOST_FOREACH(const uint64_t neighbor, *neighbors) {
				std::cerr << neighbor << ", ";
			}
			std::cerr << "\nnot equal to default neighbor list:\n";
			BOOST_FOREACH(const uint64_t neighbor, *default_neighbors) {
				std::cerr << neighbor << ", ";
			}
			std::cerr << std::endl;
			abort();
		}
	}
	}

	// try to create a neighborhood too far away
	{
	const int hood_id = 5;
	const std::vector<neigh_t> neighborhood
		= boost::assign::list_of<neigh_t>(boost::assign::list_of(-1)(2)(-3));
	if (grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood succeeded"
				<< std::endl;
			abort();
	}
	}

	// try to use an existing neighborhood id
	{
	const int hood_id = 0;
	const std::vector<neigh_t> neighborhood
		= boost::assign::list_of<neigh_t>(boost::assign::list_of(1)(2)(-2));
	if (grid.add_remote_update_neighborhood(hood_id, neighborhood)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " add_remote_update_neighborhood succeeded"
				<< std::endl;
			abort();
	}
	}

	if (comm.rank() == 0) {
		cout << "Passed" << endl;
	}

	return EXIT_SUCCESS;
}

