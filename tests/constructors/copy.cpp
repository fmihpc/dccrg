/*
Tests copy construction of dccrg.
*/

#include "boost/foreach.hpp"
#include "boost/mpi.hpp"
#include "cstdlib"
#include "ctime"
#include "iostream"
#include "mpi.h"
#include "vector"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;

/*!
Cell data in grid1.
*/
struct Cell1 {
	int data;

	Cell1() { data = -1; }

	void mpi_datatype(
		void*& address,
		int& count,
		MPI_Datatype& datatype,
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/
	) {
		address = &(this->data);
		count = 1;
		datatype = MPI_INT;
	}
};

/*!
Cell data in grid2.
*/
struct Cell2 {
	double data;

	Cell2() { data = -2; }

	void mpi_datatype(
		void*& address,
		int& count,
		MPI_Datatype& datatype,
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/
	) {
		address = &(this->data);
		count = 1;
		datatype = MPI_DOUBLE;
	}
};

int main(int argc, char* argv[])
{
	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cout << "Zoltan_Initialize failed" << endl;
		return EXIT_FAILURE;
	}

	dccrg::Dccrg<Cell1> grid1;
	if (!grid1.set_geometry(10, 10, 10, 0, 0, 0, 1, 1, 1)) {
		cerr << "Couldn't set grid geometry" << endl;
		return EXIT_FAILURE;
	}
	grid1.initialize(comm, "RCB", 1, 0);

	// check that remote neighbor update works in original grid
	// data in grid1 == process rank
	const vector<uint64_t> cells1 = grid1.get_cells();
	BOOST_FOREACH(const uint64_t cell, cells1) {
		Cell1* cell_data = grid1[cell];
		cell_data->data = comm.rank();
	}

	grid1.update_remote_neighbor_data();

	const boost::unordered_set<uint64_t>& remote_neighbors1
		= grid1.get_cells_with_remote_neighbors();

	const boost::unordered_map<uint64_t, uint64_t>& cell_process1
		= grid1.get_cell_process();

	BOOST_FOREACH(const uint64_t cell, remote_neighbors1) {
		Cell1* cell_data = grid1[cell];
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
	dccrg::Dccrg<Cell2> grid2(grid1);
	const vector<uint64_t> cells2 = grid2.get_cells();

	if (cells1.size() != cells2.size()) {
		cerr << "Rank " << comm.rank() << ": Number of cells doesn't match" << endl;
		abort();
	}

	BOOST_FOREACH(const uint64_t cell, cells2) {
		Cell2* cell_data = grid2[cell];
		if (cell_data->data != -2) {
			cerr << "Rank " << comm.rank()
				<< ": Wrong data in cell " << cell
				<< ": " << cell_data->data
				<< ", should be -2"
				<< endl;
			abort();
		}
	}

	// check that remote neighbor update works in original grid
	// data in grid2 == 2 * process rank
	BOOST_FOREACH(const uint64_t cell, cells2) {
		Cell2* cell_data = grid2[cell];
		cell_data->data = 2 * comm.rank();
	}

	grid2.update_remote_neighbor_data();

	const boost::unordered_set<uint64_t>& remote_neighbors2
		= grid2.get_cells_with_remote_neighbors();

	const boost::unordered_map<uint64_t, uint64_t>& cell_process2
		= grid2.get_cell_process();

	BOOST_FOREACH(const uint64_t cell, remote_neighbors2) {
		Cell2* cell_data = grid2[cell];
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

	if (comm.rank() == 0) {
		cout << "PASSED" << endl;
	}

	return EXIT_SUCCESS;
}

