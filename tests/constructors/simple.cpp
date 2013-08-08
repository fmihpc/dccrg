/*
Tests how well dccrg avoids unnecessary construction / destruction of cell data.
*/

#include "boost/foreach.hpp"
#include "boost/mpi.hpp"
#include "boost/tuple/tuple.hpp"
#include "cstdlib"
#include "ctime"
#include "iostream"
#include "vector"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;

struct CellData {

	double data;

	CellData()
	{
		cout << "Cell default constructed" << endl;
	}

	~CellData()
	{
		cout << "Cell default destructed" << endl;
	}

	CellData(CellData& /*given*/)
	{
		cout << "Cell copied from non-const" << endl;
	}

	CellData(const CellData& /*given*/)
	{
		cout << "Cell copied from const" << endl;
	}

	CellData& operator = (const CellData& /*given*/)
	{
		cout << "Cell assigned" << endl;
		return *this;
	}

	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/
	) {
		return boost::make_tuple(&(this->data), 1, MPI_DOUBLE);
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

	// initialize grid
	cout << "\nDccrg<CellData> grid:" << endl;
	dccrg::Dccrg<CellData> grid;

	cout << "\ngrid.initialize:" << endl;
	const boost::array<uint64_t, 3> grid_length = {{1, 1, 1}};
	grid.initialize(grid_length, comm, "RCB", 1, 0);

	cout << "\ngrid.get_cells:" << endl;
	vector<uint64_t> cells = grid.get_cells();

	cout << "\nBOOST_FOREACH(const uint64_t& cell, cells):" << endl;
	BOOST_FOREACH(const uint64_t& cell, cells) {
		cout << cell << " ";
	}
	cout << endl;

	cout << "\nexiting:" << endl;
	return EXIT_SUCCESS;
}

