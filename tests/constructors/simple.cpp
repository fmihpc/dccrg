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
		cout << "Default constructed" << endl;
	}

	~CellData()
	{
		cout << "Default destructed" << endl;
	}

	CellData(CellData& /*given*/)
	{
		cout << "Copied from non-const" << endl;
	}

	CellData(const CellData& /*given*/)
	{
		cout << "Copied from const" << endl;
	}

	CellData& operator = (const CellData& /*given*/)
	{
		cout << "Assigned" << endl;
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

	cout << "\ngrid.set_geometry:" << endl;
	if (!grid.set_geometry(1, 1, 1, 0, 0, 0, 1, 1, 1)) {
		cerr << "Couldn't set grid geometry" << endl;
		return EXIT_FAILURE;
	}

	cout << "\ngrid.initialize:" << endl;
	grid.initialize(comm, "RCB", 1, 0);

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

