/*
Tests how well dccrg avoids unnecessary construction / destruction of cell data.

Copyright 2015, 2016 Finnish Meteorological Institute
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

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) &(this->data), 1, MPI_DOUBLE);
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
		cout << "Zoltan_Initialize failed" << endl;
		return EXIT_FAILURE;
	}

	// initialize grid
	cout << "\nDccrg<CellData> grid:" << endl;
	dccrg::Dccrg<CellData> grid;

	cout << "\ngrid.initialize:" << endl;
	const std::array<uint64_t, 3> grid_length = {{1, 1, 1}};
	grid.initialize(grid_length, comm, "RCB", 1, 0);

	cout << "\ngrid.get_cells:" << endl;
	const auto cells = grid.get_cells();

	cout << "\nfor (const auto& cell: cells):" << endl;
	for (const auto& cell: cells) {
		cout << cell << " ";
	}
	cout << endl;

	cout << "\nexiting:" << endl;

	MPI_Finalize();

	return EXIT_SUCCESS;
}

