/*
Tests the Some_Reduce MPI support class of dccrg.

Copyright 2015 Finnish Meteorological Institute
Copyright 2015 Ilja Honkonen

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

#include "cstdint"
#include "cstdlib"
#include "iostream"

#include "mpi.h"

#include "dccrg_mpi_support.hpp"

using namespace std;
using namespace dccrg;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank, size;
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Comm_rank(world, &rank);
	MPI_Comm_size(world, &size);

	// previous and next processes are neighbors
	std::unordered_set<int> neighbors;
	uint64_t correct_value = 10 * (uint64_t) abs(rank);
	if (rank > 0) {
		neighbors.insert(rank - 1);
		correct_value += 10 * (uint64_t) abs(rank - 1);
	}
	if (rank < size - 1) {
		neighbors.insert(rank + 1);
		correct_value += 10 * (uint64_t) abs(rank + 1);
	}

	const uint64_t
		local_value = (uint64_t) abs(rank) * 10,
		reduced_value = dccrg::Some_Reduce()(local_value, neighbors, world);

	if (reduced_value != correct_value) {
		cout << "FAILED (Process " << rank << "): "
			<< " reduced value should be " << correct_value
			<< " but is " << reduced_value
			<< endl;
		abort();
	}

	MPI_Finalize();

	return 0;
}

