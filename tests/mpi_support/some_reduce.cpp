/*
Tests the Some_Reduce MPI support class of dccrg.
*/

#include "cstdlib"
#include "iostream"
#include "mpi.h"
#include "stdint.h"

#include "../../dccrg_mpi_support.hpp"

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
	boost::unordered_set<int> neighbors;
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
	} else if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return 0;
}

