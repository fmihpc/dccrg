/*
Tests the All_Gather MPI support class of dccrg.
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

	/*
	Each process has its rank + 1 number of values
	with values within range [rank, 2 * rank + 1]
	*/
	std::vector<uint64_t> values(rank + 1, -1);
	for (size_t i = 0; i < values.size(); i++) {
		values[i] = rank + i;
	}

	std::vector<std::vector<uint64_t> > result;
	dccrg::All_Gather()(values, result, world);

	// chech the result
	for (size_t i = 0; i < result.size(); i++) {
		for (size_t j = 0; j < result[i].size(); j++) {
			if (result[i][j] != i + j) {
				cerr << "FAILED: invalid value received from process " << i
					<< " at index " << j
					<< ": " << result[i][j]
					<< ", should be " << i + j
					<< endl;
				abort();
			}
		}
	}

	cout.flush();
	MPI_Barrier(world);

	/*
	First process doesn't send any value
	*/
	if (rank == 0) {
		values.clear();
	}

	result.clear();
	dccrg::All_Gather()(values, result, world);

	if (result[0].size() > 0) {
		cout << "FAILED: rank 0 should have no values" << endl;
		abort();
	}

	// chech the result
	for (size_t i = 0; i < result.size(); i++) {
		for (size_t j = 0; j < result[i].size(); j++) {
			if (result[i][j] != i + j) {
				cerr << "FAILED: invalid value received from process " << i
					<< " at index " << j
					<< ": " << result[i][j]
					<< ", should be " << i + j
					<< endl;
				abort();
			}
		}
	}

	cout.flush();
	MPI_Barrier(world);

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return 0;
}

