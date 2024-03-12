#include "common.hpp"
#include <gtest/gtest.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
	int mpiError = MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int ret {RUN_ALL_TESTS()};
	mpiError = MPI_Finalize();
	return ret;
}

TEST(hello, passes_on_all)
{
	int myRank {0};
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	std::cout << "Hello from " << myRank << std::endl;
	ASSERT_EQ_MPI(0, 0);
}

TEST(hello, fails_on_all)
{
	int myRank {0};
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	std::cout << "Hello from " << myRank << std::endl;
	ASSERT_EQ_MPI(myRank, -1);
}

TEST(hello, fails_on_one)
{
	int myRank {0};
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	std::cout << "Hello from " << myRank << std::endl;
	ASSERT_NE_MPI(myRank, 0);
}
