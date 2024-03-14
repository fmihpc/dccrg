#include <gtest/gtest.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
	int mpiError{MPI_SUCCESS};
	mpiError = MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int ret {RUN_ALL_TESTS()};
	mpiError = MPI_Finalize();
	return ret;
}

TEST(hello, passes_on_all)
{
	int myRank {0};
	int mpiError {MPI_SUCCESS};

	mpiError = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	ASSERT_EQ_MPI(mpiError, MPI_SUCCESS);
	ASSERT_EQ_MPI(0, 0);
}

TEST(hello, fails_on_all)
{
	int myRank {0};
	int mpiError {MPI_SUCCESS};

	mpiError = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	ASSERT_EQ_MPI(mpiError, MPI_SUCCESS);
	ASSERT_EQ_MPI(myRank, -1);
}

TEST(hello, fails_on_one)
{
	int myRank {0};
	int mpiError {MPI_SUCCESS};

	mpiError = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	ASSERT_EQ_MPI(mpiError, MPI_SUCCESS);
	ASSERT_NE_MPI(myRank, 0);
}

TEST(MPI, sendrecv)
{
	int myRank{0};
	int size{0};
	int mpiError{MPI_SUCCESS};

	mpiError = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	ASSERT_EQ_MPI(mpiError, MPI_SUCCESS);
	mpiError = MPI_Comm_size(MPI_COMM_WORLD, &size);
	ASSERT_EQ_MPI(mpiError, MPI_SUCCESS);

	int dest {myRank < size - 1 ? myRank + 1 : 0};
	int source {myRank > 0 ? myRank - 1 : size - 1};
	int send {myRank};
	int recv {-1};

	mpiError = MPI_Sendrecv(&send, 1, MPI_INT, dest, 0, &recv, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	ASSERT_EQ_MPI(mpiError, MPI_SUCCESS);
	ASSERT_EQ_MPI(source, recv);
}
