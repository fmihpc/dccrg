// Validation for MPI functionality and gtest
// hello.fails_on_all and hello.fails_on_last should fail
// And the other tests should pass

#include <gtest/gtest.h>
#include <mpi.h>
#include "../dccrg.hpp"
#include "../dccrg_cartesian_geometry.hpp"

int main(int argc, char* argv[])
{
	int mpiError{MPI_SUCCESS};
	mpiError = MPI_Init(&argc, &argv);
	::testing::InitGoogleTest(&argc, argv);
	int ret {RUN_ALL_TESTS()};
	mpiError = MPI_Finalize();
	return ret;
}

class SimpleGrid : public testing::Test {
	protected:
	SimpleGrid()
	{
		// Defaults should be fine for everything
		grid.set_initial_length({4, 4, 4}).initialize(MPI_COMM_WORLD).set_load_balancing_method("RANDOM").balance_load();

		// Simple data we can check
		for (auto cell : grid.get_cells()) {
			*grid[cell] = cell;
		}

		grid.update_copies_of_remote_neighbors();
	}

	dccrg::Dccrg<uint64_t, dccrg::Cartesian_Geometry> grid;
};

TEST_F(SimpleGrid, contents)
{
	for (auto cell : grid.get_cells()) {
		EXPECT_EQ_MPI(*grid[cell], cell);
	}
}

TEST_F(SimpleGrid, contents_after_loadbalance)
{
	grid.balance_load();
	for (auto cell : grid.get_cells()) {
		EXPECT_EQ_MPI(*grid[cell], cell);
	}
}

// TODO: currently only checks consistency between local cells
TEST_F(SimpleGrid, consistent_neighbors)
{
	auto cells = grid.get_cells();
	for (auto cell : grid.get_cells()) {
		auto* my_neighbors_of {grid.get_neighbors_of(cell)};
		EXPECT_NE_MPI(my_neighbors_of, nullptr);
		auto* my_neighbors_to {grid.get_neighbors_to(cell)};
		EXPECT_NE_MPI(my_neighbors_to, nullptr);

		for (auto [neighbor, dir] : *my_neighbors_of) {
			if (neighbor != dccrg::error_cell) {
				EXPECT_NE_MPI(grid[neighbor], nullptr);
				if (std::find(cells.begin(), cells.end(), neighbor) != cells.end()) {
					auto* other_neighbors_to {grid.get_neighbors_to(neighbor)};
					EXPECT_NE_MPI(other_neighbors_to, nullptr);
					EXPECT_NE_MPI(std::find_if(other_neighbors_to->begin(), other_neighbors_to->end(), [&cell](const std::pair<const uint64_t, std::array<int, 4>> pair){return pair.first == cell;}), other_neighbors_to->end());
				}
			}
		}

		for (auto [neighbor, dir] : *my_neighbors_to) {
			if (std::find(cells.begin(), cells.end(), neighbor) != cells.end()) {
				auto* other_neighbors_of {grid.get_neighbors_of(neighbor)};
				EXPECT_NE_MPI(other_neighbors_of, nullptr);
				EXPECT_NE_MPI(std::find_if(other_neighbors_of->begin(), other_neighbors_of->end(), [&cell](const std::pair<const uint64_t, std::array<int, 4>> pair){return pair.first == cell;}), other_neighbors_of->end());
			}
		}
	}
}
