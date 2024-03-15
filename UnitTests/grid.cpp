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

enum class GridType {
	simple,
	vlasovian,
	magnetospheric
};

class GridTest : public testing::TestWithParam<GridType> {
	protected:
	void SetUp() override
	{
		typedef dccrg::Types<3>::neighborhood_item_t neigh_t;
		std::vector<neigh_t> neighborhood;

		grid.set_initial_length({size, size, size}).set_neighborhood_length(stencil_width).set_maximum_refinement_level(reflevel).initialize(MPI_COMM_WORLD);
		switch (GetParam()) {
			case GridType::simple:
				break; // Defaults should be fine for everything
			case GridType::magnetospheric:
				for (int i = 0; i < reflevel; ++i) {
					grid.refine_completely_at({static_cast<double>(size) / 2.0, static_cast<double>(size) / 2.0, static_cast<double>(size) / 2.0});
					grid.stop_refining();
				}
				// Fallthrough
			case GridType::vlasovian:
				// Extended sysboundaries
				for (int x = -stencil_width; x <= stencil_width; ++x) {
					for (int y = -stencil_width; y <= stencil_width; ++y) {
						for (int z = -stencil_width; z <= stencil_width; ++z) {
							if (x || y || z) {
								neighborhood.push_back({x, y, z});
							}
						}
					}
				}
				grid.add_neighborhood(neighborhoods++, neighborhood);

				// Vlasov solver
				for (int d = -stencil_width; d <= stencil_width; ++d) {
					if (d) {
						neighborhood.push_back({d, 0, 0});
						neighborhood.push_back({0, d, 0});
						neighborhood.push_back({0, 0, d});
					}
				}
				grid.add_neighborhood(neighborhoods++, neighborhood);

				break;
			default:
				FAIL() << "Grid type not implemented!";
				break;
		}

		grid.set_load_balancing_method("RANDOM").balance_load();

		// Simple data we can check
		for (auto cell : grid.get_cells()) {
			*grid[cell] = cell;
		}
		grid.update_copies_of_remote_neighbors();
		for (int neighborhood = 0; neighborhood < neighborhoods; ++neighborhood) {
			grid.update_copies_of_remote_neighbors(neighborhood);
		}
	}

	dccrg::Dccrg<uint64_t, dccrg::Cartesian_Geometry> grid;
	int neighborhoods {0};
	const int stencil_width {3};
	const uint64_t size {10};
	int reflevel {2};
};

INSTANTIATE_TEST_SUITE_P(Simple, GridTest, testing::Values(GridType::simple));
INSTANTIATE_TEST_SUITE_P(Vlasovian, GridTest, testing::Values(GridType::vlasovian));
INSTANTIATE_TEST_SUITE_P(Magnetospheric, GridTest, testing::Values(GridType::magnetospheric));

TEST_P(GridTest, contents)
{
	for (auto cell : grid.get_cells()) {
		EXPECT_NE_MPI(grid[cell], nullptr);
		EXPECT_EQ_MPI(*grid[cell], cell);
	}
}

TEST_P(GridTest, contents_after_loadbalance)
{
	grid.balance_load();
	for (auto cell : grid.get_cells()) {
		EXPECT_NE_MPI(grid[cell], nullptr);
		EXPECT_EQ_MPI(*grid[cell], cell);
	}
}

TEST_P(GridTest, consistent_neighbors)
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
				std::vector<std::pair<uint64_t, std::array<int, 4>>> other_neighbors_to;
				if (std::find(cells.begin(), cells.end(), neighbor) != cells.end()) {
					auto* p {grid.get_neighbors_to(neighbor)};
					EXPECT_NE_MPI(p, nullptr);
					other_neighbors_to = *p;
				} else {
					// Warning: giga jank
					std::vector<uint64_t> found_neighbors;
					for (auto& [neigh, dir] : grid.find_neighbors_of(neighbor, grid.get_neighborhood_of(), grid.get_max_ref_lvl_diff())) {
						found_neighbors.push_back(neigh);
					}
					other_neighbors_to = grid.find_neighbors_to(neighbor, found_neighbors);
				}
				EXPECT_NE_MPI(std::find_if(other_neighbors_to.begin(), other_neighbors_to.end(), [&cell](const std::pair<const uint64_t, std::array<int, 4>> pair){return pair.first == cell;}), other_neighbors_to.end());
			}
		}

		for (auto [neighbor, dir] : *my_neighbors_to) {
			std::vector<std::pair<uint64_t, std::array<int, 4>>> other_neighbors_of;
			if (std::find(cells.begin(), cells.end(), neighbor) != cells.end()) {
				auto* p {grid.get_neighbors_of(neighbor)};
				EXPECT_NE_MPI(p, nullptr);
				other_neighbors_of = *p;
			} else {
				other_neighbors_of = grid.find_neighbors_of(neighbor, grid.get_neighborhood_of(), grid.get_max_ref_lvl_diff());
			}
			EXPECT_NE_MPI(std::find_if(other_neighbors_of.begin(), other_neighbors_of.end(), [&cell](const std::pair<const uint64_t, std::array<int, 4>> pair){return pair.first == cell;}), other_neighbors_of.end());
		}
	}
}

// TODO: cannot test consistency of remote neighbors without dccrg changes
TEST_P(GridTest, consistent_user_neighbors)
{
	auto cells = grid.get_cells();
	for (auto cell : grid.get_cells()) {
		for (int neighborhood = 0; neighborhood < neighborhoods; ++neighborhood) {
			auto* my_neighbors_of {grid.get_neighbors_of(cell, neighborhood)};
			EXPECT_NE_MPI(my_neighbors_of, nullptr);

			for (auto [neighbor, dir] : *my_neighbors_of) {
				if (neighbor != dccrg::error_cell) {
					EXPECT_NE_MPI(grid[neighbor], nullptr);
					std::vector<std::pair<uint64_t, std::array<int, 4>>> other_neighbors_to;
					if (std::find(cells.begin(), cells.end(), neighbor) != cells.end()) {
						auto* p {grid.get_neighbors_to(neighbor, neighborhood)};
						EXPECT_NE_MPI(p, nullptr);
						other_neighbors_to = *p;
						EXPECT_NE_MPI(std::find_if(other_neighbors_to.begin(), other_neighbors_to.end(), [&cell](const std::pair<const uint64_t, std::array<int, 4>> pair){return pair.first == cell;}), other_neighbors_to.end());
					}
				}
			}

			auto* my_neighbors_to {grid.get_neighbors_to(cell, neighborhood)};
			EXPECT_NE_MPI(my_neighbors_to, nullptr);

			for (auto [neighbor, dir] : *my_neighbors_to) {
				std::vector<std::pair<uint64_t, std::array<int, 4>>> other_neighbors_of;
				if (std::find(cells.begin(), cells.end(), neighbor) != cells.end()) {
					auto* p {grid.get_neighbors_of(neighbor, neighborhood)};
					EXPECT_NE_MPI(p, nullptr);
					other_neighbors_of = *p;
					EXPECT_NE_MPI(std::find_if(other_neighbors_of.begin(), other_neighbors_of.end(), [&cell](const std::pair<const uint64_t, std::array<int, 4>> pair){return pair.first == cell;}), other_neighbors_of.end());
				}
			}
		}
	}
}

// TODO test proper copies and frees of dccrg.comm
// Right now this can't be done because the getter is not a getter
TEST_P(GridTest, copy)
{
	auto other_grid = grid;

	// Local cells should be identical immediately after copy
	for (auto cell : grid.get_cells()) {
		EXPECT_NE_MPI(other_grid[cell], nullptr);
		EXPECT_EQ_MPI(*other_grid[cell], cell);
	}

	other_grid.balance_load();
	for (auto cell : other_grid.get_cells()) {
		EXPECT_NE_MPI(other_grid[cell], nullptr);
		EXPECT_EQ_MPI(*other_grid[cell], cell);
	}

	// Load balancing copy shouldn't affect original
	for (auto cell : grid.get_cells()) {
		EXPECT_NE_MPI(grid[cell], nullptr);
		EXPECT_EQ_MPI(*grid[cell], cell);
	}
}
