/*
The simplest 2 D game of life program that demonstrates the basic usage of dccrg
The program doesn't scale, takes no input and produces no output.
For good scaling see the file game_of_life.cpp, for output see the files game_of_life_with_output.cpp and dc2vtk.cpp
*/

#include "cstdlib"

#include "mpi.h"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;

// store in every cell of the grid whether the cell is alive and the number of live neighbors it has
struct game_of_life_cell {
	unsigned int
		is_alive = 0,
		live_neighbor_count = 0;

	/*
	Whenever cell data is transferred between MPI processes dccrg
	passes on this data to MPI for sending/receiving the data.
	*/
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() {
		return std::make_tuple((void*) &(this->is_alive), 1, MPI_UNSIGNED);
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
	    exit(EXIT_FAILURE);
	}

	dccrg::Dccrg<game_of_life_cell> game_grid;

	const int
		// the cells that share a vertex are considered neighbors
		neighborhood_size = 1,
		// don't allow refining cells into 8 smaller ones
		maximum_refinement_level = 0;

	// length of grid in cells of refinement level 0 (largest possible)
	const std::array<uint64_t, 3> grid_length{{10, 10, 1}};

	game_grid.initialize(
		grid_length,
		comm,
		// use the recursive coordinate bisection method for load
		// balancing (http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg_rcb.html)
		"RCB",
		neighborhood_size,
		maximum_refinement_level
	);

	// since the grid doesn't change (isn't refined / unrefined)
	// during the game, workload can be balanced just once in the beginning
	game_grid.balance_load();

	// initialize the game
	for (const auto& cell: game_grid.cells) {
		cell.data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell.data->is_alive = 1;
		} else {
			cell.data->is_alive = 0;
		}
	}

	const int turns = 100;
	for (int turn = 0; turn < turns; turn++) {

		/*
		Transfer data of cells, which are neighbors but reside
		on different processes, between processes.
		*/
		game_grid.update_copies_of_remote_neighbors();

		// get the neighbor counts for every local cell
		for (const auto& cell: game_grid.cells) {
			cell.data->live_neighbor_count = 0;

			const auto* const neighbors = game_grid.get_neighbors_of(cell.id);

			for (const auto& neighbor: *neighbors) {
				/*
				Skip neighbors that would be outside of
				the grid and are recorded as error_cell
				*/
				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = game_grid[neighbor];
				if (neighbor_data == nullptr) {
					abort();
				}

				if (neighbor_data->is_alive > 0) {
					cell.data->live_neighbor_count++;
				}
			}
		}

		// calculate the next turn
		for (const auto& cell: game_grid.cells) {
			if (cell.data->live_neighbor_count == 3) {
				cell.data->is_alive = 1;
			} else if (cell.data->live_neighbor_count != 2) {
				cell.data->is_alive = 0;
			}
		}

	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

