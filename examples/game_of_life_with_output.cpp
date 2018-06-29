/*
A simple 2 D game of life program to demonstrate the efficient parallel usage
of dccrg and shows an example of how to output dccrg grid data into a file.
*/

#include "cstddef"
#include "cstdio"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "sstream"
#include "cstdint"

#include "mpi.h"
#include "zoltan.h"

#include "../dccrg.hpp"
#include "../dccrg_cartesian_geometry.hpp"

using namespace std;

// store in every cell of the grid whether the cell is alive and the number of live neighbors it has
struct game_of_life_cell {
	unsigned int
		is_alive = 0,
		live_neighbor_count = 0;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() {
		return std::make_tuple((void*) &(this->is_alive), 1, MPI_UNSIGNED);
	}
};


/*!
Initializes game.
*/
void initialize_game(
	dccrg::Dccrg<game_of_life_cell, dccrg::Cartesian_Geometry>& grid
) {
	for (const auto& cell: grid.local_cells) {
		cell.data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell.data->is_alive = 1;
		} else {
			cell.data->is_alive = 0;
		}
	}
}


/*!
Calculates the number of live neihgbours for every cell given, all of which must be local
*/
template<class Cell_Iterator> void get_live_neighbor_counts(
	const Cell_Iterator& cells
) {
	for (const auto& cell: cells) {
		cell.data->live_neighbor_count = 0;

		for (const auto& neighbor: cell.neighbors_of) {
			if (neighbor.data->is_alive > 0) {
				cell.data->live_neighbor_count++;
			}
		}
	}
}


/*!
Applies the game of life rules to every given cell, all of which must be local
*/
void apply_rules(
	const dccrg::Dccrg<game_of_life_cell, dccrg::Cartesian_Geometry>& grid
) {
	for (const auto& cell: grid.cells) {
		if (cell.data->live_neighbor_count == 3) {
			cell.data->is_alive = 1;
		} else if (cell.data->live_neighbor_count != 2) {
			cell.data->is_alive = 0;
		}
	}
}


/*!
Writes the game state into a file named game_of_life_, postfixed with the timestep and .dc

See the file dc2vtk.cpp for an example of how to read in the saved files.
*/
bool write_game_data(
	uint64_t step,
	dccrg::Dccrg<game_of_life_cell, dccrg::Cartesian_Geometry>& grid
) {
	// get the output filename
	ostringstream basename("game_of_life_"), step_string, suffix(".dc");
	step_string.width(3);
	step_string.fill('0');
	step_string << step;

	string output_name("");
	output_name += basename.str();
	output_name += step_string.str();
	output_name += suffix.str();

	// process 0 writes the file header which consists of given time step
	std::tuple<void*, int, MPI_Datatype> file_header;
	get<0>(file_header) = (void*) &step; // header data starts from this memory address
	get<1>(file_header) = 1; // number of datatypes to write
	get<2>(file_header) = MPI_UINT64_T; // datatype to use

	if (not grid.save_grid_data(output_name, 0, file_header)) {
		abort();
	}

	return true;
}


/*!
See the comments in simple_game_of_life.cpp and game_of_life.cpp
for an explanation of the basics.
*/
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

	dccrg::Dccrg<game_of_life_cell, dccrg::Cartesian_Geometry> grid;

	const int
		neighborhood_size = 1,
		maximum_refinement_level = 0;
	const std::array<uint64_t, 3> grid_length = {{10, 10, 1}};
	grid
		.initialize(
			grid_length,
			comm,
			"RCB",
			neighborhood_size,
			maximum_refinement_level)
		.set_geometry(
			dccrg::Cartesian_Geometry_Parameters({{0, 0, 0}}, {{1, 1, 1}}))
		.balance_load();

	initialize_game(grid);

	const int turns = 10;
	for (int turn = 0; turn < turns; turn++) {

		write_game_data(turn, grid);

		grid.start_remote_neighbor_copy_updates();
		get_live_neighbor_counts(grid.inner_cells);

		grid.wait_remote_neighbor_copy_updates();
		get_live_neighbor_counts(grid.outer_cells);

		apply_rules(grid);
	}
	write_game_data(turns, grid);

	MPI_Finalize();

	return EXIT_SUCCESS;
}

