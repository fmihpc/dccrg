/*
Tests the grid with a game of life on a refined grid in 3 D with neighbors only in the ? plane

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute

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

#include "algorithm"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "sstream"
#include "unordered_set"

#include "mpi.h"
#include "zoltan.h"

#include "dccrg_stretched_cartesian_geometry.hpp"
#include "dccrg.hpp"


struct game_of_life_cell {
	unsigned int is_alive, live_neighbor_count;

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->is_alive), 1, MPI_UNSIGNED);
	}
};


using namespace std;
using namespace dccrg;

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
	if (rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	Dccrg<game_of_life_cell, Stretched_Cartesian_Geometry> game_grid;

	const std::array<uint64_t, 3> grid_length = {{5, 5, 3}};
	const double cell_length = 1.0 / grid_length[0];

	#define NEIGHBORHOOD_SIZE 1
	game_grid.initialize(grid_length, comm, "RANDOM", NEIGHBORHOOD_SIZE);

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	game_grid.set_geometry(geom_params);

	vector<uint64_t> cells = game_grid.get_cells();
	if (rank == 0) {
		cout << "Maximum refinement level of the grid: "
			<< game_grid.get_maximum_refinement_level()
			<< endl;
	}

	// refine the grid increasingly in the z direction a few times
	game_grid.balance_load();
	cells = game_grid.get_cells();
	for (const auto& cell: cells) {
		if (game_grid.geometry.get_center(cell)[2] > 1 * 1.5 * cell_length) {
			game_grid.refine_completely(cell);
		}
	}
	game_grid.stop_refining();
	game_grid.balance_load();
	cells = game_grid.get_cells();
	cout << "Process " << rank << ": number of cells after refining: " << cells.size() << endl;
	for (const auto& cell: cells) {
		if (game_grid.geometry.get_center(cell)[2] > 2 * 1.5 * cell_length) {
			game_grid.refine_completely(cell);
		}
	}
	game_grid.stop_refining();
	game_grid.balance_load();
	cells = game_grid.get_cells();
	cout << "Process " << rank << ": number of cells after refining: " << cells.size() << endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (auto& item: game_grid.cells) {
		const auto& cell = get<0>(item);
		game_of_life_cell* cell_data = get<1>(item);
		cell_data->live_neighbor_count = 0;

		const std::array<double, 3>
			cell_center = game_grid.geometry.get_center(cell),
			cell_length = game_grid.geometry.get_length(cell);

		if (fabs(0.5 + 0.1 * cell_length[1] - cell_center[1]) < 0.5 * cell_length[1]) {
			cell_data->is_alive = 1;
		} else {
			cell_data->is_alive = 0;
		}
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "tests/game_of_life/refined_" << rank << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (rank == 0) {
		visit_file.open("tests/game_of_life/refined.visit");
		visit_file << "!NBLOCKS " << comm_size << endl;
	}

	#define TIME_STEPS 25
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.update_copies_of_remote_neighbors();
		auto cells = game_grid.cells;
		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
		sort(cells.begin(), cells.end());

		if (rank == 0) {
			cout << "step: " << step << endl;
		}

		// write the game state into a file named according to the current time step
		string current_output_name("");
		ostringstream step_string;
		step_string.fill('0');
		step_string.width(5);
		step_string << step;
		current_output_name += basename.str();
		current_output_name += step_string.str();
		current_output_name += suffix.str();

		// visualize the game with visit -o game_of_life_test.visit
		if (rank == 0) {
			for (int process = 0; process < comm_size; process++) {
				visit_file << "refined_" << process << "_" << step_string.str() << suffix.str() << endl;
			}
		}


		// write the grid into a file
		game_grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// go through the grids cells and write their state into the file
		outfile << "SCALARS is_alive float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& item: cells) {
			const game_of_life_cell* const cell_data = get<1>(item);

			if (cell_data->is_alive > 0) {
				outfile << "1";
			} else {
				outfile << "0";
			}
			outfile << endl;

		}

		// write each cells live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& item: cells) {
			const game_of_life_cell* const cell_data = get<1>(item);
			outfile << cell_data->live_neighbor_count << endl;
		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& item: cells) {
			const auto& cell = get<0>(item);
			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(cell);
			outfile << neighbors->size() << endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (size_t i = 0; i < cells.size(); i++) {
			outfile << rank << endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (auto& item: cells) {
			const auto& cell = get<0>(item);
			outfile << cell << endl;
		}
		outfile.close();

		// get the neighbor counts of every cell
		// FIXME: use the (at some point common) solver from (un)refined2d and only include x and y directions in neighborhood
		for (auto& item: game_grid.cells) {
			const auto& cell = get<0>(item);
			auto* const cell_data = get<1>(item);

			cell_data->live_neighbor_count = 0;
			const auto* const neighbors = game_grid.get_neighbors_of(cell);
			if (neighbors == NULL) {
				cout << "Process " << rank
					<< ": neighbor list for cell " << cell << " not available"
					<< endl;
				exit(EXIT_FAILURE);
			}

			for (const auto& neighbor: *neighbors) {

				if (neighbor == dccrg::error_cell) {
					continue;
				}

				// only consider neighbors in the same z plane
				if (
					game_grid.geometry.get_center(cell)[2]
					!= game_grid.geometry.get_center(neighbor)[2]
				) {
					continue;
				}

				const auto* const neighbor_data = game_grid[neighbor];
				if (neighbor_data == NULL) {
					cout << "Process " << rank
						<< ": neighbor " << neighbor
						<< " data of cell " << cell << " not available"
						<< endl;
					exit(EXIT_FAILURE);
				}
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// calculate the next turn
		for (auto& item: game_grid.cells) {
			auto* const cell_data = get<1>(item);

			if (cell_data->live_neighbor_count == 3) {
				cell_data->is_alive = 1;
			} else if (cell_data->live_neighbor_count != 2) {
				cell_data->is_alive = 0;
			}
		}

	}

	if (rank == 0) {
		visit_file.close();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
