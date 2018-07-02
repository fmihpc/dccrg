/*
Tests the grid with some simple game of life patters using a hierarchical load balancing method.
Returns EXIT_SUCCESS if everything went ok.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016, 2018 Finnish Meteorological Institute

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
	    cerr << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	Dccrg<game_of_life_cell> grid;

	const std::array<uint64_t, 3> grid_length = {{34, 7, 1}};
	grid
		.set_initial_length(grid_length)
		.set_neighborhood_length(1)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("HIER")
		.initialize(comm);

	grid.add_partitioning_level(12);
	grid.add_partitioning_option(0, "LB_METHOD", "HYPERGRAPH");
	grid.add_partitioning_option(0, "IMBALANCE_TOL", "1.05");

	grid.add_partitioning_level(1);
	grid.add_partitioning_option(1, "LB_METHOD", "HYPERGRAPH");

	grid.balance_load();
	// assumes no further load balancing, AMR, etc.
	auto cells = grid.cells;
	sort(cells.begin(), cells.end());

	// initialize the game
	for (const auto& cell: cells) {
		cell.data->live_neighbor_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell.data->is_alive = 1;
		} else {
			cell.data->is_alive = 0;
		}
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "tests/game_of_life/hierarchical_test_" << rank << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (rank == 0) {
		visit_file.open("tests/game_of_life/hierarchical_test.visit");
		visit_file << "!NBLOCKS " << comm_size << endl;
		//cout << "step: ";
	}

	#define TIME_STEPS 25
	for (int step = 0; step < TIME_STEPS; step++) {

		grid.start_remote_neighbor_copy_updates();
		grid.wait_remote_neighbor_copy_updates();

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
				visit_file << "hierarchical_test_" << process
					<< "_" << step_string.str() << suffix.str()
					<< endl;
			}
		}


		// write the grid into a file
		grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// go through the grids cells and write their state into the file
		outfile << "SCALARS is_alive float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			if (cell.data->is_alive == 1) {
				outfile << "1";
			} else {
				outfile << "0";
			}
			outfile << endl;

		}

		// write each cells live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			outfile << cell.data->live_neighbor_count << endl;
		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			const auto* const neighbors = grid.get_neighbors_of(cell.id);
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
		for (const auto& cell: cells) {
			outfile << cell.id << endl;
		}
		outfile.close();


		// get the neighbor counts of every cell, starting with the cells whose neighbor data doesn't come from other processes
		vector<uint64_t> inner_cells = grid.get_local_cells_not_on_process_boundary();
		for (const auto& cell: inner_cells) {

			auto* const cell_data = grid[cell];

			cell_data->live_neighbor_count = 0;
			const auto* const neighbors = grid.get_neighbors_of(cell);
			if (neighbors == nullptr) {
				cerr << "Process " << rank
					<< ": neighbor list for cell " << cell
					<< " not available"
					<< endl;
				exit(EXIT_FAILURE);
			}

			for (const auto& neighbor_i: *neighbors) {
				const auto& neighbor = neighbor_i.first;

				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = grid[neighbor];
				if (neighbor_data == nullptr) {
					cerr << "Process " << rank
						<< ": neighbor " << neighbor
						<< " data for cell " << cell
						<< " not available"
						<< endl;
					exit(EXIT_FAILURE);
				}
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		vector<uint64_t> outer_cells = grid.get_local_cells_on_process_boundary();
		for (const auto& cell: outer_cells) {

			auto* const cell_data = grid[cell];

			cell_data->live_neighbor_count = 0;
			const auto* const neighbors = grid.get_neighbors_of(cell);
			if (neighbors == nullptr) {
				cerr << "Process " << rank
					<< ": neighbor list for cell " << cell
					<< " not available"
					<< endl;
				exit(EXIT_FAILURE);
			}

			for (const auto& neighbor_i: *neighbors) {
				const auto& neighbor = neighbor_i.first;

				if (neighbor == dccrg::error_cell) {
					continue;
				}

				const auto* const neighbor_data = grid[neighbor];
				if (neighbor_data == nullptr) {
					cerr << "Process " << rank
						<< ": neighbor " << neighbor
						<< " data for cell " << cell
						<< " not available"
						<< endl;
					exit(EXIT_FAILURE);
				}
				if (neighbor_data->is_alive) {
					cell_data->live_neighbor_count++;
				}
			}
		}

		// calculate the next turn
		for (const auto& cell: cells) {
			if (cell.data->live_neighbor_count == 3) {
				cell.data->is_alive = 1;
			} else if (cell.data->live_neighbor_count != 2) {
				cell.data->is_alive = 0;
			}
		}
	}

	if (rank == 0) {
		visit_file.close();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
