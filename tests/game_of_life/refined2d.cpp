/*
Tests the grid with a game of life on a refined grid in 2 D, emulating unrefined behaviour:
-siblings consider all of each others neighbors
-neighbor counts are considered on the level of unrefined cells

Copyright 2010, 2011, 2012, 2013, 2014,
2015 Finnish Meteorological Institute

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
#include "array"
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
	uint64_t
		// total live neighbor count for all siblings
		total_live_neighbor_count,
		// must be next to live_unrefined_neighbors
		is_alive;

	// record live neighbors of refined cells to calculate the above
	std::array<uint64_t, 3> live_unrefined_neighbors;

	// only count one sibling of an unrefined cell (their states of life should be identical)
	std::array<uint64_t, 8> child_of_processed;


	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(&(this->is_alive), 4, MPI_UINT64_T);
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

	const std::array<uint64_t, 3> grid_length = {{15, 15, 1}};

	#define NEIGHBORHOOD_SIZE 1
	game_grid.initialize(grid_length, comm, "RANDOM", NEIGHBORHOOD_SIZE);

	const double cell_length = 1.0 / grid_length[0];
	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	game_grid.set_geometry(geom_params);

	// create a blinker
	#define BLINKER_START 198
	vector<uint64_t> blinker_cells{
		BLINKER_START,
		BLINKER_START + 1,
		BLINKER_START + 2
	};
	for (const auto& cell: blinker_cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a toad
	#define TOAD_START 188
	vector<uint64_t> toad_cells{
		TOAD_START,
		TOAD_START + 1,
		TOAD_START + 2,
		TOAD_START + 1 + grid_length[0],
		TOAD_START + 2 + grid_length[0],
		TOAD_START + 3 + grid_length[0]
	};
	for (const auto& cell: toad_cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a beacon
	#define BEACON_START 137
	vector<uint64_t> beacon_cells{
		BEACON_START,
		BEACON_START + 1,
		BEACON_START - grid_length[0],
		BEACON_START + 1 - grid_length[0],
		BEACON_START + 2 - 2 * grid_length[0],
		BEACON_START + 3 - 2 * grid_length[0],
		BEACON_START + 2 - 3 * grid_length[0],
		BEACON_START + 3 - 3 * grid_length[0]
	};
	for (const auto& cell: beacon_cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a glider
	#define GLIDER_START 143
	vector<uint64_t> glider_cells{
		GLIDER_START + 1,
		GLIDER_START + 2 - grid_length[0],
		GLIDER_START - 2 * grid_length[0],
		GLIDER_START + 1 - 2 * grid_length[0],
		GLIDER_START + 2 - 2 * grid_length[0]
	};
	for (const auto& cell: glider_cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a block
	#define BLOCK_START 47
	vector<uint64_t> block_cells{
		BLOCK_START,
		BLOCK_START + 1,
		BLOCK_START - grid_length[0],
		BLOCK_START + 1 - grid_length[0]
	};
	for (const auto& cell: block_cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a beehive
	#define BEEHIVE_START 51
	vector<uint64_t> beehive_cells{
		BEEHIVE_START - grid_length[0],
		BEEHIVE_START + 1,
		BEEHIVE_START + 2,
		BEEHIVE_START + 1 - 2 * grid_length[0],
		BEEHIVE_START + 2 - 2 * grid_length[0],
		BEEHIVE_START + 3 - grid_length[0]
	};
	for (const auto& cell: beehive_cells) {
		auto* const cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// refine half of the grid randomly
	vector<uint64_t> cells = game_grid.get_cells();
	random_shuffle(cells.begin(), cells.end());
	for (int i = 0; i < int(cells.size() / 2); i++) {
		game_grid.refine_completely(cells[i]);
	}
	vector<uint64_t> new_cells = game_grid.stop_refining();
	// assign parents' state to children
	for (const auto& cell: new_cells) {
		auto
			*new_cell_data = game_grid[cell],
			*parent_data = game_grid[game_grid.get_parent(cell)];

		new_cell_data->is_alive = parent_data->is_alive;
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "tests/game_of_life/refined2d_" << rank << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (rank == 0) {
		visit_file.open("tests/game_of_life/refined2d.visit");
		visit_file << "!NBLOCKS " << comm_size << endl;
	}

	#define TIME_STEPS 25
	if (rank == 0) {
		cout << "step: ";
		cout.flush();
	}
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.update_copies_of_remote_neighbors();
		cells = game_grid.get_cells();
		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
		sort(cells.begin(), cells.end());

		if (rank == 0) {
			cout << step << " ";
			cout.flush();
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
				visit_file << "refined2d_" << process
					<< "_" << step_string.str() << suffix.str()
					<< "\n";
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
		for (const auto& cell: cells) {

			auto* const cell_data = game_grid[cell];

			if (cell_data->is_alive > 0) {
				outfile << "1";
			} else {
				outfile << "0";
			}
			outfile << endl;

		}

		// write each cells total live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {

			game_of_life_cell* cell_data = game_grid[cell];

			outfile << cell_data->total_live_neighbor_count << endl;

		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
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
		for (const auto& cell: cells) {
			outfile << cell << endl;
		}
		outfile.close();

		// get the neighbor counts of every cell
		for (const auto& cell: cells) {

			auto* const cell_data = game_grid[cell];
			if (cell_data == NULL) {
				cout << __FILE__ << ":" << __LINE__ << " no data for cell: " << cell << endl;
				exit(EXIT_FAILURE);
			}

			cell_data->total_live_neighbor_count = 0;

			for (int i = 0; i < 3; i++) {
				cell_data->live_unrefined_neighbors[i] = 0;
			}

			for (int i = 0; i < 8; i++) {
				cell_data->child_of_processed[i] = 0;
			}

			const auto* const neighbors = game_grid.get_neighbors_of(cell);
			// unrefined cells just consider neighbor counts at the level of unrefined cells
			if (game_grid.get_refinement_level(cell) == 0) {

				for (const auto& neighbor: *neighbors) {

					if (neighbor == dccrg::error_cell) {
						continue;
					}

					auto* const neighbor_data = game_grid[neighbor];
					if (neighbor_data == NULL) {
						cout << __FILE__ << ":" << __LINE__
							<< " no data for neighbor of cell " << cell
							<< ": " << neighbor
							<< endl;
						exit(EXIT_FAILURE);
					}

					if (game_grid.get_refinement_level(neighbor) == 0) {
						if (neighbor_data->is_alive) {
							cell_data->total_live_neighbor_count++;
						}
					// consider only one sibling...
					} else {

						bool sibling_processed = false;
						uint64_t parent_of_neighbor = game_grid.get_parent(neighbor);
						for (int i = 0; i < 8; i++) {
							if (cell_data->child_of_processed[i] == parent_of_neighbor) {
								sibling_processed = true;
								break;
							}
						}

						// ...by recording its parent
						if (sibling_processed) {
							continue;
						} else {
							for (int i = 0; i < 8; i++) {
								if (cell_data->child_of_processed[i] == 0) {
									cell_data->child_of_processed[i] = parent_of_neighbor;
									break;
								}
							}
						}

						if (neighbor_data->is_alive > 0) {
							cell_data->total_live_neighbor_count++;
						}
					}
				}

			// refined cells total the neighbor counts of siblings
			} else {

				for (const auto& neighbor: *neighbors) {

					if (neighbor == dccrg::error_cell) {
						continue;
					}

					auto* const neighbor_data = game_grid[neighbor];
					if (neighbor_data == NULL) {
						cout << __FILE__ << ":" << __LINE__
							<< " no data for neighbor of refined cell " << cell
							<< ": " << neighbor
							<< endl;
						exit(EXIT_FAILURE);
					}

					if (game_grid.get_refinement_level(neighbor) == 0) {

						// larger neighbors appear several times in the neighbor list
						bool neighbor_processed = false;
						for (int i = 0; i < 8; i++) {
							if (cell_data->child_of_processed[i] == neighbor) {
								neighbor_processed = true;
								break;
							}
						}

						if (neighbor_processed) {
							continue;
						} else {
							for (int i = 0; i < 8; i++) {
								if (cell_data->child_of_processed[i] == 0) {
									cell_data->child_of_processed[i] = neighbor;
									break;
								}
							}
						}

						if (neighbor_data->is_alive) {
							for (int i = 0; i < 3; i++) {
								if (cell_data->live_unrefined_neighbors[i] == 0) {
									cell_data->live_unrefined_neighbors[i] = neighbor;
									break;
								}
							}
						}

					// consider only one sibling of all parents of neighboring cells...
					} else {

						// ignore own siblings
						if (game_grid.get_parent(cell) == game_grid.get_parent(neighbor)) {
							continue;
						}

						bool sibling_processed = false;
						uint64_t parent_of_neighbor = game_grid.get_parent(neighbor);
						for (int i = 0; i < 8; i++) {
							if (cell_data->child_of_processed[i] == parent_of_neighbor) {
								sibling_processed = true;
								break;
							}
						}

						if (sibling_processed) {
							continue;
						} else {
							for (int i = 0; i < 8; i++) {
								if (cell_data->child_of_processed[i] == 0) {
									cell_data->child_of_processed[i] = parent_of_neighbor;
									break;
								}
							}
						}

						// ...by recording which parents have been considered
						if (neighbor_data->is_alive) {
							for (int i = 0; i < 3; i++) {
								if (cell_data->live_unrefined_neighbors[i] == 0) {
									cell_data->live_unrefined_neighbors[i] = parent_of_neighbor;
									break;
								}
							}
						}
					}
				}
			}

		}
		game_grid.update_copies_of_remote_neighbors();

		// get the total neighbor counts of refined cells
		for (const auto& cell: cells) {

			if (game_grid.get_refinement_level(cell) == 0) {
				continue;
			}
			auto* const cell_data = game_grid[cell];

			unordered_set<uint64_t> current_live_unrefined_neighbors;
			for (int i = 0; i < 3; i++) {
				current_live_unrefined_neighbors.insert(cell_data->live_unrefined_neighbors[i]);
			}

			const auto* const neighbors = game_grid.get_neighbors_of(cell);
			for (const auto& neighbor: *neighbors) {
				if (neighbor == dccrg::error_cell) {
					continue;
				}

				if (game_grid.get_refinement_level(neighbor) == 0) {
					continue;
				}

				// total live neighbors counts only between siblings
				if (game_grid.get_parent(cell) != game_grid.get_parent(neighbor)) {
					continue;
				}

				auto* const neighbor_data = game_grid[neighbor];
				for (int i = 0; i < 3; i++) {
					current_live_unrefined_neighbors.insert(neighbor_data->live_unrefined_neighbors[i]);
				}
			}

			current_live_unrefined_neighbors.erase(0);
			cell_data->total_live_neighbor_count += current_live_unrefined_neighbors.size();
		}

		// calculate the next turn
		for (const auto& cell: cells) {

			auto* const cell_data = game_grid[cell];

			if (cell_data->total_live_neighbor_count == 3) {
				cell_data->is_alive = 1;
			} else if (cell_data->total_live_neighbor_count != 2) {
				cell_data->is_alive = 0;
			}
		}

	}

	if (rank == 0) {
		cout << endl;
		visit_file.close();
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}
