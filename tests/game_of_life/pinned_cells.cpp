/*
As unrefined2d.cpp but pinns cells to particular processes.

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

#include "array"
#include "algorithm"
#include "cmath"
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

	Dccrg<game_of_life_cell, Stretched_Cartesian_Geometry> grid;

	const std::array<uint64_t, 3> grid_length = {{15, 15, 1}};
	const double cell_length = 1.0 / grid_length[0];

	#define NEIGHBORHOOD_SIZE 1
	grid
		.set_initial_length(grid_length)
		.set_neighborhood_length(NEIGHBORHOOD_SIZE)
		.set_maximum_refinement_level(-1)
		.set_load_balancing_method("RANDOM")
		.initialize(comm);

	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
		for (size_t i = 0; i <= grid_length[dimension]; i++) {
			geom_params.coordinates[dimension].push_back(double(i) * cell_length);
		}
	}
	grid.set_geometry(geom_params);

	// create a blinker
	#define BLINKER_START 198
	const vector<uint64_t> blinker_cells{
		BLINKER_START,
		BLINKER_START + 1,
		BLINKER_START + 2
	};
	for (const auto& cell: blinker_cells) {
		game_of_life_cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a toad
	#define TOAD_START 188
	const vector<uint64_t> toad_cells{
		TOAD_START,
		TOAD_START + 1,
		TOAD_START + 2,
		TOAD_START + 1 + grid_length[0],
		TOAD_START + 2 + grid_length[0],
		TOAD_START + 3 + grid_length[0]
	};
	for (const auto& cell: toad_cells) {
		game_of_life_cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a beacon
	#define BEACON_START 137
	const vector<uint64_t> beacon_cells{
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
		game_of_life_cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a glider
	#define GLIDER_START 143
	const vector<uint64_t> glider_cells{
		GLIDER_START + 1,
		GLIDER_START + 2 - grid_length[0],
		GLIDER_START - 2 * grid_length[0],
		GLIDER_START + 1 - 2 * grid_length[0],
		GLIDER_START + 2 - 2 * grid_length[0]
	};
	for (const auto& cell: glider_cells) {
		game_of_life_cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a block
	#define BLOCK_START 47
	const vector<uint64_t> block_cells{
		BLOCK_START,
		BLOCK_START + 1,
		BLOCK_START - grid_length[0],
		BLOCK_START + 1 - grid_length[0]
	};
	for (const auto& cell: block_cells) {
		game_of_life_cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// create a beehive
	#define BEEHIVE_START 51
	const vector<uint64_t> beehive_cells{
		BEEHIVE_START - grid_length[0],
		BEEHIVE_START + 1,
		BEEHIVE_START + 2,
		BEEHIVE_START + 1 - 2 * grid_length[0],
		BEEHIVE_START + 2 - 2 * grid_length[0],
		BEEHIVE_START + 3 - grid_length[0]
	};
	for (const auto& cell: beehive_cells) {
		game_of_life_cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->is_alive = 1;
	}

	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "tests/game_of_life/pinned_cells_" << rank << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (rank == 0) {
		visit_file.open("tests/game_of_life/pinned_cells.visit");
		visit_file << "!NBLOCKS " << comm_size << endl;
		cout << "step: ";
		cout.flush();
	}

	#define TIME_STEPS 25
	for (int step = 0; step < TIME_STEPS; step++) {

		// refine random unrefined cells and unrefine random refined cells
		// TODO merge with identical code in unrefine2d, ...
		auto cells = grid.get_cells();
		random_shuffle(cells.begin(), cells.end());

		if (step % 2 == 0) {

			for (uint64_t i = 0, refined = 0;
				i < cells.size() && refined <= grid_length[0] * grid_length[1] / (5 * comm_size);
				i++
			) {
				if (grid.get_refinement_level(cells[i]) == 0) {
					if (!grid.refine_completely(cells[i])) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Couldn't refine cell " << cells[i]
							<< std::endl;
						abort();
					}
					refined++;
				}
			}

		} else {

			for (uint64_t i = 0, unrefined = 0;
				i < cells.size() && unrefined <= grid_length[0] * grid_length[1] / (4 * comm_size);
				i++
			) {
				if (grid.get_refinement_level(cells[i]) > 0) {
					if (!grid.unrefine_completely(cells[i])) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Couldn't unrefine cell " << cells[i]
							<< std::endl;
						abort();
					}
					unrefined++;
				}
			}
		}

		vector<uint64_t> new_cells = grid.stop_refining();

		// assign parents' state to children
		for (const auto& new_cell: new_cells) {
			game_of_life_cell* new_cell_data = grid[new_cell];
			if (new_cell_data == NULL) {
				cout << __FILE__ << ":" << __LINE__
					<< " No data for created cell " << new_cell
					<< endl;
				abort();
			}
			game_of_life_cell* parent_data = grid[grid.get_parent(new_cell)];
			if (parent_data == NULL) {
				cout << __FILE__ << ":" << __LINE__
					<< " No data for parent cell " << grid.get_parent(new_cell)
					<< endl;
				abort();
			}
			new_cell_data->is_alive = parent_data->is_alive;
		}

		// "interpolate" parent cell's value from unrefined children
		auto removed_cells = grid.get_removed_cells();
		for (const auto& removed_cell: removed_cells) {
			game_of_life_cell* removed_cell_data = grid[removed_cell];
			if (removed_cell_data == NULL) {
				cout << __FILE__ << ":" << __LINE__
					<< " no data for removed cell after unrefining: " << removed_cell
					<< endl;
				abort();
			}
			game_of_life_cell* parent_data = grid[grid.mapping.get_parent(removed_cell)];
			if (parent_data == NULL) {
				cout << __FILE__ << ":" << __LINE__
					<< " no data for parent cell after unrefining: "
					<< grid.mapping.get_parent(removed_cell)
					<< endl;
				abort();
			}
			parent_data->is_alive = removed_cell_data->is_alive;
		}
		grid.clear_refined_unrefined_data();

		// pin cells to process 0 either inside or outside the circle
		for (const auto& cell: grid.get_cells()) {
			const std::array<double, 3> cell_center = grid.geometry.get_center(cell);

			const double distance
				= std::pow(cell_center[0] - 0.5, 2.0)
				+ std::pow(cell_center[1] - 0.5, 2.0);

			if (step % 2 == 0) {
				if (distance <= 0.3 * 0.3) {
					if (!grid.pin(cell, 0)) {
						cout << __FILE__ << ":" << __LINE__
							<< " Couldn't pin cell " << cell
							<< endl;
						abort();
					}
				} else {
					if (!grid.unpin(cell)) {
						cout << __FILE__ << ":" << __LINE__
							<< " Couldn't unpin cell " << cell
							<< endl;
						abort();
					}
				}
			} else {
				if (distance <= 0.3 * 0.3) {
					if (!grid.unpin(cell)) {
						cout << __FILE__ << ":" << __LINE__
							<< " Couldn't unpin cell " << cell
							<< endl;
						abort();
					}
				} else {
					if (!grid.pin(cell, 0)) {
						cout << __FILE__ << ":" << __LINE__
							<< " Couldn't pin cell " << cell
							<< endl;
						abort();
					}
				}
			}
		}

		grid.balance_load();
		grid.update_copies_of_remote_neighbors();
		cells = grid.get_cells();
		// the library writes the grid into a file in ascending
		// cell order, do the same for the grid data at every time step
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
				visit_file << "pinned_cells_" << process << "_"
					<< step_string.str() << suffix.str()
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
			const auto* const cell_data = grid[cell];
			if (cell_data->is_alive == 1) {
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
			const auto* const cell_data = grid[cell];
			outfile << cell_data->total_live_neighbor_count << endl;
		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			const auto* const neighbors = grid.get_neighbors_of(cell);
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
			auto* const cell_data = grid[cell];
			cell_data->total_live_neighbor_count = 0;

			for (int i = 0; i < 3; i++) {
				cell_data->live_unrefined_neighbors[i] = 0;
			}

			for (int i = 0; i < 8; i++) {
				cell_data->child_of_processed[i] = 0;
			}

			const auto* const neighbors = grid.get_neighbors_of(cell);
			// unrefined cells just consider neighbor counts at the level of unrefined cells
			if (grid.get_refinement_level(cell) == 0) {

				for (const auto& neighbor_i: *neighbors) {
					const auto& neighbor = neighbor_i.first;

					if (neighbor == dccrg::error_cell) {
						continue;
					}

					game_of_life_cell* neighbor_data = grid[neighbor];
					if (neighbor_data == NULL) {
						cout << __FILE__ << ":" << __LINE__
							<< " no data for neighbor of cell " << cell
							<< ": " << neighbor
							<< endl;
						exit(EXIT_FAILURE);
					}

					if (grid.get_refinement_level(neighbor) == 0) {
						if (neighbor_data->is_alive) {
							cell_data->total_live_neighbor_count++;
						}
					// consider only one sibling...
					} else {

						bool sibling_processed = false;
						uint64_t parent_of_neighbor = grid.get_parent(neighbor);
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

						if (neighbor_data->is_alive) {
							cell_data->total_live_neighbor_count++;
						}
					}
				}

			// refined cells total the neighbor counts of siblings
			} else {

				for (const auto& neighbor_i: *neighbors) {
					const auto& neighbor = neighbor_i.first;

					if (neighbor == dccrg::error_cell) {
						continue;
					}

					game_of_life_cell* neighbor_data = grid[neighbor];
					if (neighbor_data == NULL) {
						cout << __FILE__ << ":" << __LINE__
							<< " no data for neighbor of refined cell " << cell
							<< ": " << neighbor
							<< endl;
						exit(EXIT_FAILURE);
					}

					if (grid.get_refinement_level(neighbor) == 0) {
						// TODO merge solver with the one in unrefined2d, etc
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
						if (grid.get_parent(cell) == grid.get_parent(neighbor)) {
							continue;
						}

						bool sibling_processed = false;
						uint64_t parent_of_neighbor = grid.get_parent(neighbor);
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
		grid.update_copies_of_remote_neighbors();

		// get the total neighbor counts of refined cells
		for (const auto& cell: cells) {
			auto* const cell_data = grid[cell];
			if (grid.get_refinement_level(cell) == 0) {
				continue;
			}

			unordered_set<uint64_t> current_live_unrefined_neighbors;
			for (int i = 0; i < 3; i++) {
				current_live_unrefined_neighbors.insert(cell_data->live_unrefined_neighbors[i]);
			}

			const auto* const neighbors = grid.get_neighbors_of(cell);
			for (const auto& neighbor_i: *neighbors) {
				const auto& neighbor = neighbor_i.first;

				if (neighbor == dccrg::error_cell) {
					continue;
				}

				if (grid.get_refinement_level(neighbor) == 0) {
					continue;
				}

				// total live neighbors counts only between siblings
				if (grid.get_parent(cell) != grid.get_parent(neighbor)) {
					continue;
				}

				game_of_life_cell* neighbor_data = grid[neighbor];
				for (int i = 0; i < 3; i++) {
					current_live_unrefined_neighbors.insert(neighbor_data->live_unrefined_neighbors[i]);
				}
			}

			current_live_unrefined_neighbors.erase(0);
			cell_data->total_live_neighbor_count += current_live_unrefined_neighbors.size();
		}

		// calculate the next turn
		for (const auto& cell: grid.get_cells()) {
			auto* const cell_data = grid[cell];
			if (cell_data->total_live_neighbor_count == 3) {
				cell_data->is_alive = 1;
			} else if (cell_data->total_live_neighbor_count != 2) {
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

