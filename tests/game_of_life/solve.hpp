/*
A solver class for the game of life tests of dccrg.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute
Copyright 2018 Ilja Honkonen

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SOLVE_HPP
#define SOLVE_HPP

#include "cstdint"
#include "iostream"
#include "unordered_map"
#include "unordered_set"
#include "vector"

#include "dccrg.hpp"

/*!
Advances the game of life on given grid one turn.

Works only if cells within the grid have been refined <= 1 time.
*/
template<
	class Cell_Data,
	class Geometry
> void get_live_neighbors(
	dccrg::Dccrg<Cell_Data, Geometry>& grid
) {
	const auto cells = grid.get_cells();
	// get the neighbor counts of every cell
	for (const auto& cell: cells) {

		auto* const cell_data = grid[cell];
		if (cell_data == nullptr) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " no data for cell: " << cell
				<< std::endl;
			abort();
		}

		for (size_t i = 1; i < cell_data->data.size(); i++) {
			cell_data->data[i] = 0;
		}

		const auto& cell_parent_id = grid.mapping.get_level_0_parent(cell);

		const auto* const neighbors = grid.get_neighbors_of(cell);
		for (const auto& neighbor_i: *neighbors) {
			const auto& neighbor = neighbor_i.first;

			if (neighbor == dccrg::error_cell) {
				continue;
			}

			const auto& neighbor_parent_id = grid.mapping.get_level_0_parent(neighbor);
			if (neighbor_parent_id == cell_parent_id) {
				continue;
			}

			const auto* const neighbor_data = grid[neighbor];
			if (neighbor_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for neighbor of cell " << cell
					<< ": " << neighbor
					<< std::endl;
				abort();
			}

			// check that neighbor hasn't been recorded as alive
			if (neighbor_data->data[0] == 0) {
				for (size_t i = 1; i < cell_data->data.size(); i++) {
					if (cell_data->data[i] == neighbor_parent_id) {
						std::cerr << "Neighbor " << neighbor
							<< " of cell " << cell
							<< " at index " << i
							<< " should not be alive."
							<< std::endl;
						abort();
					}
				}
			} else {
				for (size_t i = 1; i < cell_data->data.size(); i++) {
					if (cell_data->data[i] == neighbor_parent_id) {
						break;
					} else if (cell_data->data[i] == dccrg::error_cell) {
						cell_data->data[i] = neighbor_parent_id;
						break;
					} else if (i == cell_data->data.size() - 1) {
						std::cerr << "No more room in live neighbor list." << std::endl;
						abort();
					}
				}
			}
		}
	}
	grid.update_copies_of_remote_neighbors();

	// spread live neighbor info between siblings
	for (const auto& cell: cells) {

		const auto cell_parent_id = grid.mapping.get_level_0_parent(cell);
		auto* const cell_data = grid[cell];

		const auto* const neighbors = grid.get_neighbors_of(cell);
		for (const auto& neighbor_i: *neighbors) {
			const auto& neighbor = neighbor_i.first;

			if (neighbor == dccrg::error_cell) {
				continue;
			}

			if (cell_parent_id != grid.mapping.get_level_0_parent(neighbor)) {
				continue;
			}

			const auto* const neighbor_data = grid[neighbor];
			for (size_t i = 1; i < neighbor_data->data.size(); i++) {
				if (neighbor_data->data[i] == dccrg::error_cell) {
					break;
				}

				for (size_t j = 1; j < cell_data->data.size(); j++) {
					if (cell_data->data[j] == dccrg::error_cell) {
						cell_data->data[j] = neighbor_data->data[i];
						break;
					} else if (cell_data->data[j] == neighbor_data->data[i]) {
						break;
					} else if (j == cell_data->data.size() - 1) {
						std::cerr << "No room in live neighbor list of cell " << cell
							<< " while trying to insert " << neighbor_data->data[i]
							<< " at index " << i << " in original list"
							<< std::endl;
						abort();
					}
				}
			}
		}
	}

	// calculate the next turn
	for (const auto& cell: cells) {
		auto* const cell_data = grid[cell];
		size_t live_neighbors = 0;
		for (size_t i = 1; i < cell_data->data.size(); i++) {
			if (cell_data->data[i] != dccrg::error_cell) {
				live_neighbors++;
			}
			cell_data->data[i] = dccrg::error_cell;
		}
		if (live_neighbors == 3) {
			cell_data->data[0] = 1;
		} else if (live_neighbors != 2) {
			cell_data->data[0] = 0;
		}
	}
}

#endif
