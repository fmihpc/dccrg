/*
Optimized version of solve.hpp, uses cached cell data pointers, etc.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute
Copyright 2016, 2018 Ilja Honkonen

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

#ifndef SOLVE_OPTIMIZED_HPP
#define SOLVE_OPTIMIZED_HPP

#include "cstdint"
#include "iostream"
#include "unordered_map"
#include "unordered_set"
#include "vector"

#include "dccrg.hpp"

#include "cell.hpp"


/*!
Advances the game of life on given grid one turn.

Works only if cells within the grid have been refined <= 1 time.
*/
template<class Cell_Data, class Geometry> void get_live_neighbors(dccrg::Dccrg<Cell_Data, Geometry>& grid)
{
	// get the neighbor counts of every cell
	for (const auto& cell: grid.local_cells()) {
		for (size_t i = 1; i < cell.data->data.size(); i++) {
			cell.data->data[i] = 0;
		}

		const auto cell_parent_id = grid.mapping.get_level_0_parent(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			const auto neighbor_parent_id = grid.mapping.get_level_0_parent(neighbor.id);
			if (neighbor_parent_id == cell_parent_id) {
				continue;
			}

			// check that neighbor hasn't been recorded as alive
			if (neighbor.data->data[0] == 0) {
				for (size_t i = 1; i < cell.data->data.size(); i++) {
					if (cell.data->data[i] == neighbor_parent_id) {
						std::cerr << "Neighbor " << neighbor.id
							<< " of cell " << cell.id
							<< " at index " << i
							<< " should not be alive."
							<< std::endl;
						abort();
					}
				}
			} else {
				for (size_t i = 1; i < cell.data->data.size(); i++) {
					if (cell.data->data[i] == neighbor_parent_id) {
						break;
					} else if (cell.data->data[i] == dccrg::error_cell) {
						cell.data->data[i] = neighbor_parent_id;
						break;
					} else if (i == cell.data->data.size() - 1) {
						std::cerr << "No more room in live neighbor list." << std::endl;
						abort();
					}
				}
			}
		}

	}
	grid.update_copies_of_remote_neighbors();

	// spread live neighbor info between siblings
	for (const auto& cell: grid.local_cells()) {
		const auto cell_parent_id = grid.mapping.get_level_0_parent(cell.id);

		for (const auto& neighbor: cell.neighbors_of) {
			if (cell_parent_id != grid.mapping.get_level_0_parent(neighbor.id)) {
				continue;
			}

			for (size_t i = 1; i < neighbor.data->data.size(); i++) {
				if (neighbor.data->data[i] == dccrg::error_cell) {
					break;
				}

				for (size_t j = 1; j < cell.data->data.size(); j++) {
					if (cell.data->data[j] == dccrg::error_cell) {
						cell.data->data[j] = neighbor.data->data[i];
						break;
					} else if (cell.data->data[j] == neighbor.data->data[i]) {
						break;
					} else if (j == cell.data->data.size() - 1) {
						std::cerr << "No room in live neighbor list of cell " << cell.id
							<< " while trying to insert " << neighbor.data->data[i]
							<< " at index " << i << " in original list"
							<< std::endl;
						abort();
					}
				}
			}
		}
	}

	// calculate the next turn
	for (const auto& cell: grid.local_cells()) {
		size_t live_neighbors = 0;
		for (size_t i = 1; i < cell.data->data.size(); i++) {
			if (cell.data->data[i] != dccrg::error_cell) {
				live_neighbors++;
			}
			cell.data->data[i] = dccrg::error_cell;
		}
		if (live_neighbors == 3) {
			cell.data->data[0] = 1;
		} else if (live_neighbors != 2) {
			cell.data->data[0] = 0;
		}
	}
}

#endif
