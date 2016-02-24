/*
Optimized version of solve.hpp, uses cached cell data pointers, etc.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute
Copyright 2016 Ilja Honkonen

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
template<class UserGeometry> class Solve
{
public:

	/*!
	Calculates live neighbors of one cell at given index in dccrg's pointer cache.

	Works with adaptive mesh refinement.
	Returns next index to process or one past end of cache.
	Does nothing if given index is outside of pointer cache.
	*/
	static size_t get_live_neighbors(dccrg::Dccrg<Cell, UserGeometry>& game_grid, size_t index)
	{
		using std::get;

		const auto& cache = game_grid.get_cell_data_pointers();
		if (index >= cache.size()) {
			return cache.size();
		}

		const auto& cell_id = get<0>(cache[index]);
		if (cell_id == dccrg::error_cell) {
			return index + 1;
		}

		const auto& cell_offset = get<2>(cache[index]);
		if (cell_offset[0] != 0 or cell_offset[1] != 0 or cell_offset[2] != 0) {
			std::cerr << "Invalid index for cell: " << index
				<< "with offset: " << cell_offset[0] << ", "
				<< cell_offset[1] << ", "
				<< cell_offset[2] << std::endl;
			abort();
		}

		const auto& cell_parent_id = game_grid.mapping.get_level_0_parent(cell_id);
		auto* const cell_data = get<1>(cache[index]);
		for (index++; index < cache.size(); index++) {

			const auto& neighbor_id = get<0>(cache[index]);
			if (neighbor_id == dccrg::error_cell) {
				index++;
				break;
			}

			const auto& offset = get<2>(cache[index]);
			if (offset[0] == 0 and offset[1] == 0 and offset[2] == 0) {
				break;
			}

			const auto& neighbor_parent_id = game_grid.mapping.get_level_0_parent(neighbor_id);
			if (neighbor_parent_id == cell_parent_id) {
				continue;
			}

			const auto* const neighbor_data = get<1>(cache[index]);

			// check that neighbor hasn't been recorded as alive
			if (neighbor_data->data[0] == 0) {
				for (size_t i = 1; i < cell_data->data.size(); i++) {
					if (cell_data->data[i] == neighbor_parent_id) {
						std::cerr << "Neighbor " << neighbor_id
							<< " of cell " << cell_id
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

		return index;
	}

	/*!
	Assigns all unseen live neighbors from siblings of cell at given index in dccrg's pointer cache.

	Returns the next index to process or one past end of cache.
	*/
	static size_t spread_live_neighbor_data(dccrg::Dccrg<Cell, UserGeometry>& game_grid, size_t index)
	{
		using std::get;

		const auto& cache = game_grid.get_cell_data_pointers();
		if (index >= cache.size()) {
			return cache.size();
		}

		const auto& cell_id = get<0>(cache[index]);
		if (cell_id == dccrg::error_cell) {
			return index + 1;
		}

		const auto& cell_offset = get<2>(cache[index]);
		if (cell_offset[0] != 0 or cell_offset[1] != 0 or cell_offset[2] != 0) {
			std::cerr << "Invalid index for cell: " << index
				<< "with offset: " << cell_offset[0] << ", "
				<< cell_offset[1] << ", "
				<< cell_offset[2] << std::endl;
			abort();
		}

		const auto& cell_parent_id = game_grid.mapping.get_level_0_parent(cell_id);

		auto* const cell_data = get<1>(cache[index]);
		for (index++; index < cache.size(); index++) {

			const auto& neighbor_id = get<0>(cache[index]);
			if (neighbor_id == dccrg::error_cell) {
				index++;
				break;
			}

			const auto& offset = get<2>(cache[index]);
			if (offset[0] == 0 and offset[1] == 0 and offset[2] == 0) {
				break;
			}

			const auto& neighbor_parent_id = game_grid.mapping.get_level_0_parent(neighbor_id);
			if (neighbor_parent_id != cell_parent_id) {
				continue;
			}

			const auto* const neighbor_data = get<1>(cache[index]);
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
						std::cerr << "No room in live neighbor list of cell " << cell_id
							<< " while trying to insert " << neighbor_data->data[i]
							<< " at index " << i << " in original list"
							<< std::endl;
						abort();
					}
				}
			}
		}

		return index;
	}

	static void get_live_neighbors(dccrg::Dccrg<Cell, UserGeometry>& game_grid)
	{
		using std::get;

		const auto& cache_size = game_grid.get_cell_data_pointers().size();


		size_t i = 0;
		while (i < cache_size) {
			i = get_live_neighbors(game_grid, i);
		}
		game_grid.update_copies_of_remote_neighbors();

		// spread live neighbor info between siblings
		size_t iterations = (1 << game_grid.get_maximum_refinement_level()) - 1;
		for (size_t iter = 0; iter < iterations; iter++) {
			i = 0;
			while (i < cache_size) {
				i = spread_live_neighbor_data(game_grid, i);
			}
			game_grid.update_copies_of_remote_neighbors();
		}

		// calculate the next turn
		for (auto& cell: game_grid.cells) {
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
};

#endif

