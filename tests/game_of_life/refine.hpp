/*
A class for refining the grid in game of life tests of dccrg.

Copyright 2011, 2012, 2013, 2014, 2015, 2016 Finnish Meteorological Institute

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

#ifndef REFINE_HPP
#define REFINE_HPP

#include "cstdint"
#include "vector"

#include "dccrg.hpp"

#include "cell.hpp"

/*!
Refines and unrefines given grid randomly.
*/
template<class UserGeometry> class Refine
{
public:

	static void refine(
		dccrg::Dccrg<Cell, UserGeometry>& grid,
		const int grid_size,
		const int step,
		const int processes
	) {
		// refine random unrefined cells and unrefine random refined cells
		std::vector<uint64_t> cells;
		for (const auto& cell: grid.local_cells) {
			cells.push_back(cell.id);
		}
		std::random_shuffle(cells.begin(), cells.end());

		if (step % 2 == 0) {

			for (int i = 0, refined = 0;
				i < int(cells.size()) && refined <= grid_size * grid_size / (5 * processes);
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

			for (int i = 0, unrefined = 0;
				i < int(cells.size()) && unrefined <= grid_size * grid_size / (4 * processes);
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

		const auto new_cells = grid.stop_refining();

		// assign parents' state to children
		for (const auto& new_cell: new_cells) {

			auto* const new_cell_data = grid[new_cell];
			if (new_cell_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for created cell " << new_cell
					<< std::endl;
				abort();
			}
			auto* const parent_data = grid[grid.get_parent(new_cell)];
			if (parent_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for parent cell " << grid.get_parent(new_cell)
					<< std::endl;
				abort();
			}
			new_cell_data->data[0] = parent_data->data[0];
		}

		// "interpolate" parent cell's value from removed children
		const auto removed_cells = grid.get_removed_cells();

		for (const auto& removed_cell: removed_cells) {

			auto* const removed_cell_data = grid[removed_cell];
			if (removed_cell_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for removed cell after unrefining: "
					<< removed_cell
					<< std::endl;
				abort();
			}

			auto* const parent_data = grid[grid.mapping.get_parent(removed_cell)];
			if (parent_data == nullptr) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for parent cell after unrefining: "
					<< grid.mapping.get_parent(removed_cell)
					<< std::endl;
				abort();
			}
			parent_data->data[0] = removed_cell_data->data[0];
		}

		grid.clear_refined_unrefined_data();
	}
};

#endif
