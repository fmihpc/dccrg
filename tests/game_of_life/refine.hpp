/*
A class for refining the grid in game of life tests of dccrg.

Copyright 2011, 2012, 2013 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef REFINE_HPP
#define REFINE_HPP

#include "boost/foreach.hpp"
#include "stdint.h"
#include "vector"

#include "../../dccrg.hpp"

#include "cell.hpp"

/*!
Initializes given grid with given size in the first game coordinate with some game of life patterns.
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
		std::vector<uint64_t> cells = grid.get_cells();
		random_shuffle(cells.begin(), cells.end());

		if (step % 2 == 0) {

			for (int i = 0, refined = 0;
				i < int(cells.size()) && refined <= grid_size * grid_size / (5 * processes);
				i++
			) {
				if (grid.get_refinement_level(cells[i]) == 0) {
					grid.refine_completely(cells[i]);
					refined++;
				}
			}

		} else {

			for (int i = 0, unrefined = 0;
				i < int(cells.size()) && unrefined <= grid_size * grid_size / (4 * processes);
				i++
			) {
				if (grid.get_refinement_level(cells[i]) > 0) {
					grid.unrefine_completely(cells[i]);
					unrefined++;
				}
			}
		}

		const std::vector<uint64_t> new_cells = grid.stop_refining();

		// assign parents' state to children
		BOOST_FOREACH(const uint64_t new_cell, new_cells) {

			Cell* new_cell_data = grid[new_cell];
			if (new_cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for created cell " << new_cell
					<< std::endl;
				abort();
			}
			Cell* parent_data = grid[grid.get_parent(new_cell)];
			if (parent_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for parent cell " << grid.get_parent(new_cell)
					<< std::endl;
				abort();
			}
			new_cell_data->data[0] = parent_data->data[0];
		}

		// "interpolate" parent cell's value from unrefined children
		const std::vector<uint64_t> removed_cells = grid.get_removed_cells();

		BOOST_FOREACH(const uint64_t removed_cell, removed_cells) {

			Cell* removed_cell_data = grid[removed_cell];
			if (removed_cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for removed cell after unrefining: "
					<< removed_cell
					<< std::endl;
				abort();
			}

			Cell* parent_data = grid[grid.get_parent_for_removed(removed_cell)];
			if (parent_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " no data for parent cell after unrefining: "
					<< grid.get_parent_for_removed(removed_cell)
					<< std::endl;
				abort();
			}
			parent_data->data[0] = removed_cell_data->data[0];
		}

		grid.clear_refined_unrefined_data();
	}
};

#endif
