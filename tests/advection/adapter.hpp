/*
Adapter class for the advection test program of dccrg.

Copyright 2012, 2013, 2014, 2015, 2016 Finnish Meteorological Institute

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

#ifndef DCCRG_ADVECTION_ADAPTER_HPP
#define DCCRG_ADVECTION_ADAPTER_HPP


#include "cmath"
#include "iostream"
#include "utility"
#include "vector"

#include "dccrg.hpp"

#include "solve.hpp"


/*!
Adapts the grid based on relative difference in variables between neighboring cells.

Refinement level target of cells (RDI = relative diff increase for a cell):
\verbatim
	   RDI               |  ref lvl
[0, 1 * diff_increase[   |     0
[1, 2 * diff_increase[   |     1
[2, 3 * diff_increase[   |     2
...
\endverbatim
*/
template<
	class CellData,
	class Geometry
> void check_for_adaptation(
	const double diff_increase,
	const double diff_threshold,
	const double unrefine_sensitivity,
	std::unordered_set<uint64_t>& cells_to_refine,
	std::unordered_set<uint64_t>& cells_not_to_unrefine,
	std::unordered_set<uint64_t>& cells_to_unrefine,
	dccrg::Dccrg<CellData, Geometry>& grid
) {
	if (grid.get_maximum_refinement_level() == 0) {
		return;
	}

	cells_to_refine.clear();
	cells_not_to_unrefine.clear();
	cells_to_unrefine.clear();

	for (const auto& cell: grid.local_cells) {
		cell.data->max_diff() = 0;
	}

	// collect maximum relative differences
	for (const auto& cell: grid.local_cells) {

		for (const auto& neighbor: cell.neighbors_of) {
			// skip non-face neighbors
			...

			const double diff
				= std::fabs(cell_data->density() - neighbor_data->density())
				/ (std::min(cell_data->density(), neighbor_data->density()) + diff_threshold);

			cell_data->max_diff() = std::max(diff, cell_data->max_diff());

			// maximize diff for local neighbor
			if (grid.is_local(neighbor)) {
				neighbor_data->max_diff() = std::max(diff, neighbor_data->max_diff());
			}
		}
	}

	// decide whether to refine or unrefine cells
	for (const auto& cell: cells) {

		const int refinement_level = grid.get_refinement_level(cell);

		// refine / unrefine if max relative difference larger / smaller than:
		const double
			refine_diff = (refinement_level + 1) * diff_increase,
			unrefine_diff = unrefine_sensitivity * refine_diff;

		const auto siblings = grid.get_all_children(grid.get_parent(cell));
		if (siblings.size() == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No siblings for cell " << cell
				<< std::endl;
			abort();
		}

		auto* const cell_data = grid[cell];

		const double diff = cell_data->max_diff();

		// refine
		if (diff > refine_diff) {

			cells_to_refine.insert(cell);

			for (const auto& sibling: siblings) {
				cells_to_unrefine.erase(sibling);
				cells_not_to_unrefine.erase(sibling);
			}

		// make sure siblings aren't unrefined
		} else if (diff >= unrefine_diff) {

			bool dont_unrefine = true;

			for (const auto& sibling: siblings) {
				if (cells_to_refine.count(sibling) > 0
				|| cells_not_to_unrefine.count(sibling) > 0) {
					dont_unrefine = false;
					break;
				}
			}

			if (dont_unrefine && grid.get_refinement_level(cell) > 0) {
				cells_not_to_unrefine.insert(cell);

				for (const auto& sibling: siblings) {
					cells_to_unrefine.erase(sibling);
				}
			}

		// unrefine
		} else {

			bool unrefine = true;

			for (const auto& sibling: siblings) {
				if (cells_to_refine.count(sibling) > 0
				|| cells_not_to_unrefine.count(sibling) > 0) {
					unrefine = false;
					break;
				}
			}

			if (unrefine && grid.get_refinement_level(cell) > 0) {
				cells_to_unrefine.insert(cell);
			}
		}
	}
}


/*!
Refines/Unrefines given cells in the grid.

Returns the number of created and removed cells.
Clears given sets of cells after executing refines.
*/
template<
	class CellData,
	class Geometry
> std::pair<uint64_t, uint64_t> adapt_grid(
	std::unordered_set<uint64_t>& cells_to_refine,
	std::unordered_set<uint64_t>& cells_not_to_unrefine,
	std::unordered_set<uint64_t>& cells_to_unrefine,
	dccrg::Dccrg<CellData, Geometry>& grid
) {
	if (grid.get_maximum_refinement_level() == 0) {
		return std::make_pair(0, 0);
	}

	// record the number of adaptivity failures
	uint64_t
		failed_refines = 0,
		failed_unrefines = 0,
		failed_dont_unrefines = 0;

	for (const auto& cell: cells_to_refine) {
		if (!grid.refine_completely(cell)) {
			failed_refines++;
		}
	}
	cells_to_refine.clear();

	for (const auto& cell: cells_not_to_unrefine) {
		if (!grid.dont_unrefine(cell)) {
			failed_dont_unrefines++;
		}
	}
	cells_not_to_unrefine.clear();

	for (const auto& cell: cells_to_unrefine) {
		if (!grid.unrefine_completely(cell)) {
			failed_unrefines++;
		}
	}
	cells_to_unrefine.clear();

	/*
	Refines
	*/

	// assign parents' state to children
	const std::vector<uint64_t> new_cells = grid.stop_refining();

	for (const auto& new_cell: new_cells) {
		auto* const new_cell_data = grid[new_cell];
		if (new_cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for created cell " << new_cell
				<< std::endl;
			abort();
		}

		auto* const parent_data = grid[grid.get_parent(new_cell)];
		if (parent_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell " << grid.get_parent(new_cell)
				<< std::endl;
			abort();
		}

		new_cell_data->density() = parent_data->density();
		new_cell_data->flux() = 0;

		const std::array<double, 3> cell_center = grid.geometry.get_center(new_cell);
		new_cell_data->vx() = get_vx(cell_center[1]);
		new_cell_data->vy() = get_vy(cell_center[0]);
		new_cell_data->vz() = 0;
	}

	/*
	Unrefines
	*/

	const std::vector<uint64_t> removed_cells = grid.get_removed_cells();

	// optimize by gathering all parents of removed cells
	std::unordered_set<uint64_t> parents;
	for (const auto& removed_cell: removed_cells) {
		parents.insert(grid.mapping.get_parent(removed_cell));
	}

	// initialize parent data
	for (const auto& parent: parents) {
		auto* const parent_data = grid[parent];
		if (parent_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell: " << parent
				<< std::endl;
			abort();
		}

		parent_data->density() = 0;
		parent_data->flux() = 0;

		const std::array<double, 3> cell_center = grid.geometry.get_center(parent);
		parent_data->vx() = get_vx(cell_center[1]);
		parent_data->vy() = get_vy(cell_center[0]);
		parent_data->vz() = 0;
	}

	// average parents' density from their children
	for (const auto& removed_cell: removed_cells) {

		auto* const removed_cell_data = grid[removed_cell];
		if (removed_cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for removed cell after unrefining: " << removed_cell
				<< std::endl;
			abort();
		}

		auto* const parent_data = grid[grid.mapping.get_parent(removed_cell)];
		if (parent_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell after unrefining: "
				<< grid.mapping.get_parent(removed_cell)
				<< std::endl;
			abort();
		}

		parent_data->density() += removed_cell_data->density() / 8;
	}

	grid.clear_refined_unrefined_data();

	return std::make_pair(new_cells.size(), removed_cells.size());
}

#endif
