/*
A solver class for the game of life tests of dccrg.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016, 2018 Finnish Meteorological Institute

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

#include "iostream"
#include "cstdint"
#include "vector"

#include "../../dccrg.hpp"

#include "cell.hpp"

/*!
Advances the game of life on given grid one turn.

Works only if cells within the grid have been refined <= 1 time.
*/
template<class Geometry> void solve(dccrg::Dccrg<Cell, Geometry>& grid) {
	// get the neighbor counts of every cell
	for (const auto& cell: grid.local_cells()) {
		for (int i = 0; i < 12; i++) {
			cell.data->data[1 + i] = 0;
		}

		// unrefined cells just consider neighbor counts at the level of unrefined cells
		if (grid.get_refinement_level(cell.id) == 0) {

			for (const auto& neighbor: cell.neighbors_of) {
				if (grid.get_refinement_level(neighbor.id) == 0) {
					if (neighbor.data->data[0]) {
						cell.data->data[1]++;
					}
				// consider only one sibling...
				} else {

					bool sibling_processed = false;
					uint64_t parent_of_neighbor = grid.mapping.get_parent(neighbor.id);
					for (int i = 0; i < 8; i++) {
						if (cell.data->data[5 + i] == parent_of_neighbor) {
							sibling_processed = true;
							break;
						}
					}

					// ...by recording its parent
					if (sibling_processed) {
						continue;
					} else {
						for (int i = 0; i < 8; i++) {
							if (cell.data->data[5 + i] == 0) {
								cell.data->data[5 + i] = parent_of_neighbor;
								break;
							}
						}
					}

					if (neighbor.data->data[0]) {
						cell.data->data[1]++;
					}
				}
			}

		// refined cells total the neighbor counts of siblings
		} else {

			for (const auto& neighbor: cell.neighbors_of) {
				if (grid.get_refinement_level(neighbor.id) == 0) {

					// larger neighbors appear several times in the neighbor list
					bool neighbor_processed = false;
					for (int i = 0; i < 8; i++) {
						if (cell.data->data[5 + i] == neighbor.id) {
							neighbor_processed = true;
							break;
						}
					}

					if (neighbor_processed) {
						continue;
					} else {
						for (int i = 0; i < 8; i++) {
							if (cell.data->data[5 + i] == 0) {
								cell.data->data[5 + i] = neighbor.id;
								break;
							}
						}
					}

					if (neighbor.data->data[0]) {
						for (int i = 0; i < 3; i++) {
							if (cell.data->data[2 + i] == 0) {
								cell.data->data[2 + i] = neighbor.id;
								break;
							}
						}
					}

				// consider only one sibling of all parents of neighboring cells...
				} else {

					// ignore own siblings
					if (grid.mapping.get_parent(cell.id) == grid.mapping.get_parent(neighbor.id)) {
						continue;
					}

					bool sibling_processed = false;
					uint64_t parent_of_neighbor = grid.mapping.get_parent(neighbor.id);
					for (int i = 0; i < 8; i++) {
						if (cell.data->data[5 + i] == parent_of_neighbor) {
							sibling_processed = true;
							break;
						}
					}

					if (sibling_processed) {
						continue;
					} else {
						for (int i = 0; i < 8; i++) {
							if (cell.data->data[5 + i] == 0) {
								cell.data->data[5 + i] = parent_of_neighbor;
								break;
							}
						}
					}

					// ...by recording which parents have been considered
					if (neighbor.data->data[0]) {
						for (int i = 0; i < 3; i++) {
							if (cell.data->data[2 + i] == 0) {
								cell.data->data[2 + i] = parent_of_neighbor;
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
	for (const auto& cell: grid.local_cells()) {
		if (grid.get_refinement_level(cell.id) == 0) {
			continue;
		}

		std::unordered_set<uint64_t> current_live_unrefined_neighbors;
		for (int i = 0; i < 3; i++) {
			current_live_unrefined_neighbors.insert(cell.data->data[2 + i]);
		}

		for (const auto& neighbor: cell.neighbors_of) {
			if (grid.get_refinement_level(neighbor.id) == 0) {
				continue;
			}

			// total live neighbors counts only between siblings
			if (grid.mapping.get_parent(cell.id) != grid.mapping.get_parent(neighbor.id)) {
				continue;
			}

			for (int i = 0; i < 3; i++) {
				current_live_unrefined_neighbors.insert(neighbor.data->data[2 + i]);
			}
		}

		current_live_unrefined_neighbors.erase(0);
		cell.data->data[1] += current_live_unrefined_neighbors.size();
	}

	// calculate the next turn
	for (const auto& cell: grid.local_cells()) {
		if (cell.data->data[1] == 3) {
			cell.data->data[0] = true;
		} else if (cell.data->data[1] != 2) {
			cell.data->data[0] = false;
		}
	}
}

#endif
