/*
A solver class for the game of life tests of dccrg.

Copyright 2010, 2011 Finnish Meteorological Institute

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

#ifndef SOLVE_HPP
#define SOLVE_HPP

#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "iostream"
#include "stdint.h"
#include "vector"

#include "../../dccrg.hpp"

#include "cell.hpp"

using namespace std;

/*!
Advances the game of life on given grid one turn.

Works only if cells within the grid have been refined <= 1 time.
*/
template<class UserGeometry> class Solve
{
public:

	static void solve(dccrg::Dccrg<Cell, UserGeometry>& game_grid)
	{
		// get the neighbour counts of every cell
		for (boost::unordered_map<uint64_t, Cell>::const_iterator
			cell_item = game_grid.begin();
			cell_item != game_grid.end();
			cell_item++
		) {
			const uint64_t cell = cell_item->first;

			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " no data for cell: " << cell
					<< std::endl;
				abort();
			}

			for (int i = 0; i < 12; i++) {
				cell_data->data[1 + i] = 0;
			}

			const vector<uint64_t>* neighbours = game_grid.get_neighbours(cell);
			// unrefined cells just consider neighbour counts at the level of unrefined cells
			if (game_grid.get_refinement_level(cell) == 0) {

				for (vector<uint64_t>::const_iterator
					neighbour = neighbours->begin();
					neighbour != neighbours->end();
					neighbour++
				) {
					if (*neighbour == 0) {
						continue;
					}

					Cell* neighbour_data = game_grid[*neighbour];
					if (neighbour_data == NULL) {
						cerr << __FILE__ << ":" << __LINE__
							<< " no data for neighbour of cell " << cell
							<< ": " << *neighbour
							<< std::endl;
						abort();
					}

					if (game_grid.get_refinement_level(*neighbour) == 0) {
						if (neighbour_data->data[0]) {
							cell_data->data[1]++;
						}
					// consider only one sibling...
					} else {

						bool sibling_processed = false;
						uint64_t parent_of_neighbour = game_grid.get_parent(*neighbour);
						for (int i = 0; i < 8; i++) {
							if (cell_data->data[5 + i] == parent_of_neighbour) {
								sibling_processed = true;
								break;
							}
						}

						// ...by recording its parent
						if (sibling_processed) {
							continue;
						} else {
							for (int i = 0; i < 8; i++) {
								if (cell_data->data[5 + i] == 0) {
									cell_data->data[5 + i] = parent_of_neighbour;
									break;
								}
							}
						}

						if (neighbour_data->data[0]) {
							cell_data->data[1]++;
						}
					}
				}

			// refined cells total the neighbour counts of siblings
			} else {

				for (vector<uint64_t>::const_iterator
					neighbour = neighbours->begin();
					neighbour != neighbours->end(); 
					neighbour++
				) {
					if (*neighbour == 0) {
						continue;
					}

					Cell* neighbour_data = game_grid[*neighbour];
					if (neighbour_data == NULL) {
						cerr << __FILE__ << ":" << __LINE__
							<< " no data for neighbour of refined cell " << cell
							<< ": " << *neighbour
							<< std::endl;
						abort();
					}

					if (game_grid.get_refinement_level(*neighbour) == 0) {

						// larger neighbours appear several times in the neighbour list
						bool neighbour_processed = false;
						for (int i = 0; i < 8; i++) {
							if (cell_data->data[5 + i] == *neighbour) {
								neighbour_processed = true;
								break;
							}
						}

						if (neighbour_processed) {
							continue;
						} else {
							for (int i = 0; i < 8; i++) {
								if (cell_data->data[5 + i] == 0) {
									cell_data->data[5 + i] = *neighbour;
									break;
								}
							}
						}

						if (neighbour_data->data[0]) {
							for (int i = 0; i < 3; i++) {
								if (cell_data->data[2 + i] == 0) {
									cell_data->data[2 + i] = *neighbour;
									break;
								}
							}
						}

					// consider only one sibling of all parents of neighbouring cells...
					} else {

						// ignore own siblings
						if (game_grid.get_parent(cell) == game_grid.get_parent(*neighbour)) {
							continue;
						}

						bool sibling_processed = false;
						uint64_t parent_of_neighbour = game_grid.get_parent(*neighbour);
						for (int i = 0; i < 8; i++) {
							if (cell_data->data[5 + i] == parent_of_neighbour) {
								sibling_processed = true;
								break;
							}
						}

						if (sibling_processed) {
							continue;
						} else {
							for (int i = 0; i < 8; i++) {
								if (cell_data->data[5 + i] == 0) {
									cell_data->data[5 + i] = parent_of_neighbour;
									break;
								}
							}
						}

						// ...by recording which parents have been considered
						if (neighbour_data->data[0]) {
							for (int i = 0; i < 3; i++) {
								if (cell_data->data[2 + i] == 0) {
									cell_data->data[2 + i] = parent_of_neighbour;
									break;
								}
							}
						}
					}
				}
			}

		}
		game_grid.update_remote_neighbour_data();

		// get the total neighbour counts of refined cells
		for (boost::unordered_map<uint64_t, Cell>::const_iterator
			cell_item = game_grid.begin();
			cell_item != game_grid.end();
			cell_item++
		) {
			const uint64_t cell = cell_item->first;

			if (game_grid.get_refinement_level(cell) == 0) {
				continue;
			}
			Cell* cell_data = game_grid[cell];

			boost::unordered_set<uint64_t> current_live_unrefined_neighbours;
			for (int i = 0; i < 3; i++) {
				current_live_unrefined_neighbours.insert(cell_data->data[2 + i]);
			}

			const vector<uint64_t>* neighbours = game_grid.get_neighbours(cell);
			for (vector<uint64_t>::const_iterator
				neighbour = neighbours->begin();
				neighbour != neighbours->end();
				neighbour++
			) {
				if (*neighbour == 0) {
					continue;
				}

				if (game_grid.get_refinement_level(*neighbour) == 0) {
					continue;
				}

				// total live neighbours counts only between siblings
				if (game_grid.get_parent(cell) != game_grid.get_parent(*neighbour)) {
					continue;
				}

				Cell* neighbour_data = game_grid[*neighbour];
				for (int i = 0; i < 3; i++) {
					current_live_unrefined_neighbours.insert(neighbour_data->data[2 + i]);
				}
			}

			current_live_unrefined_neighbours.erase(0);
			cell_data->data[1] += current_live_unrefined_neighbours.size();
		}

		// calculate the next turn
		for (boost::unordered_map<uint64_t, Cell>::const_iterator
			cell_item = game_grid.begin();
			cell_item != game_grid.end();
			cell_item++
		) {
			const uint64_t cell = cell_item->first;

			Cell* cell_data = game_grid[cell];

			if (cell_data->data[1] == 3) {
				cell_data->data[0] = true;
			} else if (cell_data->data[1] != 2) {
				cell_data->data[0] = false;
			}
		}
	}
};

#endif

