/*
A class for initializing the game of life tests of dccrg.

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

#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include "cstdint"
#include "vector"

#include "dccrg.hpp"


std::set<uint64_t> get_live_cells(const uint64_t grid_size) {
	std::set<uint64_t> live_cells;

	// create a blinker
	const uint64_t blinker_start = 198;
	live_cells.insert({
		blinker_start,
		blinker_start + 1,
		blinker_start + 2
	});

	// create a toad
	const uint64_t toad_start = 188;
	live_cells.insert({
		toad_start,
		toad_start + 1,
		toad_start + 2,
		toad_start + 1 + grid_size,
		toad_start + 2 + grid_size,
		toad_start + 3 + grid_size
	});

	// create a beacon
	const uint64_t beacon_start = 137;
	live_cells.insert({
		beacon_start,
		beacon_start + 1,
		beacon_start - grid_size,
		beacon_start + 1 - grid_size,
		beacon_start + 2 - 2 * grid_size,
		beacon_start + 3 - 2 * grid_size,
		beacon_start + 2 - 3 * grid_size,
		beacon_start + 3 - 3 * grid_size
	});

	// create a glider
	const uint64_t glider_start = 143;
	live_cells.insert({
		glider_start + 1,
		glider_start + 2 - grid_size,
		glider_start - 2 * grid_size,
		glider_start + 1 - 2 * grid_size,
		glider_start + 2 - 2 * grid_size
	});

	// create a block
	const uint64_t block_start = 47;
	live_cells.insert({
		block_start,
		block_start + 1,
		block_start - grid_size,
		block_start + 1 - grid_size
	});

	// create a beehive
	const uint64_t beehive_start = 51;
	live_cells.insert({
		beehive_start - grid_size,
		beehive_start + 1,
		beehive_start + 2,
		beehive_start + 1 - 2 * grid_size,
		beehive_start + 2 - 2 * grid_size,
		beehive_start + 3 - grid_size
	});
	return live_cells;
}

/*!
Initializes given grid with given size in the first game coordinate with some game of life patterns.
*/
template<class Cell_Data, class Geometry> void initialize(
	dccrg::Dccrg<Cell_Data, Geometry>& grid,
	const uint64_t grid_size
) {
	if (grid_size == 0) {
		abort();
	}

	const auto live_cells = get_live_cells(grid_size);

	for (const auto& cell: grid.local_cells) {
		for (auto& item: cell.data->data) {
			item = 0;
		}
		if (live_cells.count(cell.id) > 0) {
			cell.data->data[0] = 1;
		}
	}
}

#endif
