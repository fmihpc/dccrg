/*
A class for initializing the game of life tests of dccrg.

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

#ifndef INITIALIZE_HPP
#define INITIALIZE_HPP

#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "stdint.h"
#include "vector"

#include "../../dccrg.hpp"

#include "cell.hpp"

/*!
Initializes given grid with given size in the first game coordinate with some game of life patterns.
*/
template<class UserGeometry> class Initialize
{
public:

	static void initialize(dccrg::Dccrg<Cell, UserGeometry>& game_grid, const uint64_t grid_size)
	{
		if (grid_size == 0) {
			assert(false);
		}

		// create a blinker
		const uint64_t blinker_start = 198;
		const std::vector<uint64_t> blinker_cells = boost::assign::list_of
			(blinker_start)
			(blinker_start + 1)
			(blinker_start + 2);
		BOOST_FOREACH(uint64_t cell, blinker_cells) {
			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				continue;
			}
			cell_data->data[0] = 1;
		}

		// create a toad
		const uint64_t toad_start = 188;
		const std::vector<uint64_t> toad_cells = boost::assign::list_of
			(toad_start)
			(toad_start + 1)
			(toad_start + 2)
			(toad_start + 1 + grid_size)
			(toad_start + 2 + grid_size)
			(toad_start + 3 + grid_size);
		BOOST_FOREACH(uint64_t cell, toad_cells) {
			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				continue;
			}
			cell_data->data[0] = 1;
		}

		// create a beacon
		const uint64_t beacon_start = 137;
		const std::vector<uint64_t> beacon_cells = boost::assign::list_of
			(beacon_start)
			(beacon_start + 1)
			(beacon_start - grid_size)
			(beacon_start + 1 - grid_size)
			(beacon_start + 2 - 2 * grid_size)
			(beacon_start + 3 - 2 * grid_size)
			(beacon_start + 2 - 3 * grid_size)
			(beacon_start + 3 - 3 * grid_size);
		BOOST_FOREACH(uint64_t cell, beacon_cells) {
			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				continue;
			}
			cell_data->data[0] = 1;
		}

		// create a glider
		const uint64_t glider_start = 143;
		const std::vector<uint64_t> glider_cells = boost::assign::list_of
			(glider_start + 1)
			(glider_start + 2 - grid_size)
			(glider_start - 2 * grid_size)
			(glider_start + 1 - 2 * grid_size)
			(glider_start + 2 - 2 * grid_size);
		BOOST_FOREACH(uint64_t cell, glider_cells) {
			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				continue;
			}
			cell_data->data[0] = 1;
		}

		// create a block
		const uint64_t block_start = 47;
		const std::vector<uint64_t> block_cells = boost::assign::list_of
			(block_start)
			(block_start + 1)
			(block_start - grid_size)
			(block_start + 1 - grid_size);
		BOOST_FOREACH(uint64_t cell, block_cells) {
			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				continue;
			}
			cell_data->data[0] = 1;
		}

		// create a beehive
		const uint64_t beehive_start = 51;
		const std::vector<uint64_t> beehive_cells = boost::assign::list_of
			(beehive_start - grid_size)
			(beehive_start + 1)
			(beehive_start + 2)
			(beehive_start + 1 - 2 * grid_size)
			(beehive_start + 2 - 2 * grid_size)
			(beehive_start + 3 - grid_size);
		BOOST_FOREACH(uint64_t cell, beehive_cells) {
			Cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				continue;
			}
			cell_data->data[0] = 1;
		}
	}

};

#endif

