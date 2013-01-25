/*
A class for saving game of life tests of dccrg to vtk files.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute

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

#ifndef SAVE_HPP
#define SAVE_HPP

#include "algorithm"
#include "boost/foreach.hpp"
#include "iostream"
#include "stdint.h"
#include "vector"

#include "../../dccrg.hpp"

#include "cell.hpp"

/*!
Saves the current state of given game of life grid into one vtk file per process.

Visualize for example with VisIt (https://wci.llnl.gov/codes/visit/)
*/
template<class UserGeometry> class Save
{
public:

	static void save(
		const std::string& name,
		const int process,
		dccrg::Dccrg<Cell, UserGeometry>& game_grid
	) {
		std::vector<uint64_t> cells = game_grid.get_cells();
		sort(cells.begin(), cells.end());

		// write the grid into a file
		game_grid.write_vtk_file(name.c_str());
		// prepare to write the game data into the same file
		std::ofstream outfile(name.c_str(), std::ofstream::app);
		outfile << "CELL_DATA " << cells.size() << std::endl;

		// go through the grids cells and write their state into the file
		outfile << "SCALARS is_alive float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;
		BOOST_FOREACH(uint64_t cell, cells) {
			Cell* cell_data = game_grid[cell];
			outfile << cell_data->data[0] << std::endl;
		}

		// write each cells total live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;
		BOOST_FOREACH(uint64_t cell, cells) {
			Cell* cell_data = game_grid[cell];
			outfile << cell_data->data[1] << std::endl;
		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;
		BOOST_FOREACH(uint64_t cell, cells) {
			const std::vector<uint64_t>* neighbors = game_grid.get_neighbors(cell);
			outfile << neighbors->size() << std::endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;
		for (uint64_t i = 0; i < cells.size(); i++) {
			outfile << process << std::endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;
		BOOST_FOREACH(uint64_t cell, cells) {
			outfile << cell << std::endl;
		}
		outfile.close();
	}
};

#endif

