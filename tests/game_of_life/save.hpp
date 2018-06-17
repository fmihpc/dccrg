/*
A class for saving game of life tests of dccrg to vtk files.

Copyright 2010, 2011, 2012, 2013, 2014, 2015, 2016 Finnish Meteorological Institute

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

#ifndef SAVE_HPP
#define SAVE_HPP

#include "algorithm"
#include "iostream"
#include "iterator"
#include "cstdint"

#include "dccrg.hpp"

#include "cell.hpp"

/*!
Saves the current state of given game of life grid into one vtk file per process.

Visualize for example with VisIt (https://wci.llnl.gov/codes/visit/)
*/
template<class Cell_Data, class Geometry> void save(
	const std::string& name,
	const int process,
	dccrg::Dccrg<Cell_Data, Geometry>& grid
) {
	auto cells = grid.cells;
	std::sort(cells.begin(), cells.end());

	// write the grid into a file
	grid.write_vtk_file(name.c_str());
	// prepare to write the game data into the same file
	std::ofstream outfile(name.c_str(), std::ofstream::app);
	outfile << "CELL_DATA " << cells.size() << std::endl;

	// go through the grids cells and write their state into the file
	outfile << "SCALARS is_alive float 1" << std::endl;
	outfile << "LOOKUP_TABLE default" << std::endl;
	for (const auto& cell: cells) {
		outfile << cell.data->data[0] << std::endl;
	}

	// write each cells total live neighbor count
	outfile << "SCALARS live_neighbor_count float 1" << std::endl;
	outfile << "LOOKUP_TABLE default" << std::endl;
	for (const auto& cell: cells) {
		outfile << cell.data->data[1] << std::endl;
	}

	// write each cells neighbor count
	outfile << "SCALARS neighbors int 1" << std::endl;
	outfile << "LOOKUP_TABLE default" << std::endl;
	for (const auto& cell: cells) {
		const auto* const neighbors = grid.get_neighbors_of(cell.id);
		std::set<uint64_t> neighbors1, neighbors2;
		for (const auto& neighbor: *neighbors) neighbors2.insert(neighbor.first);
		neighbors2.erase(dccrg::error_cell);
		for (const auto& neighbor: cell.neighbors_of) {
			neighbors1.insert(neighbor.id);
		}
		if (neighbors1.size() != neighbors2.size()) {
			std::cerr << __FILE__ "(" << __LINE__ << "): "
				<< neighbors1.size() << " != " << neighbors2.size() << std::endl;
		}
		if (not std::equal(neighbors1.cbegin(), neighbors1.cend(), neighbors2.cbegin())) {
			std::cerr << __FILE__ "(" << __LINE__ << ")" << std::endl;
		}
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
	for (const auto& cell: cells) {
		outfile << cell.id << std::endl;
	}
	outfile.close();
}

#endif
