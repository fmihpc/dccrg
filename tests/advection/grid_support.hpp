/*
Support function for grid operations of advection test in dccrg.

Copyright 2012, 2013 Finnish Meteorological Institute

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

#ifndef DCCRG_GRID_SUPPORT_HPP
#define DCCRG_GRID_SUPPORT_HPP

#include "iostream"
#include "vector"

#include "dccrg.hpp"

namespace dccrg {

/*!
Directions used by get_neighbor_directions(...).
*/
enum direction_t {
	POS_X,
	NEG_X,
	POS_Y,
	NEG_Y,
	POS_Z,
	NEG_Z
};


/*!
Finds all face neighbors of given cell.

Also fills the neighbors' direction from given cell.
Given vectors are cleared before filling.

Assumes given dccrg instance was initialized with neighborhood
size of 0.
*/
template<class Cell_Data, class Geometry> void get_face_neighbors(
	const uint64_t cell,
	const dccrg::Dccrg<Cell_Data, Geometry>& grid,
	std::vector<uint64_t>& face_neighbors,
	std::vector<direction_t>& directions
) {
	face_neighbors.clear();
	directions.clear();

	face_neighbors.reserve(24);
	directions.reserve(24);

	const int refinement_level = grid.get_refinement_level(cell);

	const std::vector<uint64_t>* neighbors = grid.get_neighbors(cell);
	if (neighbors == NULL) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " No neighbors for cell " << cell
			<< std::endl;
		abort();
	}

	unsigned int neighbor_index = 0;

	// -z direction
	if ((*neighbors)[neighbor_index] != 0) {

		if (grid.get_refinement_level((*neighbors)[neighbor_index]) <= refinement_level) {

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Z);

		// only face neighbors
		} else {
			neighbor_index += 4;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Z);
			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Z);
			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Z);
			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Z);
		}
	}
	neighbor_index++;

	// -y direction
	if ((*neighbors)[neighbor_index] != 0) {

		if (grid.get_refinement_level((*neighbors)[neighbor_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Y);

		// solve only face neighbors
		} else {
			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Y);
			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Y);
			neighbor_index += 3;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Y);
			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_Y);
		}
	}
	neighbor_index++;

	// -x direction
	if ((*neighbors)[neighbor_index] != 0) {

		if (grid.get_refinement_level((*neighbors)[neighbor_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_X);

		// solve only face neighbors
		} else {
			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_X);
			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_X);
			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_X);
			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(NEG_X);
		}
	}
	neighbor_index++;

	// +x direction
	if ((*neighbors)[neighbor_index] != 0) {

		if (grid.get_refinement_level((*neighbors)[neighbor_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_X);

		// solve only face neighbors
		} else {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_X);

			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_X);

			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_X);

			neighbor_index += 2;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_X);

			neighbor_index++;
		}
	}
	neighbor_index++;

	// +y direction
	if ((*neighbors)[neighbor_index] != 0) {

		if (grid.get_refinement_level((*neighbors)[neighbor_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Y);

		// solve only face neighbors
		} else {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Y);

			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Y);

			neighbor_index += 3;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Y);

			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Y);

			neighbor_index += 2;
		}
	}
	neighbor_index++;

	// +z direction
	if ((*neighbors)[neighbor_index] != 0) {

		if (grid.get_refinement_level((*neighbors)[neighbor_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Z);

		// solve only face neighbors
		} else {
			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Z);

			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Z);

			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Z);

			neighbor_index++;

			face_neighbors.push_back((*neighbors)[neighbor_index]);
			directions.push_back(POS_Z);
		}
	}

	if (neighbor_index > neighbors->size()) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Added more neighbors (" << neighbor_index
			<< ") than exist: " << neighbors->size()
			<< std::endl;
		abort();
	}
}

/*
Removes from given face_neighbor cells those which are local and in negative direction.
*/
template<class Cell_Data, class Geometry> void remove_local_negative_neighbors(
	const dccrg::Dccrg<Cell_Data, Geometry>& grid,
	std::vector<uint64_t>& face_neighbors,
	std::vector<direction_t>& directions
) {
	std::vector<uint64_t> ret_neighbors;
	std::vector<direction_t> ret_directions;

	for (size_t i = 0; i < face_neighbors.size(); i++) {

		const uint64_t neighbor = face_neighbors[i];
		const direction_t direction = directions[i];

		if (
			   direction == POS_X
			|| direction == POS_Y
			|| direction == POS_Z
			|| !grid.is_local(neighbor)
		) {
			ret_neighbors.push_back(neighbor);
			ret_directions.push_back(direction);
		}
	}

	face_neighbors = ret_neighbors;
	directions = ret_directions;
}

} // namespace

#endif

