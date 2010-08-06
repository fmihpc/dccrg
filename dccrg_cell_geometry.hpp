/*
The mapping logic between a cell's id and its geometry (location, size and comparable stuff)

Copyright 2009, 2010 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_CELL_GEOMETRY_HPP
#define DCCRG_CELL_GEOMETRY_HPP


#ifdef DCCRG_ARBITRARY_STRETCH
	#ifdef DCCRG_CONSTANT_STRETCH
		#error Only one type of grid stretching can be used simultaneously
	#endif
#endif


#include "iostream"
#include "stdint.h"
#include "cassert"
#include "cstdlib"
#include "vector"


class CellGeometry
{

public:

	#ifdef DCCRG_ARBITRARY_STRETCH

	/*
	x, y and z_coordinates:
		The coordinates of unrefined cells in the respective direction
		First coordinate is the starting point of the grid, the following ith value is the endpoint of the ith unrefined cell
	*/
	CellGeometry(const std::vector<double> x_coordinates, const std::vector<double> y_coordinates, const std::vector<double> z_coordinates, const int maximum_refinement_level);

	#else

	/*
	x_start, y_start, z_start:
		the starting corner of the grid
	cell_size:
		the size of each unrefined cell in every direction
	x_length, y_length, z_length:
		the number of cells in the grid in x, y and z direction
	*/
	CellGeometry(const double x_start, const double y_start, const double z_start, const double cell_x_size, const double cell_y_size, const double cell_z_size, const uint64_t x_length, const uint64_t y_length, const uint64_t z_length, const int maximum_refinement_level);

	#endif


	/*
	Returns the refinement level of given cell (0 means unrefined)
	Returns -1 if given cell cannot exist in the current grid
	*/
	int get_refinement_level(uint64_t cell);


	/*
	The following return the length of given cell in x, y or z direction
	*/
	double get_cell_x_size(const uint64_t cell);
	double get_cell_y_size(const uint64_t cell);
	double get_cell_z_size(const uint64_t cell);


	/*
	The following return the x, y or z coordinate of the center of given cell regardless of whether it exists or has children
	*/
	double get_cell_x(const uint64_t cell);
	double get_cell_y(const uint64_t cell);
	double get_cell_z(const uint64_t cell);


	/*
	Returns the maximum refinement level of any cell in the current grid
	*/
	int get_maximum_refinement_level(void) { return this->maximum_refinement_level; }

private:

	#ifdef DCCRG_ARBITRARY_STRETCH
	/*
	The coordinates of unrefined cells in respective directions
	First value is the starting point of the grid, the following ith value is the end point of the ith unrefined cell
	*/
	std::vector<double> x_coordinates, y_coordinates, z_coordinates;
	#else
	// starting corner coordinates of the grid
	double x_start, y_start, z_start;
	// length of unrefined cells in all directions
	double cell_x_size, cell_y_size, cell_z_size;
	#endif
	// size of the grid in unrefined cells
	uint64_t x_length, y_length, z_length;
	// maximum refinemet level of any cell in the grid, 0 means unrefined
	int maximum_refinement_level;


	/*
	These return the index of the cell with given id in x, y or z direction of the grid, starting from 0.
	For cells that are larger than the smallest possible according to maximum_refinement_level, the index within the cell with minimum x, y and z value is returned

	Example with maximum refinement level 1
	index    0 1 2
		-------
	      0 | |   |
		---   |
	      1 | |   |
		-------
	 */
	uint64_t get_x_index(uint64_t cell);
	uint64_t get_y_index(uint64_t cell);
	uint64_t get_z_index(uint64_t cell);


	/*
	Returns the lengths of given cell in indices in every direction
	*/
	uint64_t get_cell_size_in_indices(const uint64_t cell);


	/*
	Returns the cell of given refinement level at given indices
	*/
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int refinement_level);


	/*
	These return the index of the starting coordinate in the coordinates vector of the given unrefined cell in the x, y and z direction respectively
	*/
	uint64_t get_unref_cell_x_coord_start_index(const uint64_t cell);
	uint64_t get_unref_cell_y_coord_start_index(const uint64_t cell);
	uint64_t get_unref_cell_z_coord_start_index(const uint64_t cell);

};


#ifdef DCCRG_ARBITRARY_STRETCH
CellGeometry::CellGeometry(const std::vector<double> x_coordinates, const std::vector<double> y_coordinates, const std::vector<double> z_coordinates, const int maximum_refinement_level = -1)
{
	if (x_coordinates.size() < 2) {
		std::cerr << "At least two coordinates are required for grid cells in the x direction" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (y_coordinates.size() < 2) {
		std::cerr << "At least two coordinates are required for grid cells in the y direction" << std::endl;
		exit(EXIT_FAILURE);
	}
	if (z_coordinates.size() < 2) {
		std::cerr << "At least two coordinates are required for grid cells in the z direction" << std::endl;
		exit(EXIT_FAILURE);
	}
	for (uint64_t i = 0; i < x_coordinates.size() - 1; i++) {
		if (x_coordinates[i] >= x_coordinates[i + 1]) {
			std::cerr << "Coordinates in the x direction must be strictly increasing" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	for (uint64_t i = 0; i < y_coordinates.size() - 1; i++) {
		if (y_coordinates[i] >= y_coordinates[i + 1]) {
			std::cerr << "Coordinates in the y direction must be strictly increasing" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	for (uint64_t i = 0; i < z_coordinates.size() - 1; i++) {
		if (z_coordinates[i] >= z_coordinates[i + 1]) {
			std::cerr << "Coordinates in the z direction must be strictly increasing" << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	this->x_coordinates.reserve(x_coordinates.size());
	this->x_coordinates.insert(this->x_coordinates.begin(), x_coordinates.begin(), x_coordinates.end());
	this->y_coordinates.reserve(y_coordinates.size());
	this->y_coordinates.insert(this->y_coordinates.begin(), y_coordinates.begin(), y_coordinates.end());
	this->z_coordinates.reserve(z_coordinates.size());
	this->z_coordinates.insert(this->z_coordinates.begin(), z_coordinates.begin(), z_coordinates.end());

	this->x_length = this->x_coordinates.size() - 1;
	this->y_length = this->y_coordinates.size() - 1;
	this->z_length = this->z_coordinates.size() - 1;


#else

CellGeometry::CellGeometry(const double x_start, const double y_start, const double z_start, const double cell_x_size, const double cell_y_size, const double cell_z_size, const uint64_t x_length, const uint64_t y_length, const uint64_t z_length, const int maximum_refinement_level = -1)
{
	this->x_start = x_start;
	this->y_start = y_start;
	this->z_start = z_start;


	if (cell_x_size <= 0) {
		std::cerr << "Cell size in x direction must be > 0" << std::endl;
		exit(EXIT_FAILURE);
	}
	this->cell_x_size = cell_x_size;

	if (cell_y_size <= 0) {
		std::cerr << "Cell size in y direction must be > 0" << std::endl;
		exit(EXIT_FAILURE);
	}
	this->cell_y_size = cell_y_size;

	if (cell_z_size <= 0) {
		std::cerr << "Cell size in z direction must be > 0" << std::endl;
		exit(EXIT_FAILURE);
	}
	this->cell_z_size = cell_z_size;


	if (x_length == 0) {
		std::cerr << "Length of the grid in unrefined cells must be > 0 in the x direction" << std::endl;
		exit(EXIT_FAILURE);
	}
	this->x_length = x_length;

	if (y_length == 0) {
		std::cerr << "Length of the grid in unrefined cells must be > 0 in the y direction" << std::endl;
		exit(EXIT_FAILURE);
	}
	this->y_length = y_length;

	if (z_length == 0) {
		std::cerr << "Length of the grid in unrefined cells must be > 0 in the z direction" << std::endl;
		exit(EXIT_FAILURE);
	}
	this->z_length = z_length;


#endif

	// get the maximum refinement level based on the size of the grid when using uint64_t for cell ids
	double max_id = uint64_t(~0), last_id = x_length * y_length * z_length;
	int refinement_level = 0;
	while (last_id / max_id < 1) {
		refinement_level++;
		last_id += double(x_length) * y_length * z_length * (uint64_t(1) << refinement_level * 3);
	}
	refinement_level--;

	// grid is too large even without refinement
	if (refinement_level < 0) {
		std::cerr << "Given grid would contain more than 2^64 - 1 unrefined cells" << std::endl;
		exit(EXIT_FAILURE);
	}


	if (maximum_refinement_level > refinement_level) {

		std::cerr << "Given maximum refinement level (" << maximum_refinement_level << ") is too large: last cell id would be " << double(x_length) * y_length * z_length * (uint64_t(1) << maximum_refinement_level * 3) << " but maximum for uint64_t is " << double(uint64_t(~0)) << std::endl;
		exit(EXIT_FAILURE);

	} else if (maximum_refinement_level < 0) {
		this->maximum_refinement_level = refinement_level;
	} else {
		this->maximum_refinement_level = maximum_refinement_level;
	}
}


int CellGeometry::get_refinement_level(uint64_t cell)
{
	if (cell == 0) {
		return -1;
	}

	int refinement_level = 0;
	uint64_t last_cell = this->x_length * this->y_length * this->z_length;

	while (last_cell < cell) {
		refinement_level++;
		last_cell += this->x_length * this->y_length * this->z_length * (uint64_t(1) << 3 * refinement_level);
	}

	if (refinement_level > this->maximum_refinement_level) {
		return -1;
	}

	return refinement_level;
}


uint64_t CellGeometry::get_x_index(uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	// substract larger cells
	int refinement_level = this->get_refinement_level(cell);
	for (int i = 0; i < refinement_level; i++) {
		cell -= this->x_length * this->y_length * this->z_length * (uint64_t(1) << i * 3);
	}

	// get the index at this cells refinement level
	cell -= 1;	// cell numbering starts from 1
	uint64_t this_level_index = cell % (this->x_length * (uint64_t(1) << refinement_level));

	assert(this_level_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level)) < this->x_length * (uint64_t(1) << this->maximum_refinement_level));
	return this_level_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level));
}

uint64_t CellGeometry::get_y_index(uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	// substract larger cells
	int refinement_level = this->get_refinement_level(cell);
	for (int i = 0; i < refinement_level; i++) {
		cell -= this->x_length * this->y_length * this->z_length * (uint64_t(1) << i * 3);
	}

	// get the index at this cells refinement level
	cell -= 1;	// cell numbering starts from 1
	uint64_t this_level_index =  (cell / (this->x_length * (uint64_t(1) << refinement_level))) % (this->y_length  * (uint64_t(1) << refinement_level));

	return this_level_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level));
}

uint64_t CellGeometry::get_z_index(uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	// substract larger cells
	int refinement_level = this->get_refinement_level(cell);
	for (int i = 0; i < refinement_level; i++) {
		cell -= this->x_length * this->y_length * this->z_length * (uint64_t(1) << i * 3);
	}

	// get the index at this cells refinement level
	cell -= 1;	// cell numbering starts from 1
	uint64_t this_level_index =  cell / (this->x_length * this->y_length * (uint64_t(1) << 2 * refinement_level));

	return this_level_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level));
}


uint64_t CellGeometry::get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int refinement_level)
{
	assert(refinement_level >= 0);
	assert(refinement_level <= this->maximum_refinement_level);
	assert(x_index < this->x_length * (uint64_t(1) << this->maximum_refinement_level));
	assert(y_index < this->y_length * (uint64_t(1) << this->maximum_refinement_level));
	assert(z_index < this->z_length * (uint64_t(1) << this->maximum_refinement_level));

	uint64_t cell = 1;

	// add larger cells
	for (int i = 0; i < refinement_level; i++) {
		cell += this->x_length * this->y_length * this->z_length * (uint64_t(1) << 3 * i);
	}

	// convert to indices of this cells refinement level
	uint64_t this_level_x_index = x_index / (uint64_t(1) << (maximum_refinement_level - refinement_level));
	uint64_t this_level_y_index = y_index / (uint64_t(1) << (maximum_refinement_level - refinement_level));
	uint64_t this_level_z_index = z_index / (uint64_t(1) << (maximum_refinement_level - refinement_level));

	// get the size of the grid in cells of this refinement level
	uint64_t this_level_x_length = x_length * (uint64_t(1) << refinement_level);
	uint64_t this_level_y_length = y_length * (uint64_t(1) << refinement_level);

	cell += this_level_x_index + this_level_y_index * this_level_x_length + this_level_z_index * this_level_x_length * this_level_y_length;

	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);
	return cell;
}


uint64_t CellGeometry::get_cell_size_in_indices(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	return uint64_t(1) << (this->maximum_refinement_level - this->get_refinement_level(cell));
}


double CellGeometry::get_cell_x_size(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

#ifdef DCCRG_ARBITRARY_STRETCH
	int refinement_level = this->get_refinement_level(cell);
	uint64_t unref_cell_x_index = this->get_unref_cell_x_coord_start_index(cell);
	return (this->x_coordinates[unref_cell_x_index + 1] - this->x_coordinates[unref_cell_x_index]) / (uint64_t(1) << refinement_level);
}
#else
	return this->cell_x_size / (uint64_t(1) << this->get_refinement_level(cell));
}
#endif

double CellGeometry::get_cell_y_size(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

#ifdef DCCRG_ARBITRARY_STRETCH
	int refinement_level = this->get_refinement_level(cell);
	uint64_t unref_cell_y_index = this->get_unref_cell_y_coord_start_index(cell);
	return (this->y_coordinates[unref_cell_y_index + 1] - this->y_coordinates[unref_cell_y_index]) / (uint64_t(1) << refinement_level);
}
#else
	return this->cell_y_size / (uint64_t(1) << this->get_refinement_level(cell));
}
#endif

double CellGeometry::get_cell_z_size(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

#ifdef DCCRG_ARBITRARY_STRETCH
	int refinement_level = this->get_refinement_level(cell);
	uint64_t unref_cell_z_index = this->get_unref_cell_z_coord_start_index(cell);
	return (this->z_coordinates[unref_cell_z_index + 1] - this->z_coordinates[unref_cell_z_index]) / (uint64_t(1) << refinement_level);
}
#else
	return this->cell_z_size / (uint64_t(1) << this->get_refinement_level(cell));
}
#endif


double CellGeometry::get_cell_x(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

#ifdef DCCRG_ARBITRARY_STRETCH
	uint64_t unref_cell_x_coord_start_index = this->get_unref_cell_x_coord_start_index(cell);
	uint64_t unref_cell_x_index = unref_cell_x_coord_start_index * (uint64_t(1) << this->maximum_refinement_level);

	double unref_cell_x_size = this->x_coordinates[unref_cell_x_coord_start_index + 1] - this->x_coordinates[unref_cell_x_coord_start_index];
	double size_of_local_index = unref_cell_x_size / (uint64_t(1) << this->maximum_refinement_level);

	return this->x_coordinates[unref_cell_x_coord_start_index] + size_of_local_index * (this->get_x_index(cell) - unref_cell_x_index) + this->get_cell_x_size(cell) / 2;
}
#else
	return this->x_start + this->get_x_index(cell) * this->cell_x_size / (uint64_t(1) << this->maximum_refinement_level) + this->get_cell_x_size(cell) / 2;
}
#endif

double CellGeometry::get_cell_y(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

#ifdef DCCRG_ARBITRARY_STRETCH
	uint64_t unref_cell_y_coord_start_index = this->get_unref_cell_y_coord_start_index(cell);
	uint64_t unref_cell_y_index = unref_cell_y_coord_start_index * (uint64_t(1) << this->maximum_refinement_level);

	double unref_cell_y_size = this->y_coordinates[unref_cell_y_coord_start_index + 1] - this->y_coordinates[unref_cell_y_coord_start_index];
	double size_of_local_index = unref_cell_y_size / (uint64_t(1) << this->maximum_refinement_level);

	return this->y_coordinates[unref_cell_y_coord_start_index] + size_of_local_index * (this->get_y_index(cell) - unref_cell_y_index) + this->get_cell_y_size(cell) / 2;
}
#else
	return this->y_start + this->get_y_index(cell) * this->cell_y_size / (uint64_t(1) << this->maximum_refinement_level) + this->get_cell_y_size(cell) / 2;
}
#endif

double CellGeometry::get_cell_z(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

#ifdef DCCRG_ARBITRARY_STRETCH
	uint64_t unref_cell_z_coord_start_index = this->get_unref_cell_z_coord_start_index(cell);
	uint64_t unref_cell_z_index = unref_cell_z_coord_start_index * (uint64_t(1) << this->maximum_refinement_level);

	double unref_cell_z_size = this->z_coordinates[unref_cell_z_coord_start_index + 1] - this->z_coordinates[unref_cell_z_coord_start_index];
	double size_of_local_index = unref_cell_z_size / (uint64_t(1) << this->maximum_refinement_level);

	return this->z_coordinates[unref_cell_z_coord_start_index] + size_of_local_index * (this->get_z_index(cell) - unref_cell_z_index) + this->get_cell_z_size(cell) / 2;
}
#else
	return this->z_start + this->get_z_index(cell) * this->cell_z_size / (uint64_t(1) << this->maximum_refinement_level) + this->get_cell_z_size(cell) / 2;
}
#endif


#ifdef DCCRG_ARBITRARY_STRETCH
uint64_t CellGeometry::get_unref_cell_x_coord_start_index(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	return this->get_x_index(cell) / (uint64_t(1) << this->maximum_refinement_level);
}

uint64_t CellGeometry::get_unref_cell_y_coord_start_index(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	return this->get_y_index(cell) / (uint64_t(1) << this->maximum_refinement_level);
}

uint64_t CellGeometry::get_unref_cell_z_coord_start_index(const uint64_t cell)
{
	assert(cell > 0);
	assert(this->get_refinement_level(cell) >= 0);
	assert(this->get_refinement_level(cell) <= this->maximum_refinement_level);

	return this->get_z_index(cell) / (uint64_t(1) << this->maximum_refinement_level);
}
#endif


#endif
