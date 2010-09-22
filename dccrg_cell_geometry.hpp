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


#include "cassert"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "stdint.h"
#include "vector"


class CellGeometry
{

public:

	/*!
	Creates a default instance of the grid, with one cell of size 1 in every direction between coordinates (0, 0, 0) and (1, 1, 1), with maximum refinement level of 0.
	Use the set_* functions below to change the defaults.
	*/
	CellGeometry();

	/*!
	Returns grid parameters to their default values.
	*/
	~CellGeometry();

	#ifdef DCCRG_ARBITRARY_STRETCH

	/*!
	Sets the coordinates of unrefined cells in the grid.
	First coordinate is the starting point of the grid and the following ith value is the endpoint of the ith unrefined cell
	At least two values must be given for each direction and all values must be strictly increasing.
	Returns true if successful, probably invalidating all previous cell geometries, e.g. their location, size, etc.
	Returns false if unsuccessful and in that case has no effect.
	Automatically maximizes maximum_refinement_level.
	*/
	bool set_coordinates(const std::vector<double> x_coordinates, const std::vector<double> y_coordinates, const std::vector<double> z_coordinates);

	#else

	/*!
	The following set the starting corner of the grid, e.g. the first face of the first unrefined cell(s) in that direction.
	At the moment these always return true.
	Automatically maximizes maximum_refinement_level.
	*/
	// FIXME: check that all cell coordinates fit into a double, and return false if some cells would be out of range
	bool set_x_start(const double x_start);
	bool set_y_start(const double y_start);
	bool set_z_start(const double z_start);

	/*!
	The following return the starting corner of the grid
	*/
	double get_x_start(void);
	double get_y_start(void);
	double get_z_start(void);


	/*!
	The following set the size of unrefined cells.
	These return true if successful, probably invalidating all previous cell areas, volumes, etc.
	Return false if unsuccessful and in that case have no effect.
	Automatically maximizes maximum_refinement_level.
	*/
	// FIXME: check that all cell coordinates fit into a double, and return false if some cells would be out of range
	bool set_cell_x_size(const double cell_x_size);
	bool set_cell_y_size(const double cell_y_size);
	bool set_cell_z_size(const double cell_z_size);

	/*!
	The following set the size of the grid in unrefined cells.
	These return true if successful, probably invalidating all previous cell information (id, geometry, etc.)
	Return false if unsuccessful and in that case have no effect.
	Automatically maximizes maximum_refinement_level.
	*/
	bool set_x_length(const uint64_t x_length);
	bool set_y_length(const uint64_t y_length);
	bool set_z_length(const uint64_t z_length);

	#endif


	/*!
	The following return the length of the grid in unrefined cells
	*/
	uint64_t get_x_length(void);
	uint64_t get_y_length(void);
	uint64_t get_z_length(void);


	/*!
	Sets the maximum refinement level of the grid (0 means unrefined), probably invalidating all previous cell indices.
	Returns true if all cells in the grid with the given maximum refinement level would fit into an uint64_t, otherwise does nothing and returns false.
	*/
	bool set_maximum_refinement_level(const int maximum_refinement_level);

	/*!
	Returns the maximum refinement level of any cell in the grid (0 means unrefined).
	*/
	int get_maximum_refinement_level(void);

	/*!
	Returns the maximum possible refinement level for a cell in the grid (0 means unrefined).
	*/
	int get_maximum_possible_refinement_level(void);


	/*!
	Returns the refinement level of given cell (0 means unrefined)
	Returns -1 if given cell cannot exist in the current grid
	*/
	int get_refinement_level(uint64_t cell);


	/*!
	The following return the length of given cell.
	*/
	double get_cell_x_size(const uint64_t cell);
	double get_cell_y_size(const uint64_t cell);
	double get_cell_z_size(const uint64_t cell);


	/*!
	The following return the centroid of given cell.
	*/
	double get_cell_x(const uint64_t cell);
	double get_cell_y(const uint64_t cell);
	double get_cell_z(const uint64_t cell);

	/*!
	The following return the centroid of a cell of given refinement level at given index
	*/
	double get_cell_x(const int refinement_level, const uint64_t x_index);
	double get_cell_y(const int refinement_level, const uint64_t y_index);
	double get_cell_z(const int refinement_level, const uint64_t z_index);


	/*!
	Returns the cell at given location, or 0 if outside of the current grid
	*/
	uint64_t get_cell(const int refinement_level, const double x, const double y, const double z);


	/*!
	These return the index of given cell, starting from 0.
	For cells that are larger than the smallest possible according to maximum_refinement_level, the index closest to the grid starting point is returned.

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

	/*!
	These return the index (starting from 0) of given location.
	TODO: check that location is within the current grid
	*/
	uint64_t get_x_index(const double x);
	uint64_t get_y_index(const double y);
	uint64_t get_z_index(const double z);



private:

	/*!
	Initializes the grid parameters to the values specified in the constructor comments.
	*/
	void init(void);


	#ifdef DCCRG_ARBITRARY_STRETCH
	/*!
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

	// BOOST serialization support
	#ifdef BOOST_VERSION
	#include "boost/serialization/serialization.hpp"
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		#ifdef DCCRG_ARBITRARY_STRETCH
		ar & x_coordinates & y_coordinates & z_coordinates;
		#else
		ar & x_start & y_start & z_start;
		ar & cell_x_size & cell_y_size & cell_z_size;
		#endif
		ar & x_length & y_length & z_length & maximum_refinement_level;
	}
	#endif


	/*!
	Returns the lengths of given cell in indices in every direction
	*/
	uint64_t get_cell_size_in_indices(const uint64_t cell);


	/*!
	Returns the cell of given refinement level at given indices
	*/
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int refinement_level);


	/*!
	These return the index in the coordinates vector of the starting coordinate of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell
	*/
	uint64_t get_unref_cell_x_coord_start_index(const uint64_t cell);
	uint64_t get_unref_cell_y_coord_start_index(const uint64_t cell);
	uint64_t get_unref_cell_z_coord_start_index(const uint64_t cell);

};


CellGeometry::CellGeometry()
{
	this->init();
}

CellGeometry::~CellGeometry()
{
	this->init();
}

void CellGeometry::init(void)
{
	#ifdef DCCRG_ARBITRARY_STRETCH

	this->x_coordinates.clear();
	this->x_coordinates.push_back(0);
	this->x_coordinates.push_back(1);

	this->y_coordinates.clear();
	this->y_coordinates.push_back(0);
	this->y_coordinates.push_back(1);

	this->z_coordinates.clear();
	this->z_coordinates.push_back(0);
	this->z_coordinates.push_back(1);

	#else

	this->x_start = 0;
	this->y_start = 0;
	this->z_start = 0;

	this->cell_x_size = 1;
	this->cell_y_size = 1;
	this->cell_z_size = 1;

	#endif

	this->x_length = 1;
	this->y_length = 1;
	this->z_length = 1;

	this->maximum_refinement_level = 0;
}


#ifdef DCCRG_ARBITRARY_STRETCH

bool CellGeometry::set_coordinates(const std::vector<double> x_coordinates, const std::vector<double> y_coordinates, const std::vector<double> z_coordinates)
{
	if (x_coordinates.size() < 2) {
		std::cerr << "At least two coordinates are required for grid cells in the x direction" << std::endl;
		return false;
	}
	if (y_coordinates.size() < 2) {
		std::cerr << "At least two coordinates are required for grid cells in the y direction" << std::endl;
		return false;
	}
	if (z_coordinates.size() < 2) {
		std::cerr << "At least two coordinates are required for grid cells in the z direction" << std::endl;
		return false;
	}

	for (uint64_t i = 0; i < x_coordinates.size() - 1; i++) {
		if (x_coordinates[i] >= x_coordinates[i + 1]) {
			std::cerr << "Coordinates in the x direction must be strictly increasing" << std::endl;
			return false;
		}
	}
	for (uint64_t i = 0; i < y_coordinates.size() - 1; i++) {
		if (y_coordinates[i] >= y_coordinates[i + 1]) {
			std::cerr << "Coordinates in the y direction must be strictly increasing" << std::endl;
			return false;
		}
	}
	for (uint64_t i = 0; i < z_coordinates.size() - 1; i++) {
		if (z_coordinates[i] >= z_coordinates[i + 1]) {
			std::cerr << "Coordinates in the z direction must be strictly increasing" << std::endl;
			return false;
		}
	}

	this->x_coordinates.clear();
	this->x_coordinates.reserve(x_coordinates.size());
	this->x_coordinates.insert(this->x_coordinates.begin(), x_coordinates.begin(), x_coordinates.end());
	this->y_coordinates.clear();
	this->y_coordinates.reserve(y_coordinates.size());
	this->y_coordinates.insert(this->y_coordinates.begin(), y_coordinates.begin(), y_coordinates.end());
	this->z_coordinates.clear();
	this->z_coordinates.reserve(z_coordinates.size());
	this->z_coordinates.insert(this->z_coordinates.begin(), z_coordinates.begin(), z_coordinates.end());

	this->x_length = this->x_coordinates.size() - 1;
	this->y_length = this->y_coordinates.size() - 1;
	this->z_length = this->z_coordinates.size() - 1;

	return true;
}

#else

bool CellGeometry::set_x_start(const double x_start)
{
	this->x_start = x_start;
	return true;
}

bool CellGeometry::set_y_start(const double y_start)
{
	this->y_start = y_start;
	return true;
}

bool CellGeometry::set_z_start(const double z_start)
{
	this->z_start = z_start;
	return true;
}


double CellGeometry::get_x_start(void)
{
	return this->x_start;
}

double CellGeometry::get_y_start(void)
{
	return this->y_start;
}

double CellGeometry::get_z_start(void)
{
	return this->z_start;
}


bool CellGeometry::set_cell_x_size(const double cell_x_size)
{
	if (cell_x_size <= 0) {
		std::cerr << "Cell size in x direction must be > 0" << std::endl;
		return false;
	}

	this->cell_x_size = cell_x_size;

	return true;
}

bool CellGeometry::set_cell_y_size(const double cell_y_size)
{
	if (cell_y_size <= 0) {
		std::cerr << "Cell size in y direction must be > 0" << std::endl;
		return false;
	}
	this->cell_y_size = cell_y_size;

	return true;
}

bool CellGeometry::set_cell_z_size(const double cell_z_size)
{
	if (cell_z_size <= 0) {
		std::cerr << "Cell size in z direction must be > 0" << std::endl;
		return false;
	}
	this->cell_z_size = cell_z_size;

	return true;
}


bool CellGeometry::set_x_length(const uint64_t x_length)
{
	if (x_length == 0) {
		std::cerr << "Length of the grid in unrefined cells must be > 0 in the x direction" << std::endl;
		return false;
	}

	uint64_t old_x_length = this->x_length;
	this->x_length = x_length;

	if (this->maximum_refinement_level > this->get_maximum_possible_refinement_level()) {
		std::cerr << "Grid could have too many cells for an uint64_t with current refinement level" << std::endl;
		this->x_length = old_x_length;
		return false;
	} else {
		return true;
	}
}

bool CellGeometry::set_y_length(const uint64_t y_length)
{
	if (y_length == 0) {
		std::cerr << "Length of the grid in unrefined cells must be > 0 in the y direction" << std::endl;
		return false;
	}

	double last_unrefined_cell = double(this->x_length) * double(y_length) * double(this->z_length);
	if (last_unrefined_cell > double(uint64_t(~0))) {
		std::cerr << "Grid would have too many unrefined cells to fit them into uint64_t (1..." << last_unrefined_cell << ")" << std::endl;
		return false;
	}

	this->y_length = y_length;
	return true;
}

bool CellGeometry::set_z_length(const uint64_t z_length)
{
	if (z_length == 0) {
		std::cerr << "Length of the grid in unrefined cells must be > 0 in the z direction" << std::endl;
		return false;
	}

	double last_unrefined_cell = double(this->x_length) * double(this->y_length) * double(z_length);
	if (last_unrefined_cell > double(uint64_t(~0))) {
		std::cerr << "Grid would have too many unrefined cells to fit them into uint64_t (1..." << last_unrefined_cell << ")" << std::endl;
		return false;
	}

	this->z_length = z_length;
	return true;
}

#endif


bool CellGeometry::set_maximum_refinement_level(const int maximum_refinement_level)
{
	if (maximum_refinement_level <= this->get_maximum_possible_refinement_level()) {
		this->maximum_refinement_level = maximum_refinement_level;
		return true;
	} else {
		return false;
	}
}


int CellGeometry::get_maximum_refinement_level(void)
{
	return this->maximum_refinement_level;
}


int CellGeometry::get_maximum_possible_refinement_level(void)
{
	double last_possible_cell = uint64_t(~0);
	double last_cell = this->x_length * this->y_length * this->z_length;
	int refinement_level = 0;

	while (last_cell <= last_possible_cell) {
		refinement_level++;
		last_cell += double(this->x_length) * double(this->y_length) * double(this->z_length) * double(uint64_t(1) << refinement_level * 3);
	}

	return refinement_level - 1;
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


uint64_t CellGeometry::get_x_index(const double x)
{
	#ifdef DCCRG_ARBITRARY_STRETCH

	assert((x >= this->x_coordinates[0]) && (x <= this->x_coordinates[this->get_x_length()]));

	uint64_t x_coord_start_index = 0;
	while (this->x_coordinates[x_coord_start_index] < x) {
		x_coord_start_index++;
	}
	x_coord_start_index--;

	double x_size_of_index = (this->x_coordinates[x_coord_start_index + 1] - this->x_coordinates[x_coord_start_index]) / this->get_cell_size_in_indices(1);

	uint64_t index_offset = 0;
	while (this->x_coordinates[x_coord_start_index] + index_offset * x_size_of_index < x) {
		index_offset++;
	}
	index_offset--;

	return x_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;

	#else

	assert((x >= this->get_x_start()) and (x <= this->get_x_start() + this->get_x_length() * this->cell_x_size));
	return uint64_t((x - this->get_x_start()) / (this->cell_x_size / (uint64_t(1) << this->maximum_refinement_level)));

	#endif
}

uint64_t CellGeometry::get_y_index(const double y)
{
	#ifdef DCCRG_ARBITRARY_STRETCH

	assert((y >= this->y_coordinates[0]) && (y <= this->y_coordinates[this->get_y_length()]));

	uint64_t y_coord_start_index = 0;
	while (this->y_coordinates[y_coord_start_index] < y) {
		y_coord_start_index++;
	}
	y_coord_start_index--;

	double y_size_of_index = (this->y_coordinates[y_coord_start_index + 1] - this->y_coordinates[y_coord_start_index]) / this->get_cell_size_in_indices(1);

	uint64_t index_offset = 0;
	while (this->y_coordinates[y_coord_start_index] + index_offset * y_size_of_index < y) {
		index_offset++;
	}
	index_offset--;

	return y_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;

	#else

	assert((y >= this->get_y_start()) and (y <= this->get_y_start() + this->get_y_length() * this->cell_y_size));
	return uint64_t((y - this->get_y_start()) / (this->cell_y_size / (uint64_t(1) << this->maximum_refinement_level)));

	#endif
}

uint64_t CellGeometry::get_z_index(const double z)
{
	#ifdef DCCRG_ARBITRARY_STRETCH

	assert((z >= this->z_coordinates[0]) && (z <= this->z_coordinates[this->get_z_length()]));

	uint64_t z_coord_start_index = 0;
	while (this->z_coordinates[z_coord_start_index] < z) {
		z_coord_start_index++;
	}
	z_coord_start_index--;

	double z_size_of_index = (this->z_coordinates[z_coord_start_index + 1] - this->z_coordinates[z_coord_start_index]) / this->get_cell_size_in_indices(1);

	uint64_t index_offset = 0;
	while (this->z_coordinates[z_coord_start_index] + index_offset * z_size_of_index < z) {
		index_offset++;
	}
	index_offset--;

	return z_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;

	#else

	assert((z >= this->get_z_start()) and (z <= this->get_z_start() + this->get_z_length() * this->cell_z_size));
	return uint64_t((z - this->get_z_start()) / (this->cell_z_size / (uint64_t(1) << this->maximum_refinement_level)));

	#endif
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


uint64_t CellGeometry::get_cell(const int refinement_level, const double x, const double y, const double z)
{
	if (refinement_level < 0 || refinement_level > this->maximum_refinement_level) {
		return 0;
	}

	return get_cell_from_indices(this->get_x_index(x), this->get_x_index(y), this->get_x_index(z), refinement_level);
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
	 if (cell == 0) {
		return std::numeric_limits<double>::quiet_NaN();
	 }

	if (this->get_refinement_level(cell) < 0 || this->get_refinement_level(cell) > this->maximum_refinement_level) {
		return std::numeric_limits<double>::quiet_NaN();
	}

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
	 if (cell == 0) {
		return std::numeric_limits<double>::quiet_NaN();
	 }

	if (this->get_refinement_level(cell) < 0 || this->get_refinement_level(cell) > this->maximum_refinement_level) {
		return std::numeric_limits<double>::quiet_NaN();
	}

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
	 if (cell == 0) {
		return std::numeric_limits<double>::quiet_NaN();
	 }

	if (this->get_refinement_level(cell) < 0 || this->get_refinement_level(cell) > this->maximum_refinement_level) {
		return std::numeric_limits<double>::quiet_NaN();
	}

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

double CellGeometry::get_cell_x(const int refinement_level, const uint64_t x_index)
{
	if (refinement_level < 0 || refinement_level > this->maximum_refinement_level) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	if (x_index > this->x_length * (uint64_t(1) << this->maximum_refinement_level)) {
		return std::numeric_limits<double>::quiet_NaN();
	}

#ifdef DCCRG_ARBITRARY_STRETCH
	uint64_t unref_cell_x_coord_start_index = x_index / (uint64_t(1) << this->maximum_refinement_level);

	double unref_cell_x_size = this->x_coordinates[unref_cell_x_coord_start_index + 1] - this->x_coordinates[unref_cell_x_coord_start_index];
	double size_of_local_index = unref_cell_x_size / (uint64_t(1) << this->maximum_refinement_level);

	return this->x_coordinates[unref_cell_x_coord_start_index] + size_of_local_index * (x_index - unref_cell_x_coord_start_index * (uint64_t(1) << this->maximum_refinement_level)) + size_of_local_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level)) / 2;
}
#else
	return this->x_start + x_index * this->cell_x_size / (uint64_t(1) << this->maximum_refinement_level) + this->cell_x_size / (uint64_t(1) << refinement_level) / 2;
}
#endif

double CellGeometry::get_cell_y(const int refinement_level, const uint64_t y_index)
{
	if (refinement_level < 0 || refinement_level > this->maximum_refinement_level) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	if (y_index > this->y_length * (uint64_t(1) << this->maximum_refinement_level)) {
		return std::numeric_limits<double>::quiet_NaN();
	}

#ifdef DCCRG_ARBITRARY_STRETCH
	uint64_t unref_cell_y_coord_start_index = y_index / (uint64_t(1) << this->maximum_refinement_level);

	double unref_cell_y_size = this->y_coordinates[unref_cell_y_coord_start_index + 1] - this->y_coordinates[unref_cell_y_coord_start_index];
	double size_of_local_index = unref_cell_y_size / (uint64_t(1) << this->maximum_refinement_level);

	return this->y_coordinates[unref_cell_y_coord_start_index] + size_of_local_index * (y_index - unref_cell_y_coord_start_index * (uint64_t(1) << this->maximum_refinement_level)) + size_of_local_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level)) / 2;
}
#else
	return this->y_start + y_index * this->cell_y_size / (uint64_t(1) << this->maximum_refinement_level) + this->cell_y_size / (uint64_t(1) << refinement_level) / 2;
}
#endif

double CellGeometry::get_cell_z(const int refinement_level, const uint64_t z_index)
{
	if (refinement_level < 0 || refinement_level > this->maximum_refinement_level) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	if (z_index > this->z_length * (uint64_t(1) << this->maximum_refinement_level)) {
		return std::numeric_limits<double>::quiet_NaN();
	}

#ifdef DCCRG_ARBITRARY_STRETCH
	uint64_t unref_cell_z_coord_start_index = z_index / (uint64_t(1) << this->maximum_refinement_level);

	double unref_cell_z_size = this->z_coordinates[unref_cell_z_coord_start_index + 1] - this->z_coordinates[unref_cell_z_coord_start_index];
	double size_of_local_index = unref_cell_z_size / (uint64_t(1) << this->maximum_refinement_level);

	return this->z_coordinates[unref_cell_z_coord_start_index] + size_of_local_index * (z_index - unref_cell_z_coord_start_index * (uint64_t(1) << this->maximum_refinement_level)) + size_of_local_index * (uint64_t(1) << (this->maximum_refinement_level - refinement_level)) / 2;
}
#else
	return this->z_start + z_index * this->cell_z_size / (uint64_t(1) << this->maximum_refinement_level) + this->cell_z_size / (uint64_t(1) << refinement_level) / 2;
}
#endif


uint64_t CellGeometry::get_x_length(void)
{
	return this->x_length;
}

uint64_t CellGeometry::get_y_length(void)
{
	return this->y_length;
}

uint64_t CellGeometry::get_z_length(void)
{
	return this->z_length;
}



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
