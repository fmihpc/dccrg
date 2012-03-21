/*
The mapping logic between cells' and their geometry (location, size and comparable stuff)
in dccrg, constant version.

Copyright 2009, 2010, 2011, 2012 Finnish Meteorological Institute

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


#ifndef DCCRG_CONSTANT_GEOMETRY_HPP
#define DCCRG_CONSTANT_GEOMETRY_HPP


#include "cassert"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "stdint.h"
#include "vector"

#include "dccrg_index.hpp"

namespace dccrg {

/*!
A geometry class in which the sizes of unrefined cells are given
by three floating points numbers.
*/
class ConstantGeometry : public Index
{

public:

	/*!
	Creates and sets the geometry of the grid to the following:
		-starting corner at (0, 0, 0)
		-size of unrefined cells in each direction: 1
	*/
	ConstantGeometry()
	{
		this->x_start = 0;
		this->y_start = 0;
		this->z_start = 0;

		this->cell_x_size = 1;
		this->cell_y_size = 1;
		this->cell_z_size = 1;
	}

	/*!
	Sets the geometry of the grid to the following:
		-starting corner at (0, 0, 0)
		-size of unrefined cells in each direction: 1
	*/
	~ConstantGeometry()
	{
		this->x_start = 0;
		this->y_start = 0;
		this->z_start = 0;

		this->cell_x_size = 1;
		this->cell_y_size = 1;
		this->cell_z_size = 1;
	}


	/*!
	Sets the grid's geometry to given values.

	x, y and z_length set the number of unrefined cells in the grid in x, y and z direction.

	x, y and z_start set the starting corner of the grid, e.g. the first face
	of the first unrefined cell(s) in x, y  and z direction.

	x, y and z_size set the size of unrefined cell in x, y and z direction.

	Returns true on success and false otherwise.
	*/
	bool set_geometry(
		const uint64_t given_x_length, const uint64_t given_y_length, const uint64_t given_z_length,
		const double given_x_start, const double given_y_start, const double given_z_start,
		const double given_cell_x_size, const double given_cell_y_size, const double given_cell_z_size
	) {
		if (!this->set_length(given_x_length, given_y_length, given_z_length)) {
			return false;
		}

		// FIXME: check that all coordinates fit into a double
		this->x_start = given_x_start;
		this->y_start = given_y_start;
		this->z_start = given_z_start;

		if (given_cell_x_size <= 0) {
			std::cerr << "Cell size in x direction must be > 0, but is " << given_cell_x_size << std::endl;
			return false;
		}
		this->cell_x_size = given_cell_x_size;

		if (given_cell_y_size <= 0) {
			std::cerr << "Cell size in y direction must be > 0, but is " << given_cell_y_size << std::endl;
			return false;
		}
		this->cell_y_size = given_cell_y_size;

		if (given_cell_z_size <= 0) {
			std::cerr << "Cell size in z direction must be > 0, but is " << given_cell_z_size << std::endl;
			return false;
		}
		this->cell_z_size = given_cell_z_size;

		return true;
	}

	/*!
	Returns the starting corner of the grid in x direction.
	*/
	double get_x_start(void) const
	{
		return this->x_start;
	}

	/*!
	Returns the starting corner of the grid in y direction.
	*/
	double get_y_start(void) const
	{
		return this->y_start;
	}

	/*!
	Returns the starting corner of the grid in z direction.
	*/
	double get_z_start(void) const
	{
		return this->z_start;
	}


	/*!
	Returns the end corner of the grid in x direction.
	*/
	double get_x_end(void) const
	{
		return this->x_start + double(this->x_length) * this->cell_x_size;
	}

	/*!
	Returns the end corner of the grid in y direction.
	*/
	double get_y_end(void) const
	{
		return this->y_start + double(this->y_length) * this->cell_y_size;
	}

	/*!
	Returns the end corner of the grid in z direction.
	*/
	double get_z_end(void) const
	{
		return this->z_start + double(this->z_length) * this->cell_z_size;
	}


	/*!
	Returns the length of unrefined cells in x direction
	*/
	double get_unrefined_cell_x_size(void) const
	{
		return this->cell_x_size;
	}

	/*!
	Returns the length of unrefined cells in y direction
	*/
	double get_unrefined_cell_y_size(void) const
	{
		return this->cell_y_size;
	}

	/*!
	Returns the length of unrefined cells in z direction
	*/
	double get_unrefined_cell_z_size(void) const
	{
		return this->cell_z_size;
	}


	/*!
	Returns the length of given cell in x direction.

	Returns a quiet nan if given an invalid cell.
	*/
	double get_cell_x_size(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->cell_x_size / double(uint64_t(1) << this->get_refinement_level(cell));
	}

	/*!
	Returns the length of given cell in y direction.
	*/
	double get_cell_y_size(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->cell_y_size / double(uint64_t(1) << this->get_refinement_level(cell));
	}

	/*!
	Returns the length of given cell in z direction.
	*/
	double get_cell_z_size(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->cell_z_size / double(uint64_t(1) << this->get_refinement_level(cell));
	}


	/*!
	Returns the center of given cell in x direction.
	*/
	double get_cell_x(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->get_indices(cell);

		return this->x_start + double(indices[0]) * this->cell_x_size / double(uint64_t(1) << this->max_refinement_level) + this->get_cell_x_size(cell) / 2;
	}

	/*!
	Returns the center of given cell in y direction.
	*/
	double get_cell_y(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->get_indices(cell);

		return this->y_start + double(indices[1]) * this->cell_y_size / double(uint64_t(1) << this->max_refinement_level) + this->get_cell_y_size(cell) / 2;
	}

	/*!
	Returns the center of given cell in z direction.
	*/
	double get_cell_z(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->get_indices(cell);

		return this->z_start + double(indices[2]) * this->cell_z_size / double(uint64_t(1) << this->max_refinement_level) + this->get_cell_z_size(cell) / 2;
	}

	/*!
	Returns the x coordinate of given cells face in negative x direction.
	*/
	double get_cell_x_min(const uint64_t cell) const
	{
		return this->get_cell_x(cell) - this->get_cell_x_size(cell) / 2;
	}

	/*!
	Returns the y coordinate of given cells face in negative y direction.
	*/
	double get_cell_y_min(const uint64_t cell) const
	{
		return this->get_cell_y(cell) - this->get_cell_y_size(cell) / 2;
	}

	/*!
	Returns the z coordinate of given cells face in negative z direction.
	*/
	double get_cell_z_min(const uint64_t cell) const
	{
		return this->get_cell_z(cell) - this->get_cell_z_size(cell) / 2;
	}


	/*!
	Returns the x coordinate of the cells face in positive x direction.
	*/
	double get_cell_x_max(const uint64_t cell) const
	{
		return this->get_cell_x(cell) + this->get_cell_x_size(cell) / 2;
	}

	/*!
	Returns the y coordinate of the cells face in positive y direction.
	*/
	double get_cell_y_max(const uint64_t cell) const
	{
		return this->get_cell_y(cell) + this->get_cell_y_size(cell) / 2;
	}

	/*!
	Returns the z coordinate of the cells face in positive z direction.
	*/
	double get_cell_z_max(const uint64_t cell) const
	{
		return this->get_cell_z(cell) + this->get_cell_z_size(cell) / 2;
	}


	/*!
	Returns the center of a cell in x direction of given refinement level at given index.
	*/
	double get_cell_x(const int refinement_level, const uint64_t x_index) const
	{
		if (refinement_level < 0 || refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (x_index > this->x_length * (uint64_t(1) << this->max_refinement_level)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->x_start
			+ double(x_index) * this->cell_x_size / double(uint64_t(1) << this->max_refinement_level)
			+ this->cell_x_size / double(uint64_t(1) << refinement_level) / 2;
	}

	/*!
	Returns the center of a cell in y direction of given refinement level at given index.
	*/
	double get_cell_y(const int refinement_level, const uint64_t y_index) const
	{
		if (refinement_level < 0 || refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (y_index > this->y_length * (uint64_t(1) << this->max_refinement_level)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->y_start
			+ double(y_index) * this->cell_y_size / double(uint64_t(1) << this->max_refinement_level)
			+ this->cell_y_size / double(uint64_t(1) << refinement_level) / 2;
	}

	/*!
	Returns the center of a cell in z direction of given refinement level at given index.
	*/
	double get_cell_z(const int refinement_level, const uint64_t z_index) const
	{
		if (refinement_level < 0 || refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (z_index > this->z_length * (uint64_t(1) << this->max_refinement_level)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->z_start
			+ double(z_index) * this->cell_z_size / double(uint64_t(1) << this->max_refinement_level)
			+ this->cell_z_size / double(uint64_t(1) << refinement_level) / 2;
	}


	/*!
	Returns the cell of given refinement level at given location, or 0 if outside of the current grid in location or refinement level
	*/
	uint64_t get_cell(const int refinement_level, const double x, const double y, const double z) const
	{
		if (refinement_level < 0 || refinement_level > this->max_refinement_level) {
			return error_cell;
		}

		return this->get_cell_from_indices(
			this->get_x_index_of_coord(x),
			this->get_y_index_of_coord(y),
			this->get_z_index_of_coord(z),
			refinement_level
		);
	}


	/*!
	Returns the real value of given x coordinate in this geometry.

	Returns given x if it is inside this geometry.
	Returns a quiet NaN if this geometry is not periodic in x direction
	and given x is outside of the geometry.
	If this geometry is periodic in x returns a value inside the
	geometry that is at the same location in the geometry as given x.
	*/
	double get_real_x(const double x) const
	{
		if (x >= this->get_x_start()
		&& x <= this->get_x_end()) {

			return x;

		} else if (!this->periodic[0]) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_x_end() - this->get_x_start();

			if (x < this->get_x_start()) {

				const double distance = this->get_x_start() - x;
				return x + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = x - this->get_x_end();
				return x - grid_size * ceil(distance/ grid_size);
			}
		}
	}

	/*!
	Returns the real value of given y coordinate in this geometry.

	Returns given y if it is inside this geometry.
	Returns a quiet NaN if this geometry is not periodic in y direction
	and given y is outside of the geometry.
	If this geometry is periodic in y returns a value inside the
	geometry that is at the same location in the geometry as given y.
	*/
	double get_real_y(const double y) const
	{
		if (y >= this->get_y_start()
		&& y <= this->get_y_end()) {

			return y;

		} else if (!this->periodic[0]) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_y_end() - this->get_y_start();

			if (y < this->get_y_start()) {

				const double distance = this->get_y_start() - y;
				return y + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = y - this->get_y_end();
				return y - grid_size * ceil(distance/ grid_size);
			}
		}
	}

	/*!
	Returns the real value of given z coordinate in this geometry.

	Returns given z if it is inside this geometry.
	Returns a quiet NaN if this geometry is not periodic in z direction
	and given z is outside of the geometry.
	If this geometry is periodic in z returns a value inside the
	geometry that is at the same location in the geometry as given z.
	*/
	double get_real_z(const double z) const
	{
		if (z >= this->get_z_start()
		&& z <= this->get_z_end()) {

			return z;

		} else if (!this->periodic[0]) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_z_end() - this->get_z_start();

			if (z < this->get_z_start()) {

				const double distance = this->get_z_start() - z;
				return z + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = z - this->get_z_end();
				return z - grid_size * ceil(distance/ grid_size);
			}
		}
	}


	/*!
	Returns the x index of given location, starting from 0.

	Returns error_index if given location is outside of the grid
	and the grid is not periodic in that direction.
	*/
	uint64_t get_x_index_of_coord(double x) const
	{
		x = this->get_real_x(x);

		if (isnan(x)
		|| x < this->get_x_start()
		|| x > this->get_x_end()) {

			return error_index;

		} else {

			return uint64_t(floor(
				(x - this->get_x_start())
				/ (this->cell_x_size / (uint64_t(1) << this->max_refinement_level))
			));
		}
	}

	/*!
	Returns the y index of given location, starting from 0.

	Returns error_index if given location is outside of the grid
	and the grid is not periodic in that direction.
	*/
	uint64_t get_y_index_of_coord(double y) const
	{
		y = this->get_real_y(y);

		if (isnan(y)
		|| y < this->get_y_start()
		|| y > this->get_y_end()) {

			return error_index;

		} else {

			return uint64_t(floor(
				(y - this->get_y_start())
				/ (this->cell_y_size / (uint64_t(1) << this->max_refinement_level))
			));
		}
	}

	/*!
	Returns the z index of given location, starting from 0.

	Returns error_index if given location is outside of the grid
	and the grid is not periodic in that direction.
	*/
	uint64_t get_z_index_of_coord(double z) const
	{
		z = this->get_real_z(z);

		if (isnan(z)
		|| z < this->get_z_start()
		|| z > this->get_z_end()) {

			return error_index;

		} else {

			return uint64_t(floor(
				(z - this->get_z_start())
				/ (this->cell_z_size / (uint64_t(1) << this->max_refinement_level))
			));
		}
	}


	// Optional BOOST serialization support
	#ifdef BOOST_SERIALIZATION_LIBRARY_VERSION
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & this->x_start
			& this->y_start
			& this->z_start
			& this->cell_x_size
			& this->cell_y_size
			& this->cell_z_size
			& this->x_length
			& this->y_length
			& this->z_length
			& this->maximum_refinement_level
			& this->periodic;
	}
	#endif


	/*!
	Sets the periodicity of the geometry.

	index = 0 == x direction.
	*/
	void set_periodicity(const size_t index, const bool value)
	{
		if (index > 2) {
			return;
		}

		this->periodic[index] = value;
	}


	/*!
	Returns whether the geometry in periodic in given direction.

	index = 0 == x direction.
	Returns false if given index > 2.
	*/
	bool is_periodic(const size_t index) const
	{
		if (index > 2) {
			return false;
		}

		return this->periodic[index];
	}


private:

	// starting corner coordinates of the grid
	double x_start, y_start, z_start;
	// length of unrefined cells in all directions
	double cell_x_size, cell_y_size, cell_z_size;

	// periodic[0] == true means that the grid wraps around in x direction
	bool periodic[3];


};	// class

}	// namespace

#endif

