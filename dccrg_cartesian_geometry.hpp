/*
Dccrg class for a cartesian geometry in which cells are cubes.

Copyright 2009, 2010, 2011, 2012, 2013 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_CARTESIAN_GEOMETRY_HPP
#define DCCRG_CARTESIAN_GEOMETRY_HPP


#include "cassert"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "stdint.h"
#include "vector"

#include "dccrg_length.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"


namespace dccrg {


/*!
\brief Geometry class for dccrg with cubic cells

A geometry class in which the sizes of cells of
refinement level 0 are given by three floating points numbers.
*/
class Cartesian_Geometry
{

public:

	/*!
	Public read-only version of the grid's length in cells of refinement level 0.

	\see Grid_Length
	*/
	const Grid_Length& length;

	/*!
	Public read-only version of the mapping of a
	cell's ids to its size and location in the grid.

	\see Mapping
	*/
	const Mapping& mapping;

	/*!
	Public read-only version of the grid's topology.

	\see Grid_Topology
	*/
	const Grid_Topology& topology;

	/*!
	Creates and sets the geometry of the grid to the following:
		- starting corner at (0, 0, 0)
		- size of unrefined cells in each direction: 1
	*/
	Cartesian_Geometry(
		const Grid_Length& given_length,
		const Mapping& given_mapping,
		const Grid_Topology& given_topology
	) :
		length(given_length),
		mapping(given_mapping),
		topology(given_topology)
	{
		this->start_x = 0;
		this->start_y = 0;
		this->start_z = 0;

		this->cell_length_x = 1;
		this->cell_length_y = 1;
		this->cell_length_z = 1;
	}

	/*!
	Sets the geometry of the grid to the following:
		- starting corner at (0, 0, 0)
		- size of unrefined cells in each direction: 1
	*/
	~Cartesian_Geometry()
	{
		this->start_x = 0;
		this->start_y = 0;
		this->start_z = 0;

		this->cell_length_x = 1;
		this->cell_length_y = 1;
		this->cell_length_z = 1;
	}


	/*!
	Sets the grid's length in cells and its geometry to given values.

	- given x, y and start_z set the starting corner of the grid, e.g.
	  the first face of the first cell(s) of refinement level 0 in
	  x, y  and z dimensions.
	- given x, y and length_z set the size of each cell of
	  refinement level 0 in x, y and z dimensions.

	Returns true on success and false otherwise.
	*/
	bool set(
		const double given_start_x,
		const double given_start_y,
		const double given_start_z,
		const double given_cell_length_x,
		const double given_cell_length_y,
		const double given_cell_length_z
	) {
		// FIXME: check that all coordinates fit into a double
		this->start_x = given_start_x;
		this->start_y = given_start_y;
		this->start_z = given_start_z;

		if (given_cell_length_x <= 0) {
			std::cerr << "Cell size in x direction must be > 0, but is " << given_cell_length_x
				<< std::endl;
			return false;
		}
		this->cell_length_x = given_cell_length_x;

		if (given_cell_length_y <= 0) {
			std::cerr << "Cell size in y direction must be > 0, but is " << given_cell_length_y
				<< std::endl;
			return false;
		}
		this->cell_length_y = given_cell_length_y;

		if (given_cell_length_z <= 0) {
			std::cerr << "Cell size in z direction must be > 0, but is " << given_cell_length_z
				<< std::endl;
			return false;
		}
		this->cell_length_z = given_cell_length_z;

		return true;
	}


	/*!
	Sets the same geometry as in the given one.
	*/
	bool set(const Cartesian_Geometry& other)
	{
		return this->set(
			other.get_start_x(),
			other.get_start_y(),
			other.get_start_z(),
			other.get_unrefined_cell_length_x(),
			other.get_unrefined_cell_length_y(),
			other.get_unrefined_cell_length_z()
		);
	}


	/*!
	Returns the starting corner of the grid in x direction.
	*/
	double get_start_x() const
	{
		return this->start_x;
	}

	/*!
	Returns the starting corner of the grid in y direction.
	*/
	double get_start_y() const
	{
		return this->start_y;
	}

	/*!
	Returns the starting corner of the grid in z direction.
	*/
	double get_start_z() const
	{
		return this->start_z;
	}


	/*!
	Returns the end corner of the grid in x direction.
	*/
	double get_end_x() const
	{
		return this->start_x + double(this->length.get()[0]) * this->cell_length_x;
	}

	/*!
	Returns the end corner of the grid in y direction.
	*/
	double get_end_y() const
	{
		return this->start_y + double(this->length.get()[1]) * this->cell_length_y;
	}

	/*!
	Returns the end corner of the grid in z direction.
	*/
	double get_end_z() const
	{
		return this->start_z + double(this->length.get()[2]) * this->cell_length_z;
	}


	/*!
	Returns the length of unrefined cells in x direction
	*/
	double get_unrefined_cell_length_x() const
	{
		return this->cell_length_x;
	}

	/*!
	Returns the length of unrefined cells in y direction
	*/
	double get_unrefined_cell_length_y() const
	{
		return this->cell_length_y;
	}

	/*!
	Returns the length of unrefined cells in z direction
	*/
	double get_unrefined_cell_length_z() const
	{
		return this->cell_length_z;
	}


	/*!
	Returns the length of given cell in x direction.

	Returns a quiet nan if given an invalid cell.
	*/
	double get_cell_length_x(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->cell_length_x / double(uint64_t(1) << this->mapping.get_refinement_level(cell));
	}

	/*!
	Returns the length of given cell in y direction.
	*/
	double get_cell_length_y(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->cell_length_y / double(uint64_t(1) << this->mapping.get_refinement_level(cell));
	}

	/*!
	Returns the length of given cell in z direction.
	*/
	double get_cell_length_z(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->cell_length_z / double(uint64_t(1) << this->mapping.get_refinement_level(cell));
	}


	/*!
	Returns the center of given cell in x direction.
	*/
	double get_cell_x(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		return this->start_x + double(indices[0]) * this->cell_length_x / double(uint64_t(1) << this->mapping.get_maximum_refinement_level()) + this->get_cell_length_x(cell) / 2;
	}

	/*!
	Returns the center of given cell in y direction.
	*/
	double get_cell_y(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		return this->start_y + double(indices[1]) * this->cell_length_y / double(uint64_t(1) << this->mapping.get_maximum_refinement_level()) + this->get_cell_length_y(cell) / 2;
	}

	/*!
	Returns the center of given cell in z direction.
	*/
	double get_cell_z(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		return this->start_z + double(indices[2]) * this->cell_length_z / double(uint64_t(1) << this->mapping.get_maximum_refinement_level()) + this->get_cell_length_z(cell) / 2;
	}

	/*!
	Returns the x coordinate of given cells face in negative x direction.
	*/
	double get_cell_x_min(const uint64_t cell) const
	{
		return this->get_cell_x(cell) - this->get_cell_length_x(cell) / 2;
	}

	/*!
	Returns the y coordinate of given cells face in negative y direction.
	*/
	double get_cell_y_min(const uint64_t cell) const
	{
		return this->get_cell_y(cell) - this->get_cell_length_y(cell) / 2;
	}

	/*!
	Returns the z coordinate of given cells face in negative z direction.
	*/
	double get_cell_z_min(const uint64_t cell) const
	{
		return this->get_cell_z(cell) - this->get_cell_length_z(cell) / 2;
	}


	/*!
	Returns the x coordinate of the cells face in positive x direction.
	*/
	double get_cell_x_max(const uint64_t cell) const
	{
		return this->get_cell_x(cell) + this->get_cell_length_x(cell) / 2;
	}

	/*!
	Returns the y coordinate of the cells face in positive y direction.
	*/
	double get_cell_y_max(const uint64_t cell) const
	{
		return this->get_cell_y(cell) + this->get_cell_length_y(cell) / 2;
	}

	/*!
	Returns the z coordinate of the cells face in positive z direction.
	*/
	double get_cell_z_max(const uint64_t cell) const
	{
		return this->get_cell_z(cell) + this->get_cell_length_z(cell) / 2;
	}


	/*!
	Returns the center of a cell in x direction of given refinement level at given index.
	*/
	double get_cell_x(const int refinement_level, const uint64_t x_index) const
	{
		if (refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (x_index > this->length.get()[0] * (uint64_t(1) << this->mapping.get_maximum_refinement_level())) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->start_x
			+ double(x_index) * this->cell_length_x / double(uint64_t(1) << this->mapping.get_maximum_refinement_level())
			+ this->cell_length_x / double(uint64_t(1) << refinement_level) / 2;
	}

	/*!
	Returns the center of a cell in y direction of given refinement level at given index.
	*/
	double get_cell_y(const int refinement_level, const uint64_t y_index) const
	{
		if (refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (y_index > this->length.get()[1] * (uint64_t(1) << this->mapping.get_maximum_refinement_level())) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->start_y
			+ double(y_index) * this->cell_length_y / double(uint64_t(1) << this->mapping.get_maximum_refinement_level())
			+ this->cell_length_y / double(uint64_t(1) << refinement_level) / 2;
	}

	/*!
	Returns the center of a cell in z direction of given refinement level at given index.
	*/
	double get_cell_z(const int refinement_level, const uint64_t z_index) const
	{
		if (refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (z_index > this->length.get()[2] * (uint64_t(1) << this->mapping.get_maximum_refinement_level())) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		return this->start_z
			+ double(z_index) * this->cell_length_z / double(uint64_t(1) << this->mapping.get_maximum_refinement_level())
			+ this->cell_length_z / double(uint64_t(1) << refinement_level) / 2;
	}


	/*!
	Returns the cell of given refinement level at given location, or 0 if outside of the current grid in location or refinement level
	*/
	uint64_t get_cell(
		const int refinement_level,
		const double x,
		const double y,
		const double z
	) const {
		if (refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return error_cell;
		}

		const Types<3>::indices_t index = {{
			this->get_x_index_of_coord(x),
			this->get_y_index_of_coord(y),
			this->get_z_index_of_coord(z)
		}};

		return this->mapping.get_cell_from_indices(index, refinement_level);
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
		if (x >= this->get_start_x()
		&& x <= this->get_end_x()) {

			return x;

		} else if (!this->topology.is_periodic(0)) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_end_x() - this->get_start_x();

			if (x < this->get_start_x()) {

				const double distance = this->get_start_x() - x;
				return x + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = x - this->get_end_x();
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
		if (y >= this->get_start_y()
		&& y <= this->get_end_y()) {

			return y;

		} else if (!this->topology.is_periodic(1)) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_end_y() - this->get_start_y();

			if (y < this->get_start_y()) {

				const double distance = this->get_start_y() - y;
				return y + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = y - this->get_end_y();
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
		if (z >= this->get_start_z()
		&& z <= this->get_end_z()) {

			return z;

		} else if (!this->topology.is_periodic(2)) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_end_z() - this->get_start_z();

			if (z < this->get_start_z()) {

				const double distance = this->get_start_z() - z;
				return z + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = z - this->get_end_z();
				return z - grid_size * ceil(distance/ grid_size);
			}
		}
	}


	/*!
	Returns the x index of given location, starting from 0.

	Returns error_index if given location is outside of the grid
	and the grid is not periodic in that direction.
	*/
	uint64_t get_x_index_of_coord(const double given_x) const
	{
		const double x = this->get_real_x(given_x);

		if (::isnan(x)
		|| x < this->get_start_x()
		|| x > this->get_end_x()) {

			return error_index;

		} else {

			return uint64_t(floor(
				(x - this->get_start_x())
				/ (this->cell_length_x / double(uint64_t(1) << this->mapping.get_maximum_refinement_level()))
			));
		}
	}

	/*!
	Returns the y index of given location, starting from 0.

	Returns error_index if given location is outside of the grid
	and the grid is not periodic in that direction.
	*/
	uint64_t get_y_index_of_coord(const double given_y) const
	{
		const double y = this->get_real_y(given_y);

		if (::isnan(y)
		|| y < this->get_start_y()
		|| y > this->get_end_y()) {

			return error_index;

		} else {

			return uint64_t(floor(
				(y - this->get_start_y())
				/ (this->cell_length_y / double(uint64_t(1) << this->mapping.get_maximum_refinement_level()))
			));
		}
	}

	/*!
	Returns the z index of given location, starting from 0.

	Returns error_index if given location is outside of the grid
	and the grid is not periodic in that direction.
	*/
	uint64_t get_z_index_of_coord(const double given_z) const
	{
		const double z = this->get_real_z(given_z);

		if (::isnan(z)
		|| z < this->get_start_z()
		|| z > this->get_end_z()) {

			return error_index;

		} else {

			return uint64_t(floor(
				(z - this->get_start_z())
				/ (this->cell_length_z / double(uint64_t(1) << this->mapping.get_maximum_refinement_level()))
			));
		}
	}


private:

	// starting corner coordinates of the grid
	double start_x, start_y, start_z;
	// length of unrefined cells in all directions
	double cell_length_x, cell_length_y, cell_length_z;


};	// class

}	// namespace

#endif

