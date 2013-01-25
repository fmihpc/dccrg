/*
The mapping logic between cells' and their geometry (location, size, etc.)
in dccrg, stretched cartesian version.

Copyright 2009, 2010, 2011, 2012, 2013 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this dccrg. If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_STRETCHED_CARTESIAN_GEOMETRY_HPP
#define DCCRG_STRETCHED_CARTESIAN_GEOMETRY_HPP


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
A geometry class in which the coordinates of unrefined cells are given
by three vectors of floating points numbers.
*/
class Stretched_Cartesian_Geometry : public Index
{

public:

	/*!
	Creates and sets the geometry of the grid to the following:
		- starting corner at (0, 0, 0)
		- size of unrefined cells in each direction: 1
	*/
	Stretched_Cartesian_Geometry()
	{
		this->x_coordinates.clear();
		this->x_coordinates.reserve(this->length_x + 1);
		for (uint64_t i = 0; i <= this->length_x; i++) {
			this->x_coordinates.push_back(i);
		}

		this->y_coordinates.clear();
		this->y_coordinates.reserve(this->length_y + 1);
		for (uint64_t i = 0; i <= this->length_y; i++) {
			this->y_coordinates.push_back(i);
		}

		this->z_coordinates.clear();
		this->z_coordinates.reserve(this->length_z + 1);
		for (uint64_t i = 0; i <= this->length_z; i++) {
			this->z_coordinates.push_back(i);
		}
	}

	/*!
	Sets the geometry of the grid to the following:
		- starting corner at (0, 0, 0)
		- size of unrefined cells in each direction: 1
	*/
	~Stretched_Cartesian_Geometry()
	{
		this->x_coordinates.clear();
		this->x_coordinates.reserve(this->length_x + 1);
		for (uint64_t i = 0; i <= this->length_x; i++) {
			this->x_coordinates.push_back(i);
		}

		this->y_coordinates.clear();
		this->y_coordinates.reserve(this->length_y + 1);
		for (uint64_t i = 0; i <= this->length_y; i++) {
			this->y_coordinates.push_back(i);
		}

		this->z_coordinates.clear();
		this->z_coordinates.reserve(this->length_z + 1);
		for (uint64_t i = 0; i <= this->length_z; i++) {
			this->z_coordinates.push_back(i);
		}
	}

	/*!
	Sets the coordinates of unrefined cells in the grid.

	First coordinate is the starting point of the grid and the following
	ith value is the endpoint of the ith unrefined cell
	At least two values must be given for each direction and all values
	must be strictly increasing.
	Returns true if successful, probably invalidating all previous cell
	geometries, e.g. their location, size, etc.
	Returns false if unsuccessful and in that case has no effect.
	Automatically maximizes max_refinement_level.
	*/
	bool set_geometry(
		const std::vector<double> x_coordinates,
		const std::vector<double> y_coordinates,
		const std::vector<double> z_coordinates
	) {
		if (x_coordinates.size() < 2) {
			std::cerr << "At least two coordinates are required for grid cells in the x direction"
				<< std::endl;
			return false;
		}
		if (y_coordinates.size() < 2) {
			std::cerr << "At least two coordinates are required for grid cells in the y direction"
				<< std::endl;
			return false;
		}
		if (z_coordinates.size() < 2) {
			std::cerr << "At least two coordinates are required for grid cells in the z direction"
				<< std::endl;
			return false;
		}

		for (uint64_t i = 0; i < x_coordinates.size() - 1; i++) {
			if (x_coordinates[i] >= x_coordinates[i + 1]) {
				std::cerr << "Coordinates in the x direction must be strictly increasing"
					<< std::endl;
				return false;
			}
		}
		for (uint64_t i = 0; i < y_coordinates.size() - 1; i++) {
			if (y_coordinates[i] >= y_coordinates[i + 1]) {
				std::cerr << "Coordinates in the y direction must be strictly increasing"
					<< std::endl;
				return false;
			}
		}
		for (uint64_t i = 0; i < z_coordinates.size() - 1; i++) {
			if (z_coordinates[i] >= z_coordinates[i + 1]) {
				std::cerr << "Coordinates in the z direction must be strictly increasing"
					<< std::endl;
				return false;
			}
		}

		this->x_coordinates.clear();
		this->x_coordinates.reserve(x_coordinates.size());
		this->x_coordinates.insert(
			this->x_coordinates.begin(),
			x_coordinates.begin(),
			x_coordinates.end()
		);

		this->y_coordinates.clear();
		this->y_coordinates.reserve(y_coordinates.size());
		this->y_coordinates.insert(
			this->y_coordinates.begin(),
			y_coordinates.begin(),
			y_coordinates.end()
		);

		this->z_coordinates.clear();
		this->z_coordinates.reserve(z_coordinates.size());
		this->z_coordinates.insert(
			this->z_coordinates.begin(),
			z_coordinates.begin(),
			z_coordinates.end()
		);

		if (!this->set_length(
			this->x_coordinates.size() - 1,
			this->y_coordinates.size() - 1,
			this->z_coordinates.size() - 1
		)) {
			return false;
		}

		return true;
	}


	/*!
	Returns the starting corner of the grid in x direction.
	*/
	double get_start_x() const
	{
		return this->x_coordinates[0];
	}

	/*!
	Returns the starting corner of the grid in y direction.
	*/
	double get_start_y() const
	{
		return this->y_coordinates[0];
	}

	/*!
	Returns the starting corner of the grid in z direction.
	*/
	double get_start_z() const
	{
		return this->z_coordinates[0];
	}


	/*!
	Returns the end corner of the grid in x direction.
	*/
	double get_x_end() const
	{
		return this->x_coordinates[this->length_x];
	}

	/*!
	Returns the end corner of the grid in y direction.
	*/
	double get_y_end() const
	{
		return this->y_coordinates[this->length_y];
	}

	/*!
	Returns the end corner of the grid in z direction.
	*/
	double get_z_end() const
	{
		return this->z_coordinates[this->length_z];
	}


	/*!
	Returns the length of given cell in x direction.
	*/
	double get_cell_x_size(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->get_refinement_level(cell) >= 0);
		assert(this->get_refinement_level(cell) <= this->max_refinement_level);

		int refinement_level = this->get_refinement_level(cell);
		uint64_t unref_cell_x_index = this->get_unref_cell_x_coord_start_index(cell);

		return
			  (this->x_coordinates[unref_cell_x_index + 1] - this->x_coordinates[unref_cell_x_index])
			/ (uint64_t(1) << refinement_level);
	}

	/*!
	Returns the length of given cell in y direction.
	*/
	double get_cell_y_size(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->get_refinement_level(cell) >= 0);
		assert(this->get_refinement_level(cell) <= this->max_refinement_level);

		int refinement_level = this->get_refinement_level(cell);
		uint64_t unref_cell_y_index = this->get_unref_cell_y_coord_start_index(cell);

		return
			  (this->y_coordinates[unref_cell_y_index + 1] - this->y_coordinates[unref_cell_y_index])
			  / (uint64_t(1) << refinement_level);
	}

	/*!
	Returns the length of given cell in z direction.
	*/
	double get_cell_z_size(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->get_refinement_level(cell) >= 0);
		assert(this->get_refinement_level(cell) <= this->max_refinement_level);

		int refinement_level = this->get_refinement_level(cell);
		uint64_t unref_cell_z_index = this->get_unref_cell_z_coord_start_index(cell);
		return
			  (this->z_coordinates[unref_cell_z_index + 1] - this->z_coordinates[unref_cell_z_index])
			  / (uint64_t(1) << refinement_level);
	}


	/*!
	Returns the center of given cell in x direction.
	*/
	double get_cell_x(const uint64_t cell) const
	{
		 if (cell == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		 }

		if (
			   this->get_refinement_level(cell) < 0
			|| this->get_refinement_level(cell) > this->max_refinement_level
		) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->get_indices(cell);

		uint64_t unref_cell_x_coord_start_index = this->get_unref_cell_x_coord_start_index(cell);
		uint64_t unref_cell_x_index
			= unref_cell_x_coord_start_index * (uint64_t(1) << this->max_refinement_level);

		double
			unref_cell_x_size
				= this->x_coordinates[unref_cell_x_coord_start_index + 1]
				- this->x_coordinates[unref_cell_x_coord_start_index],
			size_of_local_index
				= unref_cell_x_size / (uint64_t(1) << this->max_refinement_level);

		return
			  this->x_coordinates[unref_cell_x_coord_start_index]
			+ size_of_local_index * (indices[0] - unref_cell_x_index)
			+ this->get_cell_x_size(cell) / 2;
	}

	/*!
	Returns the center of given cell in y direction.
	*/
	double get_cell_y(const uint64_t cell) const
	{
		 if (cell == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		 }

		if (this->get_refinement_level(cell) < 0 || this->get_refinement_level(cell) > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->get_indices(cell);

		uint64_t unref_cell_y_coord_start_index = this->get_unref_cell_y_coord_start_index(cell);
		uint64_t unref_cell_y_index
			= unref_cell_y_coord_start_index * (uint64_t(1) << this->max_refinement_level);

		double
			unref_cell_y_size
				= this->y_coordinates[unref_cell_y_coord_start_index + 1]
				- this->y_coordinates[unref_cell_y_coord_start_index],
			size_of_local_index
				= unref_cell_y_size / (uint64_t(1) << this->max_refinement_level);

		return
			  this->y_coordinates[unref_cell_y_coord_start_index]
			+ size_of_local_index * (indices[1] - unref_cell_y_index)
			+ this->get_cell_y_size(cell) / 2;
	}

	/*!
	Returns the center of given cell in z direction.
	*/
	double get_cell_z(const uint64_t cell) const
	{
		 if (cell == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		 }

		if (this->get_refinement_level(cell) < 0 || this->get_refinement_level(cell) > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->get_indices(cell);

		uint64_t unref_cell_z_coord_start_index = this->get_unref_cell_z_coord_start_index(cell);
		uint64_t unref_cell_z_index
			= unref_cell_z_coord_start_index * (uint64_t(1) << this->max_refinement_level);

		double
			unref_cell_z_size
				= this->z_coordinates[unref_cell_z_coord_start_index + 1]
				- this->z_coordinates[unref_cell_z_coord_start_index],
			size_of_local_index
				= unref_cell_z_size / (uint64_t(1) << this->max_refinement_level);

		return
			  this->z_coordinates[unref_cell_z_coord_start_index]
			+ size_of_local_index * (indices[2] - unref_cell_z_index)
			+ this->get_cell_z_size(cell) / 2;
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

		if (x_index > this->length_x * (uint64_t(1) << this->max_refinement_level)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		uint64_t unref_cell_x_coord_start_index
			= x_index / (uint64_t(1) << this->max_refinement_level);

		double unref_cell_x_size = this->x_coordinates[unref_cell_x_coord_start_index + 1] - this->x_coordinates[unref_cell_x_coord_start_index];
		double size_of_local_index = unref_cell_x_size / (uint64_t(1) << this->max_refinement_level);

		return this->x_coordinates[unref_cell_x_coord_start_index]
			+ size_of_local_index * (x_index - unref_cell_x_coord_start_index * (uint64_t(1) << this->max_refinement_level))
			+ size_of_local_index * (uint64_t(1) << (this->max_refinement_level - refinement_level)) / 2;
	}

	/*!
	Returns the center of a cell in y direction of given refinement level at given index.
	*/
	double get_cell_y(const int refinement_level, const uint64_t y_index) const
	{
		if (refinement_level < 0 || refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (y_index > this->length_y * (uint64_t(1) << this->max_refinement_level)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		uint64_t unref_cell_y_coord_start_index
			= y_index / (uint64_t(1) << this->max_refinement_level);

		double unref_cell_y_size = this->y_coordinates[unref_cell_y_coord_start_index + 1] - this->y_coordinates[unref_cell_y_coord_start_index];
		double size_of_local_index = unref_cell_y_size / (uint64_t(1) << this->max_refinement_level);

		return this->y_coordinates[unref_cell_y_coord_start_index]
			+ size_of_local_index * (y_index - unref_cell_y_coord_start_index * (uint64_t(1) << this->max_refinement_level))
			+ size_of_local_index * (uint64_t(1) << (this->max_refinement_level - refinement_level)) / 2;
	}

	/*!
	Returns the center of a cell in z direction of given refinement level at given index.
	*/
	double get_cell_z(const int refinement_level, const uint64_t z_index) const
	{
		if (refinement_level < 0 || refinement_level > this->max_refinement_level) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (z_index > this->length_z * (uint64_t(1) << this->max_refinement_level)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		uint64_t unref_cell_z_coord_start_index
			= z_index / (uint64_t(1) << this->max_refinement_level);

		double unref_cell_z_size = this->z_coordinates[unref_cell_z_coord_start_index + 1] - this->z_coordinates[unref_cell_z_coord_start_index];
		double size_of_local_index = unref_cell_z_size / (uint64_t(1) << this->max_refinement_level);

		return this->z_coordinates[unref_cell_z_coord_start_index]
			+ size_of_local_index * (z_index - unref_cell_z_coord_start_index * (uint64_t(1) << this->max_refinement_level))
			+ size_of_local_index * (uint64_t(1) << (this->max_refinement_level - refinement_level)) / 2;
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
		if (x >= this->get_start_x()
		&& x <= this->get_x_end()) {

			return x;

		} else if (!this->periodic[0]) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_x_end() - this->get_start_x();

			if (x < this->get_start_x()) {

				const double distance = this->get_start_x() - x;
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
		if (y >= this->get_start_y()
		&& y <= this->get_y_end()) {

			return y;

		} else if (!this->periodic[0]) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_y_end() - this->get_start_y();

			if (y < this->get_start_y()) {

				const double distance = this->get_start_y() - y;
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
		if (z >= this->get_start_z()
		&& z <= this->get_z_end()) {

			return z;

		} else if (!this->periodic[0]) {

			return std::numeric_limits<double>::quiet_NaN();

		} else {
			const double grid_size = this->get_z_end() - this->get_start_z();

			if (z < this->get_start_z()) {

				const double distance = this->get_start_z() - z;
				return z + grid_size * ceil(distance/ grid_size);

			} else {

				const double distance = z - this->get_z_end();
				return z - grid_size * ceil(distance/ grid_size);
			}
		}
	}


	/*!
	Returns the x index of given location, starting from 0.

	Returns error_index if given location is outside of the grid and
	the geometry is not periodic in that direction.
	*/
	uint64_t get_x_index_of_coord(double x) const
	{
		x = this->get_real_x(x);

		if (::isnan(x)
		|| x < this->get_start_x()
		|| x > this->get_x_end()) {

			return error_index;

		}

		uint64_t x_coord_start_index = 0;
		while (this->x_coordinates[x_coord_start_index] < x) {
			x_coord_start_index++;
		}
		x_coord_start_index--;

		double x_size_of_index
			= (this->x_coordinates[x_coord_start_index + 1] - this->x_coordinates[x_coord_start_index])
			/ this->get_cell_size_in_indices(1);

		uint64_t index_offset = 0;
		while (this->x_coordinates[x_coord_start_index] + index_offset * x_size_of_index < x) {
			index_offset++;
		}
		index_offset--;

		return x_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;
	}

	/*!
	Returns the y index of given location, starting from 0.

	Returns error_index if given location is outside of the grid and
	the geometry is not periodic in that direction.
	*/
	uint64_t get_y_index_of_coord(double y) const
	{
		y = this->get_real_y(y);

		if (::isnan(y)
		|| y < this->get_start_y()
		|| y > this->get_y_end()) {

			return error_index;

		}

		uint64_t y_coord_start_index = 0;
		while (this->y_coordinates[y_coord_start_index] < y) {
			y_coord_start_index++;
		}
		y_coord_start_index--;

		double y_size_of_index
			= (this->y_coordinates[y_coord_start_index + 1] - this->y_coordinates[y_coord_start_index])
			/ this->get_cell_size_in_indices(1);

		uint64_t index_offset = 0;
		while (this->y_coordinates[y_coord_start_index] + index_offset * y_size_of_index < y) {
			index_offset++;
		}
		index_offset--;

		return y_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;
	}

	/*!
	Returns the z index of given location, starting from 0.

	Returns error_index if given location is outside of the grid and
	the geometry is not periodic in that direction.
	*/
	uint64_t get_z_index_of_coord(double z) const
	{
		z = this->get_real_z(z);

		if (::isnan(z)
		|| z < this->get_start_z()
		|| z > this->get_z_end()) {

			return error_index;

		}

		uint64_t z_coord_start_index = 0;
		while (this->z_coordinates[z_coord_start_index] < z) {
			z_coord_start_index++;
		}
		z_coord_start_index--;

		double z_size_of_index
			= (this->z_coordinates[z_coord_start_index + 1] - this->z_coordinates[z_coord_start_index])
			/ this->get_cell_size_in_indices(1);

		uint64_t index_offset = 0;
		while (this->z_coordinates[z_coord_start_index] + index_offset * z_size_of_index < z) {
			index_offset++;
		}
		index_offset--;

		return z_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;
	}


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

	// periodic[0] == true means that the grid wraps around in x direction
	bool periodic[3];

	/*!
	The coordinates of unrefined cells in respective directions
	First value is the starting point of the grid, the following
	ith value is the end point of the ith unrefined cell
	*/
	std::vector<double> x_coordinates, y_coordinates, z_coordinates;


	/*!
	Returns the x index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.
	*/
	uint64_t get_unref_cell_x_coord_start_index(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->get_refinement_level(cell) >= 0);
		assert(this->get_refinement_level(cell) <= this->max_refinement_level);

		const Types<3>::indices_t indices = this->get_indices(cell);

		return indices[0] / (uint64_t(1) << this->max_refinement_level);
	}

	/*!
	Returns the y index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.
	*/
	uint64_t get_unref_cell_y_coord_start_index(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->get_refinement_level(cell) >= 0);
		assert(this->get_refinement_level(cell) <= this->max_refinement_level);

		const Types<3>::indices_t indices = this->get_indices(cell);

		return indices[1] / (uint64_t(1) << this->max_refinement_level);
	}

	/*!
	Returns the z index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.
	*/
	uint64_t get_unref_cell_z_coord_start_index(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->get_refinement_level(cell) >= 0);
		assert(this->get_refinement_level(cell) <= this->max_refinement_level);

		const Types<3>::indices_t indices = this->get_indices(cell);

		return indices[2] / (uint64_t(1) << this->max_refinement_level);
	}

};	// class

}	// namespace

#endif

