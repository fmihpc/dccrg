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


#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_length.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"


namespace dccrg {


/*!
\brief Geometry class for dccrg with rectangular cuboid cells

A geometry class in which the coordinates of unrefined cells are given
by three vectors of floating points numbers.
The number of values in each vector must be equal to the
length of the grid + 1 in the respective dimension.
*/
class Stretched_Cartesian_Geometry
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
	Calls reset().
	*/
	Stretched_Cartesian_Geometry(
		const Grid_Length& given_length,
		const Mapping& given_mapping,
		const Grid_Topology& given_topology
	) :
		length(given_length),
		mapping(given_mapping),
		topology(given_topology)
	{
		this->reset();
	}

	/*!
	Calls reset().
	*/
	~Stretched_Cartesian_Geometry()
	{
		this->reset();
	}


	/*!
	Sets the geometry of the grid to the following:
		- starting corner at (0, 0, 0)
		- size of cells of refinement level 0 in each dimension: 1
	*/
	void reset()
	{
		this->coordinates[0].resize(this->length.get()[0] + 1);
		for (uint64_t i = 0; i <= this->length.get()[0]; i++) {
			this->coordinates[0][i] = i;
		}

		this->coordinates[1].resize(this->length.get()[1] + 1);
		for (uint64_t i = 0; i <= this->length.get()[1]; i++) {
			this->coordinates[1][i] = i;
		}

		this->coordinates[2].resize(this->length.get()[2] + 1);
		for (uint64_t i = 0; i <= this->length.get()[2]; i++) {
			this->coordinates[2][i] = i;
		}
	}


	/*!
	Returns the geometry of the grid.
	*/
	const boost::array<std::vector<double>, 3>& get() const
	{
		return this->coordinates;
	}


	/*!
	Sets the sizes of the grid's cells of refinement level 0 to given values.

	Ith value in a dimension is the starting point of the Ith cell
	in that dimension. Last value is the ending point of the last cell.
	At least two values must be given for each dimension and all values
	must be strictly increasing.
	Returns true on success.
	Returns false if unsuccessful and does nothing in that case.
	*/
	bool set(const boost::array<std::vector<double>, 3>& given_coordinates)
	{
		if (coordinates[0].size() < 2) {
			std::cerr << "At least two coordinates are required for grid cells in the x direction"
				<< std::endl;
			return false;
		}
		if (coordinates[1].size() < 2) {
			std::cerr << "At least two coordinates are required for grid cells in the y direction"
				<< std::endl;
			return false;
		}
		if (coordinates[2].size() < 2) {
			std::cerr << "At least two coordinates are required for grid cells in the z direction"
				<< std::endl;
			return false;
		}

		for (uint64_t i = 0; i < coordinates[0].size() - 1; i++) {
			if (coordinates[0][i] >= coordinates[0][i + 1]) {
				std::cerr << "Coordinates in the x direction must be strictly increasing"
					<< std::endl;
				return false;
			}
		}
		for (uint64_t i = 0; i < coordinates[1].size() - 1; i++) {
			if (coordinates[1][i] >= coordinates[1][i + 1]) {
				std::cerr << "Coordinates in the y direction must be strictly increasing"
					<< std::endl;
				return false;
			}
		}
		for (uint64_t i = 0; i < coordinates[2].size() - 1; i++) {
			if (coordinates[2][i] >= coordinates[2][i + 1]) {
				std::cerr << "Coordinates in the z direction must be strictly increasing"
					<< std::endl;
				return false;
			}
		}

		for (size_t i = 0; i < this->length.get().size(); i++) {
			if (coordinates[i].size() != this->length.get()[i] + 1) {
				std::cerr << "Number of values in dimension " << i
					<< " must be length of the grid + 1 (" << this->length.get()[i] + 1
					<< ") but is " << coordinates[i].size()
					<< std::endl;
				return false;
			}
		}

		this->coordinates[0].clear();
		this->coordinates[0].reserve(given_coordinates[0].size());
		this->coordinates[0].insert(
			this->coordinates[0].begin(),
			given_coordinates[0].begin(),
			given_coordinates[0].end()
		);

		this->coordinates[1].clear();
		this->coordinates[1].reserve(given_coordinates[1].size());
		this->coordinates[1].insert(
			this->coordinates[1].begin(),
			given_coordinates[1].begin(),
			given_coordinates[1].end()
		);

		this->coordinates[2].clear();
		this->coordinates[2].reserve(given_coordinates[2].size());
		this->coordinates[2].insert(
			this->coordinates[2].begin(),
			given_coordinates[2].begin(),
			given_coordinates[2].end()
		);

		return true;
	}

	/*!
	Sets the same geometry as in the given one.
	*/
	bool set(const Stretched_Cartesian_Geometry& other)
	{
		for (size_t dimension = 0; dimension < this->length.get().size(); dimension++) {
			if (!this->length.get()[dimension] == other.length.get()[dimension]) {
				std::cerr << "Geometries must have the same grid length but have "
					<< this->length.get()[dimension] << " and "
					<< other.length.get()[dimension] << " in dimension " << dimension
					<< std::endl;
				return false;
			}
		}

		return this->set(other.get());
	}

	/*!
	Sets the same geometry as in the given one.
	*/
	bool set(const Cartesian_Geometry& other)
	{
		for (size_t dimension = 0; dimension < this->length.get().size(); dimension++) {
			if (!this->length.get()[dimension] == other.length.get()[dimension]) {
				std::cerr << "Geometries must have the same grid length but have "
					<< this->length.get()[dimension] << " and "
					<< other.length.get()[dimension] << " in dimension " << dimension
					<< std::endl;
				return false;
			}
		}

		const boost::array<uint64_t, 3> grid_length = {{
			other.length.get()[0],
			other.length.get()[1],
			other.length.get()[2]
		}};

		const boost::array<double, 3>
			start = {{
				other.get_start_x(),
				other.get_start_y(),
				other.get_start_z()
			}},
			cell_length = {{
				other.get_unrefined_cell_length_x(),
				other.get_unrefined_cell_length_y(),
				other.get_unrefined_cell_length_z()
			}};

		boost::array<std::vector<double>, 3> coordinates;
		for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
			for (uint64_t i = 0; i <= grid_length[dimension]; i++) {
				coordinates[dimension].push_back(
					start[dimension] + double(i) * cell_length[dimension]
				);
			}
		}

		return this->set(coordinates);
	}


	/*!
	Returns the starting corner of the grid in x direction.
	*/
	double get_start_x() const
	{
		return this->coordinates[0][0];
	}

	/*!
	Returns the starting corner of the grid in y direction.
	*/
	double get_start_y() const
	{
		return this->coordinates[1][0];
	}

	/*!
	Returns the starting corner of the grid in z direction.
	*/
	double get_start_z() const
	{
		return this->coordinates[2][0];
	}


	/*!
	Returns the end corner of the grid in x direction.
	*/
	double get_end_x() const
	{
		return this->coordinates[0][this->length.get()[0]];
	}

	/*!
	Returns the end corner of the grid in y direction.
	*/
	double get_end_y() const
	{
		return this->coordinates[1][this->length.get()[1]];
	}

	/*!
	Returns the end corner of the grid in z direction.
	*/
	double get_end_z() const
	{
		return this->coordinates[2][this->length.get()[2]];
	}


	/*!
	Returns the length of given cell in x direction.
	*/
	double get_cell_length_x(const uint64_t cell) const
	{
		assert(cell != error_cell);
		assert(this->mapping.get_refinement_level(cell) >= 0);
		assert(this->mapping.get_refinement_level(cell) <= this->mapping.get_maximum_refinement_level());

		int refinement_level = this->mapping.get_refinement_level(cell);
		uint64_t unref_cell_x_index = this->get_unref_cell_x_coord_start_index(cell);

		return
			  (this->coordinates[0][unref_cell_x_index + 1] - this->coordinates[0][unref_cell_x_index])
			/ (uint64_t(1) << refinement_level);
	}

	/*!
	Returns the length of given cell in y direction.
	*/
	double get_cell_length_y(const uint64_t cell) const
	{
		assert(cell != error_cell);
		assert(this->mapping.get_refinement_level(cell) >= 0);
		assert(this->mapping.get_refinement_level(cell) <= this->mapping.get_maximum_refinement_level());

		int refinement_level = this->mapping.get_refinement_level(cell);
		uint64_t unref_cell_y_index = this->get_unref_cell_y_coord_start_index(cell);

		return
			  (this->coordinates[1][unref_cell_y_index + 1] - this->coordinates[1][unref_cell_y_index])
			  / (uint64_t(1) << refinement_level);
	}

	/*!
	Returns the length of given cell in z direction.
	*/
	double get_cell_length_z(const uint64_t cell) const
	{
		assert(cell != error_cell);
		assert(this->mapping.get_refinement_level(cell) >= 0);
		assert(this->mapping.get_refinement_level(cell) <= this->mapping.get_maximum_refinement_level());

		int refinement_level = this->mapping.get_refinement_level(cell);
		uint64_t unref_cell_z_index = this->get_unref_cell_z_coord_start_index(cell);
		return
			  (this->coordinates[2][unref_cell_z_index + 1] - this->coordinates[2][unref_cell_z_index])
			  / (uint64_t(1) << refinement_level);
	}


	/*!
	Returns the center of given cell in x direction.
	*/
	double get_cell_x(const uint64_t cell) const
	{
		 if (cell == error_cell) {
			return std::numeric_limits<double>::quiet_NaN();
		 }

		if (
			   this->mapping.get_refinement_level(cell) < 0
			|| this->mapping.get_refinement_level(cell) > this->mapping.get_maximum_refinement_level()
		) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		const uint64_t
			unref_cell_x_coord_start_index
				= this->get_unref_cell_x_coord_start_index(cell),
			unref_cell_x_index
				= unref_cell_x_coord_start_index
				* (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		const double
			unref_cell_length_x
				= this->coordinates[0][unref_cell_x_coord_start_index + 1]
				- this->coordinates[0][unref_cell_x_coord_start_index],
			size_of_local_index
				= unref_cell_length_x
				/ (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		return
			  this->coordinates[0][unref_cell_x_coord_start_index]
			+ size_of_local_index * (indices[0] - unref_cell_x_index)
			+ this->get_cell_length_x(cell) / 2;
	}

	/*!
	Returns the center of given cell in y direction.
	*/
	double get_cell_y(const uint64_t cell) const
	{
		 if (cell == error_cell) {
			return std::numeric_limits<double>::quiet_NaN();
		 }

		if (this->mapping.get_refinement_level(cell) < 0
		|| this->mapping.get_refinement_level(cell) > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		const uint64_t
			unref_cell_y_coord_start_index
				= this->get_unref_cell_y_coord_start_index(cell),
			unref_cell_y_index
				= unref_cell_y_coord_start_index
				* (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		const double
			unref_cell_length_y
				= this->coordinates[1][unref_cell_y_coord_start_index + 1]
				- this->coordinates[1][unref_cell_y_coord_start_index],
			size_of_local_index
				= unref_cell_length_y
				/ (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		return
			  this->coordinates[1][unref_cell_y_coord_start_index]
			+ size_of_local_index * (indices[1] - unref_cell_y_index)
			+ this->get_cell_length_y(cell) / 2;
	}

	/*!
	Returns the center of given cell in z direction.
	*/
	double get_cell_z(const uint64_t cell) const
	{
		 if (cell == error_cell) {
			return std::numeric_limits<double>::quiet_NaN();
		 }

		if (this->mapping.get_refinement_level(cell) < 0
		|| this->mapping.get_refinement_level(cell) > this->mapping.get_maximum_refinement_level()) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		const uint64_t
			unref_cell_z_coord_start_index
				= this->get_unref_cell_z_coord_start_index(cell),
			unref_cell_z_index
				= unref_cell_z_coord_start_index
				* (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		const double
			unref_cell_length_z
				= this->coordinates[2][unref_cell_z_coord_start_index + 1]
				- this->coordinates[2][unref_cell_z_coord_start_index],
			size_of_local_index
				= unref_cell_length_z / (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		return
			  this->coordinates[2][unref_cell_z_coord_start_index]
			+ size_of_local_index * (indices[2] - unref_cell_z_index)
			+ this->get_cell_length_z(cell) / 2;
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

		uint64_t unref_cell_x_coord_start_index
			= x_index / (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		double
			unref_cell_length_x
				= this->coordinates[0][unref_cell_x_coord_start_index + 1]
				- this->coordinates[0][unref_cell_x_coord_start_index],
			size_of_local_index
				= unref_cell_length_x
				/ (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		return
			this->coordinates[0][unref_cell_x_coord_start_index]
			+ size_of_local_index * (x_index - unref_cell_x_coord_start_index * (uint64_t(1) << this->mapping.get_maximum_refinement_level()))
			+ size_of_local_index * (uint64_t(1) << (this->mapping.get_maximum_refinement_level() - refinement_level)) / 2;
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

		uint64_t unref_cell_y_coord_start_index
			= y_index / (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		double
			unref_cell_length_y
				= this->coordinates[1][unref_cell_y_coord_start_index + 1]
				- this->coordinates[1][unref_cell_y_coord_start_index],
			size_of_local_index
				= unref_cell_length_y
				/ (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		return
			this->coordinates[1][unref_cell_y_coord_start_index]
			+ size_of_local_index * (y_index - unref_cell_y_coord_start_index * (uint64_t(1) << this->mapping.get_maximum_refinement_level()))
			+ size_of_local_index * (uint64_t(1) << (this->mapping.get_maximum_refinement_level() - refinement_level)) / 2;
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

		uint64_t unref_cell_z_coord_start_index
			= z_index / (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		double
			unref_cell_length_z
				= this->coordinates[2][unref_cell_z_coord_start_index + 1]
				- this->coordinates[2][unref_cell_z_coord_start_index],
			size_of_local_index
				= unref_cell_length_z
				/ (uint64_t(1) << this->mapping.get_maximum_refinement_level());

		return
			this->coordinates[2][unref_cell_z_coord_start_index]
			+ size_of_local_index * (z_index - unref_cell_z_coord_start_index * (uint64_t(1) << this->mapping.get_maximum_refinement_level()))
			+ size_of_local_index * (uint64_t(1) << (this->mapping.get_maximum_refinement_level() - refinement_level)) / 2;
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

		} else if (!this->topology.is_periodic(0)) {

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

		} else if (!this->topology.is_periodic(0)) {

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

	Returns error_index if given location is outside of the grid and
	the geometry is not periodic in that direction.
	*/
	uint64_t get_x_index_of_coord(double x) const
	{
		x = this->get_real_x(x);

		if (::isnan(x)
		|| x < this->get_start_x()
		|| x > this->get_end_x()) {

			return error_index;

		}

		uint64_t x_coord_start_index = 0;
		while (this->coordinates[0][x_coord_start_index] < x) {
			x_coord_start_index++;
		}
		x_coord_start_index--;

		double length_x_of_index
			= (this->coordinates[0][x_coord_start_index + 1] - this->coordinates[0][x_coord_start_index])
			/ this->mapping.get_cell_length_in_indices(1);

		uint64_t index_offset = 0;
		while (this->coordinates[0][x_coord_start_index] + index_offset * length_x_of_index < x) {
			index_offset++;
		}
		index_offset--;

		return x_coord_start_index * this->mapping.get_cell_length_in_indices(1) + index_offset;
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
		|| y > this->get_end_y()) {

			return error_index;

		}

		uint64_t y_coord_start_index = 0;
		while (this->coordinates[1][y_coord_start_index] < y) {
			y_coord_start_index++;
		}
		y_coord_start_index--;

		double length_y_of_index
			= (this->coordinates[1][y_coord_start_index + 1] - this->coordinates[1][y_coord_start_index])
			/ this->mapping.get_cell_length_in_indices(1);

		uint64_t index_offset = 0;
		while (this->coordinates[1][y_coord_start_index] + index_offset * length_y_of_index < y) {
			index_offset++;
		}
		index_offset--;

		return y_coord_start_index * this->mapping.get_cell_length_in_indices(1) + index_offset;
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
		|| z > this->get_end_z()) {

			return error_index;

		}

		uint64_t z_coord_start_index = 0;
		while (this->coordinates[2][z_coord_start_index] < z) {
			z_coord_start_index++;
		}
		z_coord_start_index--;

		double length_z_of_index
			= (this->coordinates[2][z_coord_start_index + 1] - this->coordinates[2][z_coord_start_index])
			/ this->mapping.get_cell_length_in_indices(1);

		uint64_t index_offset = 0;
		while (this->coordinates[2][z_coord_start_index] + index_offset * length_z_of_index < z) {
			index_offset++;
		}
		index_offset--;

		return z_coord_start_index * this->mapping.get_cell_length_in_indices(1) + index_offset;
	}


	/*!
	Returns the x index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.
	*/
	uint64_t get_unref_cell_x_coord_start_index(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->mapping.get_refinement_level(cell) >= 0);
		assert(this->mapping.get_refinement_level(cell) <= this->mapping.get_maximum_refinement_level());

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		return indices[0] / (uint64_t(1) << this->mapping.get_maximum_refinement_level());
	}

	/*!
	Returns the y index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.
	*/
	uint64_t get_unref_cell_y_coord_start_index(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->mapping.get_refinement_level(cell) >= 0);
		assert(this->mapping.get_refinement_level(cell) <= this->mapping.get_maximum_refinement_level());

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		return indices[1] / (uint64_t(1) << this->mapping.get_maximum_refinement_level());
	}

	/*!
	Returns the z index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.
	*/
	uint64_t get_unref_cell_z_coord_start_index(const uint64_t cell) const
	{
		assert(cell > 0);
		assert(this->mapping.get_refinement_level(cell) >= 0);
		assert(this->mapping.get_refinement_level(cell) <= this->mapping.get_maximum_refinement_level());

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		return indices[2] / (uint64_t(1) << this->mapping.get_maximum_refinement_level());
	}



private:

	/*!
	The coordinates of unrefined cells in respective directions
	First value is the starting point of the grid, the following
	ith value is the end point of the ith unrefined cell
	*/
	boost::array<std::vector<double>, 3> coordinates;

};	// class

}	// namespace

#endif

