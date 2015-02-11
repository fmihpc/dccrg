/*
Dccrg class for a cartesian geometry in which cells are cubes.

Copyright 2009, 2010, 2011, 2012, 2013,
2014, 2015 Finnish Meteorological Institute

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
#include "mpi.h"
#include "stdint.h"
#include "vector"

#include "dccrg_length.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_mpi_support.hpp"
#include "dccrg_topology.hpp"
#include "dccrg_types.hpp"


namespace dccrg {


/*!
\brief Parameters required for the Cartesian_Geometry class of dccrg.

Is given to the Cartesian_Geometry::set() function.
*/
class Cartesian_Geometry_Parameters
{
public:
	std::array<double, 3>
		//! starting coordinate of the grid
		start,
		//! length of cells of refinement level 0
		level_0_cell_length;

	Cartesian_Geometry_Parameters()
	{
		this->start[0] =
		this->start[1] =
		this->start[2] = 0;

		this->level_0_cell_length[0] =
		this->level_0_cell_length[1] =
		this->level_0_cell_length[2] = 1;
	}

	Cartesian_Geometry_Parameters(
		const std::array<double, 3>& given_start,
		const std::array<double, 3>& given_level_0_cell_length
	) {
		this->start = given_start;

		for (size_t dimension = 0; dimension < this->level_0_cell_length.size(); dimension++) {
			if (given_level_0_cell_length[dimension] <= 0) {
				std::cerr << "Cell length in dimension " << dimension
					<< " must be > 0, but is "
					<< given_level_0_cell_length[dimension]
					<< std::endl;
				abort();
			}
		}
		this->level_0_cell_length = given_level_0_cell_length;
	}
};


/*!
\brief Geometry class for dccrg with cubic cells

A geometry class in which the sizes of cells of
refinement level 0 are given by three floating points numbers.
*/
class Cartesian_Geometry
{

public:

	/*!
	Unique identifier of this geometry class, used when
	storing the geometry to a file.
	*/
	static const int geometry_id = 1;

	/*!
	Parameter type that is defined by every geometry class
	and used to refer to their own parameter type.
	*/
	typedef Cartesian_Geometry_Parameters Parameters;

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

	\see
	Grid_Length Mapping Grid_Topology
	*/
	Cartesian_Geometry(
		const Grid_Length& given_length,
		const Mapping& given_mapping,
		const Grid_Topology& given_topology
	) :
		length(given_length),
		mapping(given_mapping),
		topology(given_topology)
	{}

	/*!
	Sets the geometry of the grid to the following:
		- starting corner at (0, 0, 0)
		- size of unrefined cells in each direction: 1
	*/
	~Cartesian_Geometry()
	{
		this->parameters.start[0] = 0;
		this->parameters.start[1] = 0;
		this->parameters.start[2] = 0;

		this->parameters.level_0_cell_length[0] = 1;
		this->parameters.level_0_cell_length[1] = 1;
		this->parameters.level_0_cell_length[2] = 1;
	}


	/*!
	Returns the parameters of the grid's geometry.
	*/
	const Parameters& get() const
	{
		return this->parameters;
	}


	/*!
	Sets the geometry of the grid to given values.

	Returns true on success. On failure returns false
	and has no effect.
	*/
	bool set(const Parameters& given_parameters)
	{
		for (size_t
			dimension = 0;
			dimension < this->parameters.level_0_cell_length.size();
			dimension++
		) {
			if (given_parameters.level_0_cell_length[dimension] <= 0) {
				std::cerr << "Cell length in dimension " << dimension
					<< " must be > 0, but is "
					<< given_parameters.level_0_cell_length[dimension]
					<< std::endl;
				return false;
			}
		}

		// TODO: check that all cell coordinates fit into a double

		this->parameters = given_parameters;

		return true;
	}


	/*!
	Sets the same geometry as in the given one.
	*/
	bool set(const Cartesian_Geometry& other)
	{
		return this->set(other.get());
	}


	/*!
	Returns the starting corner of the grid.

	Starting corner is defined as the corner
	with minimum value of the coordinate in each
	dimension.
	*/
	std::array<double, 3> get_start() const
	{
		const std::array<double, 3> ret_val = {{
			this->parameters.start[0],
			this->parameters.start[1],
			this->parameters.start[2]
		}};
		return ret_val;
	}


	/*!
	Returns the end corner of the grid.

	End corner is defined as the corner with
	maximum value of the coordinate in each
	dimension.
	*/
	std::array<double, 3> get_end() const
	{
		// in cells of refinement level 0
		const std::array<uint64_t, 3> length = this->length.get();

		const std::array<double, 3>
			grid_start =  this->get_start(),
			ret_val = {{
				grid_start[0]
				+ double(length[0])
					* this->parameters.level_0_cell_length[0],

				grid_start[1]
				+ double(length[1])
					* this->parameters.level_0_cell_length[1],

				grid_start[2]
				+ double(length[2])
					* this->parameters.level_0_cell_length[2],
			}};

		return ret_val;
	}


	/*!
	Returns the length of cells of refinement level 0.
	*/
	std::array<double, 3> get_level_0_cell_length() const
	{
		return this->parameters.level_0_cell_length;
	}


	/*!
	Returns the length of given cell.

	A quiet NaN is returned if the given error_cell,
	or the given cell cannot exist in the grid.
	*/
	std::array<double, 3> get_length(const uint64_t cell) const
	{
		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {

			const std::array<double, 3> error_val = {{
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN()
			}};

			return error_val;
		}

		const double scaling_factor = 1.0 / double(uint64_t(1) << refinement_level);
		const std::array<double, 3> ret_val = {{
			this->parameters.level_0_cell_length[0] * scaling_factor,
			this->parameters.level_0_cell_length[1] * scaling_factor,
			this->parameters.level_0_cell_length[2] * scaling_factor
		}};

		return ret_val;
	}


	/*!
	Returns the center of given cell.

	A quiet NaN is returned if given error_cell,
	or the given cell cannot exist in the grid.
	*/
	std::array<double, 3> get_center(const uint64_t cell) const
	{

		const int refinement_level = this->mapping.get_refinement_level(cell);

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {

			const std::array<double, 3> error_val = {{
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN()
			}};

			return error_val;
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);
		const int max_ref_lvl = this->mapping.get_maximum_refinement_level();

		const std::array<double, 3>
			grid_start = this->get_start(),
			cell_length = this->get_length(cell),
			level_0_cell_length = this->get_level_0_cell_length(),
			ret_val = {{
				grid_start[0]
				+ double(indices[0])
					* level_0_cell_length[0]
					/ double(uint64_t(1) << max_ref_lvl)
				+ cell_length[0] / 2,

				grid_start[1]
				+ double(indices[1])
					* level_0_cell_length[1]
					/ double(uint64_t(1) << max_ref_lvl)
				+ cell_length[1] / 2,

				grid_start[2]
				+ double(indices[2])
					* level_0_cell_length[2]
					/ double(uint64_t(1) << max_ref_lvl)
				+ cell_length[2] / 2
			}};

		return ret_val;
	}


	/*!
	Returns the cell's corner closest to starting corner of the grid.

	In other words if the cell occupies the range:
	\verbatim
	[x1..x2],
	[y1..y2],
	...
	\endverbatim
	returns (x1, y1, ...).
	*/
	std::array<double, 3> get_min(const uint64_t cell) const
	{
		const std::array<double, 3>
			center = this->get_center(cell),
			length = this->get_length(cell),
			ret_val = {{
				center[0] - length[0] / 2,
				center[1] - length[1] / 2,
				center[2] - length[2] / 2
			}};

		return ret_val;
	}


	/*!
	Returns the cell's corner furthest from the starting corner of the grid.

	In other words if the cell occupies the range:
	\verbatim
	[x1..x2],
	[y1..y2],
	...
	\endverbatim
	returns (x2, y2, ...).
	*/
	std::array<double, 3> get_max(const uint64_t cell) const
	{
		const std::array<double, 3>
			center = this->get_center(cell),
			length = this->get_length(cell),
			ret_val = {{
				center[0] + length[0] / 2,
				center[1] + length[1] / 2,
				center[2] + length[2] / 2
			}};

		return ret_val;
	}


	/*!
	Returns the center of a cell of given refinement level at given index.

	A quiet NaN is returned if given error_cell,
	or the given cell cannot exist in the grid.
	*/
	std::array<double, 3> get_center(
		const int refinement_level,
		const Types<3>::indices_t index
	) const {

		const std::array<double, 3> error_val = {{
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		}};

		if (refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return error_val;
		}

		const uint64_t index_scaling_factor
			= uint64_t(1) << this->mapping.get_maximum_refinement_level();

		const Types<3>::indices_t max_index = {{
			this->length.get()[0] * index_scaling_factor,
			this->length.get()[1] * index_scaling_factor,
			this->length.get()[2] * index_scaling_factor
		}};

		if (
			index[0] > max_index[0]
			|| index[1] > max_index[1]
			|| index[2] > max_index[2]
		) {
			return error_val;
		}

		const double
			coordinate_scaling_factor = 1.0 / double(index_scaling_factor),
			cell_offset_scaling_factor = 1.0 / double(uint64_t(1) << refinement_level) / 2;

		const std::array<double, 3>
			grid_start =  this->get_start(),
			ret_val = {{
				grid_start[0]
				+ double(index[0])
					* this->parameters.level_0_cell_length[0]
					* coordinate_scaling_factor
				+ this->parameters.level_0_cell_length[0]
					* cell_offset_scaling_factor,

				grid_start[1]
				+ double(index[1])
					* this->parameters.level_0_cell_length[1]
					* coordinate_scaling_factor
				+ this->parameters.level_0_cell_length[1]
					* cell_offset_scaling_factor,

				grid_start[2]
				+ double(index[2])
					* this->parameters.level_0_cell_length[2]
					* coordinate_scaling_factor
				+ this->parameters.level_0_cell_length[2]
					* cell_offset_scaling_factor
			}};

		return ret_val;
	}


	/*!
	Returns a cell of given refinement level at given location.

	Returns error_cell if given a location outside of the current
	grid either in coordinate or refinement level.
	*/
	uint64_t get_cell(
		const int refinement_level,
		const std::array<double, 3>& coordinate
	) const {
		if (refinement_level < 0
		|| refinement_level > this->mapping.get_maximum_refinement_level()) {
			return error_cell;
		}

		const Types<3>::indices_t indices = this->get_indices(coordinate);

		return this->mapping.get_cell_from_indices(indices, refinement_level);
	}


	/*!
	Returns the real value of given coordinate in this geometry.

	For each dimension:
	Returns given coordinate if it is inside this geometry.
	Returns a quiet NaN if this geometry is not periodic in
	the dimension and given coordinate is outside of the geometry.
	If this geometry is periodic in the dimension returns a value
	inside the geometry that is at the same location in the
	geometry as given coordinate.
	*/
	std::array<double, 3> get_real_coordinate(
		const std::array<double, 3>& given_coordinate
	) const {
		const std::array<double, 3>
			start = this->get_start(),
			end = this->get_end();

		std::array<double, 3> ret_val = {{
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN(),
			std::numeric_limits<double>::quiet_NaN()
		}};

		for (size_t dimension = 0; dimension < given_coordinate.size(); dimension++) {

			if (given_coordinate[dimension] >= start[dimension]
			&& given_coordinate[dimension] <= end[dimension]) {

				ret_val[dimension] = given_coordinate[dimension];

			} else if (this->topology.is_periodic(dimension)) {

				const double grid_length = end[dimension] - start[dimension];

				if (given_coordinate[dimension] < start[dimension]) {

					const double distance = start[dimension] - given_coordinate[dimension];

					ret_val[dimension]
						= given_coordinate[dimension]
						+ grid_length * ceil(distance / grid_length);

				} else {

					const double distance = given_coordinate[dimension] - end[dimension];

					ret_val[dimension]
						= given_coordinate[dimension]
						- grid_length * ceil(distance / grid_length);
				}
			}
		}

		return ret_val;
	}


	/*!
	Returns the index (starting from 0) of given coordinate in each dimension.

	Returns error_index for coordinates outside of the grid in a particular
	dimension if the grid is not periodic in that dimension.
	*/
	Types<3>::indices_t get_indices(const std::array<double, 3>& given_coordinate) const
	{
		Types<3>::indices_t ret_val = {{
			error_index,
			error_index,
			error_index
		}};

		const std::array<double, 3>
			grid_start = this->get_start(),
			grid_end = this->get_end(),
			coordinate = this->get_real_coordinate(given_coordinate),
			level_0_cell_length = this->get_level_0_cell_length();


		for (size_t dimension = 0; dimension < given_coordinate.size(); dimension++) {

			if (coordinate[dimension] >= grid_start[dimension]
			&& coordinate[dimension] <= grid_end[dimension]) {

				ret_val[dimension] = uint64_t(
					floor(
						(coordinate[dimension] - grid_start[dimension])
						/ (level_0_cell_length[dimension]
							/ double(uint64_t(1) << this->mapping.get_maximum_refinement_level())
						)
					)
				);

			}
		}

		return ret_val;
	}


	/*!
	Writes the geometry into given open file starting at given offset.

	Returns true on success, false otherwise.

	The number of bytes written by this function can be obtained
	from geometry_data_size().
	*/
	bool write(MPI_File file, MPI_Offset offset) const
	{
		int ret_val = -1;

		const int temp_id = Cartesian_Geometry::geometry_id;
		ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) &temp_id,
			1,
			MPI_INT,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write geometry id to file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += sizeof(int);

		ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) this->parameters.start.data(),
			3,
			MPI_DOUBLE,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write geometry start to file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += 3 * sizeof(double);

		ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) this->parameters.level_0_cell_length.data(),
			3,
			MPI_DOUBLE,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write level 0 cell length to file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		return true;
	}


	/*!
	Reads the geometry from given open file starting at given offset.

	Returns true on success, false otherwise.
	Unless you know what you're doing must be called by all processes
	with identical arguments.
	*/
	bool read(MPI_File file, MPI_Offset offset)
	{
		int
			read_geometry_id = Cartesian_Geometry::geometry_id + 1,
			ret_val = -1;

		ret_val = MPI_File_read_at(
			file,
			offset,
			(void*) &read_geometry_id,
			1,
			MPI_INT,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read geometry id from file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += sizeof(int);

		// TODO: don't error out if given No_Geometry
		if (read_geometry_id != Cartesian_Geometry::geometry_id) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Wrong geometry: " << read_geometry_id
				<< ", should be " << Cartesian_Geometry::geometry_id
				<< std::endl;
			return false;
		}


		Parameters read_parameters;

		ret_val = MPI_File_read_at(
			file,
			offset,
			(void*) read_parameters.start.data(),
			3,
			MPI_DOUBLE,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read geometry start from file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += 3 * sizeof(double);

		ret_val = MPI_File_read_at(
			file,
			offset,
			(void*) read_parameters.level_0_cell_length.data(),
			3,
			MPI_DOUBLE,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read level 0 cell length from file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		if (!this->set(read_parameters)) {
			return false;
		}

		return true;
	}


	/*!
	Returns the number of bytes that will be required / was required for geometry data.
	*/
	size_t data_size() const
	{
		return sizeof(int) + 6 * sizeof(double);
	}



private:

	Parameters parameters;

};	// class

}	// namespace

#endif

