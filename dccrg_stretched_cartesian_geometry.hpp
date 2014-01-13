/*
The mapping logic between cells' and their geometry (location, size, etc.)
in dccrg, stretched cartesian version.

Copyright 2009, 2010, 2011, 2012, 2013,
2014 Finnish Meteorological Institute

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


#include "cmath"
#include "cstdint"
#include "cstdlib"
#include "iostream"
#include "limits"
#include "vector"


#include "dccrg_cartesian_geometry.hpp"
#include "dccrg_length.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_topology.hpp"


namespace dccrg {


/*!
\brief Parameters required for the Stretched_Cartesian_Geometry class of dccrg.

Is given to the Stretched_Cartesian_Geometry::set() function.
*/
class Stretched_Cartesian_Geometry_Parameters
{
public:
	/*!
	The coordinates of unrefined cells in respective dimensions.

	First value is the starting point of the grid, the following
	ith value is the end point of the ith cell of refinement level 0.
	*/
	std::array<std::vector<double>, 3> coordinates;
};


/*!
\brief Geometry class for dccrg with rectangular cuboid cells

A geometry class in which the coordinates of cells of
refinement level 0 are given by three vectors of floating
points numbers. The number of values in each vector must
be equal to the length of the grid + 1 in the respective dimension.
*/
class Stretched_Cartesian_Geometry
{

public:

	/*!
	Unique identifier of this geometry class, used when
	storing the geometry to a file.
	*/
	static const int geometry_id = 2;

	/*!
	Parameter type that is defined by every geometry class
	and used to refer to their own parameter type.
	*/
	typedef Stretched_Cartesian_Geometry_Parameters Parameters;

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
		for (size_t
			dimension = 0;
			dimension < this->parameters.coordinates.size();
			dimension++
		) {
			this->parameters.coordinates[dimension].resize(this->length.get()[dimension] + 1);
			for (uint64_t i = 0; i <= this->length.get()[dimension]; i++) {
				this->parameters.coordinates[dimension][i] = i;
			}
		}
	}


	/*!
	Returns the geometry of the grid.
	*/
	const Parameters& get() const
	{
		return this->parameters;
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
	bool set(const Parameters& given_parameters)
	{
		for (size_t
			dimension = 0;
			dimension < given_parameters.coordinates.size();
			dimension++
		) {
			if (given_parameters.coordinates[dimension].size() < 2) {
				std::cerr << "At least two coordinates are required for grid cells in the "
					<< dimension << " dimension"
					<< std::endl;
				return false;
			}

			if (given_parameters.coordinates[dimension].size() != this->length.get()[dimension] + 1) {
				std::cerr << "Number of values in dimension " << dimension
					<< " must be length of the grid + 1 (" << this->length.get()[dimension] + 1
					<< ") but is " << given_parameters.coordinates[dimension].size()
					<< std::endl;
				return false;
			}

			for (uint64_t i = 0; i < given_parameters.coordinates[dimension].size() - 1; i++) {
				if (
					given_parameters.coordinates[dimension][i]
					>= given_parameters.coordinates[dimension][i + 1]
				) {
					std::cerr << "Coordinates in the " << dimension
						<< " dimension must be strictly increasing"
						<< std::endl;
					return false;
				}
			}
		}

		this->parameters = given_parameters;

		return true;
	}

	/*!
	Sets the same geometry as in the given one.
	*/
	bool set(const Stretched_Cartesian_Geometry& other)
	{
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

		const std::array<uint64_t, 3> grid_length = other.length.get();

		const std::array<double, 3>
			start = other.get_start(),
			cell_length = other.get_level_0_cell_length();

		Parameters parameters;
		for (size_t dimension = 0; dimension < grid_length.size(); dimension++) {
			for (uint64_t i = 0; i <= grid_length[dimension]; i++) {
				parameters.coordinates[dimension].push_back(
					start[dimension] + double(i) * cell_length[dimension]
				);
			}
		}

		return this->set(parameters);
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
			this->parameters.coordinates[0][0],
			this->parameters.coordinates[1][0],
			this->parameters.coordinates[2][0]
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
		const std::array<uint64_t, 3> length = this->length.get();
		const std::array<double, 3> ret_val = {{
			this->parameters.coordinates[0][length[0]],
			this->parameters.coordinates[1][length[1]],
			this->parameters.coordinates[2][length[2]]
		}};

		return ret_val;
	}


	/*!
	Returns the length of given cell.

	A quiet NaN is returned if given error_cell,
	or the given cell cannot exist in the grid.
	*/
	std::array<double, 3> get_length(const uint64_t cell) const
	{
		const int
			refinement_level = this->mapping.get_refinement_level(cell),
			max_ref_lvl = this->mapping.get_maximum_refinement_level();

		if (cell == error_cell
		|| refinement_level < 0
		|| refinement_level > max_ref_lvl) {

			const std::array<double, 3> error_val = {{
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN()
			}};

			return error_val;
		}

		const std::array<uint64_t, 3> coord_start_indices
			= this->get_level_0_cell_coord_start_index(cell);

		const uint64_t length_in_indices = uint64_t(1) << refinement_level;

		const std::array<double, 3> ret_val = {{
			(this->parameters.coordinates[0][coord_start_indices[0] + 1]
				- this->parameters.coordinates[0][coord_start_indices[0]])
			/ length_in_indices,

			(this->parameters.coordinates[1][coord_start_indices[1] + 1]
				- this->parameters.coordinates[1][coord_start_indices[1]])
			/ length_in_indices,

			(this->parameters.coordinates[2][coord_start_indices[2] + 1]
				- this->parameters.coordinates[2][coord_start_indices[2]])
			/ length_in_indices
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

		const int
			ref_lvl = this->mapping.get_refinement_level(cell),
			max_ref_lvl = this->mapping.get_maximum_refinement_level();

		if (cell == error_cell
		|| ref_lvl < 0
		|| ref_lvl > max_ref_lvl) {

			const std::array<double, 3> error_val = {{
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN(),
				std::numeric_limits<double>::quiet_NaN()
			}};

			return error_val;
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);

		const std::array<uint64_t, 3> coord_start_indices
			= this->get_level_0_cell_coord_start_index(cell);

		const uint64_t level_0_length_in_indices = uint64_t(1) << max_ref_lvl;

		const std::array<double, 3>
			length_of_index = {{
				(this->parameters.coordinates[0][coord_start_indices[0] + 1]
					- this->parameters.coordinates[0][coord_start_indices[0]])
				/ level_0_length_in_indices,

				(this->parameters.coordinates[1][coord_start_indices[1] + 1]
					- this->parameters.coordinates[1][coord_start_indices[1]])
				/ level_0_length_in_indices,

				(this->parameters.coordinates[2][coord_start_indices[2] + 1]
					- this->parameters.coordinates[2][coord_start_indices[2]])
				/ level_0_length_in_indices,
			}},
			cell_length = this->get_length(cell),
			ret_val = {{
				this->parameters.coordinates[0][coord_start_indices[0]]
				+ length_of_index[0] * (indices[0] - coord_start_indices[0] * level_0_length_in_indices)
				+ cell_length[0] / 2,

				this->parameters.coordinates[1][coord_start_indices[1]]
				+ length_of_index[1] * (indices[1] - coord_start_indices[1] * level_0_length_in_indices)
				+ cell_length[1] / 2,

				this->parameters.coordinates[2][coord_start_indices[2]]
				+ length_of_index[2] * (indices[2] - coord_start_indices[2] * level_0_length_in_indices)
				+ cell_length[2] / 2,
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

	A quiet NaN is returned if given an invalid
	refinement level or index.
	*/
	std::array<double, 3> get_center(
		const Types<3>::indices_t indices,
		const int refinement_level
	) const {
		return this->get_center(this->mapping.get_cell_from_indices(indices, refinement_level));
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
	Types<3>::indices_t get_indices(
		const std::array<double, 3>& coordinate
	) const {
		Types<3>::indices_t ret_val = {{
			error_index,
			error_index,
			error_index
		}};

		const std::array<double, 3>
			grid_start = this->get_start(),
			grid_end = this->get_end();

		const int max_ref_lvl = this->mapping.get_maximum_refinement_level();
		const uint64_t level_0_length_in_indices = uint64_t(1) << max_ref_lvl;

		for (size_t dimension = 0; dimension < coordinate.size(); dimension++) {

			if (coordinate[dimension] >= grid_start[dimension]
			&& coordinate[dimension] <= grid_end[dimension]) {

				const std::vector<double>& current_coords = this->parameters.coordinates[dimension];

				// find where given coord starts in the coord vector
				uint64_t coord_start_index = 0;
				while (current_coords[coord_start_index] < coordinate[dimension]) {
					coord_start_index++;
				}
				coord_start_index--;

				// length of an index in this level 0 cell
				const double length_of_index
					= (current_coords[coord_start_index + 1] - current_coords[coord_start_index])
					/ double(level_0_length_in_indices);

				// find the offset in indices of given coord inside this level 0 cell
				uint64_t index_offset = 0;
				while (
					current_coords[coord_start_index] + index_offset * length_of_index
					< coordinate[dimension]
				) {
					index_offset++;
				}
				index_offset--;

				ret_val[dimension] = coord_start_index * level_0_length_in_indices + index_offset;
			}
		}

		return ret_val;
	}


	/*!
	Returns the x index in the coordinates vector of the starting coordinate
	of an unrefined cell that is the (grand, grandgrand, ...)parent of given cell.

	Returns error_index if given an invalid cell.
	*/
	std::array<uint64_t, 3> get_level_0_cell_coord_start_index(const uint64_t cell) const
	{
		const int
			ref_lvl = this->mapping.get_refinement_level(cell),
			max_ref_lvl = this->mapping.get_maximum_refinement_level();

		if (cell == error_cell || ref_lvl < 0 || ref_lvl > max_ref_lvl) {
			const std::array<uint64_t, 3> error_val = {{
				error_index,
				error_index,
				error_index
			}};

			return error_val;
		}

		const Types<3>::indices_t indices = this->mapping.get_indices(cell);
		const std::array<uint64_t, 3> ret_val = {{
			indices[0] / (uint64_t(1) << max_ref_lvl),
			indices[1] / (uint64_t(1) << max_ref_lvl),
			indices[2] / (uint64_t(1) << max_ref_lvl)
		}};

		return ret_val;

	}


	/*!
	Writes the geometry into given open file starting at given offset.

	Returns true on success, false otherwise.

	The number of bytes written by this function can be obtained
	from data_size().
	*/
	bool write(MPI_File file, MPI_Offset offset) const
	{
		int ret_val = -1;

		const int temp_id = Stretched_Cartesian_Geometry::geometry_id;
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

		// write number of coordinates in each dimension
		std::array<uint64_t, 3> number_of_coordinates = {{
			this->parameters.coordinates[0].size(),
			this->parameters.coordinates[1].size(),
			this->parameters.coordinates[2].size()
		}};
		ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) number_of_coordinates.data(),
			3,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write geometry start from file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += 3 * sizeof(uint64_t);

		for (size_t dimension = 0; dimension < this->parameters.coordinates.size(); dimension++) {
			ret_val = MPI_File_write_at(
				file,
				offset,
				(void*) &(this->parameters.coordinates[dimension][0]),
				(int) this->parameters.coordinates[dimension].size(),
				MPI_DOUBLE,
				MPI_STATUS_IGNORE
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Couldn't read coordinates in dimension " << dimension
					<< ": " << Error_String()(ret_val)
					<< std::endl;
				return false;
			}
			offset += this->parameters.coordinates[dimension].size() * sizeof(double);
		}

		return true;
	}


	/*!
	Reads the geometry from given open file starting at given offset.

	Returns true on success, false otherwise.
	*/
	bool read(MPI_File file, MPI_Offset offset)
	{
		int
			read_geometry_id = Stretched_Cartesian_Geometry::geometry_id + 1,
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

		// TODO: don't error out if given No_Geometry of Cartesian_Geometry
		if (read_geometry_id != Stretched_Cartesian_Geometry::geometry_id) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Wrong geometry: " << read_geometry_id
				<< ", should be " << Stretched_Cartesian_Geometry::geometry_id
				<< std::endl;
			return false;
		}

		// read number of coordinates in each dimension
		std::array<uint64_t, 3> number_of_coordinates = {{0, 0, 0}};
		ret_val = MPI_File_read_at(
			file,
			offset,
			(void*) number_of_coordinates.data(),
			3,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read geometry start from file: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += 3 * sizeof(uint64_t);


		Parameters read_parameters;
		for (size_t dimension = 0; dimension < this->parameters.coordinates.size(); dimension++) {
			read_parameters.coordinates[dimension].resize(number_of_coordinates[dimension]);

			ret_val = MPI_File_read_at(
				file,
				offset,
				(void*) &(read_parameters.coordinates[dimension][0]),
				(int) number_of_coordinates[dimension],
				MPI_DOUBLE,
				MPI_STATUS_IGNORE
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Couldn't read coordinates in dimension " << dimension
					<< ": " << Error_String()(ret_val)
					<< std::endl;
				return false;
			}
			offset += number_of_coordinates[dimension] * sizeof(double);
		}

		if (!this->set(read_parameters)) {
			return false;
		}

		return true;
	}


	/*!
	Returns the number of bytes that will be required / was required for geometry data.

	Returns the correct value only after set() or read() have been called successfully.
	*/
	size_t data_size() const
	{
		size_t ret_val = 3 * sizeof(uint64_t);

		for (size_t dimension = 0; dimension < this->parameters.coordinates.size(); dimension++) {
			ret_val += this->parameters.coordinates[dimension].size() * sizeof(double);
		}

		return ret_val;
	}



private:

	Parameters parameters;

};	// class

}	// namespace

#endif

