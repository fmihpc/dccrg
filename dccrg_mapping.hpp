/*
Dccrg class for mapping cell ids to their size and location.

Copyright 2009, 2010, 2011, 2012, 2013,
2014, 2015, 2016 Finnish Meteorological Institute

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

#ifndef DCCRG_MAPPING_HPP
#define DCCRG_MAPPING_HPP


#include "array"
#include "cmath"
#include "cstdint"
#include "iostream"

#include "dccrg_length.hpp"
#include "dccrg_mpi_support.hpp"
#include "dccrg_types.hpp"


namespace dccrg {

//! Indicates a non-existing cell or an error when dealing with cells
static const uint64_t error_cell = 0;

//! Indicates a non-existing index or an error when dealing with indices
static const uint64_t error_index = 0xFFFFFFFFFFFFFFFF;


/*!
\brief Maps cells ids to their size and location.

Also handles calculations related to cell indices.

\see
Figures 2, 3 and the related text in
http://dx.doi.org/10.1016/j.cpc.2012.12.017
or
http://arxiv.org/abs/1212.3496.
*/
class Mapping
{

private:

	//! read-write version of Grid_Length for internal use
	Grid_Length length_rw;


public:


	//! see Grid_Length
	const Grid_Length& length;


	/*!
	Creates a grid with size of 1 cell and maximum refinement level 0.

	Length of the grid is 1 cell in each dimension.
	*/
	Mapping() : length(length_rw)
	{
		this->max_refinement_level = 0;
		this->last_cell = 1;
	}

	/*!
	Creates a grid with size of 1 cell and maximum refinement level 0.

	Length of the grid is 1 cell in each dimension.
	*/
	Mapping(
		const std::array<uint64_t, 3>& given_length
	) : length(length_rw)
	{
		if (!this->length_rw.set(given_length)) {
			abort();
		}
		this->max_refinement_level = 0;
		this->update_last_cell();
	}


	/*!
	Sets the grid to a default constructed state.
	*/
	~Mapping()
	{
		this->max_refinement_level = 0;
		this->last_cell = 1;
	}


	//! see Grid_Length::set_length()
	bool set_length(const std::array<uint64_t, 3>& given_length)
	{
		if (!this->length_rw.set(given_length)) {
			return false;
		}

		this->update_last_cell();

		return true;
	}


	/*!
	Returns the maximum refinement level of any cell in the grid (0 means unrefined).
	*/
	const int& get_maximum_refinement_level() const
	{
		return this->max_refinement_level;
	}

	/*!
	Sets the maximum refinement level of the grid.

	0 means unrefined.
	Probably invalidating all previous cell indices.
	Returns true if successful, otherwise does nothing and returns false.
	*/
	bool set_maximum_refinement_level(const int given_refinement_level)
	{
		if (given_refinement_level <= this->get_maximum_possible_refinement_level()) {
			this->max_refinement_level = given_refinement_level;
			this->update_last_cell();
			return true;
		} else {
			return false;
		}
	}


	/*!
	Returns the cell of given refinement level at given indices.

	Returns error_cell if any of given indices is invalid.
	*/
	uint64_t get_cell_from_indices(
		const Types<3>::indices_t& indices,
		const int refinement_level
	) const {
		if (indices[0] >= this->length.get()[0] * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (indices[1] >= this->length.get()[1] * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (indices[2] >= this->length.get()[2] * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (refinement_level < 0) {
			return error_cell;
		}

		if (refinement_level > this->max_refinement_level) {
			return error_cell;
		}

		// cell numbering starts at 1
		uint64_t cell = 1;

		// add ids of larger cells
		for (int i = 0; i < refinement_level; i++) {
			cell +=
				  this->length.get()[0]
				* this->length.get()[1]
				* this->length.get()[2]
				* (uint64_t(1) << (i * 3));
		}

		// convert to indices of this cell's refinement level
		const Types<3>::indices_t this_level_indices = {{
			indices[0] / (uint64_t(1) << (this->max_refinement_level - refinement_level)),
			indices[1] / (uint64_t(1) << (this->max_refinement_level - refinement_level)),
			indices[2] / (uint64_t(1) << (this->max_refinement_level - refinement_level))
		}};

		// get the length of the grid in terms of cells of this refinement level
		const std::array<uint64_t, 2> this_level_length = {{
			this->length.get()[0] * (uint64_t(1) << refinement_level),
			this->length.get()[1] * (uint64_t(1) << refinement_level)
		}};

		cell
			+= this_level_indices[0]
			+ this_level_indices[1] * this_level_length[0]
			+ this_level_indices[2] * this_level_length[0] * this_level_length[1];

		return cell;
	}


	/*!
	Returns the indices of given cell.

	See the definition of indices_t for an explanation about them.
	Returned indices are invalid if given a cell outside the valid range.
	*/
	Types<3>::indices_t get_indices(uint64_t cell) const
	{
		if (cell == error_cell || cell > this->last_cell) {
			const Types<3>::indices_t error_indices = {{error_index, error_index, error_index}};
			return error_indices;
		}

		// substract ids of larger cells
		const int refinement_level = this->get_refinement_level(cell);
		for (int i = 0; i < refinement_level; i++) {
			cell -=
				  this->length.get()[0]
				* this->length.get()[1]
				* this->length.get()[2]
				* (uint64_t(1) << (i * 3));
		}

		cell -= 1;	// cell numbering starts from 1
		const Types<3>::indices_t indices = {{

			  (cell % (this->length.get()[0] * (uint64_t(1) << refinement_level)))
			* (uint64_t(1) << (max_refinement_level - refinement_level)),

			((cell / (this->length.get()[0] * (uint64_t(1) << refinement_level)))
				% (this->length.get()[1] * (uint64_t(1) << refinement_level)))
			* (uint64_t(1) << (max_refinement_level - refinement_level)),

			(cell / (
				this->length.get()[0]
				* this->length.get()[1]
				* (uint64_t(1) << (2 * refinement_level))
			))
			* (uint64_t(1) << (max_refinement_level - refinement_level))
		}};

		return indices;
	}


	/*!
	Returns the refinement level of given cell (0 means unrefined).

	Returns -1 if given an invalid cell.
	*/
	int get_refinement_level(const uint64_t cell) const
	{
		if (cell == error_cell || cell > this->last_cell) {
			return -1;
		}

		int refinement_level = 0;
		uint64_t current_last = 0;

		while (refinement_level <= this->max_refinement_level) {
			current_last +=
				  this->length.get()[0]
				* this->length.get()[1]
				* this->length.get()[2]
				* (uint64_t(1) << 3 * refinement_level);

			if (cell <= current_last) {
				break;
			}

			refinement_level++;
		}

		if (refinement_level > this->max_refinement_level) {
			return -1;
		}

		return refinement_level;
	}


	/*!
	Returns the lengths of given cell in indices in every direction.

	Returns error_index if given an invalid cell.
	*/
	uint64_t get_cell_length_in_indices(const uint64_t cell) const
	{
		if (cell == error_cell) {
			return error_index;
		}

		const int refinement_level = this->get_refinement_level(cell);

		if (refinement_level < 0) {
			return error_index;
		}

		return uint64_t(1) << (this->max_refinement_level - refinement_level);
	}


	/*!
	Returns the maximum possible refinement level for a cell in the grid (0 means unrefined).
	*/
	int get_maximum_possible_refinement_level() const
	{
		const uint64_t grid_length
			= this->length.get()[0] * this->length.get()[1] * this->length.get()[2];
		int refinement_level = 0;
		double current_last = 0;
		while (current_last <= double(~uint64_t(0))) {
			// TODO: don't assume 3 dimensions
			current_last += double(grid_length) * pow(double(8), double(refinement_level));
			refinement_level++;
		}

		return refinement_level - 2;
	}


	/*!
	Returns the parent of given cell.

	Returns the given cell if its refinement level == 0.
	Returns error_cell if given cell's refinement level > maximum refinement level or < 0.

	\see get_refinement_level()
	*/
	uint64_t get_parent(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (
			refinement_level < 0
			or refinement_level > this->max_refinement_level
		) {
			return error_cell;
		}

		if (refinement_level == 0) {
			return cell;
		}

		return this->get_cell_from_indices(this->get_indices(cell), refinement_level - 1);
	}


	/*!
	Returns all children of given cell.

	Returns error_cells if childrens' refinement level would exceed max_refinement_level.
	 */
	std::array<uint64_t, 8> get_all_children(const uint64_t cell) const
	{
		std::array<uint64_t, 8> children{
			error_cell, error_cell, error_cell, error_cell,
			error_cell, error_cell, error_cell, error_cell
		};

		if (cell == error_cell) {
			return children;
		}

		// given cell cannot have children
		int refinement_level = this->get_refinement_level(cell);
		if (refinement_level >= this->get_maximum_refinement_level()) {
			return children;
		}

		Types<3>::indices_t indices = this->get_indices(cell);

		// get indices of next refinement level within this cell
		refinement_level++;
		const uint64_t index_offset
			= (uint64_t(1) << (this->get_maximum_refinement_level() - refinement_level));

		size_t i = 0;
		for (uint64_t
			z_index_offset = 0;
			z_index_offset < 2 * index_offset;
			z_index_offset += index_offset
		)
		for (uint64_t
			y_index_offset = 0;
			y_index_offset < 2 * index_offset;
			y_index_offset += index_offset
		)
		for (uint64_t
			x_index_offset = 0;
			x_index_offset < 2 * index_offset;
			x_index_offset += index_offset
		) {
			const Types<3>::indices_t index{
				indices[0] + x_index_offset,
				indices[1] + y_index_offset,
				indices[2] + z_index_offset
			};

			children[i++] = this->get_cell_from_indices(index, refinement_level);
		}

		return children;
	}


	/*!
	Returns given cell and its siblings or just given cell if refinement level 0.

	Returns error_cells if given cell with too large refinement level.
	*/
	std::array<uint64_t, 8> get_siblings(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);
		if (
			refinement_level < 0
			or refinement_level > this->get_maximum_refinement_level()
		) {
			return {
				error_cell, error_cell, error_cell, error_cell,
				error_cell, error_cell, error_cell, error_cell
			};
		}

		if (refinement_level == 0) {
			return {
				cell, error_cell, error_cell, error_cell,
				error_cell, error_cell, error_cell, error_cell
			};
		}

		return this->get_all_children(this->get_parent(cell));
	}


	/*!
	Returns the refinement level 0 parent of given cell.

	If given a refinement level 0 cell returns the same cell.
	Returns dccrg::error_cell otherwise.
	*/
	uint64_t get_level_0_parent(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);

		if (refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return error_cell;
		}

		if (refinement_level == 0) {
			return cell;
		}

		return this->get_cell_from_indices(this->get_indices(cell), 0);
	}


	/*!
	Returns the last valid cell id in the grid.
	*/
	uint64_t get_last_cell() const
	{
		return this->last_cell;
	}


	//! Format in which mapping data is stored into a file.
	typedef std::array<uint8_t, 3> mapping_file_data_t;


	/*!
	Reads the mapping from given open file starting at given offset.

	Returns true on success, false otherwise.
	Reads grid length in level 0 cells and cells' maximum refinement level.
	*/
	bool read(MPI_File file, MPI_Offset offset)
	{
		int
			ret_val = -1,
			max_ref_lvl = -1;

		Grid_Length::type length = {{0, 0, 0}};
		ret_val = MPI_File_read_at_all(
			file,
			offset,
			(void*) length.data(),
			3,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read length data: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += sizeof(Grid_Length::type);

		if (!this->set_length(length)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't set length of grid"
				<< std::endl;
			return false;
		}

		ret_val = MPI_File_read_at_all(
			file,
			offset,
			(void*) &max_ref_lvl,
			1,
			MPI_INT,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read maximum refinement level: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		if (!this->set_maximum_refinement_level(max_ref_lvl)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Couldn't set maximum refinement level to " << max_ref_lvl
				<< std::endl;
			return false;
		}

		return true;
	}


	/*!
	Writes the mapping into given open file starting at given offset.

	Returns true on success, false otherwise.
	*/
	bool write(MPI_File file, MPI_Offset offset) const
	{
		int ret_val = -1;

		Grid_Length::type length = this->length.get();
		ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) length.data(),
			3,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write length data: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}
		offset += sizeof(Grid_Length::type);

		ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) &(this->max_refinement_level),
			1,
			MPI_INT,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write maximum refinement level: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		return true;
	}


	/*!
	Returns the number of bytes required for mapping data.
	*/
	size_t data_size() const
	{
		return sizeof(Grid_Length::type) + sizeof(int);
	}



private:

	//! maximum refinemet level of any cell in the grid, 0 means unrefined
	int max_refinement_level;

	/*!
	Last valid cell id based on grid lengths and maximum
	refinement level of cells in the grid
	*/
	uint64_t last_cell;

	/*!
	Set the value of last_cell based on current grid lengths and max_refinement_level.
	*/
	void update_last_cell()
	{
		const uint64_t grid_length
			= this->length.get()[0] * this->length.get()[1] * this->length.get()[2];
		this->last_cell = 0;
		for (int i = 0; i <= this->max_refinement_level; i++) {
			this->last_cell += grid_length * (uint64_t(1) << (i * 3));
		}
	}


};	// class

}	// namespace

#endif

