/*
Dccrg class for storing the grid length and related operations.

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

#ifndef DCCRG_LENGTH_HPP
#define DCCRG_LENGTH_HPP


#include "array"
#include "iostream"
#include "stdint.h"


namespace dccrg {

/*!
\brief Length of the grid in cells of refinement level 0 and related operations.
*/
class Grid_Length
{

public:

	/*!
	Represents the length of the grid in cells of
	refinement level 0 in each dimension.
	*/
	typedef std::array<uint64_t, 3> type;


	/*!
	Sets the size of the grid in cells of refinement level 0 to 1.
	*/
	Grid_Length()
	{
		this->length[0] =
		this->length[1] =
		this->length[2] = 1;
	}

	/*!
	Creates a grid of given length in cells of refinement level 0.
	*/
	Grid_Length(const Grid_Length::type& given_length)
	{
		if (!this->set(given_length)) {
			abort();
		}
	}


	/*!
	Sets the grid to a default constructed state.
	*/
	~Grid_Length()
	{
		this->length[0] =
		this->length[1] =
		this->length[2] = 1;
	}


	/*!
	The length of the grid in cells of refinement level 0.
	*/
	const Grid_Length::type& get() const
	{
		return this->length;
	}


	/*!
	Sets the length of the grid in cells of refinement level 0.

	Returns true if successful, probably invalidating all previous cell
	information (cell numbers, indices, etc.)
	Returns false if unsuccessful and in that case has no effect.
	Automatically maximizes the maximum refinement level of the grid.
	*/
	bool set(const Grid_Length::type& given_length)
	{
		if (
			   given_length[0] == 0
			|| given_length[1] == 0
			|| given_length[2] == 0
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "All lengths given must be > 0 but are "
				<< given_length[0] << " "
				<< given_length[1] << " "
				<< given_length[2]
				<< std::endl;
			return false;
		}

		const Grid_Length::type old_length = this->length;

		this->length = given_length;

		const double max_number_of_initial_cells
			= double(this->length[0])
			+ double(this->length[1])
			+ double(this->length[2]);

		if (max_number_of_initial_cells > double(~uint64_t(0))) {
			std::cerr
				<< "Grid would have too many cells of refinement level 0 for uint64_t: "
				<< max_number_of_initial_cells
				<< " (length_x, length_y, length_z: "
				<< this->length[0] << " " << this->length[1] << " " << this->length[2]
				<< ")"
				<< std::endl;

			this->length = old_length;
			return false;
		}

		return true;
	}



private:

	Grid_Length::type length;

};

}	// namespace

#endif

