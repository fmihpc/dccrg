/*
A class representing one cell in the game of life tests of dccrg.

Copyright 2011, 2012, 2013 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RESTART_CELL_HPP
#define RESTART_CELL_HPP

#include "boost/array.hpp"
#include "mpi.h"
#include "stdint.h"

/*!
Game of life cell with 8 * 13 bytes of data.
*/
class Cell
{
public:

	/*!
	data[0] > 0 if cell is alive
	data[1] total number of live neighbors for all siblings
	The following are used only with cell refinement
	data[2..4] live neighbors of refined cells to calculate the above
	data[5..12] process only one sibling (assume identical state),
	            these record the parents of processed siblings
	*/
	boost::array<uint64_t, 13> data;

	// == true when saving and restarting
	static bool transfer_only_life;

	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

	#error DCCRG_TRANSFER_USING_BOOST_MPI cannot be defined when compiling this program

	#else // ifdef DCCRG_TRANSFER_USING_BOOST_MPI

	void mpi_datatype(
		void*& address,
		int& count,
		MPI_Datatype& datatype,
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/
	) {
		address = &(this->data);

		if (Cell::transfer_only_life) {
			count = 1;
		} else {
			count = 13;
		}

		datatype = MPI_UINT64_T;
	}

	#endif // ifdef DCCRG_TRANSFER_USING_BOOST_MPI
};

bool Cell::transfer_only_life = false;

#endif

