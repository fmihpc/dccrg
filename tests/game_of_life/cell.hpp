/*
A class representing one cell in the game of life tests of dccrg.

Copyright 2011 Finnish Meteorological Institute

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

#ifndef CELL_HPP
#define CELL_HPP

#include "boost/array.hpp"
#include "stdint.h"

/*!
Game of life cell with 8 * (13 + EXTRA_DATA) bytes of data.
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

	// use boost::mpi for data transfers over MPI
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER

	template<typename Archiver> void serialize(
		Archiver& ar,
		const unsigned int /*version*/
	) {
		ar & data;
	}

	// use MPI directly for data transfers
	#else

	void* at(void)
	{
		return this;
	}

	#ifdef DCCRG_USER_MPI_DATA_TYPE
	MPI_Datatype mpi_datatype(void) const
	{
		MPI_Datatype type;
		MPI_Type_contiguous(sizeof(this->data), MPI_BYTE, &type);
		return type;
	}
	#else
	static size_t size(void)
	{
		// processes don't need other processes' live neighbour info
		return 13 * sizeof(uint64_t);
	}
	#endif

	#endif	// ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
};

#endif

