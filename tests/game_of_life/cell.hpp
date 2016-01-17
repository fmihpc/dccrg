/*
A class representing one cell in the game of life tests of dccrg.

Copyright 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CELL_HPP
#define CELL_HPP

#include "array"
#include "cstdint"
#include "tuple"

#include "mpi.h"

/*!
Game of life cell with 8 * 13 bytes of data.
*/
struct Cell
{
	/*!
	data[0] > 0 if cell is alive
	data[1..8] record which neighbors are alive,
	           with adaptive mesh refinement record
	           live (grand...)parents of refinement
	           level 0
	*/
	std::array<uint64_t, 9> data;

	std::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype(
		// test detail::get_cell_datatype
		#ifdef SEND_SINGLE_CELLS
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
		#endif
	) const {
		return std::make_tuple((void*) &(this->data), this->data.size(), MPI_UINT64_T);
	}
};

#endif

