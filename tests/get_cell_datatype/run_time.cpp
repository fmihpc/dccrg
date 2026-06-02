/*
Test for dccrg's get_cell_mpi_datatype logic.

Copyright 2014, 2015, 2016, 2018 Ilja Honkonen

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

#include "iostream"
#include "tuple"

#include "mpi.h"

#include "dccrg_get_cell_datatype.hpp"

struct Cell3 {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) {
		return std::make_tuple((void*) NULL, 3, MPI_DATATYPE_NULL);
	}
};

struct Cell4 {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) const {
		return std::make_tuple((void*) NULL, 4, MPI_DATATYPE_NULL);
	}
};

struct Cell5 {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) const {
		return std::make_tuple((void*) NULL, 5, MPI_DATATYPE_NULL);
	}

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) {
		return std::make_tuple((void*) NULL, 6, MPI_DATATYPE_NULL);
	}
};

#define CHECK_DATATYPE_COUNT(cell, returned_count) \
std::tie( \
	address, \
	count, \
	datatype \
) = dccrg::detail::get_cell_mpi_datatype( \
	cell, \
	0, \
	0, \
	0, \
	false, \
	0 \
); \
if (count != returned_count) { \
	std::cerr << __FILE__ << ":" << __LINE__ << std::endl; \
	abort(); \
}

int main(int /*argc*/, char** /*argv*/)
{
	void* address = NULL;
	int count = -1;
	MPI_Datatype datatype = MPI_DATATYPE_NULL;

	Cell3 c3;
	CHECK_DATATYPE_COUNT(c3, 3)

	const Cell4 c4;
	CHECK_DATATYPE_COUNT(c4, 4)

	const Cell5 c5_1;
	CHECK_DATATYPE_COUNT(c5_1, 5)

	Cell5 c5_2;
	CHECK_DATATYPE_COUNT(c5_2, 6)

	return EXIT_SUCCESS;
}
