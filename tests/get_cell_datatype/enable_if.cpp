/*
Tests for enable_if machinery for dccrg.

Copyright 2014, 2015 Ilja Honkonen

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/


#include "cstdlib"
#include "iostream"

#include "boost/function_types/property_tags.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/tti/has_member_function.hpp"
#include "mpi.h"

struct Cell1 {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) const {
		return std::make_tuple((void*) NULL, 1, MPI_DATATYPE_NULL);
	}

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) {
		return std::make_tuple((void*) NULL, 2, MPI_DATATYPE_NULL);
	}

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype() const
	{
		return std::make_tuple((void*) NULL, 3, MPI_DATATYPE_NULL);
	}

	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) NULL, 4, MPI_DATATYPE_NULL);
	}
};

struct Cell2 {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple((void*) NULL, 5, MPI_DATATYPE_NULL);
	}
};

struct Cell3 {};


BOOST_TTI_HAS_MEMBER_FUNCTION(get_mpi_datatype)


int main(int /*argc*/, char** /*argv*/)
{
	void* address = NULL;
	int count = -1;
	MPI_Datatype datatype = MPI_DATATYPE_NULL;

	static_assert(
		has_member_function_get_mpi_datatype<
			Cell1,
			std::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<
				const uint64_t,
				const int,
				const int,
				const bool,
				const int
			>,
			boost::function_types::const_qualified
		>::value,
		"Error"
	);
	const Cell1 c1_1;
	std::tie(
		address,
		count,
		datatype
	) = c1_1.get_mpi_datatype(0, 0, 0, false, 0);
	if (count != 1) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	static_assert(
		has_member_function_get_mpi_datatype<
			Cell1,
			std::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<
				const uint64_t,
				const int,
				const int,
				const bool,
				const int
			>
		>::value,
		"Error"
	);
	Cell1 c1_2;
	std::tie(
		address,
		count,
		datatype
	) = c1_2.get_mpi_datatype(0, 0, 0, false, 0);
	if (count != 2) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}


	static_assert(
		has_member_function_get_mpi_datatype<
			Cell1,
			std::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>,
			boost::function_types::const_qualified
		>::value,
		"Error"
	);
	const Cell1 c1_3;
	std::tie(
		address,
		count,
		datatype
	) = c1_3.get_mpi_datatype();
	if (count != 3) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	static_assert(
		has_member_function_get_mpi_datatype<
			Cell1,
			std::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>
		>::value,
		"Error"
	);
	Cell1 c1_4;
	std::tie(
		address,
		count,
		datatype
	) = c1_4.get_mpi_datatype();
	if (count != 4) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}


	static_assert(
		has_member_function_get_mpi_datatype<
			Cell2,
			std::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>
		>::value,
		"Error"
	);
	Cell2 c2_1;
	std::tie(
		address,
		count,
		datatype
	) = c2_1.get_mpi_datatype();
	if (count != 5) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}


	static_assert(
		not has_member_function_get_mpi_datatype<
			Cell3,
			std::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>
		>::value,
		"Error"
	);


	return EXIT_SUCCESS;
}

