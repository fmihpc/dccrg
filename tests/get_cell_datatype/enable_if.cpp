/*
Tests for enable_if machinery for dccrg.

Copyright 2014 Ilja Honkonen

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


#include "boost/function_types/property_tags.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/tti/has_member_function.hpp"
#include "boost/static_assert.hpp"
#include "boost/tuple/tuple.hpp"
#include "cstdlib"
#include "iostream"
#include "mpi.h"

struct Cell1 {
	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) const {
		return boost::make_tuple((void*) NULL, 1, MPI_DATATYPE_NULL);
	}

	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/,
		const int /*neighborhood_id*/
	) {
		return boost::make_tuple((void*) NULL, 2, MPI_DATATYPE_NULL);
	}

	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype() const
	{
		return boost::make_tuple((void*) NULL, 3, MPI_DATATYPE_NULL);
	}

	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype()
	{
		return boost::make_tuple((void*) NULL, 4, MPI_DATATYPE_NULL);
	}
};

struct Cell2 {
	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype()
	{
		return boost::make_tuple((void*) NULL, 5, MPI_DATATYPE_NULL);
	}
};

struct Cell3 {};


BOOST_TTI_HAS_MEMBER_FUNCTION(get_mpi_datatype)


int main(int /*argc*/, char** /*argv*/)
{
	void* address = NULL;
	int count = -1;
	MPI_Datatype datatype = MPI_DATATYPE_NULL;

	BOOST_STATIC_ASSERT((
		has_member_function_get_mpi_datatype<
			Cell1,
			boost::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<
				const uint64_t,
				const int,
				const int,
				const bool,
				const int
			>,
			boost::function_types::const_qualified
		>::value
	));
	const Cell1 c1_1;
	boost::tie(
		address,
		count,
		datatype
	) = c1_1.get_mpi_datatype(0, 0, 0, false, 0);
	if (count != 1) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	BOOST_STATIC_ASSERT((
		has_member_function_get_mpi_datatype<
			Cell1,
			boost::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<
				const uint64_t,
				const int,
				const int,
				const bool,
				const int
			>
		>::value
	));
	Cell1 c1_2;
	boost::tie(
		address,
		count,
		datatype
	) = c1_2.get_mpi_datatype(0, 0, 0, false, 0);
	if (count != 2) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}


	BOOST_STATIC_ASSERT((
		has_member_function_get_mpi_datatype<
			Cell1,
			boost::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>,
			boost::function_types::const_qualified
		>::value
	));
	const Cell1 c1_3;
	boost::tie(
		address,
		count,
		datatype
	) = c1_3.get_mpi_datatype();
	if (count != 3) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	BOOST_STATIC_ASSERT((
		has_member_function_get_mpi_datatype<
			Cell1,
			boost::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>
		>::value
	));
	Cell1 c1_4;
	boost::tie(
		address,
		count,
		datatype
	) = c1_4.get_mpi_datatype();
	if (count != 4) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}


	BOOST_STATIC_ASSERT((
		has_member_function_get_mpi_datatype<
			Cell2,
			boost::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>
		>::value
	));
	Cell2 c2_1;
	boost::tie(
		address,
		count,
		datatype
	) = c2_1.get_mpi_datatype();
	if (count != 5) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}


	BOOST_STATIC_ASSERT((
		not has_member_function_get_mpi_datatype<
			Cell3,
			boost::tuple<void*, int, MPI_Datatype>,
			boost::mpl::vector<>
		>::value
	));


	return EXIT_SUCCESS;
}

