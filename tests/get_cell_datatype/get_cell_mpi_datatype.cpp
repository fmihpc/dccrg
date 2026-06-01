/*
Tests get_cell_mpi_datatype function of dccrg.

Copyright 2014, 2015, 2016, 2018 Ilja Honkonen
Copyright 2018 Finnish Meteorological Institute

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


#include "cstdlib"
#include "iostream"
#include "tuple"

#include "mpi.h"

#include "dccrg_get_cell_datatype.hpp"

int main(int /*argc*/, char** /*argv*/) {
	void* address = NULL;
	int count = -1;
	MPI_Datatype datatype = MPI_DATATYPE_NULL;

	signed char char_;
	std::tie(address, count, datatype) = dccrg::detail::get_cell_mpi_datatype(char_, 0, 0, 0, false, 0);
	if (count != -1) { // using get_mpi_datatype_basic
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}
	if (address != nullptr) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}
	if (datatype != MPI_DATATYPE_NULL) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	return EXIT_SUCCESS;
}
