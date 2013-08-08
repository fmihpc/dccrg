/*
A class for game of life restart test of dccrg.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute

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

#ifndef IO_HPP
#define IO_HPP

#include "algorithm"
#include "boost/foreach.hpp"
#include "cstdio"
#include "iostream"
#include "mpi.h"
#include "stdint.h"
#include "vector"

#include "../../dccrg.hpp"

#include "cell.hpp"

template<class UserGeometry> class IO
{
public:

	/*!
	Saves the current state of given game of life grid into the given dc file.

	Header format:
	uint64_t time step
	double   start_x
	double   start_y
	double   start_z
	double   cell_length_x
	double   cell_length_y
	double   cell_length_z
	uint64_t length_x in cells
	uint64_t length_y in cells
	uint64_t length_z in cells
	int      maximum_refinement_level
	*/
	static void save(
		const std::string& name,
		const uint64_t step,
		MPI_Comm& comm,
		dccrg::Dccrg<Cell, UserGeometry>& grid
	) {
		int rank = 0, comm_size = 0;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &comm_size);

		Cell::transfer_only_life = true;
		std::remove(name.c_str());

		MPI_Barrier(comm);

		// MPI_File_open wants a non-constant string
		char* name_c_string = new char [name.size() + 1];
		strncpy(name_c_string, name.c_str(), name.size() + 1);

		MPI_File outfile;

		int result = MPI_File_open(
			comm,
			name_c_string,
			MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL,
			&outfile
		);

		delete [] name_c_string;

		if (result != MPI_SUCCESS) {
			char mpi_error_string[MPI_MAX_ERROR_STRING + 1];
			int string_length;
			MPI_Error_string(result, mpi_error_string, &string_length);
			mpi_error_string[string_length + 1] = '\0';
			std::cerr << "Couldn't open file " << name_c_string
				<< ": " << mpi_error_string
				<< std::endl;
			abort();
		}

		// process 0 writes the header
		MPI_Offset header_size = sizeof(int) + 4 * sizeof(uint64_t) + 6 * sizeof(double);
		if (rank == 0) {
			uint8_t* buffer = new uint8_t [header_size];
			assert(buffer != NULL);

			size_t offset = 0;
			memcpy(buffer + offset, &step, sizeof(uint64_t));
			offset += sizeof(uint64_t);

			const boost::array<double, 3> grid_start = grid.geometry.get_start();
			memcpy(buffer + offset, grid_start.data(), 3 * sizeof(double));
			offset += 3 * sizeof(double);

			const boost::array<double, 3> cell_length = grid.geometry.get_length(1);
			memcpy(buffer + offset, cell_length.data(), 3 * sizeof(double));
			offset += 3 * sizeof(double);

			const boost::array<uint64_t, 3> grid_length = grid.length.get();
			memcpy(buffer + offset, grid_length.data(), 3 * sizeof(uint64_t));
			offset += 3 * sizeof(uint64_t);

			const int max_ref_lvl = grid.get_maximum_refinement_level();
			memcpy(buffer + offset, &max_ref_lvl, sizeof(int));
			offset += sizeof(int);

			result = MPI_File_write_at_all(
				outfile,
				0,
				buffer,
				header_size,
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);

			delete [] buffer;

			if (result != MPI_SUCCESS) {
				char mpi_error_string[MPI_MAX_ERROR_STRING + 1];
				int string_length;
				MPI_Error_string(result, mpi_error_string, &string_length);
				mpi_error_string[string_length + 1] = '\0';
				std::cerr << "Process " << rank
					<< " Couldn't write cell list to file " << name
					<< ": " << mpi_error_string
					<< std::endl;
				abort();
			}

		} else {
			result = MPI_File_write_at_all(
				outfile,
				0,
				(void*) &step,
				0,
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);
		}

		MPI_File_close(&outfile);

		if (!grid.save_grid_data(name, header_size)) {
			std::cerr << "Process " << rank
				<< " Writing grid to file " << name << " failed"
				<< std::endl;
			abort();
		}
		Cell::transfer_only_life = false;
	}


	static void load(
		const std::string& name,
		uint64_t& step,
		MPI_Comm& comm,
		dccrg::Dccrg<Cell, UserGeometry>& grid
	) {
		Cell::transfer_only_life = true;

		// MPI_File_open wants a non-constant string
		char* name_c_string = new char [name.size() + 1];
		strncpy(name_c_string, name.c_str(), name.size() + 1);

		MPI_File infile;

		int result = MPI_File_open(
			comm,
			name_c_string,
			MPI_MODE_RDONLY,
			MPI_INFO_NULL,
			&infile
		);

		delete [] name_c_string;

		if (result != MPI_SUCCESS) {
			std::cerr << "Couldn't open file " << name_c_string
				<< ": " << dccrg::Error_String()(result)
				<< std::endl;
			abort();
		}

		// read only time step from header, other grid parameters are known
		result = MPI_File_read_at_all(
			infile,
			0,
			&step,
			sizeof(uint64_t),
			MPI_BYTE,
			MPI_STATUS_IGNORE
		);

		if (result != MPI_SUCCESS) {
			std::cerr << "Couldn't read time step: "
				<< dccrg::Error_String()(result)
				<< std::endl;
			abort();
		}

		MPI_File_close(&infile);

		MPI_Offset header_size = sizeof(int) + 4 * sizeof(uint64_t) + 6 * sizeof(double);
		grid.load_grid_data(name, header_size);

		Cell::transfer_only_life = false;
	}
};

#endif

