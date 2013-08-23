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
#include "../../dccrg_mpi_support.hpp"

#include "cell.hpp"

template<class UserGeometry> class IO
{
public:

	/*!
	Saves the current state of given game of life grid into the given file.

	The game of life header consists of the given time step.
	*/
	static void save(
		const std::string& name,
		uint64_t step,
		dccrg::Dccrg<Cell, UserGeometry>& grid
	) {
		int rank = 0, comm_size = 0;
		MPI_Comm comm = grid.get_communicator();
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
			std::cerr << "Couldn't open file " << name_c_string
				<< ": " << dccrg::Error_String()(result)
				<< std::endl;
			abort();
		}

		// write time step
		if (rank == 0) {
			result = MPI_File_write_at(
				outfile,
				0,
				(void*) &step,
				1,
				MPI_UINT64_T,
				MPI_STATUS_IGNORE
			);
			if (result != MPI_SUCCESS) {
				std::cerr << "Process " << rank
					<< " Couldn't write cell list to file " << name
					<< ": " << dccrg::Error_String()(result)
					<< std::endl;
				abort();
			}
		}

		MPI_File_close(&outfile);
		MPI_Barrier(comm);

		if (!grid.save_grid_data(name, sizeof(uint64_t))) {
			std::cerr << "Process " << rank
				<< " Writing grid to file " << name << " failed"
				<< std::endl;
			abort();
		}
		Cell::transfer_only_life = false;
	}


	/*!
	Loads a game of life from given file into given grid.

	The grid will be initialized with the given communicator.
	Returns the time step of the game read from the given file.
	*/
	static uint64_t load(
		MPI_Comm& comm,
		const std::string& name,
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
		uint64_t ret_val = 0;
		result = MPI_File_read_at_all(
			infile,
			0,
			&ret_val,
			1,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);

		if (result != MPI_SUCCESS) {
			std::cerr << "Couldn't read time step: "
				<< dccrg::Error_String()(result)
				<< std::endl;
			abort();
		}

		MPI_File_close(&infile);
		MPI_Barrier(comm);

		if (!grid.load_grid_data(
			name,
			sizeof(uint64_t),
			comm,
			"RANDOM"
		)) {
			std::cerr << "Couldn't load grid data"
				<< std::endl;
			abort();
		}

		Cell::transfer_only_life = false;

		return ret_val;
	}
};

#endif

