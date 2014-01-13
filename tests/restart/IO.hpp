/*
A class for game of life restart test of dccrg.

Copyright 2010, 2011, 2012, 2013, 2014 Finnish Meteorological Institute

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
#include "cstdint"
#include "cstdio"
#include "iostream"
#include "mpi.h"
#include "tuple"
#include "vector"

#include "../../dccrg.hpp"
#include "../../dccrg_cartesian_geometry.hpp"
#include "../../dccrg_mpi_support.hpp"

#include "cell.hpp"

template<class UserGeometry> class IO
{
public:

	/*!
	Saves the current state of given game of life grid into the given file.

	The game of life header consists of the given time step.
	Must be called by all processes with identical arguments,
	an existing file with given name is removed.
	*/
	static void save(
		const std::string& name,
		uint64_t step,
		dccrg::Dccrg<Cell, UserGeometry>& grid
	) {
		int rank = 0;
		MPI_Comm comm = grid.get_communicator();
		MPI_Comm_rank(comm, &rank);
		if (rank == 0) {
			std::remove(name.c_str());
		}
		MPI_Barrier(comm);
		MPI_Comm_free(&comm);

		Cell::transfer_only_life = true;

		std::tuple<void*, int, MPI_Datatype> header;
		if (rank == 0) {
			std::get<0>(header) = &step;
			std::get<1>(header) = 1;
			std::get<2>(header) = MPI_UINT64_T;
		} else {
			// give a valid address just in case
			std::get<0>(header) = &(std::get<1>(header));
			std::get<1>(header) = 0;
			std::get<2>(header) = MPI_BYTE;
		}

		if (!grid.save_grid_data(name, 0, header)) {
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

		// read only time step from header, other grid parameters are known
		uint64_t ret_val = 0;

		std::tuple<void*, int, MPI_Datatype> header;
		std::get<0>(header) = &ret_val;
		std::get<1>(header) = 1;
		std::get<2>(header) = MPI_UINT64_T;

		if (!grid.load_grid_data(
			name,
			0,
			header,
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

