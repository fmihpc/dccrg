/*
A class representing one cell in the advection test of dccrg.

Copyright 2012, 2013, 2014 Finnish Meteorological Institute

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
#include "mpi.h"
#include "stdint.h"

class Cell
{
public:

	/*
	data[0]: density of the stuff being advected
	data[1]: vx in the center of the cell
	data[2]: vy
	data[3]: vz
	data[4] > 0: flux of density into this cell
	data[5]: maximum relative difference in density between this cell and and its neighbors
	*/
	std::array<double, 6> data;

	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

	template<typename Archiver> void serialize(
		Archiver& ar,
		const unsigned int /*version*/
	) {
		// transfer cell density and velocities to other processes
		// TODO: only transfer velocities after they have changed
		ar & data[0] & data[1] & data[2] & data[3];
	}

	#else

	// returns MPI_Datatype corresponding to cell data to transfer
	std::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype() const
	{
		// transfer cell density and velocities to other processes
		// TODO: only transfer velocities after they have changed
		return std::make_tuple((void*) &(this->data), 4, MPI_DOUBLE);
	}

	#endif

	/*
	References to variables stored in the cell
	*/

	// operator versions
	double& operator [] (const std::size_t i)
	{
		return this->data[i];
	}

	const double& operator [] (const std::size_t i) const
	{
		return this->data[i];
	}

	// density
	double& density()
	{
		return this->data[0];
	}


	const double& density() const
	{
		return this->data[0];
	}

	// vx
	double& vx()
	{
		return this->data[1];
	}


	const double& vx() const
	{
		return this->data[1];
	}

	// vy
	double& vy()
	{
		return this->data[2];
	}


	const double& vy() const
	{
		return this->data[2];
	}

	// vz
	double& vz()
	{
		return this->data[3];
	}


	const double& vz() const
	{
		return this->data[3];
	}

	// flux
	double& flux()
	{
		return this->data[4];
	}


	const double& flux() const
	{
		return this->data[4];
	}

	// max difference
	double& max_diff()
	{
		return this->data[5];
	}


	const double& max_diff() const
	{
		return this->data[5];
	}
};

#endif

