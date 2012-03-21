/*
A class representing one cell in the particle test of dccrg.

Copyright 2012 Finnish Meteorological Institute

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
#include "vector"
#include "iostream"


/*!
Class storing an arbitrary number of particles.
*/
class Cell
{
public:

	int number_of_particles;

	// coordinates of particles in this cell
	std::vector<boost::array<double, 3> > particles;

	/*
	The number of particles is transferred over MPI if
	this is false, otherwise the particle coordinates
	are transferred.
	*/
	static bool transfer_particles;

	// returns the starting address of data to send
	void* at(void)
	{
		if (!Cell::transfer_particles) {
			return &(this->number_of_particles);
		} else {
			if (this->particles.size() > 0) {
				return &(this->particles[0]);
			} else {
				// return a sane address just in case
				return &(this->number_of_particles);
			}
		}
	}

	// returns the length in bytes to transfer between processes for dccrg
	MPI_Datatype mpi_datatype()
	{
		MPI_Datatype datatype;

		if (!Cell::transfer_particles) {
			MPI_Type_contiguous(1, MPI_INT, &datatype);
		} else {
			MPI_Type_contiguous(
				this->particles.size() * sizeof(boost::array<double, 3>),
				MPI_BYTE,
				&datatype
			);
		}

		return datatype;
	}


	// reserves space for particle data coming over MPI.
	void prepare_to_receive_particles()
	{
		this->particles.resize(this->number_of_particles);
	}


	Cell()
	{
		this->number_of_particles = 0;
	}

	~Cell()
	{
		this->number_of_particles = 0;
	}

};

#endif

