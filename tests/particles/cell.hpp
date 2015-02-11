/*
A class representing one cell in the particle test of dccrg.

Copyright 2012, 2013, 2014, 2015 Finnish Meteorological Institute

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


#include "vector"
#include "iostream"

#include "mpi.h"


/*!
Class storing an arbitrary number of particles.
*/
class Cell
{
public:

	unsigned int number_of_particles = 0;

	// coordinates of particles in this cell
	std::vector<std::array<double, 3> > particles;

	/*
	The number of particles is transferred over MPI if
	this is false, otherwise the particle coordinates
	are transferred.
	*/
	static bool transfer_particles;


	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		void* address = NULL;
		int count = -1;
		MPI_Datatype datatype = MPI_DATATYPE_NULL;

		if (Cell::transfer_particles) {

			if (this->particles.size() > 0) {
				address = &(this->particles[0]);
			} else {
				// return a sane address just in case
				address = &(this->number_of_particles);
			}

			count = this->particles.size() * 3;
			datatype = MPI_DOUBLE;

		} else {

			address = &(this->number_of_particles);
			count = 1;
			datatype = MPI_UNSIGNED;

		}

		return std::make_tuple(address, count, datatype);
	}


	// reserves space for particle data coming over MPI.
	void resize()
	{
		this->particles.resize(this->number_of_particles);
	}
};

#endif

