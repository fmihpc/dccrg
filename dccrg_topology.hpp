/*
Dccrg class for storing the topology of the grid.

Copyright 2009, 2010, 2011, 2012, 2013 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_TOPOLOGY_HPP
#define DCCRG_TOPOLOGY_HPP


#include "boost/array.hpp"


namespace dccrg {

/*!
\brief Stores information about the topology of the grid.

Stores whether the grid wraps around in each dimension.
*/
class Grid_Topology
{

public:

	/*!
	Creates a topology in which the grid does not wrap around.
	*/
	Grid_Topology()
	{
		this->periodic[0] =
		this->periodic[1] =
		this->periodic[2] = false;
	}

	/*!
	Creates the given topology.

	true means that the grid wraps around in the dimension,
	false means that it does not.
	*/
	Grid_Topology(
		const boost::array<bool, 3>& given_periodicity
	) :
		periodic(given_periodicity)
	{}


	/*!
	Sets the periodicity of the geometry.

	index = 0 == x direction.
	Returns true on success and false otherwise.
	*/
	bool set_periodicity(const size_t index, const bool value)
	{
		if (index > 2) {
			return false;
		}

		this->periodic[index] = value;
		return true;
	}


	/*!
	Returns whether the topology is periodic in given direction.

	index = 0 == x direction.
	Returns false if given index > 2.
	*/
	bool is_periodic(const size_t index) const
	{
		if (index > 2) {
			return false;
		}

		return this->periodic[index];
	}


private:

	// true means grid wraps around in that dimension
	boost::array<bool, 3> periodic;

};

}	// namespace

#endif

