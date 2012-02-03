/*
A class representing one cell in the advection test of dccrg.

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

class Cell
{
public:

	/*
	data[0]: value, amount of the stuff being advected in this cell
	data[1] > 0: flux into this cell
	data[2]: max relative difference in value between this cell and and its neighbors
	*/
	boost::array<double, 3> data;

	void* at(void)
	{
		return &(this->data[0]);
	}

	static size_t size(void)
	{
		// don't transfer flux or difference
		return sizeof(double);
	}
};

#endif

