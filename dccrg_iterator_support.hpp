/*
Supporting routines for dccrg iterators.

Copyright 2013 Finnish Meteorological Institute

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


#ifndef DCCRG_ITERATOR_SUPPORT_HPP
#define DCCRG_ITERATOR_SUPPORT_HPP


namespace dccrg {


template <class Cell_Data, class Geometry> class Dccrg;


/*!
\brief Returns true if given cell is not on the process boundary, false otherwise.
*/
template<class Cell_Data, class Geometry> class Is_Inner_Cell
{
friend class Dccrg<Cell_Data, Geometry>;

public:
	//! Does all the work
	bool operator()(const typename Dccrg<Cell_Data, Geometry>::cell_and_data_pair_t& item)
	{
		if (this->grid.get_local_cells_on_process_boundary_internal().count(item.first)) {
			return false;
		} else {
			return true;
		}
	}

protected:
	/*!
	Given grid is used to query whether a cell is local or not.
	*/
	Is_Inner_Cell(const Dccrg<Cell_Data, Geometry>& given_grid) : grid(given_grid) {}

private:
	const Dccrg<Cell_Data, Geometry>& grid;
};


/*!
\brief Returns true if given cell is on the process boundary, false otherwise.
*/
template<class Cell_Data, class Geometry> class Is_Outer_Cell
{
friend class Dccrg<Cell_Data, Geometry>;

public:
	//! Does all the work
	bool operator()(const typename Dccrg<Cell_Data, Geometry>::cell_and_data_pair_t& item)
	{
		if (this->grid.get_local_cells_on_process_boundary_internal().count(item.first)) {
			return true;
		} else {
			return false;
		}
	}

protected:
	/*!
	Given grid is used to query whether a cell is local or not.
	*/
	Is_Outer_Cell(const Dccrg<Cell_Data, Geometry>& given_grid) : grid(given_grid) {}

private:
	const Dccrg<Cell_Data, Geometry>& grid;
};

} // namespace

#endif

