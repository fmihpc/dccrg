/*
A class representing one cell in the advection test of dccrg.

Copyright 2012, 2013, 2014, 2015, 2016,
2018 Finnish Meteorological Institute

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

#ifndef CELL_HPP
#define CELL_HPP

#include "cstdint"

#include "mpi.h"

class Cell
{
public:

	static bool transfer_all_data;

	/*
	data[0]: density of the stuff being advected
	data[1]: vx in the center of the cell
	data[2]: vy
	data[3]: vz
	data[4] > 0: flux of density into this cell
	data[5]: maximum relative difference in density between this cell and its neighbors
	data[6]: length of cell in x dimension
	data[7]: length in y
	data[8]: length in z
	*/
	std::array<double, 9> data;

	// returns MPI_Datatype corresponding to cell data to transfer
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		int nr_to_transfer = 1;
		if (transfer_all_data) {
			nr_to_transfer = this->data.size();
		}
		// transfer cell density or all data to other processes
		return std::make_tuple((void*) this->data.data(), nr_to_transfer, MPI_DOUBLE);
	}

	/*
	References to variables stored in the cell
	*/

	// density
	const double& density() const
	{
		return this->data[0];
	}
	// https://stackoverflow.com/a/123995
	double& density()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[0]);
	}

	// vx
	const double& vx() const
	{
		return this->data[1];
	}
	double& vx()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[1]);
	}

	// vy
	const double& vy() const
	{
		return this->data[2];
	}
	double& vy()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[2]);
	}

	// vz
	const double& vz() const
	{
		return this->data[3];
	}
	double& vz()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[3]);
	}

	// flux
	const double& flux() const
	{
		return this->data[4];
	}
	double& flux()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[4]);
	}

	// max difference
	const double& max_diff() const
	{
		return this->data[5];
	}
	double& max_diff()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[5]);
	}

	// x length
	const double& length_x() const
	{
		return this->data[6];
	}
	double& length_x()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[6]);
	}

	// y length
	const double& length_y() const
	{
		return this->data[7];
	}
	double& length_y()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[7]);
	}

	// z length
	const double& length_z() const
	{
		return this->data[8];
	}
	double& length_z()
	{
		return const_cast<double&>(static_cast<const Cell&>(*this).data[8]);
	}
};

struct Is_Local {
	bool is_local = false;
	template<
		class Grid, class Cell_Item, class Neighbor_Item
	> void update(
		const Grid& grid, const Cell_Item&, const Neighbor_Item& neighbor, const int&, const Is_Local&
	) {
		is_local = grid.is_local(neighbor.id);
	}
};

struct Center {
	std::array<double, 3> center;
	template<
		class Grid, class Cell_Item
	> void update(
		const Grid& grid, const Cell_Item& cell, const Center&
	) {
		center = grid.geometry.get_center(cell.id);
	}
};

#endif
