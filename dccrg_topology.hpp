/*
Dccrg class for storing the topology of the grid.

Copyright 2009, 2010, 2011, 2012, 2013,
2014, 2015, 2016 Finnish Meteorological Institute

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


#include "array"

#include "dccrg_mpi_support.hpp"


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
		const std::array<bool, 3>& given_periodicity
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


	//! Format in which topology data is stored into a file.
	typedef std::array<uint8_t, 3> topology_file_data_t;


	/*!
	Reads the topology from given open file starting at given offset.

	Returns true on success, false otherwise.
	*/
	bool read(MPI_File file, MPI_Offset offset)
	{
		topology_file_data_t data = {{0, 0, 0}};
		const int ret_val = MPI_File_read_at_all(
			file,
			offset,
			(void*) data.data(),
			3,
			MPI_UINT8_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't read topology data: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		for (size_t dimension = 0; dimension < this->periodic.size(); dimension++) {
			const bool is_periodic = (data[dimension] > 0) ? true : false;
			if (!this->set_periodicity(dimension, is_periodic)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Couldn't set periodicity in dimension: " << dimension
					<< std::endl;
				return false;
			}
		}

		return true;
	}


	/*!
	Writes the topology into given open file starting at given offset.

	Returns true on success, false otherwise.
	*/
	bool write(MPI_File file, MPI_Offset offset) const
	{
		topology_file_data_t data = {{0, 0, 0}};
		if (this->is_periodic(0)) {
			data[0] = uint8_t(1);
		}
		if (this->is_periodic(1)) {
			data[1] = uint8_t(1);
		}
		if (this->is_periodic(2)) {
			data[2] = uint8_t(1);
		}

		const int ret_val = MPI_File_write_at(
			file,
			offset,
			(void*) data.data(),
			3,
			MPI_UINT8_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't write topology data: " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		return true;
	}


	/*!
	Returns the number of bytes required for topology data.
	*/
	size_t data_size() const
	{
		return sizeof(topology_file_data_t);
	}



private:

	// true means grid wraps around in that dimension
	std::array<bool, 3> periodic;

};

}	// namespace

#endif

