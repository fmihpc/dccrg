/*
Saver class for the advection test program of dccrg.

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

#ifndef DCCRG_ADVECTION_SAVE_HPP
#define DCCRG_ADVECTION_SAVE_HPP


#include "boost/array.hpp"
#include "boost/mpi.hpp"
#include "cstring"
#include "string"
#include "vector"

#include "dccrg.hpp"


/*!
Saves the given simulation as a .dc file of the given name.

Returns the number of bytes written by this process.
Must be called simultaneously by all processes.
*/
class Save
{

public:

	template<
		class CellData,
		class Geometry
	> size_t operator()(
		const std::string& filename,
		MPI_Comm& comm,
		const dccrg::Dccrg<CellData, Geometry>& grid
	) {
		int rank = 0, comm_size = 0;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &comm_size);

		const std::string header(
			"2d advection test file\n\n"
			"Data after end of header and a line break:\n"
			"1 uint64_t 0x1234567890abcdef for checking endiannes of data\n"
			"1 double   grid start coordinate in x direction\n"
			"1 double   grid start coordinate in y direction\n"
			"1 double   grid start coordinate in z direction\n"
			"1 double   x size of unrefined spatial cells\n"
			"1 double   y size of unrefined spatial cells\n"
			"1 double   z size of unrefined spatial cells\n"
			"1 uint64_t length of the grid in unrefined cells in x direction\n"
			"1 uint64_t length of the grid in unrefined cells in y direction\n"
			"1 uint64_t length of the grid in unrefined cells in z direction\n"
			"1 uint8_t  maximum refinement level of the grid\n"
			"1 uint64_t cell id\n"
			"1 uint32_t cell process number\n"
			"1 double   density\n"
			"1 double   max relative difference in density between this cell and its neighbors\n"
			"1 double   vx\n"
			"1 double   vy\n"
			"1 double   vz\n"
			"1 uint64_t cell id\n"
			"...\n"
			"end of header\n"
		);

		// set output filename
		std::string output_name(filename);
		output_name += ".dc";

		MPI_File outfile;

		// MPI_File_open wants a non-constant string
		char* output_name_c_string = new char [output_name.size() + 1];
		output_name.copy(output_name_c_string, output_name.size());
		output_name_c_string[output_name.size()] = '\0';

		/*
		Contrary to what http://www.open-mpi.org/doc/v1.4/man3/MPI_File_open.3.php writes,
		MPI_File_open doesn't truncate the file with OpenMPI 1.4.1 on Ubuntu, so use a
		fopen call first (http://www.opengroup.org/onlinepubs/009695399/functions/fopen.html)
		*/
		if (rank == 0) {
			FILE* i = fopen(output_name_c_string, "w");
			fflush(i);
			fclose(i);
		}
		MPI_Barrier(comm);

		int result = MPI_File_open(
			comm,
			output_name_c_string,
			MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL,
			&outfile
		);

		if (result != MPI_SUCCESS) {
			char mpi_error_string[MPI_MAX_ERROR_STRING + 1];
			int mpi_error_string_length;
			MPI_Error_string(result, mpi_error_string, &mpi_error_string_length);
			mpi_error_string[mpi_error_string_length + 1] = '\0';
			std::cerr << "Couldn't open file " << output_name_c_string
				<< ": " << mpi_error_string
				<< std::endl;
			// TODO throw an exception instead
			abort();
		}

		delete [] output_name_c_string;

		// figure out how many bytes each process will write and where
		uint64_t bytes = 0, offset = 0;

		// collect data from this process into one buffer for writing
		uint8_t* buffer = NULL;

		const std::vector<uint64_t> cells = grid.get_cells();

		// header
		if (rank == 0) {
			bytes += header.size() * sizeof(char)
				+ 6 * sizeof(double)
				+ 4 * sizeof(uint64_t)
				+ sizeof(uint8_t);
		}

		// bytes of cell data
		bytes += cells.size() * (sizeof(uint64_t) + sizeof(uint32_t) + 5 * sizeof(double));

		buffer = new uint8_t [bytes];

		// header
		if (rank == 0) {
			{
			memcpy(buffer + offset, header.c_str(), header.size() * sizeof(char));
			offset += header.size() * sizeof(char);
			}

			const uint64_t endiannes = 0x1234567890abcdef;
			memcpy(buffer + offset, &endiannes, sizeof(uint64_t));
			offset += sizeof(uint64_t);

			const boost::array<double, 3> grid_start = grid.geometry.get_start();
			memcpy(buffer + offset, grid_start.data(), 3 * sizeof(double));
			offset += 3 * sizeof(double);

			const boost::array<double, 3> cell_length = grid.geometry.get_length(1);
			memcpy(buffer + offset, cell_length.data(), 3 * sizeof(double));
			offset += 3 * sizeof(double);

			const boost::array<uint64_t, 3> grid_length = grid.length.get();
			memcpy(buffer + offset, grid_length.data(), 3 * sizeof(uint64_t));
			offset += 3 * sizeof(uint64_t);

			const uint8_t max_ref_lvl = uint8_t(grid.get_maximum_refinement_level());
			memcpy(buffer + offset, &max_ref_lvl, sizeof(uint8_t));
			offset += sizeof(uint8_t);
		}

		// save cell data
		for (uint64_t i = 0; i < cells.size(); i++) {

			// cell id
			const uint64_t cell = cells[i];
			memcpy(buffer + offset, &cell, sizeof(uint64_t));
			offset += sizeof(uint64_t);

			// process number
			const uint32_t process = grid.get_process(cell);
			memcpy(buffer + offset, &process, sizeof(uint32_t));
			offset += sizeof(uint32_t);

			const CellData* const data = grid[cell];

			// density
			const double density = data->density();
			memcpy(buffer + offset, &density, sizeof(double));
			offset += sizeof(double);

			// max relative difference in density
			const double diff = data->max_diff();
			memcpy(buffer + offset, &diff, sizeof(double));
			offset += sizeof(double);

			// vx
			const double vx = data->vx();
			memcpy(buffer + offset, &vx, sizeof(double));
			offset += sizeof(double);

			// vy
			const double vy = data->vy();
			memcpy(buffer + offset, &vy, sizeof(double));
			offset += sizeof(double);

			// vz
			const double vz = data->vz();
			memcpy(buffer + offset, &vz, sizeof(double));
			offset += sizeof(double);
		}

		std::vector<uint64_t> all_bytes(comm_size, 0);
		if (
			MPI_Allgather(
				&bytes,
				1,
				MPI_UINT64_T,
				&(all_bytes[0]),
				1,
				MPI_UINT64_T,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__ << "MPI_Allgather failed." << std::endl;
			abort();
		}

		// calculate offset of this process in the file
		MPI_Offset mpi_offset = 0;
		for (int i = 0; i < rank; i++) {
			mpi_offset += all_bytes[i];
		}

		MPI_Status status;
		MPI_File_write_at_all(outfile, mpi_offset, (void*)buffer, bytes, MPI_BYTE, &status);
		//if (status...

		delete [] buffer;

		MPI_File_close(&outfile);

		return offset;
	}

};	// class Save

#endif

