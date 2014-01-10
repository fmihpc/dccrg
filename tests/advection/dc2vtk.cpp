/*
Turns .dc files produced by advection test into .vtk files, visualize for
example with VisIt (https://wci.llnl.gov/codes/visit/)

Copyright 2012, 2013, 2014 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "algorithm"
#include "boost/unordered_map.hpp"
#include "cmath"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "iostream"
#include "stdint.h"

#include "../../dccrg_cartesian_geometry.hpp"

#include "cell.hpp"

using namespace std;
using namespace dccrg;

int main(int argc, char* argv[])
{
	int result;

	if (argc == 1) {
		cerr << "Give at least one .dc file as an argument" << endl;
		exit(EXIT_FAILURE);
	}

	for (int arg = 1; arg < argc; arg++) {

		FILE* infile = fopen(argv[arg], "r");
		if (infile == NULL) {
			cerr << "Couldn't open file " << argv[arg] << endl;
			continue;
		}

		// find end of header
		const string eoh("end of header\n");
		string line;
		bool header_found = false;

		char c = fgetc(infile);
		while (c != EOF) {
			line += c;

			if (c == '\n') {
				if (line == eoh) {
					header_found = true;
					break;
				} else {
					line.clear();
				}
			}

			c = fgetc(infile);
		}

		if (!header_found) {
			cerr << "End of header not found in file " << argv[arg] << endl;
			fclose(infile);
			continue;
		}

		uint64_t endiannes;
		result = fread(&endiannes, sizeof(uint64_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read endiannes" << endl;
			exit(EXIT_FAILURE);
		}

		if (endiannes != 0x1234567890abcdef) {
			cerr << "Unsupported endiannes in file " << argv[arg]
				<< " read " << hex << endiannes
				<< " but should've read 0x1234567890abcdef"
				<< endl;
			fclose(infile);
			continue;
		}

		// read grid data
		double x_start;
		result = fread(&x_start, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read x_start" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "x_start: " << x_start << endl;

		double y_start;
		result = fread(&y_start, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read y_start" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "y_start: " << y_start  << endl;

		double z_start;
		result = fread(&z_start, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read z_start" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "z_start: " << z_start  << endl;

		double cell_x_size;
		result = fread(&cell_x_size, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell_x_size" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "cell_x_size: " << cell_x_size << endl;

		double cell_y_size;
		result = fread(&cell_y_size, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell_y_size" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "cell_y_size: " << cell_y_size << endl;

		double cell_z_size;
		result = fread(&cell_z_size, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell_z_size" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "cell_z_size: " << cell_z_size << endl;

		uint64_t x_length;
		result = fread(&x_length, sizeof(uint64_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read x_length" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "x_length: " << x_length << endl;

		uint64_t y_length;
		result = fread(&y_length, sizeof(uint64_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read y_length" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "y_length: " << y_length << endl;

		uint64_t z_length;
		result = fread(&z_length, sizeof(uint64_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read z_length" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "z_length: " << z_length << endl;

		uint8_t max_ref_level;
		result = fread(&max_ref_level, sizeof(uint8_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read maximum refinemen level" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "max_ref_level: " << uint16_t(max_ref_level) << endl;	// uint8_t == unsigned char so formatting fails

		// read in cell data
		boost::unordered_map<uint64_t, Cell> data;
		boost::unordered_map<uint64_t, uint16_t> cell_process;
		boost::unordered_map<uint64_t, boost::array<double, 3> > velocity;
		do {
			uint64_t cell;
			result = fread(&cell, sizeof(uint64_t), 1, infile);
			if (result != 1) {
				break;
			}

			if (cell == 0) {
				cerr << "Got an invalid cell from file" << endl;
				exit(EXIT_FAILURE);
			}

			cell_process[cell];
			data[cell];

			// cell process
			uint32_t process;
			result = fread(&process, sizeof(uint32_t), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read process of cell " << cell << endl;
				exit(EXIT_FAILURE);
			}
			cell_process[cell] = process;

			double value;

			// density
			result = fread(&value, sizeof(double), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read density in cell " << cell << endl;
				exit(EXIT_FAILURE);
			}
			data[cell].data[0] = value;

			// max relative diff
			result = fread(&value, sizeof(double), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read max relative difference in cell " << cell << endl;
				exit(EXIT_FAILURE);
			}
			data[cell].data[2] = value;

			// vx
			result = fread(&value, sizeof(double), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read vx in cell " << cell << endl;
				exit(EXIT_FAILURE);
			}
			velocity[cell][0] = value;

			// vy
			result = fread(&value, sizeof(double), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read vy in cell " << cell << endl;
				exit(EXIT_FAILURE);
			}
			velocity[cell][1] = value;

			// vz
			result = fread(&value, sizeof(double), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read vz in cell " << cell << endl;
				exit(EXIT_FAILURE);
			}
			velocity[cell][2] = value;

		} while (result == 1);

		const Grid_Topology topology;

		Mapping mapping;
		const boost::array<uint64_t, 3> grid_length = {{x_length, y_length, z_length}};
		if (!mapping.set_length(grid_length)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't set grid length to "
				<< grid_length[0] << ", "
				<< grid_length[1] << ", "
				<< grid_length[2] << " cells of refinement level 0"
				<< std::endl;
			abort();
		}
		if (!mapping.set_maximum_refinement_level(max_ref_level)) {
			std::cerr << "Couldn't set maximum refinement level of grid." << std::endl;
			exit(EXIT_FAILURE);
		}

		Cartesian_Geometry geometry(mapping.length, mapping, topology);

		Cartesian_Geometry::Parameters parameters;
		parameters.start[0] = x_start;
		parameters.start[1] = y_start;
		parameters.start[2] = z_start;
		parameters.level_0_cell_length[0] = cell_x_size;
		parameters.level_0_cell_length[1] = cell_y_size;
		parameters.level_0_cell_length[2] = cell_z_size;

		if (!geometry.set(parameters)) {
			std::cerr << "Couldn't set grid geometry." << std::endl;
			exit(EXIT_FAILURE);
		}

		// write grid data to a .vtk file
		const string input_name(argv[arg]),
			current_output_name(input_name.substr(0, input_name.size() - 2) + "vtk");

		std::ofstream outfile(current_output_name.c_str());
		if (!outfile.is_open()) {
			std::cerr << "Couldn't open file " << current_output_name << std::endl;
			exit(EXIT_FAILURE);
		}

		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << "GUMICS data" << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

		// write cells in a known order
		vector<uint64_t> cells;
		cells.reserve(data.size());
		for (boost::unordered_map<uint64_t, Cell>::const_iterator
			item = data.begin();
			item != data.end();
			item++
		) {
			cells.push_back(item->first);
		}
		sort(cells.begin(), cells.end());

		// write separate points for every cells' corners
		outfile << "POINTS " << cells.size() * 8 << " float" << std::endl;
		for (unsigned int i = 0; i < cells.size(); i++) {
			const boost::array<double, 3>
				cell_min = geometry.get_min(cells[i]),
				cell_max = geometry.get_max(cells[i]);

			outfile
				<< cell_min[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
				<< cell_max[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
				<< cell_min[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
				<< cell_max[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
				<< cell_min[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
				<< cell_max[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
				<< cell_min[0] << " " << cell_max[1] << " " << cell_max[2] << "\n"
				<< cell_max[0] << " " << cell_max[1] << " " << cell_max[2] << "\n";
		}

		// map cells to written points
		outfile << "CELLS " << cells.size() << " " << cells.size() * 9 << std::endl;
		for (unsigned int j = 0; j < cells.size(); j++) {
			outfile << "8 ";
			for (int i = 0; i < 8; i++) {
				 outfile << j * 8 + i << " ";
			}
			outfile << std::endl;
		}

		// cell types
		outfile << "CELL_TYPES " << cells.size() << std::endl;
		for (unsigned int i = 0; i < cells.size(); i++) {
			outfile << 11 << std::endl;
		}

		/*
		Write simulation data
		*/
		outfile << "CELL_DATA " << cells.size() << endl;

		// cells' process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << cell_process.at(*cell) << "\n";
		}

		// cells' refinement level
		outfile << "SCALARS refinement_level int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << mapping.get_refinement_level(*cell) << "\n";
		}

		// value
		outfile << "SCALARS density float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			if (fabs(data.at(*cell).data[0]) > 1e-37) {
				outfile << data.at(*cell).data[0] << "\n";
			} else {
				outfile << "0\n";
			}
		}

		// max relative difference
		outfile << "SCALARS max_diff float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			if (fabs(data.at(*cell).data[2]) > 1e-37) {
				outfile << data.at(*cell).data[2] << "\n";
			} else {
				outfile << "0\n";
			}
		}

		// velocity
		outfile << "VECTORS velocity float" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			if (fabs(velocity.at(*cell)[0]) > 1e-37) {
				outfile << velocity.at(*cell)[0] << " ";
			} else {
				outfile << "0 ";
			}

			if (fabs(velocity.at(*cell)[1]) > 1e-37) {
				outfile << velocity.at(*cell)[1] << " ";
			} else {
				outfile << "0 ";
			}

			if (fabs(velocity.at(*cell)[2]) > 1e-37) {
				outfile << velocity.at(*cell)[2] << " ";
			} else {
				outfile << "0 ";
			}

			outfile << "\n";
		}

		outfile.close();
		fclose(infile);
	}

	return EXIT_SUCCESS;
}
