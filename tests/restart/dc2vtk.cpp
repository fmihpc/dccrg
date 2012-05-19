/*
Turns the files saved from game_of_life_with_output.cpp into .vtk files, visualize for example with VisIt (https://wci.llnl.gov/codes/visit/)
*/

#include "algorithm"
#include "boost/unordered_map.hpp"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "iostream"
#include "stdint.h"

#include "../../dccrg_constant_geometry.hpp"

using namespace std;
using namespace dccrg;

int main(int argc, char* argv[])
{
	int result;

	for (int arg = 1; arg < argc; arg++) {

		FILE* infile = fopen(argv[arg], "r");
		if (infile == NULL) {
			cerr << "Couldn't open file " << argv[arg] << endl;
			continue;
		}

		uint64_t time_step;
		result = fread(&time_step, sizeof(uint64_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read time step" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "time step: " << time_step << endl;

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

		int max_ref_level;
		result = fread(&max_ref_level, sizeof(int), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read maximum refinement level" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "max_ref_level: " << max_ref_level << endl;

		uint64_t number_of_cells;
		result = fread(&number_of_cells, sizeof(uint64_t), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read number of cells" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "number of cells: " << number_of_cells << endl;

		// read in game data (cell, is_alive)
		std::vector<std::pair<uint64_t, uint64_t> > cell_data(number_of_cells);

		// read in cell list
		for (uint64_t i = 0; i < number_of_cells; i++) {
			result = fread(&(cell_data[i].first), sizeof(uint64_t), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read id of " << i + 1 << "th cell" << endl;
				exit(EXIT_FAILURE);
			}

			// skip data offset
			uint64_t temp;
			result = fread(&temp, sizeof(uint64_t), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read data offset of cell " << cell_data[i].first << endl;
				exit(EXIT_FAILURE);
			}
		}

		// read in cell data
		for (uint64_t i = 0; i < number_of_cells; i++) {
			result = fread(&(cell_data[i].second), sizeof(uint64_t), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read is_alive of cell " << cell_data[i].first << endl;
				exit(EXIT_FAILURE);
			}
		}

		ConstantGeometry geometry;
		geometry.set_geometry(
			x_length, y_length, z_length,
			x_start, y_start, z_start,
			cell_x_size, cell_y_size, cell_z_size
		);
		geometry.set_maximum_refinement_level(max_ref_level);

		// write the game data to a .vtk file
		const string input_name(argv[arg]),
			current_output_name(input_name.substr(0, input_name.size() - 2) + "vtk");

		std::ofstream outfile(current_output_name.c_str());
		if (!outfile.is_open()) {
			std::cerr << "Couldn't open file " << current_output_name << std::endl;
			exit(EXIT_FAILURE);
		}

		outfile
			<< "# vtk DataFile Version 2.0\n"
			<< "Game of Life data\n"
			<< "ASCII\n"
			<< "DATASET UNSTRUCTURED_GRID"
			<< std::endl;

		// write separate points for every cells' corners
		outfile << "POINTS " << cell_data.size() * 8 << " float" << std::endl;
		for (uint64_t i = 0; i < cell_data.size(); i++) {
			const uint64_t cell = cell_data[i].first;
			outfile << geometry.get_cell_x_min(cell) << " " << geometry.get_cell_y_min(cell) << " " << geometry.get_cell_z_min(cell) << std::endl;
			outfile << geometry.get_cell_x_max(cell) << " " << geometry.get_cell_y_min(cell) << " " << geometry.get_cell_z_min(cell) << std::endl;
			outfile << geometry.get_cell_x_min(cell) << " " << geometry.get_cell_y_max(cell) << " " << geometry.get_cell_z_min(cell) << std::endl;
			outfile << geometry.get_cell_x_max(cell) << " " << geometry.get_cell_y_max(cell) << " " << geometry.get_cell_z_min(cell) << std::endl;
			outfile << geometry.get_cell_x_min(cell) << " " << geometry.get_cell_y_min(cell) << " " << geometry.get_cell_z_max(cell) << std::endl;
			outfile << geometry.get_cell_x_max(cell) << " " << geometry.get_cell_y_min(cell) << " " << geometry.get_cell_z_max(cell) << std::endl;
			outfile << geometry.get_cell_x_min(cell) << " " << geometry.get_cell_y_max(cell) << " " << geometry.get_cell_z_max(cell) << std::endl;
			outfile << geometry.get_cell_x_max(cell) << " " << geometry.get_cell_y_max(cell) << " " << geometry.get_cell_z_max(cell) << std::endl;
		}

		// map cells to written points
		outfile << "CELLS " << cell_data.size() << " " << cell_data.size() * 9 << "\n";
		for (unsigned int j = 0; j < cell_data.size(); j++) {
			outfile << "8 ";
			for (int i = 0; i < 8; i++) {
				 outfile << j * 8 + i << " ";
			}
			outfile << "\n";
		}

		// cell types
		outfile << "CELL_TYPES " << cell_data.size() << "\n";
		for (unsigned int i = 0; i < cell_data.size(); i++) {
			outfile << 11 << "\n";
		}

		outfile << "CELL_DATA " << cell_data.size() << "\n";
		outfile << "SCALARS is_alive int 1" << "\n";
		outfile << "LOOKUP_TABLE default" << "\n";
		for (uint64_t i = 0; i < cell_data.size(); i++) {
			outfile << cell_data[i].second << "\n";
		}

		fclose(infile);
	}

	return EXIT_SUCCESS;
}
