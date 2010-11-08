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

#include "../dccrg_cell_geometry.hpp"

using namespace std;

int main(int argc, char* argv[])
{
	int result;

	for (int arg = 1; arg < argc; arg++) {

		FILE* infile = fopen(argv[arg], "r");
		if (infile == NULL) {
			cerr << "Couldn't open file " << argv[arg] << endl;
			continue;
		}

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

		double cell_size;
		result = fread(&cell_size, sizeof(double), 1, infile);
		if (result != 1) {
			cerr << "Couldn't read cell_size" << endl;
			exit(EXIT_FAILURE);
		}
		//cout << "cell_size: " << cell_size << endl;

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

		// read in game data
		boost::unordered_map<uint64_t, int> game_data;
		do {
			uint64_t cell;
			result = fread(&cell, sizeof(uint64_t), 1, infile);
			if (result != 1) {
				break;
			}

			uint64_t is_alive;
			result = fread(&is_alive, sizeof(uint64_t), 1, infile);
			if (result != 1) {
				cerr << "Couldn't read is_alive for cell " << cell << endl;
				exit(EXIT_FAILURE);
			}

			//cout << "cell " << cell << ", is_alive " << is_alive << endl;
			game_data[cell] = int(is_alive);
		} while (result == 1);

		CellGeometry geometry;
		geometry.set_x_start(x_start);
		geometry.set_y_start(y_start);
		geometry.set_z_start(z_start);
		geometry.set_cell_x_size(cell_size);
		geometry.set_cell_y_size(cell_size);
		geometry.set_cell_z_size(cell_size);
		geometry.set_x_length(x_length);
		geometry.set_y_length(y_length);
		geometry.set_z_length(z_length);
		geometry.set_maximum_refinement_level(0);

		// write the game data to a .vtk file
		string current_output_name(argv[arg]), suffix(".vtk");
		current_output_name += suffix;

		std::ofstream outfile(current_output_name);
		if (!outfile.is_open()) {
			std::cerr << "Couldn't open file " << current_output_name << std::endl;
			exit(EXIT_FAILURE);
		}

		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << "Game of Life data" << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

		// write cells in a known order
		vector<uint64_t> cells;
		cells.reserve(game_data.size());
		for (boost::unordered_map<uint64_t, int>::const_iterator data = game_data.begin(); data != game_data.end(); data++) {
			cells.push_back(data->first);
		}
		sort(cells.begin(), cells.end());

		// write separate points for every cells' corners
		outfile << "POINTS " << cells.size() * 8 << " float" << std::endl;
		for (unsigned int i = 0; i < cells.size(); i++) {
			outfile << geometry.get_cell_min_x(cells[i]) << " " << geometry.get_cell_min_y(cells[i]) << " " << geometry.get_cell_min_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_max_x(cells[i]) << " " << geometry.get_cell_min_y(cells[i]) << " " << geometry.get_cell_min_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_min_x(cells[i]) << " " << geometry.get_cell_max_y(cells[i]) << " " << geometry.get_cell_min_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_max_x(cells[i]) << " " << geometry.get_cell_max_y(cells[i]) << " " << geometry.get_cell_min_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_min_x(cells[i]) << " " << geometry.get_cell_min_y(cells[i]) << " " << geometry.get_cell_max_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_max_x(cells[i]) << " " << geometry.get_cell_min_y(cells[i]) << " " << geometry.get_cell_max_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_min_x(cells[i]) << " " << geometry.get_cell_max_y(cells[i]) << " " << geometry.get_cell_max_z(cells[i]) << std::endl;
			outfile << geometry.get_cell_max_x(cells[i]) << " " << geometry.get_cell_max_y(cells[i]) << " " << geometry.get_cell_max_z(cells[i]) << std::endl;
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

		outfile << "CELL_DATA " << cells.size() << endl;
		outfile << "SCALARS is_alive int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << game_data[*cell] << endl;
		}

		fclose(infile);
	}

	return EXIT_SUCCESS;
}
