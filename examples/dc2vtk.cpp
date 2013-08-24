/*
Turns the files saved from game_of_life_with_output.cpp into .vtk files,
visualize for example with VisIt (https://wci.llnl.gov/codes/visit/)
*/

#include "algorithm"
#include "boost/unordered_map.hpp"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "iostream"
#include "stdint.h"

#include "../dccrg_length.hpp"
#include "../dccrg_mapping.hpp"
#include "../dccrg_no_geometry.hpp"
#include "../dccrg_topology.hpp"

using namespace std;
using namespace dccrg;

int main(int argc, char* argv[])
{
	size_t result;

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

		// read in game data
		boost::unordered_map<uint64_t, uint64_t> game_data;
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

			game_data[cell] = is_alive;
		} while (result == 1);

		// use default topology where the grid isn't periodic
		const Grid_Topology topology;

		// mapping of cell ids to their (logical) size and location
		Mapping mapping;
		const boost::array<uint64_t, 3> grid_length = {{10, 10, 1}};
		mapping.set_length(grid_length);
		mapping.set_maximum_refinement_level(max_ref_level);

		// No_Geometry doesn't change after being constructed
		No_Geometry geometry(mapping.length, mapping, topology);

		// write the game data to a .vtk file
		const string input_name(argv[arg]),
			current_output_name(input_name.substr(0, input_name.size() - 2) + "vtk");

		std::ofstream outfile(current_output_name.c_str());
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
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			data = game_data.begin();
			data != game_data.end();
			data++
		) {
			cells.push_back(data->first);
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
