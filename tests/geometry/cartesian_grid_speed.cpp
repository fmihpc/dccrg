/*
Tests the speed of geometry operations with a constant cell size grid
*/

#include "ctime"
#include "iostream"
#include "stdint.h"

#include "../../dccrg_cartesian_geometry.hpp"
#include "../../dccrg_length.hpp"
#include "../../dccrg_mapping.hpp"
#include "../../dccrg_topology.hpp"

using namespace std;
using namespace dccrg;

int main(void)
{
	clock_t before, after;

	cout << "\nCartesian grid:" << endl;

	const Grid_Topology topology;

	Mapping mapping;
	const boost::array<uint64_t, 3> grid_length = {{100, 200, 300}};
	if (!mapping.set_length(grid_length)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set grid length to "
			<< grid_length[0] << ", "
			<< grid_length[1] << ", "
			<< grid_length[2] << " cells of refinement level 0"
			<< std::endl;
		abort();
	}
	mapping.set_maximum_refinement_level(mapping.get_maximum_possible_refinement_level());
	cout << "\tMaximum refinement level: " << mapping.get_maximum_refinement_level() << endl;

	Cartesian_Geometry geometry(mapping.length, mapping, topology);
	if (!geometry.set(0, 0, 0, 1, 1.1, 1.2)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set grid geometry"
			<< std::endl;
		abort();
	}

	const uint64_t cells = 100000000;
	double avg_size = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		avg_size += geometry.get_cell_length_x(cell);
	}
	after = clock();
	avg_size /= cells;
	cout << "\tAverage cell x size: " << avg_size;
	cout << ", time for " << double(cells) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	avg_size = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		avg_size += geometry.get_cell_length_y(cell);
	}
	after = clock();
	avg_size /= cells;
	cout << "\tAverage cell y size: " << avg_size;
	cout << ", time for " << double(cells) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	avg_size = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		avg_size += geometry.get_cell_length_z(cell);
	}
	after = clock();
	avg_size /= cells;
	cout << "\tAverage cell z size: " << avg_size;
	cout << ", time for " << double(cells) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;



	double avg_pos = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		avg_pos += geometry.get_cell_x(cell);
	}
	after = clock();
	avg_pos /= cells;
	cout << "\tAverage cell x position: " << avg_pos;
	cout << ", time for " << double(cells) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	avg_pos = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		avg_pos += geometry.get_cell_y(cell);
	}
	after = clock();
	avg_pos /= cells;
	cout << "\tAverage cell y position: " << avg_pos;
	cout << ", time for " << double(cells) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	avg_pos = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		avg_pos += geometry.get_cell_z(cell);
	}
	after = clock();
	avg_pos /= cells;
	cout << "\tAverage cell z position: " << avg_pos;
	cout << ", time for " << double(cells) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	return 0;
}
