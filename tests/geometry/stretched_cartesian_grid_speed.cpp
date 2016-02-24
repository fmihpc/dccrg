/*
Tests the speed of geometry operations with a constant cell size grid
*/

#include "cstdint"
#include "ctime"
#include "iostream"
#include "vector"

#include "dccrg_length.hpp"
#include "dccrg_mapping.hpp"
#include "dccrg_stretched_cartesian_geometry.hpp"
#include "dccrg_topology.hpp"

using namespace std;
using namespace dccrg;

int main()
{
	clock_t before, after;

	cout << "\nArbitrarily streched grid:" << endl;

	const Grid_Topology topology;

	const std::array<uint64_t, 3> grid_length = {{100, 200, 300}};

	Mapping mapping;
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

	Stretched_Cartesian_Geometry geometry(mapping.length, mapping, topology);

	Stretched_Cartesian_Geometry::Parameters parameters;
	for (size_t i = 0; i < grid_length.size(); i++) {
		parameters.coordinates[i].reserve(grid_length[i]);
		for (double coordinate = 0; coordinate <= double(grid_length[i]); coordinate++) {
			parameters.coordinates[i].push_back((1.0 + double(i) / 10) * coordinate);
		}
	}

	if (!geometry.set(parameters)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set grid geometry"
			<< std::endl;
		abort();
	}

	const uint64_t cells = 100000000;
	double
		avg_size_x = 0,
		avg_size_y = 0,
		avg_size_z = 0;

	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		const std::array<double, 3> length = geometry.get_length(cell);
		avg_size_x += length[0];
		avg_size_y += length[1];
		avg_size_z += length[2];
	}
	after = clock();
	avg_size_x /= cells;
	avg_size_y /= cells;
	avg_size_z /= cells;
	cout << "\tAverage cell x, y, z size: "
		<< avg_size_x << " "
		<< avg_size_y << " "
		<< avg_size_z
		<< ", time for " << double(cells) << " cells: "
		<< double(after - before) / CLOCKS_PER_SEC << " s"
		<< endl;


	double
		avg_pos_x = 0,
		avg_pos_y = 0,
		avg_pos_z = 0;

	before = clock();
	for (uint64_t cell = 1; cell <= cells; cell++) {
		const std::array<double, 3> center = geometry.get_center(cell);
		avg_pos_x += center[0];
		avg_pos_y += center[1];
		avg_pos_z += center[2];
	}
	after = clock();
	avg_pos_x /= cells;
	avg_pos_y /= cells;
	avg_pos_z /= cells;
	cout << "\tAverage cell x, y, z position: "
		<< avg_pos_x << " "
		<< avg_pos_y << " "
		<< avg_pos_z
		<< ", time for " << double(cells) << " cells: "
		<< double(after - before) / CLOCKS_PER_SEC << " s"
		<< endl;

	return 0;
}

