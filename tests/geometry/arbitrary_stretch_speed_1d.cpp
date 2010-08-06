/*
Tests the speed of geometry operations with a constant cell size grid
*/

#include "ctime"
#include "iostream"
#include "stdint.h"
#include "vector"

#define DCCRG_ARBITRARY_STRETCH
#include "../../dccrg_cell_geometry.hpp"

using namespace std;

int main(void)
{
	clock_t before, after;

	#define CELLS 100000000
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	y_coordinates.reserve(CELLS);
	for (int i = 0; i <= CELLS / (1 + 8 + 64 + 512) + 1; i++) {
		x_coordinates.push_back(i);
	}
	y_coordinates.push_back(0);
	y_coordinates.push_back(1);
	z_coordinates.push_back(0);
	z_coordinates.push_back(1);

	CellGeometry geometry(x_coordinates, y_coordinates, z_coordinates);
	cout << "Maximum refinement level: " << geometry.get_maximum_refinement_level() << endl;

	double avg_size = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= CELLS; cell++) {
		avg_size += geometry.get_cell_x_size(cell);
	}
	after = clock();
	avg_size /= CELLS;
	cout << "Arbitrarily stretched grid's average cell x size " << avg_size;
	cout << ", time for " << double(CELLS) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	double avg_x_pos = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= CELLS; cell++) {
		avg_x_pos += geometry.get_cell_x(cell);
	}
	after = clock();
	avg_x_pos /= CELLS;
	cout << "Arbitrarily stretched grid's average cell x position " << avg_x_pos;
	cout << ", time for " << double(CELLS) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	return 0;
}
