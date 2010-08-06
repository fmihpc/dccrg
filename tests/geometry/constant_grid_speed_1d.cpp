/*
Tests the speed of geometry operations with a constant cell size grid
*/

#include "ctime"
#include "iostream"
#include "stdint.h"

#include "../../dccrg_cell_geometry.hpp"

using namespace std;

int main(void)
{
	clock_t before, after;

	#define CELLS 100000000
	CellGeometry geometry(0, 0, 0, 1, 1, 1, CELLS / (1 + 8 + 64 + 512) + 1, 1, 1);
	cout << "Maximum refinement level: " << geometry.get_maximum_refinement_level() << endl;

	double avg_size = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= CELLS; cell++) {
		avg_size += geometry.get_cell_x_size(cell);
	}
	after = clock();
	avg_size /= CELLS;
	cout << "Constant grid average cell x size " << avg_size;
	cout << ", time for " << double(CELLS) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	double avg_x_pos = 0;
	before = clock();
	for (uint64_t cell = 1; cell <= CELLS; cell++) {
		avg_x_pos += geometry.get_cell_x(cell);
	}
	after = clock();
	avg_x_pos /= CELLS;
	cout << "Constant grid average cell x position " << avg_x_pos;
	cout << ", time for " << double(CELLS) << " cells: " << double(after - before) / CLOCKS_PER_SEC << " s" << endl;

	return 0;
}
