/*
Advection equation solver program for testing dccrg.

Copyright 2012 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


TODO:
Initial condition is the one in (figures 9.4 - 9.9 of):
LeVeque, R. J., High-resolution conservative algorithms for advection in
incompressible flow, SIAM J. Numer. Anal., 33, 627-665, 1996
but the used solver is probably the simplest possible.
*/

#include "algorithm"
#include "boost/array.hpp"
#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "iomanip"
#include "iostream"
#include "string"
#include "zoltan.h"

#define DCCRG_CELL_DATA_SIZE_FROM_USER
#include "../../dccrg.hpp"

#include "cell.hpp"

using namespace std;
using namespace boost::mpi;
using namespace dccrg;


double get_vx(const double y)
{
	return -y + 0.5;
	//return 0.1;
}

double get_vy(const double x)
{
	return +x - 0.5;
	//return 0.1;
}

double get_vz(void)
{
	return 0;
}


/*!
Initializes the simulation's local cells and the copies of remote neighbors.

Only z direction supported at the moment.
*/
template<class CellData> void initial_condition(Dccrg<CellData>& grid) {
	
	// initialize own cells
	for (typename boost::unordered_map<uint64_t, CellData>::const_iterator
		item = grid.begin();
		item != grid.end();
		item++
	) {
		const uint64_t cell_id = item->first;

		CellData* cell = grid[cell_id];
		if (cell == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "No data for cell " << cell_id
				<< std::endl;
			abort();
		}

		const double x = grid.get_cell_x(cell_id),
			y = grid.get_cell_y(cell_id);

		for (unsigned int i = 0; i < cell->data.size(); i++) {
			cell->data[i] = 0;
		}

		/*if (x > 0.1 && x < 0.3
		&& y > 0.1 && y < 0.3) {
			cell->data[0] = 1;
		}*/

		// initial condition
		const double radius = 0.15;

		// smooth hump
		const double hump_x0 = 0.25, hump_y0 = 0.5,
			hump_r = min(sqrt(pow(x - hump_x0, 2.0) + pow(y - hump_y0, 2.0)), radius) / radius,
			hump_value = 0.25 * (1 + cos(M_PI * hump_r));

		// slotted disk
		//const double disk_x0 = 0.5, disk_y0 = 0.75;
		// rotating cone
		//const double cone_x0 = 0.5, cone_y0 = 0.25;

		cell->data[0] = hump_value;
	}

	grid.update_remote_neighbour_data();
}


string get_output_file_name(const int time_step, const string& basename)
{
	ostringstream step_string;
	step_string << setw(7) << setfill('0') << time_step;
	return basename + step_string.str();
}


/*!
Directions used in the advection solver.
*/
enum direction_t {
	POS_X,
	NEG_X,
	POS_Y,
	NEG_Y,
	POS_Z,
	NEG_Z
};


/*!
Finds all neighbors of given cell into which stuff might advect.

Also fills neighbors' direction from given cell.
Given vectors are cleared before filling.

Assumes flux is solved only in positive directions between local
cells and that negative direction is solved only if neg. dir. cell
isn't local, e.g. fluxes on process boundaries are solved by
both processes.
*/
template<class CellData> void get_neighbor_directions(
	vector<uint64_t>& face_neighbors,
	vector<direction_t>& directions,
	const uint64_t cell,
	const Dccrg<CellData>& grid
) {
	face_neighbors.clear();
	directions.clear();

	face_neighbors.reserve(24);
	directions.reserve(24);

	const int refinement_level = grid.get_refinement_level(cell);

	const vector<uint64_t>* neighbours = grid.get_neighbours(cell);
	if (neighbours == NULL) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " No neighbors for cell " << cell
			<< std::endl;
		abort();
	}

	unsigned int neighbour_index = 0;

	// -z direction
	if ((*neighbours)[neighbour_index] != 0) {

		if (grid.get_refinement_level((*neighbours)[neighbour_index]) <= refinement_level) {

			// flux between local cells is calculated in positive direction
			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Z);
			}

		// only face neighbours
		} else {
			neighbour_index += 4;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Z);
			}
			neighbour_index++;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Z);
			}
			neighbour_index++;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Z);
			}
			neighbour_index++;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Z);
			}
		}
	}
	neighbour_index++;

	// -y direction
	if ((*neighbours)[neighbour_index] != 0) {

		if (grid.get_refinement_level((*neighbours)[neighbour_index]) <= refinement_level) {
			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Y);
			}

		// solve only face neighbours
		} else {
			neighbour_index += 2;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Y);
			}
			neighbour_index++;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Y);
			}
			neighbour_index += 3;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Y);
			}
			neighbour_index++;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_Y);
			}
		}
	}
	neighbour_index++;

	// -x direction
	if ((*neighbours)[neighbour_index] != 0) {

		if (grid.get_refinement_level((*neighbours)[neighbour_index]) <= refinement_level) {
			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_X);
			}

		// solve only face neighbours
		} else {
			neighbour_index++;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_X);
			}
			neighbour_index += 2;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_X);
			}
			neighbour_index += 2;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_X);
			}
			neighbour_index += 2;

			if (!grid.is_local((*neighbours)[neighbour_index])) {
				face_neighbors.push_back((*neighbours)[neighbour_index]);
				directions.push_back(NEG_X);
			}
		}
	}
	neighbour_index++;

	// +x direction
	if ((*neighbours)[neighbour_index] != 0) {

		if (grid.get_refinement_level((*neighbours)[neighbour_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_X);

		// solve only face neighbours
		} else {
			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_X);

			neighbour_index += 2;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_X);

			neighbour_index += 2;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_X);

			neighbour_index += 2;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_X);

			neighbour_index++;
		}
	}
	neighbour_index++;

	// +y direction
	if ((*neighbours)[neighbour_index] != 0) {

		if (grid.get_refinement_level((*neighbours)[neighbour_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Y);

		// solve only face neighbours
		} else {
			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Y);

			neighbour_index++;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Y);

			neighbour_index += 3;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Y);

			neighbour_index++;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Y);

			neighbour_index += 2;
		}
	}
	neighbour_index++;

	// +z direction
	if ((*neighbours)[neighbour_index] != 0) {

		if (grid.get_refinement_level((*neighbours)[neighbour_index]) <= refinement_level) {
			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Z);

		// solve only face neighbours
		} else {
			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Z);

			neighbour_index++;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Z);

			neighbour_index++;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Z);

			neighbour_index++;

			face_neighbors.push_back((*neighbours)[neighbour_index]);
			directions.push_back(POS_Z);
		}
	}

	if (neighbour_index > neighbours->size()) {
		cerr << __FILE__ << ":" << __LINE__
			<< " Added more neighbors (" << neighbour_index
			<< ") than exist: " << neighbours->size()
			<< endl;
		abort();
	}
}


/*!
Calculates fluxes into and out of given local cells.

The total flux to copies of remote neihghbors will be incorrect.
FIXME: AMR case
*/
template<class CellData> void solve(
	const double dt,
	const std::vector<uint64_t>& cells,
	Dccrg<CellData>& grid
) {
	BOOST_FOREACH(uint64_t cell, cells) {

		if (!grid.is_local(cell)) {
			cerr << __FILE__ << ":" << __LINE__
				<< " Cell " << cell
				<< " isn't local"
				<< endl;
			abort();
		}

		CellData* cell_data = grid[cell];
		if (cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
			abort();
		}

		const double cell_value = cell_data->data[0],
			cell_x_size = grid.get_cell_x_size(cell),
			cell_y_size = grid.get_cell_y_size(cell),
			cell_z_size = grid.get_cell_z_size(cell),
			cell_volume = cell_x_size * cell_y_size * cell_z_size;

		vector<uint64_t> neighbors_to_solve;
		vector<direction_t> directions;

		get_neighbor_directions<CellData>(neighbors_to_solve, directions, cell, grid);

		for (uint64_t i = 0; i < neighbors_to_solve.size(); i++) {

			const uint64_t neighbor = neighbors_to_solve[i];
			// direction indicates where neighbor is located from cell
			const direction_t direction = directions[i];

			CellData* neighbor_data = grid[neighbor];
			if (neighbor_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << neighbor << endl;
				abort();
			}

			const double neighbor_value = neighbor_data->data[0],
				neighbor_x_size = grid.get_cell_x_size(neighbor),
				neighbor_y_size = grid.get_cell_y_size(neighbor),
				neighbor_z_size = grid.get_cell_z_size(neighbor),
				neighbor_volume = neighbor_x_size * neighbor_y_size * neighbor_z_size;

			const double volume = min(cell_volume, neighbor_volume);

			// get area shared between cell and current neighbor
			double area = -1;
			switch (direction) {
			case POS_X:
			case NEG_X:
				area = min(cell_y_size * cell_z_size, neighbor_y_size * neighbor_z_size);
				break;

			case POS_Y:
			case NEG_Y:
				area = min(cell_x_size * cell_z_size, neighbor_x_size * neighbor_z_size);
				break;

			case POS_Z:
			case NEG_Z:
				area = min(cell_x_size * cell_y_size, neighbor_x_size * neighbor_y_size);
				break;
			}

			/*
			Solve flux
			*/

			// positive flux through a face goes into positive direction
			double flux = 0;

			// location of face between cell and its current neighbor
			double x = 0, y = 0;

			double value = 0, vx = 0, vy = 0;
			switch (direction) {
			case POS_X:
				y = grid.get_cell_y(cell);
				vx = get_vx(y);
				value = (vx >= 0) ? cell_value : neighbor_value;
				flux = value * (dt * vx * area / volume);
				break;

			case POS_Y:
				x = grid.get_cell_x(cell);
				vy = get_vy(x);
				value = (vy >= 0) ? cell_value : neighbor_value;
				flux = value * (dt * vy * area / volume);
				break;

			case POS_Z:
				// no flux in z direction
				break;

			case NEG_X:
				y = grid.get_cell_y(cell);
				vx = get_vx(y);
				value = (vx >= 0) ? neighbor_value : cell_value;
				flux = value * (dt * vx * area / volume);
				break;

			case NEG_Y:
				x = grid.get_cell_x(cell);
				vy = get_vy(x);
				value = (vy >= 0) ? neighbor_value : cell_value;
				flux = value * (dt * vy * area / volume);
				break;

			case NEG_Z:
				// no flux in z direction
				break;
			}

			// save flux
			switch (direction) {
			case POS_X:
			case POS_Y:
			case POS_Z:
				cell_data->data[1] -= flux;
				neighbor_data->data[1] += flux;
				break;

			case NEG_X:
			case NEG_Y:
			case NEG_Z:
				cell_data->data[1] += flux;
				neighbor_data->data[1] -= flux;
				break;
			}
		}
	}
}


/*!
Applies fluxes to local cells and zeroes the fluxes afterwards.
*/
template<class CellData> void apply_fluxes(Dccrg<CellData>& grid)
{
	const vector<uint64_t> cells = grid.get_cells();

	BOOST_FOREACH(uint64_t cell, cells) {
		CellData* cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "No data for cell " << cell
				<< std::endl;
			abort();
		}

		cell_data->data[0] += cell_data->data[1];
		cell_data->data[1] = 0;
	}
}


/*!
Adapts the grid based on relative difference in variables between neighboring cells.

Refinement level target of cells (RDI = relative diff increase for a cell):
\verbatim
       RDI               |  ref lvl
[0, 1 * diff_increase[   |     0
[1, 2 * diff_increase[   |     1
[2, 3 * diff_increase[   |     2
...
\endverbatim
*/
template<class CellData> void check_for_adaptation(
	const double diff_increase,
	const double diff_threshold,
	const double unrefine_sensitivity,
	boost::unordered_set<uint64_t>& cells_to_refine,
	boost::unordered_set<uint64_t>& cells_not_to_unrefine,
	boost::unordered_set<uint64_t>& cells_to_unrefine,
	Dccrg<CellData>& grid
) {
	if (grid.get_maximum_refinement_level() == 0) {
		return;
	}

	cells_to_refine.clear();
	cells_not_to_unrefine.clear();
	cells_to_unrefine.clear();

	const vector<uint64_t> cells = grid.get_cells();

	// maximum relative difference between local cells and their neighbors
	boost::unordered_map<uint64_t, double> max_diffs;

	BOOST_FOREACH(uint64_t cell, cells) {
		max_diffs[cell] = 0;
	}

	// collect maximum relative differences
	BOOST_FOREACH(uint64_t cell, cells) {

		CellData* cell_data = grid[cell];
		if (cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< endl;
			abort();
		}

		// get neighbors with which to compare
		vector<uint64_t> neighbors_to_compare;
		vector<direction_t> directions;
		get_neighbor_directions<CellData>(neighbors_to_compare, directions, cell, grid);

		BOOST_FOREACH(uint64_t neighbor, neighbors_to_compare) {

			CellData* neighbor_data = grid[neighbor];
			if (neighbor_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " No data for neighbor " << neighbor
					<< endl;
				abort();
			}

			const double diff = fabs(cell_data->data[0] - neighbor_data->data[0]) / (min(cell_data->data[0], neighbor_data->data[0]) + diff_threshold);
			max_diffs.at(cell) = std::max(diff, max_diffs.at(cell));

			// maximize diff for local neighbor
			if (max_diffs.count(neighbor) > 0) {
				max_diffs.at(neighbor) = std::max(diff, max_diffs.at(neighbor));
			}
		}
	}

	// decide whether to refine or unrefine cells
	BOOST_FOREACH(uint64_t cell, cells) {

		const int refinement_level = grid.get_refinement_level(cell);

		// refine / unrefine if max relative difference larger / smaller than:
		const double refine_diff = (refinement_level + 1) * diff_increase;
		const double unrefine_diff = unrefine_sensitivity * refine_diff;

		const vector<uint64_t> siblings = grid.get_all_children(grid.get_parent(cell));
		if (siblings.size() == 0) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No siblings for cell " << cell
				<< endl;
			abort();
		}

		const double diff = max_diffs.at(cell);

		// refine
		if (diff > refine_diff) {

			cells_to_refine.insert(cell);

			BOOST_FOREACH(uint64_t sibling, siblings) {
				cells_to_unrefine.erase(sibling);
				cells_not_to_unrefine.erase(sibling);
			}

		// make sure siblings aren't unrefined
		} else if (diff >= unrefine_diff) {

			bool dont_unrefine = true;

			BOOST_FOREACH(uint64_t sibling, siblings) {
				if (cells_to_refine.count(sibling) > 0
				|| cells_not_to_unrefine.count(sibling) > 0) {
					dont_unrefine = false;
					break;
				}
			}

			if (dont_unrefine && grid.get_refinement_level(cell) > 0) {
				cells_not_to_unrefine.insert(cell);

				BOOST_FOREACH(uint64_t sibling, siblings) {
					cells_to_unrefine.erase(sibling);
				}
			}

		// unrefine
		} else {

			bool unrefine = true;

			BOOST_FOREACH(uint64_t sibling, siblings) {
				if (cells_to_refine.count(sibling) > 0
				|| cells_not_to_unrefine.count(sibling) > 0) {
					unrefine = false;
					break;
				}
			}

			if (unrefine && grid.get_refinement_level(cell) > 0) {
				cells_to_unrefine.insert(cell);
			}			
		}
	}
}


/*!
Refines/Unrefines given cells in the grid.

Returns the number of created and removed cells.
Clears given sets of cells after executing refines.
*/
template<class CellData> std::pair<uint64_t, uint64_t>
adapt_grid(
	boost::unordered_set<uint64_t>& cells_to_refine,
	boost::unordered_set<uint64_t>& cells_not_to_unrefine,
	boost::unordered_set<uint64_t>& cells_to_unrefine,
	Dccrg<CellData>& grid
) {
	BOOST_FOREACH(uint64_t cell, cells_to_refine) {
		grid.refine_completely(cell);
	}
	cells_to_refine.clear();

	BOOST_FOREACH(uint64_t cell, cells_not_to_unrefine) {
		grid.dont_unrefine(cell);
	}
	cells_not_to_unrefine.clear();

	BOOST_FOREACH(uint64_t cell, cells_to_unrefine) {
		grid.unrefine_completely(cell);
	}
	cells_to_unrefine.clear();

	// assign parents' state to children
	const vector<uint64_t> new_cells = grid.stop_refining();

	BOOST_FOREACH(uint64_t new_cell, new_cells) {
		CellData* new_cell_data = grid[new_cell];
		if (new_cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for created cell " << new_cell
				<< std::endl;
			abort();
		}

		CellData* parent_data = grid[grid.get_parent(new_cell)];
		if (parent_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell " << grid.get_parent(new_cell)
				<< std::endl;
			abort();
		}

		new_cell_data->data[0] = parent_data->data[0];
		new_cell_data->data[1] = 0;
	}

	// average parent cell's value from unrefined children
	const vector<uint64_t> removed_cells = grid.get_removed_cells();

	boost::unordered_set<uint64_t> parents;
	BOOST_FOREACH(uint64_t removed_cell, removed_cells) {
		parents.insert(grid.get_parent_for_removed(removed_cell));
	}

	BOOST_FOREACH(uint64_t parent, parents) {
		CellData* parent_data = grid[parent];
		if (parent_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell: " << parent
				<< std::endl;
			abort();
		}

		parent_data->data[0] = 0;
	}

	BOOST_FOREACH(uint64_t removed_cell, removed_cells) {

		CellData* removed_cell_data = grid[removed_cell];
		if (removed_cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for removed cell after unrefining: " << removed_cell
				<< std::endl;
			abort();
		}

		CellData* parent_data = grid[grid.get_parent_for_removed(removed_cell)];
		if (parent_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell after unrefining: " << grid.get_parent_for_removed(removed_cell)
				<< std::endl;
			abort();
		}

		parent_data->data[0] += removed_cell_data->data[0] / 8;
	}

	grid.clear_refined_unrefined_data();

	return std::make_pair(new_cells.size(), removed_cells.size());
}


/*!
Saves the given simulation as a .dc file of the given name.

Returns the number of bytes written by this process.
Must be called simultaneously by all processes.

Data is saved in parallel using one call to MPI_File_write_at_all.
*/
template<class CellData> size_t save(
	const string& filename,
	communicator comm,
	const Dccrg<CellData>& grid
) {
	string header;
	header += "2d advection test file\n\n";
	header += "Data after end of header and a line break:\n";
	header += "1 uint64_t 0x1234567890abcdef for checking endiannes of data\n";
	header += "1 double   grid start coordinate in x direction\n";
	header += "1 double   grid start coordinate in y direction\n";
	header += "1 double   grid start coordinate in z direction\n";
	header += "1 double   x size of unrefined spatial cells\n";
	header += "1 double   y size of unrefined spatial cells\n";
	header += "1 double   z size of unrefined spatial cells\n";
	header += "1 uint64_t length of the grid in unrefined cells in x direction\n";
	header += "1 uint64_t length of the grid in unrefined cells in y direction\n";
	header += "1 uint64_t length of the grid in unrefined cells in z direction\n";
	header += "1 uint8_t  maximum refinement level of the grid\n";
	header += "1 uint64_t cell id\n";
	header += "1 uint32_t cell process number\n";
	header += "1 double   value\n";
	header += "1 double   max relative difference in value between this cell and its neighbors\n";
	header += "1 double   vx\n";
	header += "1 double   vy\n";
	header += "1 double   vz\n";
	header += "1 uint64_t cell id\n";
	header += "...\n";
	header += "end of header\n";

	// store max relative change for both neighboring cells if they are local
	boost::unordered_map<uint64_t, double> max_diff;

	// set output filename
	string output_name(filename);
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
	if (comm.rank() == 0) {
		FILE* i = fopen(output_name_c_string, "w");
		fflush(i);
		fclose(i);
	}
	comm.barrier();

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
		cerr << "Couldn't open file " << output_name_c_string
			<< ": " << mpi_error_string
			<< endl;
		// TODO throw an exception instead
		abort();
	}

	delete [] output_name_c_string;

	// figure out how many bytes each process will write and where
	size_t bytes = 0, offset = 0;

	// collect data from this process into one buffer for writing
	uint8_t* buffer = NULL;

	const vector<uint64_t> cells = grid.get_cells();

	// header
	if (comm.rank() == 0) {
		bytes += header.size() * sizeof(char)
			+ 6 * sizeof(double)
			+ 4 * sizeof(uint64_t)
			+ sizeof(uint8_t);
	}

	// bytes of cell data
	bytes += cells.size() * (sizeof(uint64_t) + sizeof(uint32_t) + 5 * sizeof(double));

	buffer = new uint8_t [bytes];

	// header
	if (comm.rank() == 0) {
		{
		memcpy(buffer + offset, header.c_str(), header.size() * sizeof(char));
		offset += header.size() * sizeof(char);
		}

		const uint64_t endiannes = 0x1234567890abcdef;
		memcpy(buffer + offset, &endiannes, sizeof(uint64_t));
		offset += sizeof(uint64_t);

		const double x_start = grid.get_x_start();
		memcpy(buffer + offset, &x_start, sizeof(double));
		offset += sizeof(double);

		const double y_start = grid.get_y_start();
		memcpy(buffer + offset, &y_start, sizeof(double));
		offset += sizeof(double);

		const double z_start = grid.get_z_start();
		memcpy(buffer + offset, &z_start, sizeof(double));
		offset += sizeof(double);

		const double cell_x_size = grid.get_cell_x_size(1);
		memcpy(buffer + offset, &cell_x_size, sizeof(double));
		offset += sizeof(double);

		const double cell_y_size = grid.get_cell_y_size(1);
		memcpy(buffer + offset, &cell_y_size, sizeof(double));
		offset += sizeof(double);

		const double cell_z_size = grid.get_cell_z_size(1);
		memcpy(buffer + offset, &cell_z_size, sizeof(double));
		offset += sizeof(double);

		const uint64_t x_length = grid.get_x_length();
		memcpy(buffer + offset, &x_length, sizeof(uint64_t));
		offset += sizeof(uint64_t);

		const uint64_t y_length = grid.get_y_length();
		memcpy(buffer + offset, &y_length, sizeof(uint64_t));
		offset += sizeof(uint64_t);

		const uint64_t z_length = grid.get_z_length();
		memcpy(buffer + offset, &z_length, sizeof(uint64_t));
		offset += sizeof(uint64_t);

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

		// value
		const double value = data->data[0];
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);

		// max relative difference in value
		const double diff = data->data[2];
		memcpy(buffer + offset, &diff, sizeof(double));
		offset += sizeof(double);

		// vx
		const double y = grid.get_cell_y(cell),
			vx = get_vx(y);
		memcpy(buffer + offset, &vx, sizeof(double));
		offset += sizeof(double);

		// vy
		const double x = grid.get_cell_y(cell),
			vy = get_vy(x);
		memcpy(buffer + offset, &vy, sizeof(double));
		offset += sizeof(double);

		// vz
		const double vz = get_vz();
		memcpy(buffer + offset, &vz, sizeof(double));
		offset += sizeof(double);
	}

	vector<size_t> all_bytes;
	all_gather(comm, bytes, all_bytes);

	// calculate offset of this process in the file
	MPI_Offset mpi_offset = 0;
	for (int i = 0; i < comm.rank(); i++) {
		mpi_offset += all_bytes[i];
	}

	MPI_Status status;
	MPI_File_write_at_all(outfile, mpi_offset, (void*)buffer, bytes, MPI_BYTE, &status);
	//if (status...

	delete [] buffer;

	MPI_File_close(&outfile);

	return offset;
}


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	/*
	Options
	*/
	//char direction;
	bool verbose = false;
	char direction;
	unsigned int cells, tmax;
	int max_ref_lvl = 0, save_n, balance_n, adapt_n;
	string load_balancing_method;
	double relative_diff, unrefine_sensitivity, diff_threshold;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("cells",
			boost::program_options::value<unsigned int>(&cells)->default_value(10000),
			"Total number of unrefined cells at the start of the simulation")
		/*("max_ref_lvl",
			boost::program_options::value<int>(&max_ref_lvl)->default_value(1),
			"Maximum refinement level of cells in the grid (0 means unrefined)")*/
		("relative-diff",
			boost::program_options::value<double>(&relative_diff)->default_value(0.1),
			"Maximum relative difference in variables for a cell which to keep at maximum refinement level")
		("diff-threshold",
			boost::program_options::value<double>(&diff_threshold)->default_value(0.01),
			"TODO")
		("unrefine-sensitivity",
			boost::program_options::value<double>(&unrefine_sensitivity)->default_value(0.5),
			"TODO")
		("save_n",
			boost::program_options::value<int>(&save_n)->default_value(0),
			"Save results every arg'th time step (0 = only save initial and final result, -1 = never save)")
		("tmax",
			boost::program_options::value<unsigned int>(&tmax)->default_value(5000),
			"Duration of run in time steps")
		("load_balancing_method",
			boost::program_options::value<string>(&load_balancing_method)->default_value("HYPERGRAPH"),
			"Use arg as load balancing method")
		("balance_n",
			boost::program_options::value<int>(&balance_n)->default_value(0),
			"Balance computational load every argth time step (-1 == never balance load)")
		("adapt_n",
			boost::program_options::value<int>(&adapt_n)->default_value(1),
			"Check for grid adaptation every argth timestep")
		("direction",
			boost::program_options::value<char>(&direction)->default_value('z'),
			"Create a 2d grid with normal into direction arg (x, y or z)")
		("verbose", "Print information during the simulation");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (comm.rank() == 0) {
			cout << options << endl;
		}
		comm.barrier();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	// check simulation parameters
	if (save_n < -1) {
		cerr << "save_n must be >= -1" << endl;
		return EXIT_FAILURE;
	}

	if (balance_n < -1) {
		cerr << "balance_n must be >= -1" << endl;
		return EXIT_FAILURE;
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}
	if (verbose && comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	// transform user-given parameters to internal units
	cells = (unsigned int) round(sqrt(double(cells)));

	// initialize grid
	Dccrg<Cell> grid;

	switch (direction) {
	case 'x':
		if (!grid.set_geometry(
			1, cells, cells,
			0, 0, 0,
			1.0 / cells, 1.0 / cells, 1.0 / cells
		)) {
			cerr << __FILE__ << ":" << __LINE__ << ": Couldn't set grid geometry" << endl;
			abort();
		}

		grid.initialize(
			comm,
			load_balancing_method.c_str(),
			// only cells sharing a face are considered neighbors
			0,
			max_ref_lvl,
			false, true, true
		);
		break;

	case 'y':
		if (!grid.set_geometry(
			cells, 1, cells,
			0, 0, 0,
			1.0 / cells, 1.0 / cells, 1.0 / cells
		)) {
			cerr << __FILE__ << ":" << __LINE__ << ": Couldn't set grid geometry" << endl;
			abort();
		}

		grid.initialize(
			comm,
			load_balancing_method.c_str(),
			// only cells sharing a face are considered neighbors
			0,
			max_ref_lvl,
			true, false, true
		);
		break;

	case 'z':
		if (!grid.set_geometry(
			cells, cells, 1,
			0, 0, 0,
			1.0 / cells, 1.0 / cells, 1.0 / cells
		)) {
			cerr << __FILE__ << ":" << __LINE__ << ": Couldn't set grid geometry" << endl;
			abort();
		}

		grid.initialize(
			comm,
			load_balancing_method.c_str(),
			// only cells sharing a face are considered neighbors
			0,
			max_ref_lvl,
			true, true, false
		);
		break;

	default:
		cerr << "Unsupported direction given: " << direction << endl;
		break;
	}

	initial_condition<Cell>(grid);

	// save initial state
	string base_output_name("2d_");
	unsigned int files_saved = 0;
	if (save_n > -1) {
		if (verbose && comm.rank() == 0) {
			cout << "Saving initial state of simulation" << endl;
		}
		save<Cell>(get_output_file_name(0, base_output_name), comm, grid);
		files_saved++;
	}

	if (verbose && comm.rank() == 0) {
		cout << "Starting simulation" << endl;
	}

	// prerefine up to maximum refinement level
	/*for (int ref_lvl = 0; ref_lvl < max_ref_lvl; ref_lvl++) {
		check_for_adaptation<Cell>(
			relative_diff / grid.get_maximum_refinement_level(),
			diff_threshold,
			unrefine_sensitivity,
			cells_to_refine,
			cells_not_to_unrefine,
			cells_to_unrefine,
			grid
		);

		const std::pair<uint64_t, uint64_t> adapted_cells = adapt_grid<Cell>(
			cells_to_refine,
			cells_not_to_unrefine,
			cells_to_unrefine,
			grid
		);
	}*/

	vector<uint64_t> inner_cells = grid.get_cells_with_local_neighbours();
	vector<uint64_t> outer_cells = grid.get_cells_with_remote_neighbour();

	// record solution time for inner cells and amount of neighbor data received
	double inner_solve_time = 0, outer_solve_time = 0, neighbor_receive_size = 0;

	for (unsigned int step = 0; step < tmax; step++) {

		grid.start_remote_neighbour_data_update();

		// solve inner cells
		const double inner_solve_start = MPI_Wtime();
		solve<Cell>(0.001, inner_cells, grid);
		inner_solve_time += MPI_Wtime() - inner_solve_start;

		// wait for remote neighbor data
		grid.wait_neighbour_data_update_receives();

		// solve outer cells
		const double outer_solve_start = MPI_Wtime();
		solve<Cell>(0.001, outer_cells, grid);
		outer_solve_time += MPI_Wtime() - outer_solve_start;

		// wait until local data has been sent
		grid.wait_neighbour_data_update_sends();

		neighbor_receive_size += Cell::size() * (grid.get_number_of_update_receive_cells() + grid.get_number_of_update_send_cells());

		/*
		Starting from this point local cells and copies of remote cells have
		data from the same timestep (variables not fluxes which aren't transferred anyway).
		*/

		// check where to adapt the grid
		boost::unordered_set<uint64_t> cells_to_refine, cells_not_to_unrefine, cells_to_unrefine;
		if (adapt_n > 0 && step % adapt_n == 0) {

			if (verbose && comm.rank() == 0) {
				cout << "Checking which cells to adapt in the grid" << endl;
			}

			check_for_adaptation<Cell>(
				relative_diff / grid.get_maximum_refinement_level(),
				diff_threshold,
				unrefine_sensitivity,
				cells_to_refine,
				cells_not_to_unrefine,
				cells_to_unrefine,
				grid
			);
		}

		// save simulation state
		if (save_n > 0 && step % save_n == 0) {
			if (verbose && comm.rank() == 0) {
				cout << "Saving simulation at step " << step << endl;
			}
			save<Cell>(get_output_file_name(step, base_output_name), comm, grid);
			files_saved++;
		}

		/*
		Up to this point local cells and copies of remote cells have
		data from the same timestep (variables not fluxes which aren't transferred anyway).
		*/

		// apply fluxes
		apply_fluxes<Cell>(grid);

		// adapt the grid
		if (adapt_n > 0 && step % adapt_n == 0) {

			if (verbose && comm.rank() == 0) {
				cout << "Adapting grid" << endl;
			}

			const std::pair<uint64_t, uint64_t> adapted_cells = adapt_grid<Cell>(
				cells_to_refine,
				cells_not_to_unrefine,
				cells_to_unrefine,
				grid
			);
			inner_cells = grid.get_cells_with_local_neighbours();
			outer_cells = grid.get_cells_with_remote_neighbour();
		}

		// balance load
		if (balance_n > 0 && step % balance_n == 0) {

			if (verbose && comm.rank() == 0) {
				cout << "Balancing load" << endl;
			}

			grid.balance_load();

			inner_cells = grid.get_cells_with_local_neighbours();
			outer_cells = grid.get_cells_with_remote_neighbour();
		}
	}

	if (save_n > -1) {
		if (verbose && comm.rank() == 0) {
			cout << "Saving final state of simulation" << endl;
		}

		save<Cell>(get_output_file_name(tmax, base_output_name), comm, grid);
		files_saved++;
	}

	// gather statistics about solving time and tranferred data
	double min_inner_solve_time = 0, max_inner_solve_time = 0, total_inner_solve_time = 0,
		min_outer_solve_time = 0, max_outer_solve_time = 0, total_outer_solve_time = 0,
		min_receive_size = 0, max_receive_size = 0, total_receive_size = 0,
		// fractions of the above
		min_fraction = 0, max_fraction = 0, total_fraction = 0;

	if (comm.rank() == 0) {
		reduce(comm, inner_solve_time, min_inner_solve_time, boost::mpi::minimum<double>(), 0);
		reduce(comm, inner_solve_time, max_inner_solve_time, boost::mpi::maximum<double>(), 0);
		reduce(comm, inner_solve_time, total_inner_solve_time, std::plus<double>(), 0);
		reduce(comm, outer_solve_time, min_outer_solve_time, boost::mpi::minimum<double>(), 0);
		reduce(comm, outer_solve_time, max_outer_solve_time, boost::mpi::maximum<double>(), 0);
		reduce(comm, outer_solve_time, total_outer_solve_time, std::plus<double>(), 0);
		reduce(comm, neighbor_receive_size, min_receive_size, boost::mpi::minimum<double>(), 0);
		reduce(comm, neighbor_receive_size, max_receive_size, boost::mpi::maximum<double>(), 0);
		reduce(comm, neighbor_receive_size, total_receive_size, std::plus<double>(), 0);
		reduce(comm, neighbor_receive_size / inner_solve_time, min_fraction, boost::mpi::minimum<double>(), 0);
		reduce(comm, neighbor_receive_size / inner_solve_time, max_fraction, boost::mpi::maximum<double>(), 0);
		reduce(comm, neighbor_receive_size / inner_solve_time, total_fraction, std::plus<double>(), 0);
	} else {
		reduce(comm, inner_solve_time, boost::mpi::minimum<double>(), 0);
		reduce(comm, inner_solve_time, boost::mpi::maximum<double>(), 0);
		reduce(comm, inner_solve_time, std::plus<double>(), 0);
		reduce(comm, outer_solve_time, boost::mpi::minimum<double>(), 0);
		reduce(comm, outer_solve_time, boost::mpi::maximum<double>(), 0);
		reduce(comm, outer_solve_time, std::plus<double>(), 0);
		reduce(comm, neighbor_receive_size, boost::mpi::minimum<double>(), 0);
		reduce(comm, neighbor_receive_size, boost::mpi::maximum<double>(), 0);
		reduce(comm, neighbor_receive_size, std::plus<double>(), 0);
		reduce(comm, neighbor_receive_size / inner_solve_time, boost::mpi::minimum<double>(), 0);
		reduce(comm, neighbor_receive_size / inner_solve_time, boost::mpi::maximum<double>(), 0);
		reduce(comm, neighbor_receive_size / inner_solve_time, std::plus<double>(), 0);
	}

	if (comm.rank() == 0) {
		cout << endl;
		cout << "Processes used: " << comm.size() << endl;
		cout << "Initial grid size: " << cells * cells << endl;
		cout << "Total timesteps calculated: " << tmax << endl;
		cout << "Total files saved: " << files_saved << endl;
		cout << "Inner cell solution time / step (s, avg, max, min):          "
			<< total_inner_solve_time / comm.size() / tmax << "\t"
			<< max_inner_solve_time / tmax << "\t"
			<< min_inner_solve_time / tmax << endl;
		cout << "Outer cell solution time / step (s, avg, max, min):          "
			<< total_outer_solve_time / comm.size() / tmax << "\t"
			<< max_outer_solve_time / tmax << "\t"
			<< min_outer_solve_time / tmax << endl;
		cout << "Remote neighbor data receive size / step (B, avg, max, min): "
			<< total_receive_size / comm.size() / tmax << "\t"
			<< max_receive_size / tmax << "\t"
			<< min_receive_size / tmax << endl;
		cout << "Per process fractions of the above (B / s, avg, max, min):   "
			<< total_fraction / comm.size() << "\t"
			<< max_fraction << "\t"
			<< min_fraction << endl;
	}

	return EXIT_SUCCESS;
}

