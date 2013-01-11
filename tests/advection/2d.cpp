/*
Advection equation solver program for testing dccrg.

Copyright 2012 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg.  If not, see <http://www.gnu.org/licenses/>.


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
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "iomanip"
#include "iostream"
#include "mpi.h"
#include "string"
#include "zoltan.h"

#define DCCRG_CELL_DATA_SIZE_FROM_USER
#include "../../dccrg.hpp"

#include "cell.hpp"
#include "grid_support.hpp"
#include "save.hpp"

using namespace std;
using namespace boost::mpi;
using namespace dccrg;

#define GRID_START_X 0
#define GRID_START_Y 0
#define GRID_START_Z 0

#define GRID_END_X 1.0
#define GRID_END_Y 1.0
#define GRID_END_Z 1.0


double get_vx(const double y)
{
	return -y + 0.5;
}

double get_vy(const double x)
{
	return +x - 0.5;
}

double get_vz(const double /*a*/)
{
	return 0;
}


string get_output_file_name(const double time_step, const string& basename)
{
	ostringstream step_string;
	step_string << setw(7) << setfill('0') << int(time_step * 1000) << "_ms";
	return basename + step_string.str();
}


/*!
Initializes the simulation's local cells and the copies of remote neighbors.

Only z direction supported at the moment.
*/
template<class CellData> void initial_condition(Dccrg<CellData>& grid)
{
	const vector<uint64_t> cells = grid.get_cells();
	// initialize own cells
	BOOST_FOREACH(const uint64_t& cell_id, cells) {

		CellData* cell = grid[cell_id];
		if (cell == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "No data for cell " << cell_id
				<< std::endl;
			abort();
		}

		for (unsigned int i = 0; i < cell->data.size(); i++) {
			cell->data[i] = 0;
		}

		const double x = grid.get_cell_x(cell_id),
			y = grid.get_cell_y(cell_id),
			//z = grid.get_cell_z(cell_id),
			radius = 0.15;

		// velocities
		cell->vx() = get_vx(y);
		cell->vy() = get_vy(x);
		cell->vz() = get_vz(0);

		/*
		Densities
		*/

		// smooth hump
		const double hump_x0 = 0.25, hump_y0 = 0.5,
			hump_r = min(sqrt(pow(x - hump_x0, 2.0) + pow(y - hump_y0, 2.0)), radius) / radius,
			hump_density = 0.25 * (1 + cos(M_PI * hump_r));

		// TODO: slotted disk
		//const double disk_x0 = 0.5, disk_y0 = 0.75;
		// TODO: rotating cone
		//const double cone_x0 = 0.5, cone_y0 = 0.25;

		cell->density() = hump_density;
	}

	grid.update_remote_neighbor_data();
}


/*!
Calculates fluxes into and out of given local cells.

The total flux to copies of remote neihghbors will be incorrect.
*/
template<class CellData> void solve(
	const double dt,
	const std::vector<uint64_t>& cells,
	Dccrg<CellData>& grid
) {
	BOOST_FOREACH(const uint64_t& cell, cells) {

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

		const double cell_density = cell_data->density(),
			cell_x_size = grid.get_cell_x_size(cell),
			cell_y_size = grid.get_cell_y_size(cell),
			cell_z_size = grid.get_cell_z_size(cell),
			cell_volume = cell_x_size * cell_y_size * cell_z_size;

		vector<uint64_t> neighbors_to_solve;
		vector<direction_t> directions;

		get_face_neighbors<CellData>(cell, grid, neighbors_to_solve, directions);
		// flux between two local cells is solved only in positive direction
		remove_local_negative_neighbors(grid, neighbors_to_solve, directions);

		for (uint64_t i = 0; i < neighbors_to_solve.size(); i++) {

			const uint64_t neighbor = neighbors_to_solve[i];
			// direction indicates where neighbor is located from cell
			const direction_t direction = directions[i];

			CellData* neighbor_data = grid[neighbor];
			if (neighbor_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << neighbor << endl;
				abort();
			}

			const double neighbor_density = neighbor_data->density(),
				neighbor_x_size = grid.get_cell_x_size(neighbor),
				neighbor_y_size = grid.get_cell_y_size(neighbor),
				neighbor_z_size = grid.get_cell_z_size(neighbor),
				neighbor_volume = neighbor_x_size * neighbor_y_size * neighbor_z_size;

			// get area shared between cell and current neighbor
			double min_area = -1;
			switch (direction) {
			case POS_X:
			case NEG_X:
				min_area = min(cell_y_size * cell_z_size, neighbor_y_size * neighbor_z_size);
				break;

			case POS_Y:
			case NEG_Y:
				min_area = min(cell_x_size * cell_z_size, neighbor_x_size * neighbor_z_size);
				break;

			case POS_Z:
			case NEG_Z:
				min_area = min(cell_x_size * cell_y_size, neighbor_x_size * neighbor_y_size);
				break;
			}

			/*
			Solve flux
			*/

			// positive flux through a face goes into positive direction
			double flux = 0;

			// velocity interpolated to shared face
			const double vx = (cell_x_size * neighbor_data->vx() + neighbor_x_size * cell_data->vx()) / (cell_x_size + neighbor_x_size),
				vy = (cell_y_size * neighbor_data->vy() + neighbor_y_size * cell_data->vy()) / (cell_y_size + neighbor_y_size),
				vz = (cell_z_size * neighbor_data->vz() + neighbor_z_size * cell_data->vz()) / (cell_z_size + neighbor_z_size);

			switch (direction) {
			case POS_X:
				if (vx >= 0) {
					flux = cell_density * dt * vx * min_area;
				} else {
					flux = neighbor_density * dt * vx * min_area;
				}
				break;

			case POS_Y:
				if (vy >= 0) {
					flux = cell_density * dt * vy * min_area;
				} else {
					flux = neighbor_density * dt * vy * min_area;
				}
				break;

			case POS_Z:
				if (vz >= 0) {
					flux = cell_density * dt * vz * min_area;
				} else {
					flux = neighbor_density * dt * vz * min_area;
				}
				break;

			case NEG_X:
				if (vx >= 0) {
					flux = neighbor_density * dt * vx * min_area;
				} else {
					flux = cell_density * dt * vx * min_area;
				}
				break;

			case NEG_Y:
				if (vy >= 0) {
					flux = neighbor_density * dt * vy * min_area;
				} else {
					flux = cell_density * dt * vy * min_area;
				}
				break;

			case NEG_Z:
				if (vz >= 0) {
					flux = neighbor_density * dt * vz * min_area;
				} else {
					flux = cell_density * dt * vz * min_area;
				}
				break;
			}

			// save flux
			switch (direction) {
			case POS_X:
			case POS_Y:
			case POS_Z:
				cell_data->flux() -= flux / cell_volume;
				neighbor_data->flux() += flux / neighbor_volume;
				break;

			case NEG_X:
			case NEG_Y:
			case NEG_Z:
				cell_data->flux() += flux / cell_volume;
				neighbor_data->flux() -= flux / neighbor_volume;
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

	BOOST_FOREACH(const uint64_t& cell, cells) {
		CellData* cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "No data for cell " << cell
				<< std::endl;
			abort();
		}

		cell_data->density() += cell_data->flux();
		cell_data->flux() = 0;
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

	BOOST_FOREACH(const uint64_t& cell, cells) {
		CellData* cell_data = grid[cell];
		if (cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< endl;
			abort();
		}

		cell_data->max_diff() = 0;
	}

	// collect maximum relative differences
	BOOST_FOREACH(const uint64_t& cell, cells) {

		CellData* cell_data = grid[cell];

		// get neighbors with which to compare
		vector<uint64_t> neighbors_to_compare;
		vector<direction_t> directions;
		get_face_neighbors<CellData>(cell, grid, neighbors_to_compare, directions);

		BOOST_FOREACH(const uint64_t& neighbor, neighbors_to_compare) {

			CellData* neighbor_data = grid[neighbor];
			if (neighbor_data == NULL) {
				cerr << __FILE__ << ":" << __LINE__
					<< " No data for neighbor " << neighbor
					<< endl;
				abort();
			}

			const double diff = fabs(cell_data->density() - neighbor_data->density())
				/ (min(cell_data->density(), neighbor_data->density()) + diff_threshold);
			cell_data->max_diff() = std::max(diff, cell_data->max_diff());

			// maximize diff for local neighbor
			if (grid.is_local(neighbor)) {
				neighbor_data->max_diff() = std::max(diff, neighbor_data->max_diff());
			}
		}
	}

	// decide whether to refine or unrefine cells
	BOOST_FOREACH(const uint64_t& cell, cells) {

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

		CellData* cell_data = grid[cell];

		const double diff = cell_data->max_diff();

		// refine
		if (diff > refine_diff) {

			cells_to_refine.insert(cell);

			BOOST_FOREACH(const uint64_t& sibling, siblings) {
				cells_to_unrefine.erase(sibling);
				cells_not_to_unrefine.erase(sibling);
			}

		// make sure siblings aren't unrefined
		} else if (diff >= unrefine_diff) {

			bool dont_unrefine = true;

			BOOST_FOREACH(const uint64_t& sibling, siblings) {
				if (cells_to_refine.count(sibling) > 0
				|| cells_not_to_unrefine.count(sibling) > 0) {
					dont_unrefine = false;
					break;
				}
			}

			if (dont_unrefine && grid.get_refinement_level(cell) > 0) {
				cells_not_to_unrefine.insert(cell);

				BOOST_FOREACH(const uint64_t& sibling, siblings) {
					cells_to_unrefine.erase(sibling);
				}
			}

		// unrefine
		} else {

			bool unrefine = true;

			BOOST_FOREACH(const uint64_t& sibling, siblings) {
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
	if (grid.get_maximum_refinement_level() == 0) {
		return std::make_pair(0, 0);
	}

	BOOST_FOREACH(const uint64_t& cell, cells_to_refine) {
		grid.refine_completely(cell);
	}
	cells_to_refine.clear();

	BOOST_FOREACH(const uint64_t& cell, cells_not_to_unrefine) {
		grid.dont_unrefine(cell);
	}
	cells_not_to_unrefine.clear();

	BOOST_FOREACH(const uint64_t& cell, cells_to_unrefine) {
		grid.unrefine_completely(cell);
	}
	cells_to_unrefine.clear();

	/*
	Refines
	*/

	// assign parents' state to children
	const vector<uint64_t> new_cells = grid.stop_refining();

	BOOST_FOREACH(const uint64_t& new_cell, new_cells) {
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

		new_cell_data->density() = parent_data->density();
		new_cell_data->flux() = 0;
		new_cell_data->vx() = get_vx(grid.get_cell_y(new_cell));
		new_cell_data->vy() = get_vy(grid.get_cell_x(new_cell));
		new_cell_data->vz() = 0;
	}

	/*
	Unrefines
	*/

	const vector<uint64_t> removed_cells = grid.get_removed_cells();

	// optimize by gathering all parents of removed cells
	boost::unordered_set<uint64_t> parents;
	BOOST_FOREACH(const uint64_t& removed_cell, removed_cells) {
		parents.insert(grid.get_parent_for_removed(removed_cell));
	}

	// initialize parent data
	BOOST_FOREACH(const uint64_t& parent, parents) {
		CellData* parent_data = grid[parent];
		if (parent_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for parent cell: " << parent
				<< std::endl;
			abort();
		}

		parent_data->density() = 0;
		parent_data->flux() = 0;
		parent_data->vx() = get_vx(grid.get_cell_y(parent));
		parent_data->vy() = get_vy(grid.get_cell_x(parent));
		parent_data->vz() = 0;
	}

	// average parents' density from their children
	BOOST_FOREACH(const uint64_t& removed_cell, removed_cells) {

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

		parent_data->density() += removed_cell_data->density() / 8;
	}

	grid.clear_refined_unrefined_data();

	return std::make_pair(new_cells.size(), removed_cells.size());
}


/*!
Returns the largest allowed global time step.

Must be called simultaneously on all processes.
Assumes that cells of same refinement level have the same size
per dimension.
*/
template<class CellData> double get_max_time_step(
	MPI_Comm& comm,
	const Dccrg<CellData>& grid
) {
	const vector<uint64_t> cells = grid.get_cells();

	double min_step = std::numeric_limits<double>::max();

	BOOST_FOREACH(const uint64_t& cell_id, cells) {

		CellData* cell = grid[cell_id];
		if (cell == NULL) {
			cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell_id
				<< endl;
			abort();
		}

		const double min_step_x = grid.get_cell_x_size(cell_id) / fabs(cell->vx()),
			min_step_y = grid.get_cell_y_size(cell_id) / fabs(cell->vy()),
			min_step_z = grid.get_cell_z_size(cell_id) / fabs(cell->vz()),
			current_min_step = min(min_step_x, min(min_step_y, min_step_z));

		if (min_step > current_min_step) {
			min_step = current_min_step;
		}
	}

	double result = 0;
	if (
		MPI_Allreduce(
			&min_step,
			&result,
			1,
			MPI_DOUBLE,
			MPI_MIN,
			comm
		) != MPI_SUCCESS
	) {
		std::cerr << __FILE__ << ":" << __LINE__ << "MPI_Allreduce failed." << std::endl;
		abort();
	}

	return result;
}


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	/*
	Options
	*/
	//char direction;
	bool verbose = false;
	char direction = 'z';
	unsigned int cells;
	int max_ref_lvl, save_n, balance_n, adapt_n;
	string load_balancing_method;
	double tmax, relative_diff, unrefine_sensitivity, diff_threshold, cfl;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("cells",
			boost::program_options::value<unsigned int>(&cells)->default_value(400),
			"Total number of unrefined cells at the start of the simulation")
		("max-ref-lvl",
			boost::program_options::value<int>(&max_ref_lvl)->default_value(2),
			"Maximum refinement level of cells in the grid (0 means unrefined)")
		("relative-diff",
			boost::program_options::value<double>(&relative_diff)->default_value(0.025),
			"Maximum relative difference in variables for a cell which to keep at maximum refinement level")
		("diff-threshold",
			boost::program_options::value<double>(&diff_threshold)->default_value(0.25),
			"TODO")
		("unrefine-sensitivity",
			boost::program_options::value<double>(&unrefine_sensitivity)->default_value(0.5),
			"TODO")
		("save-n",
			boost::program_options::value<int>(&save_n)->default_value(0),
			"Save results every arg'th time step (0 = only save initial and final result,"
			" -1 = never save)")
		("tmax",
			boost::program_options::value<double>(&tmax)->default_value(25.5),
			"Duration of run in seconds")
		("load-balancing-method",
			boost::program_options::value<string>(&load_balancing_method)->default_value("RCB"),
			"Use arg as load balancing method")
		("balance-n",
			boost::program_options::value<int>(&balance_n)->default_value(25),
			"Balance computational load every argth time step (-1 == never balance load)")
		("adapt-n",
			boost::program_options::value<int>(&adapt_n)->default_value(1),
			"Check for grid adaptation every argth timestep")
		/*("direction",
			boost::program_options::value<char>(&direction)->default_value('z'),
			"Create a 2d grid with normal into direction arg (x, y or z)")*/
		("cfl",
			boost::program_options::value<double>(&cfl)->default_value(0.5),
			"Fraction of ... to use (0..1)")
		("verbose", "Print information during the simulation");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(boost::program_options::parse_command_line(argc, argv, options), option_variables);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Barrier(comm);
		MPI_Finalize();
		return EXIT_SUCCESS;
	}

	if (option_variables.count("verbose") > 0) {
		verbose = true;
	}

	// check simulation parameters
	if (save_n < -1) {
		cerr << "save_n must be >= -1" << endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	if (balance_n < -1) {
		cerr << "balance_n must be >= -1" << endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	if (cfl < 0 || cfl > 1) {
		cerr << "cfl must be >= 0 and <= 1" << endl;
		MPI_Finalize();
		return EXIT_FAILURE;
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}
	if (verbose && rank == 0) {
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
			GRID_START_X, GRID_START_Y, GRID_START_Z,
			GRID_END_X / cells, GRID_END_Y / cells, GRID_END_Z / cells
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
			GRID_START_X, GRID_START_Y, GRID_START_Z,
			GRID_END_X / cells, GRID_END_Y / cells, GRID_END_Z / cells
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
			GRID_START_X, GRID_START_Y, GRID_START_Z,
			GRID_END_X / cells, GRID_END_Y / cells, GRID_END_Z / cells
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

	if (balance_n > -1) {
		grid.balance_load();
	}

	// apply initial condition 1st time for prerefining the grid
	initial_condition<Cell>(grid);

	boost::unordered_set<uint64_t> cells_to_refine, cells_not_to_unrefine, cells_to_unrefine;

	// prerefine up to maximum refinement level
	for (int ref_lvl = 0; ref_lvl < max_ref_lvl; ref_lvl++) {
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

		// apply initial condition on a finer grid
		initial_condition<Cell>(grid);

		grid.update_remote_neighbor_data();
	}

	double dt = get_max_time_step(comm, grid);
	if (verbose && rank == 0) {
		cout << "Initial timestep: " << dt << endl;
	}

	// save initial state
	#ifndef DEBUG
	const string base_output_name("2d_");
	#else
	const string base_output_name("2d_debug_");
	#endif
	unsigned int files_saved = 0;
	if (save_n > -1) {
		if (verbose && rank == 0) {
			cout << "Saving initial state of simulation" << endl;
		}
		Save<Cell>::save(get_output_file_name(0, base_output_name), comm, grid);
		files_saved++;
	}

	if (verbose && rank == 0) {
		cout << "Starting simulation" << endl;
	}

	vector<uint64_t> inner_cells = grid.get_cells_with_local_neighbors();
	vector<uint64_t> outer_cells = grid.get_cells_with_remote_neighbor();

	// record solution time for inner cells and amount of neighbor data received
	double inner_solve_time = 0, outer_solve_time = 0, neighbor_receive_size = 0;

	double time = 0;
	unsigned int step = 0;
	while (time < tmax) {

		if (verbose && rank == 0) {
			cout << "Simulation time: " << time << endl;
		}

		grid.start_remote_neighbor_data_update();

		// solve inner cells
		const double inner_solve_start = MPI_Wtime();
		solve<Cell>(cfl * dt, inner_cells, grid);
		inner_solve_time += MPI_Wtime() - inner_solve_start;

		// wait for remote neighbor data
		grid.wait_neighbor_data_update_receives();

		// solve outer cells
		const double outer_solve_start = MPI_Wtime();
		solve<Cell>(cfl * dt, outer_cells, grid);
		outer_solve_time += MPI_Wtime() - outer_solve_start;

		// wait until local data has been sent
		grid.wait_neighbor_data_update_sends();

		neighbor_receive_size +=
			4 * sizeof(double)
			* (grid.get_number_of_update_receive_cells() + grid.get_number_of_update_send_cells());

		/*
		Starting from this point local cells and copies of remote cells have
		data from the same timestep (flux and max_diff isn't transferred).
		*/

		// check where to adapt the grid
		if (adapt_n > 0 && step % adapt_n == 0) {

			if (verbose && rank == 0) {
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
			if (verbose && rank == 0) {
				cout << "Saving simulation at " << time << endl;
			}
			Save<Cell>::save(get_output_file_name(time, base_output_name), comm, grid);
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

			if (verbose && rank == 0) {
				cout << "Adapting grid" << endl;
			}

			const std::pair<uint64_t, uint64_t> adapted_cells = adapt_grid<Cell>(
				cells_to_refine,
				cells_not_to_unrefine,
				cells_to_unrefine,
				grid
			);
			inner_cells = grid.get_cells_with_local_neighbors();
			outer_cells = grid.get_cells_with_remote_neighbor();

			// update maximum allowed time step
			dt = get_max_time_step(comm, grid);
			if (verbose && rank == 0) {
				cout << "New timestep: " << dt << endl;
			}
		}

		// balance load
		if (balance_n > 0 && step % balance_n == 0) {

			if (verbose && rank == 0) {
				cout << "Balancing load" << endl;
			}

			grid.balance_load();

			inner_cells = grid.get_cells_with_local_neighbors();
			outer_cells = grid.get_cells_with_remote_neighbor();
		}

		step++;
		time += dt;
	}

	if (save_n > -1) {
		if (verbose && rank == 0) {
			cout << "Saving final state of simulation" << endl;
		}

		Save<Cell>::save(get_output_file_name(tmax, base_output_name), comm, grid);
		files_saved++;
	}

	// gather statistics about solving time and tranferred data
	double min_inner_solve_time = 0, max_inner_solve_time = 0, total_inner_solve_time = 0,
		min_outer_solve_time = 0, max_outer_solve_time = 0, total_outer_solve_time = 0,
		min_receive_size = 0, max_receive_size = 0, total_receive_size = 0,
		// fractions of the above
		min_fraction = 0, max_fraction = 0, total_fraction = 0;

	MPI_Reduce(&inner_solve_time, &min_inner_solve_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&inner_solve_time, &max_inner_solve_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&inner_solve_time, &total_inner_solve_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	MPI_Reduce(&outer_solve_time, &min_outer_solve_time, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&outer_solve_time, &max_outer_solve_time, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&outer_solve_time, &total_outer_solve_time, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	MPI_Reduce(&neighbor_receive_size, &min_receive_size, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&neighbor_receive_size, &max_receive_size, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&neighbor_receive_size, &total_receive_size, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	double fraction = neighbor_receive_size / inner_solve_time;
	MPI_Reduce(&fraction, &min_fraction, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
	MPI_Reduce(&fraction, &max_fraction, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	MPI_Reduce(&fraction, &total_fraction, 1, MPI_DOUBLE, MPI_SUM, 0, comm);

	if (rank == 0) {
		cout << endl;
		cout << "Processes used: " << comm_size << endl;
		cout << "Initial grid size: " << cells * cells << endl;
		cout << "Total timesteps calculated: " << tmax << endl;
		cout << "Total files saved: " << files_saved << endl;
		cout << "Inner cell solution time / step (s, avg, max, min):          "
			<< total_inner_solve_time / comm_size / tmax << "\t"
			<< max_inner_solve_time / tmax << "\t"
			<< min_inner_solve_time / tmax << endl;
		cout << "Outer cell solution time / step (s, avg, max, min):          "
			<< total_outer_solve_time / comm_size / tmax << "\t"
			<< max_outer_solve_time / tmax << "\t"
			<< min_outer_solve_time / tmax << endl;
		cout << "Remote neighbor data receive size / step (B, avg, max, min): "
			<< total_receive_size / comm_size / tmax << "\t"
			<< max_receive_size / tmax << "\t"
			<< min_receive_size / tmax << endl;
		cout << "Per process fractions of the above (B / s, avg, max, min):   "
			<< total_fraction / comm_size << "\t"
			<< max_fraction << "\t"
			<< min_fraction << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

