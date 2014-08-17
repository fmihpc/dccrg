/*
Tests the parallel Poisson solver implemented on top of dccrg in 3d.

Copyright 2012, 2013, 2014 Finnish Meteorological Institute
Copyright 2014 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "boost/foreach.hpp"
#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/static_assert.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "stdint.h"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

#include "poisson_solve.hpp"


using namespace boost::mpi;
using namespace dccrg;
using namespace std;


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


/*
Returns the analytic solution to the test Poisson's equation.
*/
double get_solution_value(const double x, const double y, const double z)
{
	return sin(x) * cos(2 * y) * sin(z/4);
}

double get_rhs_value(const double x, const double y, const double z)
{
	return -(81.0 / 16.0) * get_solution_value(x, y, z);
}


/*
Returns the p-norm of the difference of solution from exact.

Solution in refined cells is averaged to refinement level 0
cells before calculating norm.
*/
template<class Geometry> double get_p_norm(
	const std::vector<uint64_t>& cells,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid,
	const double p_of_norm
) {
	boost::unordered_map<uint64_t, double> avg_solution;

	BOOST_FOREACH(const uint64_t cell, cells) {
		const Poisson_Cell* const data = grid[cell];
		if (data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for last cell " << cell
				<< std::endl;
			abort();
		}

		const uint64_t lvl_0_parent = grid.mapping.get_level_0_parent(cell);
		if (lvl_0_parent == dccrg::error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		if (avg_solution.count(lvl_0_parent) == 0) {
			avg_solution[lvl_0_parent] = 0;
		}

		const int ref_lvl = grid.get_refinement_level(cell);
		if (ref_lvl < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		avg_solution.at(lvl_0_parent) += data->solution / std::pow(8.0, double(ref_lvl));
	}

	double local = 0, global = 0;
	for (boost::unordered_map<uint64_t, double>::const_iterator
		item = avg_solution.cbegin();
		item != avg_solution.cend();
		item++
	) {
		const uint64_t cell = item->first;
		const double solution = item->second;

		const boost::array<double, 3> cell_center = grid.geometry.get_center(cell);
		const double analytic_solution = get_solution_value(
			cell_center[0],
			cell_center[1],
			cell_center[2]
		);

		local += std::pow(
			std::fabs(solution - analytic_solution),
			p_of_norm
		);
	}
	MPI_Comm comm = grid.get_communicator();
	MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Comm_free(&comm);
	global = std::pow(global, 1.0 / p_of_norm);

	return global;
}


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cerr << "Zoltan_Initialize failed" << endl;
	    return EXIT_FAILURE;
	}

	const uint64_t
		cells_x = 8,
		cells_y = 8,
		cells_z = 8;

	const double
		cell_length_x = 2 * M_PI / cells_x,
		cell_length_y = 1 * M_PI / cells_y,
		cell_length_z = 8 * M_PI / cells_z;

	Dccrg<Poisson_Cell, Cartesian_Geometry> grid;

	const boost::array<uint64_t, 3> grid_size = {{cells_x, cells_y, cells_z}};
	grid.initialize(grid_size, comm, "RANDOM", 0, -1, true, true, true);

	Cartesian_Geometry::Parameters geom_params;
	geom_params.start[0] =
	geom_params.start[1] =
	geom_params.start[2] = 0;
	geom_params.level_0_cell_length[0] = cell_length_x;
	geom_params.level_0_cell_length[1] = cell_length_y;
	geom_params.level_0_cell_length[2] = cell_length_z;
	grid.set_geometry(geom_params);

	// refine after balance so all children stay on the same process
	grid.balance_load();

	// refine middle cells a few times
	for (size_t count = 0; count < 2; count++) {
		const std::vector<uint64_t> temp_cells = grid.get_cells();

		BOOST_FOREACH(const uint64_t cell, temp_cells) {
			const boost::array<double, 3>
				cell_start = grid.geometry.get_min(cell),
				cell_end = grid.geometry.get_max(cell);

			if (
				(cell_start[0] < 1.01 * M_PI and cell_end[0] > 0.99 * M_PI)
				and (cell_start[1] < 0.51 * M_PI and cell_end[1] > 0.49 * M_PI)
				and (cell_start[2] < 4.01 * M_PI and cell_end[2] > 3.99 * M_PI)
			) {
				grid.refine_completely(cell);
			}
		}

		grid.stop_refining();
	}

	const std::vector<uint64_t> cells = grid.get_cells();

	// initialize
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* const cell_data = grid[cell];

		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< std::endl;
			abort();
		}

		const boost::array<double, 3> cell_center = grid.geometry.get_center(cell);

		cell_data->solution = 0;
		cell_data->rhs = get_rhs_value(
			cell_center[0],
			cell_center[1],
			cell_center[2]
		);
	}

	Poisson_Solve solver;
	solver.solve(cells, grid);

	const double
		p_of_norm = 2,
		norm  = get_p_norm(cells, grid, p_of_norm);

	if (comm.rank() == 0) {
		//cout << "Norm: " << norm << endl;
	}

	if (norm < 0.35) {
		if (comm.rank() == 0) {
			cout << "PASSED" << endl;
		}
		return EXIT_SUCCESS;
	} else {
		if (comm.rank() == 0) {
			cout << "FAILED" << endl;
		}
		return EXIT_FAILURE;
	}
}

