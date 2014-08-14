/*
Tests the parallel Poisson solver implemented on top of dccrg in 2d.

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
double get_solution_value(const double x, const double y)
{
	return sin(x) * cos(2 * y);
}

double get_rhs_value(const double x, const double y)
{
	return -5 * get_solution_value(x, y);
}


/*
Offsets solution in given grid so that average is equal to analytic solution.
*/
template<class Geometry> void normalize_solution(
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Poisson_Cell, Geometry>& grid
) {
	if (cells.size() == 0) {
		return;
	}

	double avg_solved = 0, avg_analytic = 0;
	BOOST_FOREACH(const uint64_t cell, cells) {

		const boost::array<double, 3> cell_center = grid.geometry.get_center(cell);
		if (grid.length.get()[0] == 1) {
			avg_analytic += get_solution_value(
				cell_center[1],
				cell_center[2]
			);
		} else if (grid.length.get()[1] == 1) {
			avg_analytic += get_solution_value(
				cell_center[0],
				cell_center[2]
			);
		} else if (grid.length.get()[2] == 1) {
			avg_analytic += get_solution_value(
				cell_center[0],
				cell_center[1]
			);
		}

		Poisson_Cell* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for last cell " << cell
				<< std::endl;
			abort();
		}

		avg_solved += cell_data->solution;
	}
	avg_analytic /= cells.size();
	avg_solved /= cells.size();

	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for last cell " << cell
				<< std::endl;
			abort();
		}

		cell_data->solution -= avg_solved - avg_analytic;
	}
}


/*
Returns the p-norm of the difference of solution from exact.
*/
template<class Geometry> double get_p_norm(
	const std::vector<uint64_t>& cells,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid,
	const double p_of_norm
) {
	int dimensions = 0;
	if (grid.length.get()[0] > 1) {
		dimensions++;
	}
	if (grid.length.get()[1] > 1) {
		dimensions++;
	}
	if (grid.length.get()[2] > 1) {
		dimensions++;
	}
	if (dimensions != 2) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Invalid number of dimensions: " << dimensions
			<< std::endl;
		abort();
	}

	double local = 0, global = 0;
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* data = grid[cell];
		if (data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for last cell " << cell
				<< std::endl;
			abort();
		}

		double analytic_solution = std::numeric_limits<double>::max();
		const boost::array<double, 3> cell_center = grid.geometry.get_center(cell);

		if (grid.length.get()[0] == 1) {
			analytic_solution = get_solution_value(
				cell_center[1],
				cell_center[2]
			);
		} else if (grid.length.get()[1] == 1) {
			analytic_solution = get_solution_value(
				cell_center[0],
				cell_center[2]
			);
		} else if (grid.length.get()[2] == 1) {
			analytic_solution = get_solution_value(
				cell_center[0],
				cell_center[1]
			);
		}

		local += std::pow(
			fabs(data->solution - analytic_solution),
			p_of_norm
		);
	}
	MPI_Comm temp = grid.get_communicator();
	MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, temp);
	MPI_Comm_free(&temp);
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

	int success = 0, global_success = 0; // > 0 means failure

	double
		old_norm_n_n_x = std::numeric_limits<double>::max(),
		old_norm_n_n_y = std::numeric_limits<double>::max(),
		old_norm_n_n_z = std::numeric_limits<double>::max(),
		old_norm_2n_n_x = std::numeric_limits<double>::max(),
		old_norm_2n_n_y = std::numeric_limits<double>::max(),
		old_norm_2n_n_z = std::numeric_limits<double>::max(),
		old_norm_n_2n_x = std::numeric_limits<double>::max(),
		old_norm_n_2n_y = std::numeric_limits<double>::max(),
		old_norm_n_2n_z = std::numeric_limits<double>::max();

	// check with number of cells: 2, 2; 4, 2; 2, 4; 4, 4; 8, 4; 4, 8; 8, 8; ...
	const size_t max_number_of_cells = 128; // per dimension
	size_t current_cells_min = 4, current_cells_max = 2 * current_cells_min;

	for (size_t cells_y = current_cells_min; cells_y <= current_cells_max; cells_y *= 2)
	for (size_t cells_x = current_cells_min; cells_x <= current_cells_max; cells_x *= 2) {

		if (cells_y == 2 * cells_x) {
			current_cells_min = current_cells_max;
			if (current_cells_max < max_number_of_cells) {
				current_cells_max *= 2;
			}
		}

		const double
			cell_length_x = 2 * M_PI / cells_x,
			cell_length_y = 1 * M_PI / cells_y;

		/*
		Parallel 2d solution in each dimension
		*/

		Poisson_Solve solver;
		//solver.set_verbosity(true);
		Dccrg<Poisson_Cell, Cartesian_Geometry> grid_x, grid_y, grid_z;

		const boost::array<uint64_t, 3>
			grid_length_x = {{1, cells_x, cells_y}},
			grid_length_y = {{cells_x, 1, cells_y}},
			grid_length_z = {{cells_x, cells_y, 1}};

		grid_x.initialize(grid_length_x, comm, "RCB", 0, 0, true, true, true);
		grid_y.initialize(grid_length_y, comm, "RCB", 0, 0, true, true, true);
		grid_z.initialize(grid_length_z, comm, "RCB", 0, 0, true, true, true);

		Cartesian_Geometry::Parameters geom_params;
		geom_params.start[0] =
		geom_params.start[1] =
		geom_params.start[2] = 0;
		geom_params.level_0_cell_length[0] =
		geom_params.level_0_cell_length[1] =
		geom_params.level_0_cell_length[2] = 1;

		geom_params.level_0_cell_length[1] = cell_length_x;
		geom_params.level_0_cell_length[2] = cell_length_y;
		grid_x.set_geometry(geom_params);
		geom_params.level_0_cell_length[1] = 1;
		geom_params.level_0_cell_length[2] = 1;

		geom_params.level_0_cell_length[0] = cell_length_x;
		geom_params.level_0_cell_length[2] = cell_length_y;
		grid_y.set_geometry(geom_params);
		geom_params.level_0_cell_length[0] = 1;
		geom_params.level_0_cell_length[2] = 1;

		geom_params.level_0_cell_length[0] = cell_length_x;
		geom_params.level_0_cell_length[1] = cell_length_y;
		grid_z.set_geometry(geom_params);
		geom_params.level_0_cell_length[0] = 1;
		geom_params.level_0_cell_length[1] = 1;

		const std::vector<uint64_t> initial_cells = grid_x.get_cells();

		// emulate RANDOM load balance but make local cells identical in grid_x, y and z
		BOOST_FOREACH(const uint64_t cell, initial_cells) {
			const int target_process = cell % comm.size();
			grid_x.pin(cell, target_process);
			grid_y.pin(cell, target_process);
			grid_z.pin(cell, target_process);
		}
		grid_x.balance_load(false);
		grid_y.balance_load(false);
		grid_z.balance_load(false);
		grid_x.unpin_all_cells();
		grid_y.unpin_all_cells();
		grid_z.unpin_all_cells();

		const std::vector<uint64_t> cells = grid_x.get_cells();

		// initialize
		BOOST_FOREACH(const uint64_t cell, cells) {
			Poisson_Cell
				*data_x = grid_x[cell],
				*data_y = grid_y[cell],
				*data_z = grid_z[cell];

			if (data_x == NULL || data_y == NULL || data_z == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			const boost::array<double, 3>
				cell_center_x = grid_x.geometry.get_center(cell),
				cell_center_y = grid_y.geometry.get_center(cell),
				cell_center_z = grid_z.geometry.get_center(cell);

			data_x->rhs = get_rhs_value(
				cell_center_x[1],
				cell_center_x[2]
			);
			data_y->rhs = get_rhs_value(
				cell_center_y[0],
				cell_center_y[2]
			);
			data_z->rhs = get_rhs_value(
				cell_center_z[0],
				cell_center_z[1]
			);

			data_x->solution =
			data_y->solution =
			data_z->solution = 0;
		}

		solver.solve(cells, grid_x);
		solver.solve(cells, grid_y);
		solver.solve(cells, grid_z);

		normalize_solution(cells, grid_x);
		normalize_solution(cells, grid_y);
		normalize_solution(cells, grid_z);

		// check that parallel solutions with more cells have smaller norms
		const double
			p_of_norm = 2,
			norm_x  = get_p_norm(cells, grid_x, p_of_norm),
			norm_y  = get_p_norm(cells, grid_y, p_of_norm),
			norm_z  = get_p_norm(cells, grid_z, p_of_norm);

		if (cells_x == cells_y) {

			if (norm_x > old_norm_n_n_x) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between x and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_x
						<< ", should be <= " << old_norm_n_n_x
						<< std::endl;
				}
			}
			old_norm_n_n_x = norm_x;

			if (norm_y > old_norm_n_n_y) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between y and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_y
						<< ", should be <= " << old_norm_n_n_y
						<< std::endl;
				}
			}
			old_norm_n_n_y = norm_y;

			if (norm_z > old_norm_n_n_z) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between z and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_z
						<< ", should be <= " << old_norm_n_n_z
						<< std::endl;
				}
			}
			old_norm_n_n_z = norm_z;

		}

		if (cells_x == 2 * cells_y) {

			if (norm_x > old_norm_n_2n_x) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between x and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_x
						<< ", should be <= " << old_norm_n_2n_x
						<< std::endl;
				}
			}
			old_norm_n_2n_x = norm_x;

			if (norm_y > old_norm_n_2n_y) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between y and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_y
						<< ", should be <= " << old_norm_n_2n_y
						<< std::endl;
				}
			}
			old_norm_n_2n_y = norm_y;

			if (norm_z > old_norm_n_2n_z) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between z and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_z
						<< ", should be <= " << old_norm_n_2n_z
						<< std::endl;
				}
			}
			old_norm_n_2n_z = norm_z;

		}

		if (2 * cells_x == cells_y) {

			if (norm_x > old_norm_2n_n_x) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between x and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_x
						<< ", should be <= " << old_norm_2n_n_x
						<< std::endl;
				}
			}
			old_norm_2n_n_x = norm_x;

			if (norm_y > old_norm_2n_n_y) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between y and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_y
						<< ", should be <= " << old_norm_2n_n_y
						<< std::endl;
				}
			}
			old_norm_2n_n_y = norm_y;

			if (norm_z > old_norm_2n_n_z) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< "-norm between z and analytic is too large with "
						<< cells_x << "x" << cells_y
						<< " cells: " << norm_z
						<< ", should be <= " << old_norm_2n_n_z
						<< std::endl;
				}
			}
			old_norm_2n_n_z = norm_z;
		}

		MPI_Allreduce(&success, &global_success, 1, MPI_INT, MPI_SUM, comm);
		if (global_success > 0) {
			break;
		}
	}

	if (success == 0) {
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

