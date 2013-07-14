/*
Adaptive mesh refinement version of poisson1d.

Copyright 2013 Finnish Meteorological Institute

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

#include "poisson_solve.hpp"

using namespace boost::mpi;
using namespace dccrg;
using namespace std;


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;

double get_rhs_value(const double x)
{
	if (x < 0 || x > 2 * M_PI) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " x must be in the range [0, 2 * pi]" << x
			<< std::endl;
		abort();
	}

	return sin(x);
}

/*
Returns the analytic solution to the test Poisson's equation.
*/
double get_solution_value(const double x)
{
	if (x < 0 || x > 2 * M_PI) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " x must be in the range [0, 2 * pi]: " << x
			<< std::endl;
		abort();
	}

	return -sin(x);
}

/*
Returns the p-norm of the difference of solution from exact.

Given offset is added to the exact solution before calculating the norm.
*/
template<class Geometry> double get_p_norm(
	const std::vector<uint64_t>& cells,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid,
	const double p_of_norm
) {
	double local = 0, global = 0;

	BOOST_FOREACH(const uint64_t cell, cells) {
		// assumes grid is 1d
		double coord = -1;
		if (grid.length.get()[0] > 1) {
			coord = grid.geometry.get_cell_x(cell);
		}
		if (grid.length.get()[1] > 1) {
			coord = grid.geometry.get_cell_y(cell);
		}
		if (grid.length.get()[2] > 1) {
			coord = grid.geometry.get_cell_z(cell);
		}

		Poisson_Cell* data = grid[cell];
		local += std::pow(
			fabs(data->solution - get_solution_value(coord)),
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

	double old_norm = std::numeric_limits<double>::max();
	size_t old_number_of_cells = 1;

	const size_t max_number_of_cells = 4096;
	for (size_t
		number_of_cells = 2;
		number_of_cells <= max_number_of_cells;
		number_of_cells *= 2
	) {
		const double cell_length = 2 * M_PI / number_of_cells;

		Poisson_Solve solver;
		dccrg::Dccrg<Poisson_Cell> grid_x, grid_y, grid_z, grid_reference;

		grid_x.geometry.set(0, 0, 0, cell_length, 1, 1);
		grid_y.geometry.set(0, 0, 0, 1, cell_length, 1);
		grid_z.geometry.set(0, 0, 0, 1, 1, cell_length);
		grid_reference.geometry.set(0, 0, 0, cell_length, 1, 1);

		const boost::array<uint64_t, 3>
			grid_length_x = {{number_of_cells, 1, 1}},
			grid_length_y = {{1, number_of_cells, 1}},
			grid_length_z = {{1, 1, number_of_cells}};

		grid_x.initialize(grid_length_x, comm, "RCB", 0, -1, true, true, true);
		grid_y.initialize(grid_length_y, comm, "RCB", 0, -1, true, true, true);
		grid_z.initialize(grid_length_z, comm, "RCB", 0, -1, true, true, true);
		grid_reference.initialize(grid_length_x, MPI_COMM_SELF, "RCB", 0, 0, true, true, true);

		// refine every other cell once
		const std::vector<uint64_t> initial_cells = grid_x.get_cells();
		BOOST_FOREACH(const uint64_t cell, initial_cells) {
			if (cell % 2 == 0) {
				grid_x.refine_completely(cell);
				grid_y.refine_completely(cell);
				grid_z.refine_completely(cell);
			}
		}
		grid_x.stop_refining();
		grid_y.stop_refining();
		grid_z.stop_refining();

		// emulate a RANDOM but predictable load balance
		const std::vector<uint64_t>
			refined_cells_x = grid_x.get_cells(),
			refined_cells_y = grid_y.get_cells(),
			refined_cells_z = grid_z.get_cells();

		BOOST_FOREACH(const uint64_t cell, refined_cells_x) {
			const int target_process = cell % comm.size();
			grid_x.pin(cell, target_process);
		}
		BOOST_FOREACH(const uint64_t cell, refined_cells_y) {
			const int target_process = cell % comm.size();
			grid_y.pin(cell, target_process);
		}
		BOOST_FOREACH(const uint64_t cell, refined_cells_z) {
			const int target_process = cell % comm.size();
			grid_z.pin(cell, target_process);
		}
		grid_x.balance_load(false);
		grid_y.balance_load(false);
		grid_z.balance_load(false);
		grid_x.unpin_all_cells();
		grid_y.unpin_all_cells();
		grid_z.unpin_all_cells();

		const std::vector<uint64_t>
			cells_x = grid_x.get_cells(),
			cells_y = grid_y.get_cells(),
			cells_z = grid_z.get_cells(),
			cells_reference = grid_reference.get_cells();

		// initialize
		BOOST_FOREACH(const uint64_t cell, cells_x) {
			Poisson_Cell *data_x = grid_x[cell];

			if (data_x == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_x.geometry.get_cell_x(cell);
			data_x->rhs = get_rhs_value(coord);
			data_x->solution = 0;
		}
		BOOST_FOREACH(const uint64_t cell, cells_y) {
			Poisson_Cell *data_y = grid_y[cell];

			if (data_y == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_y.geometry.get_cell_y(cell);
			data_y->rhs = get_rhs_value(coord);
			data_y->solution = 0;
		}
		BOOST_FOREACH(const uint64_t cell, cells_z) {
			Poisson_Cell *data_z = grid_z[cell];

			if (data_z == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_z.geometry.get_cell_z(cell);
			data_z->rhs = get_rhs_value(coord);
			data_z->solution = 0;
		}
		BOOST_FOREACH(const uint64_t cell, cells_reference) {
			Poisson_Cell *data_reference = grid_reference[cell];

			if (data_reference == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_reference.geometry.get_cell_x(cell);
			data_reference->rhs = get_rhs_value(coord);
			data_reference->solution = 0;
		}

		solver.solve(cells_x, grid_x);
		solver.solve(cells_y, grid_y);
		solver.solve(cells_z, grid_z);
		solver.solve(cells_reference, grid_reference);

		// check that parallel solutions are close to analytic
		const double
			p_of_norm = 2,
			norm_x = get_p_norm(cells_x, grid_x, p_of_norm),
			norm_y = get_p_norm(cells_y, grid_y, p_of_norm),
			norm_z = get_p_norm(cells_z, grid_z, p_of_norm),
			norm_reference = get_p_norm(cells_reference, grid_reference, p_of_norm);

		if (norm_x > old_norm) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of x solution larger from exact with " << number_of_cells
					<< " cells (" << norm_x << ") than with " << old_number_of_cells
					<< " cells (" << old_norm << ")"
					<< std::endl;
			}
		}
		if (norm_y > old_norm) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of y solution larger from exact with " << number_of_cells
					<< " cells (" << norm_y << ") than with " << old_number_of_cells
					<< " cells"
					<< std::endl;
			}
		}
		if (norm_z > old_norm) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of z solution larger from exact with " << number_of_cells
					<< " cells (" << norm_z << ") than with " << old_number_of_cells
					<< " cells"
					<< std::endl;
			}
		}

		// check that AMR solution is better (in some cases) than reference
		if (number_of_cells == 2
		|| number_of_cells == 32
		|| number_of_cells == 512) {
			if (norm_x > norm_reference) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< " norm of x solution (" << norm_x
						<< ") larger than reference (" << norm_reference
						<< ") with " << number_of_cells
						<< " cells"
						<< std::endl;
				}
			}
			if (norm_y > norm_reference) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< " norm of y solution (" << norm_y
						<< ") larger than reference (" << norm_reference
						<< ") with " << number_of_cells
						<< " cells"
						<< std::endl;
				}
			}
			if (norm_z > norm_reference) {
				success = 1;
				if (comm.rank() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " " << p_of_norm
						<< " norm of z solution (" << norm_z
						<< ") larger than reference (" << norm_reference
						<< ") with " << number_of_cells
						<< " cells"
						<< std::endl;
				}
			}
		}

		old_norm = norm_x;
		old_number_of_cells = number_of_cells;

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

