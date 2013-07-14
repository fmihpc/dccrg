/*
Version of poisson1d that skip certain cells as defined by user.

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
		const double coord = grid.geometry.get_cell_x(cell);
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

	const size_t max_number_of_cells = 32768;
	for (size_t
		number_of_cells = 2;
		number_of_cells <= max_number_of_cells;
		number_of_cells *= 2
	) {
		const double cell_length = 2 * M_PI / number_of_cells;

		Poisson_Solve solver(10, 1e-15, 2, 10);
		dccrg::Dccrg<Poisson_Cell> grid, grid_reference;

		// create one layer of cells to skip around x dimension
		grid          .geometry.set(0, 0, 0, cell_length, 1, 1);
		grid_reference.geometry.set(0, 0, 0, cell_length, 1, 1);

		const boost::array<uint64_t, 3>
			grid_length = {{number_of_cells, 3, 3}},
			grid_reference_length = {{number_of_cells, 1, 1}};

		grid          .initialize(grid_length, comm, "RCB", 0, 0, true, false, false);
		grid_reference.initialize(grid_reference_length, comm, "RCB", 0, 0, true, false, false);

		// get cells in which to solve
		const std::vector<uint64_t>
			initial_cells = grid.get_cells(),
			initial_cells_reference = grid_reference.get_cells();

		// emulate RANDOM but predictable load balance
		BOOST_FOREACH(const uint64_t cell, initial_cells) {
			const dccrg::Types<3>::indices_t indices = grid.mapping.get_indices(cell);
			// move cells in middle of grid to same process as in reference grid
			if (indices[1] == 1 && indices[2] == 1) {
				const int target_process = indices[0] % comm.size();
				grid.pin(cell, target_process);
			}
		}
		BOOST_FOREACH(const uint64_t cell, initial_cells_reference) {
			const dccrg::Types<3>::indices_t indices = grid_reference.mapping.get_indices(cell);
			const int target_process = indices[0] % comm.size();
			grid_reference.pin(cell, target_process);
		}
		grid          .balance_load(false);
		grid_reference.balance_load(false);
		grid          .unpin_all_cells();
		grid_reference.unpin_all_cells();

		// get cells in which to solve and boundary cells to skip
		const std::vector<uint64_t>
			cells_reference = grid_reference.get_cells(),
			all_cells = grid.get_cells();

		boost::unordered_set<uint64_t> cells_to_skip;

		std::vector<uint64_t> cells;
		BOOST_FOREACH(const uint64_t cell, all_cells) {
			const dccrg::Types<3>::indices_t indices = grid.mapping.get_indices(cell);
			if (indices[1] == 1 && indices[2] == 1) {
				cells.push_back(cell);
			} else {
				cells_to_skip.insert(cell);
			}
		}

		if (cells.size() != cells_reference.size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Length of cell lists not equal: " << cells.size() << " (";
			BOOST_FOREACH(const uint64_t cell, cells) {
				cout << cell << ",";
			}
			cout << "); should be " << cells_reference.size() << " (";
			BOOST_FOREACH(const uint64_t cell, cells_reference) {
				cout << cell << ",";
			}
			cout << ")" << std::endl;
			abort();
		}

		// also skip some remote neighbors
		const std::vector<uint64_t> remote_neighbors = grid.get_remote_cells_on_process_boundary();
		BOOST_FOREACH(const uint64_t cell, remote_neighbors) {
			const dccrg::Types<3>::indices_t indices = grid.mapping.get_indices(cell);
			if (indices[1] != 1 ||  indices[2] != 1) {
				cells_to_skip.insert(cell);
			}
		}

		// initialize
		BOOST_FOREACH(const uint64_t cell, cells) {
			Poisson_Cell *data = grid[cell];
			if (data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			const double coord = grid.geometry.get_cell_x(cell);
			data->rhs = get_rhs_value(coord);
			data->solution = 0;
		}
		BOOST_FOREACH(const uint64_t cell, cells_reference) {
			Poisson_Cell *data = grid_reference[cell];
			if (data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			const double coord = grid_reference.geometry.get_cell_x(cell);
			data->rhs = get_rhs_value(coord);
			data->solution = 0;
		}

		solver.solve(cells, grid, cells_to_skip);
		solver.solve(cells_reference, grid_reference);

		// check that parallel solutions are close to analytic
		const double
			p_of_norm = 2,
			norm = get_p_norm(cells, grid, p_of_norm),
			norm_reference = get_p_norm(cells_reference, grid_reference, p_of_norm);

		if (norm > old_norm) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of solution larger from exact with " << number_of_cells
					<< " cells (" << norm << ") than with " << old_number_of_cells
					<< " cells (" << old_norm << ")"
					<< std::endl;
			}
		}
		if (norm_reference > old_norm) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of reference solution larger from exact with " << number_of_cells
					<< " cells (" << norm_reference << ") than with " << old_number_of_cells
					<< " cells"
					<< std::endl;
			}
		}

		old_norm = norm;
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

