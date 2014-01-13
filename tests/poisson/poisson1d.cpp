/*
Tests the parallel Poisson solver implemented on top of dccrg in 1d.

Copyright 2012, 2013, 2014 Finnish Meteorological Institute

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
#include "reference_poisson_solve.hpp"

using namespace boost::mpi;
using namespace dccrg;
using namespace std;


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


/*
Returns the p-norm of calculated and reference solutions to Poisson's equation.

Given offset is added to the exact solution before calculating the norm.
*/
template<class Geometry> double get_p_norm(
	const std::vector<uint64_t>& cells,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid,
	const Reference_Poisson_Solve& reference,
	const double p_of_norm
) {
	double local = 0, global = 0;
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* data = grid[cell];
		local += std::pow(fabs(data->solution - reference.get_solution(cell - 1)), p_of_norm);
	}
	MPI_Comm temp = grid.get_communicator();
	MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, temp);
	MPI_Comm_free(&temp);
	global = std::pow(global, 1.0 / p_of_norm);

	return global;
}

/*!
Returns the p-norm between given solutions to Poisson's equation.

Given cells must be local in both grids and both grids must also
have identical communicators and structure except for their orientation.
*/
template<class Geometry> double get_p_norm(
	const std::vector<uint64_t>& cells,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid1,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid2,
	const double p_of_norm
) {
	double local = 0, global = 0;
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell
			*data1 = grid1[cell],
			*data2 = grid2[cell];

		local += std::pow(fabs(data1->solution - data2->solution), p_of_norm);
	}
	MPI_Comm temp = grid1.get_communicator();
	MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, temp);
	MPI_Comm_free(&temp);
	global = std::pow(global, 1.0 / p_of_norm);

	return global;
}


/*!
Offsets the given parallel solution so that it is 0 in the last cell.
*/
template<class Geometry> void offset_solution(
	const uint64_t last_cell,
	const std::vector<uint64_t>& cells,
	const dccrg::Dccrg<Poisson_Cell, Geometry>& grid
) {
	MPI_Comm comm = grid.get_communicator();

	// globally get process with the last cell
	int proc_with_last = 0, proc_local = 0;
	if (grid.is_local(last_cell)) {
		proc_local = grid.get_rank();
	}
	MPI_Allreduce(&proc_local, &proc_with_last, 1, MPI_INT, MPI_SUM, comm);

	// proc with last cell tells others the solution in last cell
	double solution = 0;
	if (grid.is_local(last_cell)) {
		Poisson_Cell* data = grid[last_cell];
		if (data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for last cell " << last_cell
				<< std::endl;
			abort();
		}
		solution = data->solution;
	}
	MPI_Bcast(&solution, 1, MPI_DOUBLE, proc_with_last, comm);

	MPI_Comm_free(&comm);

	// offset solutions
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell *data = grid[cell];
		data->solution -= solution;
	}
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

	const size_t max_number_of_cells = 32768;
	for (size_t
		number_of_cells = 2;
		number_of_cells <= max_number_of_cells;
		number_of_cells *= 2
	) {
		const double cell_length = 2 * M_PI / number_of_cells;

		// reference solution of rhs = sin(x) where x is in range [0, 2 pi]
		Reference_Poisson_Solve reference_solver(number_of_cells, cell_length);
		for (size_t i = 0; i < number_of_cells; i++) {
			reference_solver.get_rhs(i) = sin((i + 0.5) * cell_length);
		}

		reference_solver.solve();

		/*
		Parallel 1d solution in each dimension
		*/

		// using more iterations doesn't help
		Poisson_Solve solver(10, 0, 1e-7, 2, 10, false);
		dccrg::Dccrg<
			Poisson_Cell,
			dccrg::Cartesian_Geometry
		> grid_x, grid_y, grid_z, grid_serial;

		const std::array<uint64_t, 3>
			grid_length_x = {{number_of_cells, 1, 1}},
			grid_length_y = {{1, number_of_cells, 1}},
			grid_length_z = {{1, 1, number_of_cells}};

		grid_x.initialize(grid_length_x, comm, "RCB", 0, 0, true, true, true);
		grid_y.initialize(grid_length_y, comm, "RCB", 0, 0, true, true, true);
		grid_z.initialize(grid_length_z, comm, "RCB", 0, 0, true, true, true);
		grid_serial.initialize(grid_length_x, MPI_COMM_SELF, "RCB", 0, 0, true, true, true);

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start[0] =
		geom_params.start[1] =
		geom_params.start[2] = 0;
		geom_params.level_0_cell_length[0] =
		geom_params.level_0_cell_length[1] =
		geom_params.level_0_cell_length[2] = 1;

		geom_params.level_0_cell_length[0] = cell_length;
		grid_x.set_geometry(geom_params);
		grid_serial.set_geometry(geom_params);
		geom_params.level_0_cell_length[0] = 1;

		geom_params.level_0_cell_length[1] = cell_length;
		grid_y.set_geometry(geom_params);
		geom_params.level_0_cell_length[1] = 1;

		geom_params.level_0_cell_length[2] = cell_length;
		grid_z.set_geometry(geom_params);
		geom_params.level_0_cell_length[2] = 1;

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

		// initialize parallel
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

			data_x->rhs =
			data_y->rhs =
			data_z->rhs = reference_solver.get_rhs(cell - 1);

			data_x->solution =
			data_y->solution =
			data_z->solution = 0;
		}

		// initialize serial
		const std::vector<uint64_t> cells_serial = grid_serial.get_cells();
		BOOST_FOREACH(const uint64_t cell, cells_serial) {
			Poisson_Cell* data = grid_serial[cell];
			if (data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}
			data->rhs = reference_solver.get_rhs(cell - 1);
			data->solution = 0;
		}

		solver.solve(cells, grid_x);
		solver.solve(cells, grid_y);
		solver.solve(cells, grid_z);
		solver.solve(cells_serial, grid_serial);

		offset_solution(number_of_cells, cells, grid_x);
		offset_solution(number_of_cells, cells, grid_y);
		offset_solution(number_of_cells, cells, grid_z);
		offset_solution(number_of_cells, cells_serial, grid_serial);

		// check that serial solution close to parallel
		const double
			p_of_norm = 2,
			// comm of first grid used in Allreduce so serial must be second
			norm_serial = get_p_norm(cells, grid_x, grid_serial, p_of_norm),
			// parallel solver doesn't get closer regardless of stopping residual
			norm_threshold = 3e-7;

		if (norm_serial > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << comm.rank()
					<< ": " << p_of_norm
					<< "-norm between x and serial is too large with " << number_of_cells
					<< " cells: " << norm_serial
					<< std::endl;
			}
		}

		// check that parallel solutions are close to reference and each other
		const double
			norm_x  = get_p_norm(cells, grid_x, reference_solver, p_of_norm),
			norm_y  = get_p_norm(cells, grid_y, reference_solver, p_of_norm),
			norm_z  = get_p_norm(cells, grid_z, reference_solver, p_of_norm),
			norm_xy = get_p_norm(cells, grid_x, grid_y, p_of_norm),
			norm_xz = get_p_norm(cells, grid_x, grid_z, p_of_norm),
			norm_yz = get_p_norm(cells, grid_y, grid_z, p_of_norm);

		if (norm_x > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << comm.rank()
					<< ": " << p_of_norm
					<< "-norm between x and reference is too large with " << number_of_cells
					<< " cells: " << norm_x
					<< std::endl;
			}
		}
		if (norm_y > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << comm.rank()
					<< ": " << p_of_norm
					<< "-norm between y and reference is too large with " << number_of_cells
					<< " cells: " << norm_y
					<< std::endl;
			}
		}
		if (norm_z > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << comm.rank()
					<< ": " << p_of_norm
					<< "-norm between z and reference is too large with " << number_of_cells
					<< " cells: " << norm_z
					<< std::endl;
			}
		}

		if (norm_xy > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< "-norm between x and y is too large: " << norm_xy
					<< std::endl;
			}
		}
		if (norm_xz > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< "-norm between x and z differs too much: " << norm_xz
					<< std::endl;
			}
		}
		if (norm_yz > norm_threshold) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< "-norm between y and z differs too much: " << norm_yz
					<< std::endl;
			}
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

