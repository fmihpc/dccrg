/*
Adaptive mesh refinement version of poisson1d.

Copyright 2013, 2014, 2015 Finnish Meteorological Institute
Copyright 2014, 2015 Ilja Honkonen

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

#include "boost/program_options.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "stdint.h"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

#include "poisson_solve.hpp"

using namespace dccrg;
using namespace std;


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


double get_solution_value(const double x)
{
	return std::pow(sin(x/2), 2.0);
}

double get_rhs_value(const double x)
{
	return 0.5 * (std::pow(cos(x/2), 2.0) - std::pow(sin(x/2), 2.0));
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

	double avg_solved = 0, avg_analytic = 0, divisor = 0;
	for (const auto& cell: cells) {

		const int ref_lvl = grid.get_refinement_level(cell);
		if (ref_lvl < 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid refinement level for cell " << cell
				<< std::endl;
			abort();
		}
		const double value_factor = 1.0 / std::pow(double(8), double(ref_lvl));
		divisor += value_factor;

		const std::array<double, 3> cell_center = grid.geometry.get_center(cell);
		if (grid.length.get()[0] > 1) {
			avg_analytic += value_factor * get_solution_value(cell_center[0]);
		} else if (grid.length.get()[1] > 1) {
			avg_analytic += value_factor * get_solution_value(cell_center[1]);
		} else if (grid.length.get()[2] > 1) {
			avg_analytic += value_factor * get_solution_value(cell_center[2]);
		}

		Poisson_Cell* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< std::endl;
			abort();
		}

		avg_solved += value_factor * cell_data->solution;
	}
	avg_analytic /= divisor;
	avg_solved /= divisor;

	for (const auto& cell: cells) {
		Poisson_Cell* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< std::endl;
			abort();
		}

		cell_data->solution -= avg_solved - avg_analytic;
	}
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

	for (const auto& cell: cells) {
		const std::array<double, 3> cell_center = grid.geometry.get_center(cell);
		double coord = -1;
		if (grid.length.get()[0] > 1) {
			coord = cell_center[0];
		}
		if (grid.length.get()[1] > 1) {
			coord = cell_center[1];
		}
		if (grid.length.get()[2] > 1) {
			coord = cell_center[2];
		}

		const Poisson_Cell* const data = grid[cell];

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
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

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
		number_of_cells = 32;
		number_of_cells <= max_number_of_cells;
		number_of_cells *= 2
	) {
		const double cell_length = 2 * M_PI / number_of_cells;

		Poisson_Solve solver;
		dccrg::Dccrg<
			Poisson_Cell,
			dccrg::Cartesian_Geometry
		> grid_x, grid_y, grid_z, grid_reference;

		const std::array<uint64_t, 3>
			grid_length_x = {{number_of_cells, 1, 1}},
			grid_length_y = {{1, number_of_cells, 1}},
			grid_length_z = {{1, 1, number_of_cells}};

		grid_x.initialize(grid_length_x, comm, "RCB", 0, -1, true, true, true);
		grid_y.initialize(grid_length_y, comm, "RCB", 0, -1, true, true, true);
		grid_z.initialize(grid_length_z, comm, "RCB", 0, -1, true, true, true);
		grid_reference.initialize(grid_length_x, MPI_COMM_SELF, "RCB", 0, 0, true, true, true);

		dccrg::Cartesian_Geometry::Parameters geom_params;
		geom_params.start[0] =
		geom_params.start[1] =
		geom_params.start[2] = 0;
		geom_params.level_0_cell_length[0] =
		geom_params.level_0_cell_length[1] =
		geom_params.level_0_cell_length[2] = 1;

		geom_params.level_0_cell_length[0] = cell_length;
		grid_x.set_geometry(geom_params);
		grid_reference.set_geometry(geom_params);
		geom_params.level_0_cell_length[0] = 1;

		geom_params.level_0_cell_length[1] = cell_length;
		grid_y.set_geometry(geom_params);
		geom_params.level_0_cell_length[1] = 1;

		geom_params.level_0_cell_length[2] = cell_length;
		grid_z.set_geometry(geom_params);
		geom_params.level_0_cell_length[2] = 1;

		// refine every other cell once
		const std::vector<uint64_t> initial_cells = grid_x.get_cells();
		for (const auto& cell: initial_cells) {
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

		for (const auto& cell: refined_cells_x) {
			const int target_process = cell % comm_size;
			grid_x.pin(cell, target_process);
		}
		for (const auto& cell: refined_cells_y) {
			const int target_process = cell % comm_size;
			grid_y.pin(cell, target_process);
		}
		for (const auto& cell: refined_cells_z) {
			const int target_process = cell % comm_size;
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
		for (const auto& cell: cells_x) {
			Poisson_Cell *data_x = grid_x[cell];

			if (data_x == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_x.geometry.get_center(cell)[0];
			data_x->rhs = get_rhs_value(coord);
			data_x->solution = 0;
		}
		for (const auto& cell: cells_y) {
			Poisson_Cell *data_y = grid_y[cell];

			if (data_y == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_y.geometry.get_center(cell)[1];
			data_y->rhs = get_rhs_value(coord);
			data_y->solution = 0;
		}
		for (const auto& cell: cells_z) {
			Poisson_Cell *data_z = grid_z[cell];

			if (data_z == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_z.geometry.get_center(cell)[2];
			data_z->rhs = get_rhs_value(coord);
			data_z->solution = 0;
		}
		for (const auto& cell: cells_reference) {
			Poisson_Cell *data_reference = grid_reference[cell];

			if (data_reference == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			const double coord = grid_reference.geometry.get_center(cell)[0];
			data_reference->rhs = get_rhs_value(coord);
			data_reference->solution = 0;
		}

		solver.solve(cells_x, grid_x);
		solver.solve(cells_y, grid_y);
		solver.solve(cells_z, grid_z);
		solver.solve(cells_reference, grid_reference);
		normalize_solution(cells_x, grid_x);
		normalize_solution(cells_y, grid_y);
		normalize_solution(cells_z, grid_z);
		normalize_solution(cells_reference, grid_reference);

		// check that parallel solutions are close to analytic
		const double
			p_of_norm = 2,
			norm_x = get_p_norm(cells_x, grid_x, p_of_norm),
			norm_y = get_p_norm(cells_y, grid_y, p_of_norm),
			norm_z = get_p_norm(cells_z, grid_z, p_of_norm),
			norm_reference = get_p_norm(cells_reference, grid_reference, p_of_norm);

		if (norm_x > old_norm) {
			success = 1;
			if (rank == 0) {
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
			if (rank == 0) {
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
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of z solution larger from exact with " << number_of_cells
					<< " cells (" << norm_z << ") than with " << old_number_of_cells
					<< " cells"
					<< std::endl;
			}
		}

		// check that AMR solution is better than reference
		if (norm_x > norm_reference) {
			success = 1;
			if (rank == 0) {
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
			if (rank == 0) {
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
			if (rank == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< " norm of z solution (" << norm_z
					<< ") larger than reference (" << norm_reference
					<< ") with " << number_of_cells
					<< " cells"
					<< std::endl;
			}
		}

		old_norm = norm_x;
		old_number_of_cells = number_of_cells;

		MPI_Allreduce(&success, &global_success, 1, MPI_INT, MPI_SUM, comm);
		if (global_success > 0) {
			break;
		}
	}

	MPI_Finalize();

	if (success == 0) {
		if (rank == 0) {
			cout << "PASSED" << endl;
		}
		return EXIT_SUCCESS;
	} else {
		if (rank == 0) {
			cout << "FAILED" << endl;
		}
		return EXIT_FAILURE;
	}
}

