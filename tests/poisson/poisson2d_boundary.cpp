/*
Version of poisson2d.cpp with boundaries.

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

#include "array"
#include "boost/mpi.hpp"
#include "cmath"
#include "cstdint"
#include "cstdlib"
#include "iostream"
#include "vector"

#include "dccrg.hpp"
#include "dccrg_cartesian_geometry.hpp"

#include "poisson_solve.hpp"


using namespace boost::mpi;
using namespace dccrg;
using namespace std;


int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


double get_solution_value(const double x, const double y)
{
	return x*x * std::exp(-y*y);
}

double get_rhs_value(const double x, const double y)
{
	return 2 * (1 + 2 * x*x * y*y - x*x) * std::exp(-y*y);
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
	for(const auto& cell: cells) {

		const std::array<double, 3> cell_center = grid.geometry.get_center(cell);
		avg_analytic += get_solution_value(cell_center[0], cell_center[1]);

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

	for(const auto& cell: cells) {
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
	double local = 0, global = 0;
	for(const auto& cell: cells) {
		Poisson_Cell* data = grid[cell];
		if (data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for last cell " << cell
				<< std::endl;
			abort();
		}

		double analytic_solution = std::numeric_limits<double>::max();
		const std::array<double, 3> cell_center = grid.geometry.get_center(cell);

		analytic_solution = get_solution_value(cell_center[0], cell_center[1]);

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

	int success = 0; // > 0 means failure

	const size_t max_number_of_cells = 64; // per dimension
	for (size_t nr_of_cells = 8; nr_of_cells <= max_number_of_cells; nr_of_cells *= 2) {

		const double
			cell_length_x = 4.0 / nr_of_cells,
			cell_length_y = 6.0 / nr_of_cells;

		Dccrg<Poisson_Cell, Cartesian_Geometry> grid;

		const std::array<uint64_t, 3> grid_length = {{nr_of_cells + 4, nr_of_cells + 4, 1}};

		grid.initialize(grid_length, comm, "RANDOM", 0, 0, false, false, false);

		Cartesian_Geometry::Parameters geom_params;
		geom_params.start[0] = -2 - 2 * cell_length_x;
		geom_params.start[1] = -3 -2 * cell_length_y;
		geom_params.start[2] = 0;
		geom_params.level_0_cell_length[0] = cell_length_x;
		geom_params.level_0_cell_length[1] = cell_length_y;
		geom_params.level_0_cell_length[2] = 1;
		grid.set_geometry(geom_params);

		grid.balance_load();
		const std::vector<uint64_t> all_cells = grid.get_cells();
		std::vector<uint64_t> solve_cells, boundary_cells, cells_to_skip;

		// classify local cells
		for(const auto& cell: all_cells) {
			dccrg::Types<3>::indices_t cell_indices = grid.mapping.get_indices(cell);

			if (
				cell_indices[0] == 0
				or cell_indices[0] == grid_length[0] - 1
				or cell_indices[1] == 0
				or cell_indices[1] == grid_length[1] - 1
			) {
				cells_to_skip.push_back(cell);
			} else if (
				cell_indices[0] >= 2
				and cell_indices[0] <= grid_length[0] - 3
				and cell_indices[1] >= 2
				and cell_indices[1] <= grid_length[1] - 3
			) {
				solve_cells.push_back(cell);
			} else {
				boundary_cells.push_back(cell);
			}
		}

		// initialize
		for(const auto& cell: solve_cells) {
			Poisson_Cell* const data = grid[cell];

			if (data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			const std::array<double, 3> cell_center = grid.geometry.get_center(cell);

			data->rhs = get_rhs_value(
				cell_center[0],
				cell_center[1]
			);

			data->solution = 0;
		}
		for(const auto& cell: boundary_cells) {
			Poisson_Cell* const data = grid[cell];

			if (data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			const std::array<double, 3> cell_center = grid.geometry.get_center(cell);

			data->rhs = get_rhs_value(
				cell_center[0],
				cell_center[1]
			);

			data->solution = get_solution_value(
				cell_center[0],
				cell_center[1]
			);
		}

		Poisson_Solve solver;
		solver.solve(solve_cells, grid, cells_to_skip);
		normalize_solution(solve_cells, grid);

		const double
			p_of_norm = 2,
			norm  = get_p_norm(solve_cells, grid, p_of_norm);

		if (
			(nr_of_cells == 8 and norm > 0.51)
			or (nr_of_cells == 16 and norm > 0.65)
			or (nr_of_cells == 32 and norm > 0.47)
			or (nr_of_cells == 64 and norm > 0.51)
		) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " " << p_of_norm
					<< "-norm between solved and analytic is too large with "
					<< nr_of_cells << " cells: " << norm
					<< std::endl;
			}
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

