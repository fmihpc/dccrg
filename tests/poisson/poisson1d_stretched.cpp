/*
A stretched grid version of poisson1d test.

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

#include "boost/array.hpp"
#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/static_assert.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "stdint.h"
#include "vector"

#include "dccrg_stretched_cartesian_geometry.hpp"
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
	const double p_of_norm,
	const double offset
) {
	double local = 0, global = 0;

	BOOST_FOREACH(const uint64_t cell, cells) {
		// assumes grid is 1d
		const boost::array<double, 3> cell_center = grid.geometry.get_center(cell);
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

		Poisson_Cell* data = grid[cell];
		local += std::pow(
			fabs(data->solution - (get_solution_value(coord) + offset)),
			p_of_norm
		);
	}
	MPI_Comm temp = grid.get_communicator();
	MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, temp);
	MPI_Comm_free(&temp);
	global = std::pow(global, 1.0 / p_of_norm);

	return global;
}

/*!
Returns the solution in the last cell.
*/
template<class Geometry> double get_last_cell_solution(
	const uint64_t last_cell,
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

	return solution;
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

	double old_norm_reference = std::numeric_limits<double>::max();
	size_t old_number_of_cells = 1;

	const size_t max_number_of_cells = 32768;
	for (size_t
		number_of_cells = 2;
		number_of_cells <= max_number_of_cells;
		number_of_cells *= 2
	) {
		Dccrg<Poisson_Cell, Stretched_Cartesian_Geometry> grid_stretched;
		Dccrg<Poisson_Cell, Cartesian_Geometry> grid_reference;

		const boost::array<uint64_t, 3> grid_length = {{number_of_cells, 1, 1}};

		grid_stretched.initialize(grid_length, comm, "RCB", 0, 0, true, true, true);
		grid_reference.initialize(grid_length, comm, "RCB", 0, 0, true, true, true);

		Cartesian_Geometry::Parameters geom_params;
		geom_params.start[0] =
		geom_params.start[1] =
		geom_params.start[2] = 0;
		geom_params.level_0_cell_length[0] = 2 * M_PI / number_of_cells;
		geom_params.level_0_cell_length[1] = 1;
		geom_params.level_0_cell_length[2] = 1;
		grid_reference.set_geometry(geom_params);

		Stretched_Cartesian_Geometry::Parameters stretched_geom_params;
		boost::array<vector<double>, 3>& coordinates = stretched_geom_params.coordinates;
		coordinates[0].push_back(0);
		coordinates[0].push_back(2 * M_PI);
		coordinates[1].push_back(0);
		coordinates[1].push_back(1);
		coordinates[2].push_back(0);
		coordinates[2].push_back(1);

		/*
		Divide cells using given ratio with smaller cells
		towared higher gradient of solution.
		*/
		const double ratio = 3;

		// divide equally into first four cells
		if (number_of_cells >= 2) {
			coordinates[0].insert(coordinates[0].begin() + 1, M_PI);
		}
		if (number_of_cells >= 4) {
			coordinates[0].insert(coordinates[0].begin() + 1, M_PI / 2);
			coordinates[0].insert(coordinates[0].begin() + 3, 3 * M_PI / 2);
		}

		while (coordinates[0].size() < number_of_cells) {

			const size_t new_size = 2 * (coordinates[0].size() - 1) + 1;
			for (size_t i = 1; i < new_size; i += 2) {

				double new_coord = -1;
				const double middle = (coordinates[0][i] + coordinates[0][i - 1]) / 2;

				// smaller cell towards 0
				if (middle < M_PI / 2
				|| (middle >= M_PI && middle < 3 * M_PI / 2)) {

					new_coord
						= coordinates[0][i]
						- (coordinates[0][i] - coordinates[0][i - 1]) * ratio / (ratio + 1);

				// smaller cell towards 2 pi
				} else {

					new_coord
						= coordinates[0][i - 1]
						+ (coordinates[0][i] - coordinates[0][i - 1]) * ratio / (ratio + 1);

				}
				coordinates[0].insert(coordinates[0].begin() + i, new_coord);
			}
		}

		grid_stretched.set_geometry(stretched_geom_params);

		Poisson_Solve solver(10, 1e-5, 2, 10);

		const std::vector<uint64_t> initial_cells = grid_stretched.get_cells();

		// emulate RANDOM load balance but make local cells identical in grid_x, y and z
		BOOST_FOREACH(const uint64_t cell, initial_cells) {
			const int target_process = cell % comm.size();
			grid_stretched.pin(cell, target_process);
			grid_reference.pin(cell, target_process);
		}
		grid_stretched.balance_load(false);
		grid_reference.balance_load(false);
		grid_stretched.unpin_all_cells();
		grid_reference.unpin_all_cells();

		const std::vector<uint64_t> cells = grid_stretched.get_cells();

		// initialize
		BOOST_FOREACH(const uint64_t cell, cells) {
			Poisson_Cell
				*data_stretched = grid_stretched[cell],
				*data_reference = grid_reference[cell];

			if (data_stretched == NULL || data_reference == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			const double
				stretched_coord = grid_stretched.geometry.get_center(cell)[0],
				reference_coord = grid_reference.geometry.get_center(cell)[0];

			data_stretched->rhs = get_rhs_value(stretched_coord);
			data_reference->rhs = get_rhs_value(reference_coord);

			data_stretched->solution =
			data_reference->solution = 0;
		}

		solver.solve(cells, grid_stretched);
		solver.solve(cells, grid_reference);

		// check that parallel solutions are close to analytic
		const double
			p_of_norm = 2,
			norm_stretched = get_p_norm(cells, grid_stretched, p_of_norm, 0),
			norm_reference = get_p_norm(cells, grid_reference, p_of_norm, 0);

		if (number_of_cells == 8
		&& norm_stretched > norm_reference) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " 2 norm of stretched solution from exact (" << norm_stretched
					<< ") larger than 2 norm of reference from exact (" << norm_reference << ")"
					<< std::endl;
			}
		}

		// the stretched version isn't supposed to work better with more cells,
		// at least not with a ratio of 3 in length between neighboring cells
		if (norm_reference > old_norm_reference) {
			success = 1;
			if (comm.rank() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " 2 norm of reference solution from exact larger with " << number_of_cells
					<< " cells (" << norm_reference << ") than with " << old_number_of_cells
					<< " cells (" << old_norm_reference << ")"
					<< std::endl;
			}
		}

		old_norm_reference = norm_reference;
		old_number_of_cells = number_of_cells;

		MPI_Allreduce(&success, &global_success, 1, MPI_INT, MPI_SUM, comm);
		if (global_success > 0) {
			break;
		}

		// write solutions into a gnuplottable file
		if (grid_stretched.get_comm_size() == 1
		&& number_of_cells <= 32) {

			const std::string name_prefix(
				"poisson1d_stretched_"
				+ boost::lexical_cast<std::string>(number_of_cells)
			);

			std::ofstream plot_file(name_prefix + ".gnuplot");

			plot_file <<
				"set term svg\n"
				"set output '" + name_prefix + ".svg\n"
				"set title 'Solution with "
				+ boost::lexical_cast<std::string>(number_of_cells) + " cells\n"
				"set xrange [0 : 2 * pi]\n"
				"set key outside horizontal bottom\n"
				"plot -sin(x) lw 2 lt -1 t 'analytic', "
				"'-' using 1:2 pt 5 lt 1 t 'stretched', "
				"'-' using 1:2 pt 7 lt 3 t 'reference'\n";

			// stretched
			BOOST_FOREACH(const uint64_t cell, cells) {
				Poisson_Cell* data_stretched = grid_stretched[cell];

				if (data_stretched == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for cell " << cell
						<< std::endl;
					abort();
				}

				const double stretched_coord = grid_stretched.geometry.get_center(cell)[0];

				plot_file << stretched_coord << " " << data_stretched->solution << "\n";
			}
			plot_file << "end\n";

			// reference
			BOOST_FOREACH(const uint64_t cell, cells) {
				Poisson_Cell* data_reference = grid_reference[cell];

				if (data_reference == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for cell " << cell
						<< std::endl;
					abort();
				}

				const double reference_coord = grid_reference.geometry.get_center(cell)[0];

				plot_file << reference_coord << " " << data_reference->solution << "\n";
			}
			plot_file << "end" << endl;
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

