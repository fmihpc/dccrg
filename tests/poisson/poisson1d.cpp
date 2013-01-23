/*
Tests the parallel Poisson solver implemented on top of dccrg.

Copyright 2012, 2013 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
#include "reference_poisson_solve.hpp"

using namespace boost::mpi;
using namespace dccrg;
using namespace std;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cerr << "Zoltan_Initialize failed" << endl;
	    return EXIT_FAILURE;
	}

	const uint64_t number_of_cells = 10;

	/*
	Reference solution
	*/
	Reference_Poisson_Solve reference_solver(number_of_cells);
	reference_solver.get_rhs(number_of_cells / 2) = 1;

	reference_solver.solve();

	/*
	Parallel 1d solution in each dimension
	*/
	dccrg::Dccrg<Poisson_Cell> grid_x, grid_y, grid_z;

	grid_x.set_geometry(number_of_cells, 1, 1, 0, 0, 0, 1, 1, 1);
	grid_y.set_geometry(1, number_of_cells, 1, 0, 0, 0, 1, 1, 1);
	grid_z.set_geometry(1, 1, number_of_cells, 0, 0, 0, 1, 1, 1);

	grid_x.initialize(comm, "RANDOM", 0, 0, true, false, false);
	grid_y.initialize(comm, "RANDOM", 0, 0, false, true, false);
	grid_z.initialize(comm, "RANDOM", 0, 0, false, false, true);

	grid_x.balance_load();
	grid_y.balance_load();
	grid_z.balance_load();

	grid_x.update_remote_neighbor_data();
	grid_y.update_remote_neighbor_data();
	grid_z.update_remote_neighbor_data();

	// initialize
	for (uint64_t cell = 1; cell <= number_of_cells; cell++) { // assume no AMR, etc.
		Poisson_Cell
			*data_x = grid_x[cell],
			*data_y = grid_y[cell],
			*data_z = grid_z[cell];

		if (data_x == NULL || data_y == NULL || data_z == NULL) {
			abort();
		}

		data_x->rhs = data_y->rhs = data_z->rhs = -reference_solver.get_rhs(cell - 1);
		data_x->solution = data_y->solution = data_z->solution = 0;
	}

	std::vector<uint64_t> cells;
	for (uint64_t cell = 1; cell <= number_of_cells; cell++) {
		cells.push_back(cell);
	}
	Poisson_Solve(cells, grid_x);
	Poisson_Solve(cells, grid_y);
	Poisson_Solve(cells, grid_z);

	cout << "\nReference rhs:\n";
	BOOST_FOREACH(const double rhs, reference_solver.get_complete_rhs()) {
		cout << rhs << " ";
	}

	cout << "\nrhs:\n";
	for (uint64_t cell = 1; cell <= number_of_cells; cell++) {
		Poisson_Cell
			*data_x = grid_x[cell],
			*data_y = grid_y[cell],
			*data_z = grid_z[cell];
		cout << -data_x->rhs << " ";
		if (data_x->rhs != data_y->rhs
		|| data_x->rhs != data_z->rhs) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " rhs not identical between dimensions :"
				<< data_x->rhs << " " << data_y->rhs << " " << data_z->rhs
				<< std::endl;
			abort();
		}
	}

	cout << "\n\nReference solution:\n";
	BOOST_FOREACH(const double solution, reference_solver.get_complete_solution()) {
		cout << solution << " ";
	}

	/*
	Scale parallel solution so that it is 0 in the last cell
	*/

	// get process that has the last cell
	int
		proc_with_last_x = 0,
		proc_with_last_y = 0,
		proc_with_last_z = 0,
		proc_local = 0;

	// x dimension
	if (grid_x.is_local(number_of_cells)) {
		proc_local = grid_x.get_rank();
	} else {
		proc_local = 0;
	}
	MPI_Allreduce(&proc_local, &proc_with_last_x, 1, MPI_INT, MPI_SUM, grid_x.get_comm());

	// y dimension
	if (grid_y.is_local(number_of_cells)) {
		proc_local = grid_y.get_rank();
	} else {
		proc_local = 0;
	}
	MPI_Allreduce(&proc_local, &proc_with_last_y, 1, MPI_INT, MPI_SUM, grid_y.get_comm());

	// z dimension
	if (grid_z.is_local(number_of_cells)) {
		proc_local = grid_z.get_rank();
	} else {
		proc_local = 0;
	}
	MPI_Allreduce(&proc_local, &proc_with_last_z, 1, MPI_INT, MPI_SUM, grid_z.get_comm());

	// scatter solution from last cell
	double
		solution_in_last_x = 0,
		solution_in_last_y = 0,
		solution_in_last_z = 0,
		solution_local = 0;

	// x
	if (grid_x.is_local(number_of_cells)) {
		solution_local = grid_x[number_of_cells]->solution;
	} else {
		solution_local = 0;
	}
	MPI_Scatter(
		&solution_local, 1, MPI_DOUBLE,
		&solution_in_last_x, 1, MPI_DOUBLE,
		proc_with_last_x,
		grid_x.get_comm()
	);

	// y
	if (grid_y.is_local(number_of_cells)) {
		solution_local = grid_y[number_of_cells]->solution;
	} else {
		solution_local = 0;
	}
	MPI_Scatter(
		&solution_local, 1, MPI_DOUBLE,
		&solution_in_last_y, 1, MPI_DOUBLE,
		proc_with_last_y,
		grid_y.get_comm()
	);

	// z
	if (grid_z.is_local(number_of_cells)) {
		solution_local = grid_z[number_of_cells]->solution;
	} else {
		solution_local = 0;
	}
	MPI_Scatter(
		&solution_local, 1, MPI_DOUBLE,
		&solution_in_last_z, 1, MPI_DOUBLE,
		proc_with_last_z,
		grid_z.get_comm()
	);

	// scale solutions
	for (uint64_t cell = 1; cell <= number_of_cells; cell++) {
		Poisson_Cell
			*data_x = grid_x[cell],
			*data_y = grid_y[cell],
			*data_z = grid_z[cell];

		data_x->solution -= solution_in_last_x;
		data_y->solution -= solution_in_last_y;
		data_z->solution -= solution_in_last_z;
	}

	cout << "\nSolution:\n";
	for (uint64_t cell = 1; cell <= number_of_cells; cell++) {

		Poisson_Cell
			*data_x = grid_x[cell],
			*data_y = grid_y[cell],
			*data_z = grid_z[cell];

		if (data_x->solution != data_y->solution
		|| data_x->solution != data_z->solution) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " solution not identical between dimensions :"
				<< data_x->solution << " " << data_y->solution << " " << data_z->solution
				<< std::endl;
			abort();
		}

		const double
			solution = data_x->solution,
			reference = reference_solver.get_solution(cell - 1);

		if (fabs(solution) > 1e-10
		&& fabs(solution - reference) / solution > 1e-11) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " solution too different from reference ("
				<< "solution: " << solution
				<< ", reference " << reference << ")"
				<< std::endl;
			abort();
		}

		cout << data_x->solution << " ";
	}
	cout << endl;

	return EXIT_SUCCESS;
}

