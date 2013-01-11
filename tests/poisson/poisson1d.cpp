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

	const size_t number_of_cells = 5;

	/*
	Reference solution
	*/
	Reference_Poisson_Solve reference_solver(number_of_cells);
	//reference_solver.get_rhs(number_of_cells / 2 - 1) = 1;
	reference_solver.get_rhs(number_of_cells / 2) = 1;

	reference_solver.solve();

	/*
	Parallel solution
	*/
	dccrg::Dccrg<Poisson_Cell> grid;
	grid.set_geometry(number_of_cells,  1,  1, 0, 0, 0, 1, 1, 1);
	grid.initialize(comm, "RANDOM", 0, 0, false, false, false);
	grid.balance_load();
	grid.update_remote_neighbor_data();

	std::vector<uint64_t> cells = grid.get_cells();
	std::sort(cells.begin(), cells.end());

	// initialize
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* data = grid[cell];
		if (data == NULL) {
			abort();
		}
		data->rhs = reference_solver.get_rhs(cell - 1);
		data->solution = 0;
	}

	Poisson_Solve(cells, grid);

	cout << "\nReference rhs:\n";
	BOOST_FOREACH(const double rhs, reference_solver.get_complete_rhs()) {
		cout << rhs << " ";
	}

	cout << "\nrhs:\n";
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* data = grid[cell];
		cout << data->rhs << " ";
	}

	cout << "\n\nReference solution:\n";
	BOOST_FOREACH(const double solution, reference_solver.get_complete_solution()) {
		cout << solution << " ";
	}

	/*
	Scale parallel solution so that it is 0 in the last cell
	*/

	// which process has the last cell
	int proc_with_last = 0, proc_local = 0;
	if (grid.is_local(number_of_cells)) {
		proc_local = grid.get_rank();
	}
	MPI_Allreduce(&proc_local, &proc_with_last, 1, MPI_INT, MPI_SUM, grid.get_comm());

	// scatter solution from last cell
	double solution_in_last = 0, solution_local = 0;
	if (grid.is_local(number_of_cells)) {
		solution_local = grid[number_of_cells]->solution;
	}
	MPI_Scatter(
		&solution_local, 1, MPI_DOUBLE,
		&solution_in_last, 1, MPI_DOUBLE,
		proc_with_last,
		grid.get_comm()
	);

	// scale solutions
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* data = grid[cell];
		data->solution -= solution_in_last;
	}

	cout << "\nSolution:\n";
	BOOST_FOREACH(const uint64_t cell, cells) {
		Poisson_Cell* data = grid[cell];
		cout << data->solution / data->scaling_factor << " ";
	}
	cout << endl;

	return EXIT_SUCCESS;
}

