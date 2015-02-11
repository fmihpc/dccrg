/*
Test the reference Poisson solver.

Copyright 2013, 2014, 2015 Finnish Meteorological Institute

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

#include "cmath"
#include "cstdlib"
#include "iostream"
#include "limits"

#include "reference_poisson_solve.hpp"

using namespace std;

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
If p < 0 calcualtes the infinite norm and returns 0 if p == 0.
*/
double get_p_norm(
	const Reference_Poisson_Solve& solver,
	const double p,
	const double offset
) {
	const double dx = solver.get_dx();
	const std::vector<double>& solution = solver.get_complete_solution();

	double norm = 0;
	for (size_t i = 0; i < solution.size(); i++) {
		norm += pow(
			fabs(solution[i] - (get_solution_value((i + 0.5) * dx) + offset)),
			p
		);
	}
	norm = pow(norm, 1.0 / p);

	return norm;
}

int main()
{
	double old_norm = std::numeric_limits<double>::max();
	size_t old_number_of_cells = 1;

	const size_t max_number_of_cells = 32768;
	for (size_t
		number_of_cells = 2;
		number_of_cells <= max_number_of_cells;
		number_of_cells *= 2
	) {
		const double dx = 2 * M_PI / number_of_cells;
		Reference_Poisson_Solve solver(number_of_cells, dx);

		// fill rhs
		for (size_t i = 0; i < number_of_cells; i++) {
			solver.get_rhs(i) = get_rhs_value((i + 0.5) * dx);
		}

		solver.solve();

		// offset of exact solution from solved in last cell
		const double offset
			= solver.get_solution(number_of_cells - 1)
			- get_solution_value((number_of_cells - 0.5) * dx);

		const double norm = get_p_norm(solver, 2, offset);

		if (norm > old_norm) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " 2 norm of solution from exact larger with " << number_of_cells
				<< " cells (" << norm << ") than with " << old_number_of_cells
				<< " cells (" << old_norm << ")"
				<< std::endl;
			cout << "FAILED" << endl;
			return EXIT_FAILURE;
		}
		old_norm = norm;
		old_number_of_cells = number_of_cells;
	}

	cout << "PASSED" << endl;

	return EXIT_SUCCESS;
}

