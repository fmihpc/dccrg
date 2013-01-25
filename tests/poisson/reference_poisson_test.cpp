/*
Test for the reference Poisson solver.
*/

#include "boost/foreach.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"

#include "reference_poisson_solve.hpp"

using namespace std;

double get_rhs_value(const double x)
{
	if (x < 0 || x > 2 * M_PI) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " x must be in the range [0, 2 * pi]" << x
			<< std::endl;
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
	}

	return -sin(x);
}

int main()
{
	const size_t number_of_cells = 10;
	const double dx = 2 * M_PI / number_of_cells;
	Reference_Poisson_Solve solver(number_of_cells);

	// fill rhs
	for (size_t i = 0; i < number_of_cells; i++) {
		solver.get_rhs(i) = get_rhs_value((i + 0.5) * dx);
	}

	solver.solve();

	cout << "RHS:\n";
	BOOST_FOREACH(const double rhs, solver.get_complete_rhs()) {
		cout << rhs << " ";
	}
	cout << endl;

	cout << "Solution:\n";
	BOOST_FOREACH(const double solution, solver.get_complete_solution()) {
		cout << solution << " ";
	}
	cout << endl;

	cout << "Exact scaled solution:\n";
	const double last_exact_value = get_solution_value((number_of_cells - 0.5) * dx);
	for (size_t i = 0; i < number_of_cells; i++) {
		cout << get_solution_value((i + 0.5) * dx) - last_exact_value << " ";
	}
	cout << endl;

	return 0;
}

