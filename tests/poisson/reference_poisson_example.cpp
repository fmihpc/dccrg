/*
Example of using the reference Poisson solver.
*/

#include "cstdlib"
#include "iostream"

#include "reference_poisson_solve.hpp"

using namespace std;

int main()
{
	const size_t number_of_cells = 10;
	Reference_Poisson_Solve solver(number_of_cells, 1);

	solver.get_rhs(number_of_cells / 2 - 1) = 1;
	solver.get_rhs(number_of_cells / 2) =  -1;

	solver.solve();

	cout << "Charge:\n";
	for (const auto& charge: solver.get_complete_rhs()) {
		cout << charge << " ";
	}
	cout << endl;

	cout << "Potential:\n";
	for (const auto& potential: solver.get_complete_solution()) {
		cout << potential << " ";
	}
	cout << endl;

	return 0;
}

