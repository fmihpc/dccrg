/*
Serial solver of Poisson's equation in one dimension.

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

#ifndef REFERENCE_POISSON_SOLVE_HPP
#define REFERENCE_POISSON_SOLVE_HPP

#include "boost/foreach.hpp"
#include "cstdlib"
#include "iostream"
#include "vector"

/*!
Solves the equation nabla^2 f = a in one dimension.

f is referred to as solution and a as rhs.
Assumes that total value of rhs == 0 and f == 0 in the last cell.

Cells in the grid are just indices into std::vector so start from 0, etc.

Implements the algorithm from Computer Simulation Using Particles
by Roger W. Hockney, J. James W. Eastwood.
*/
class Reference_Poisson_Solve
{

public:

	Reference_Poisson_Solve()
	{
		this->initialized = false;
	}

	Reference_Poisson_Solve(const size_t number_of_cells)
	{
		this->initialize(number_of_cells);
	}

	~Reference_Poisson_Solve()
	{
		this->initialized = false;
	}


	/*!
	Initializes the solver to a grid with given number of cells.

	Clears previous solution and rhs.
	*/
	void initialize(const size_t number_of_cells)
	{
		this->solution.clear();
		this->solution.resize(number_of_cells);
		this->rhs.clear();
		this->rhs.resize(number_of_cells);
		this->initialized = true;
	}

	/*!
	Returns the rhs in given cell.
	*/
	double& get_rhs(const size_t index)
	{
		return this->rhs.at(index);
	}

	/*!
	const version of get_rhs().
	*/
	const double& get_rhs(const size_t index) const
	{
		return this->rhs.at(index);
	}

	/*!
	Returns a reference to rhs in all cells.
	*/
	const std::vector<double>& get_complete_rhs() const
	{
		return this->rhs;
	}


	/*!
	Returns solution in given cell.
	*/
	const double& get_solution(const size_t index) const
	{
		return this->solution.at(index);
	}

	/*!
	Returns a reference to the solution in all cells.
	*/
	const std::vector<double>& get_complete_solution() const
	{
		return this->solution;
	}


	/*!
	Solves the equation.

	rhs must have been set before calling this function.
	rhs is scaled in this function so that total rhs == 0.
	*/
	void solve()
	{
		if (!this->initialized) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " solve() called before initialize()"
				<< std::endl;
			abort();
		}

		if (this->rhs.size() == 0) {
			return;
		}

		// make total rhs == 0
		double avg_charge = 0;
		for (size_t i = 0; i < this->rhs.size(); i++) {
			avg_charge += this->rhs[i];
		}

		avg_charge /= this->rhs.size();
		for (size_t i = 0; i < this->rhs.size(); i++) {
			this->rhs[i] -= avg_charge;
		}

		// solution in last cell is defined thus
		this->solution[this->solution.size() - 1] = 0;

		if (this->rhs.size() == 1) {
			return;
		}

		/*
		The algorithm in the book assumes rhs is charge density
		which is negated in the Poisson's equation for electric
		potential. To solve for general rhs negate the rhs in the
		book.
		*/
		this->solution[0] = 0;
		for (size_t i = 0; i < this->rhs.size(); i++) {
			this->solution[0] -= (i + 1) * this->rhs[i];
		}
		this->solution[0] /= this->rhs.size();

		this->solution[1] = -this->rhs[0] + 2 * this->solution[0];

		for (size_t i = 2; i < this->rhs.size(); i++) {
			this->solution[i]
				= -this->rhs[i - 1]
				+ 2 * this->solution[i - 1]
				- this->solution[i - 2];
		}
	}


private:

	bool initialized;

	std::vector<double> solution, rhs;
};

#endif

