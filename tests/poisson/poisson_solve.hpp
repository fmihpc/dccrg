/*
Parallel Poisson solver build on top of dccrg.

Copyright 2012, 2013, 2014 Finnish Meteorological Institute

Dccrg free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DCCRG_POISSON_SOLVE_HPP
#define DCCRG_POISSON_SOLVE_HPP

#include "boost/foreach.hpp"
#include "boost/format.hpp"
#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/static_assert.hpp"
#include "boost/tuple/tuple.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "mpi.h"
#include "stdint.h"
#include "vector"

#include "dccrg.hpp"


/*
Remember to allocate space for the static variable of
Poisson_Cell by placing the following inside main(...),
for example:
int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;
*/

/*!
Dccrg cell data class used by the parallel Poisson solver class.
*/
class Poisson_Cell
{
public:

	Poisson_Cell()
	{
		this->rhs            =
		this->solution       =
		this->best_solution  =
		this->p0             =
		this->p1             =
		this->r0             =
		this->r1             =
		this->A_dot_p0       =
		this->scaling_factor =
		this->f_x_pos        =
		this->f_x_neg        =
		this->f_y_pos        =
		this->f_y_neg        =
		this->f_z_pos        =
		this->f_z_neg        = 0;
	}

	double
		// rhs is b and solution is x in the equation A . x = b
		rhs, solution,

		// saved solution whenever minimum residual obtained
		best_solution,

		/*
		Same variables are used as in 2.7.6 of numerical recipes
		with the modification that bold r is r0 and bold r with
		dash on top is r1; same for p.
		*/
		p0, p1, r0, r1,

		/*
		The recurrence in numerical recipes does A . p0 two times,
		this variable caches the result from first to save CPU time.
		*/
		A_dot_p0,

		/*
		Factor by which solution variables in this cell are scaled
		so that diagonal of matrix A has only ones
		*/
		scaling_factor,

		/*
		Factors that data in neighboring cells is multiplied with
		when calculating A . i for each local cell where i = p0, p1.
		The neighbors' factors are used when calculating traspose(A) . i.
		One factor per direction.
		*/
		f_x_pos, f_x_neg, f_y_pos, f_y_neg, f_z_pos, f_z_neg;


	/*
	Decides which data to transfer over MPI.
	*/
	static int transfer_switch;
	static const int
		// data related to solving the equation
		SOLVING   = 0,
		// geometry factors
		GEOMETRY  = 1,
		// solution when initializing
		INIT      = 2;

	// tells dccrg what to transfer, assumes no padding between variables
	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype()
	{
		void* address = NULL;
		int count = -1;
		MPI_Datatype datatype = MPI_DATATYPE_NULL;

		switch(Poisson_Cell::transfer_switch) {
		case Poisson_Cell::SOLVING:
			address = &(this->p0);
			count = 2;
			datatype = MPI_DOUBLE;
			break;
		case Poisson_Cell::GEOMETRY:
			address = &(this->scaling_factor);
			count = 7;
			datatype = MPI_DOUBLE;
			break;
		case Poisson_Cell::INIT:
			address = &(this->solution);
			count = 1;
			datatype = MPI_DOUBLE;
			break;
		default:
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid transfer switch: " << Poisson_Cell::transfer_switch
				<< std::endl;
			abort();
			break;
		}

		return boost::make_tuple(address, count, datatype);
	}
};


/*!
Parallel solver for the Poisson equation built on top of dccrg.
*/
class Poisson_Solve
{

public:

	/*!
	Creates a solver with default parameters.

	Use the other constructor to customize parameters.
	*/
	Poisson_Solve()
	{
		this->max_iterations = 1000;
		this->min_iterations = 0;
		this->stop_residual = 1e-15;
		this->p_of_norm = 2;
		this->stop_after_residual_increase = 10;
		this->verbose = false;
	};

	/*!
	Creates a solver with given parameters.

	When solving no more than given_max_iterations iterations will be done.
	Solving will stop if:
		- more than given_max_iterations iterations have been calculated
		- p-norm of the solution is below given_stop_residual
		- residual has increased by a factor stop_after_residual_increase
		  from its minimum encountered so far
	but not before given_min_iterations number if iterations have been done.
	*/
	Poisson_Solve(
		const unsigned int given_max_iterations,
		const unsigned int given_min_iterations,
		const double given_stop_residual,
		const double given_p_of_norm,
		const double given_stop_after_residual_increase,
		const bool given_verbosity
	) {
		this->max_iterations = given_max_iterations;
		this->min_iterations = given_min_iterations;
		this->stop_residual = given_stop_residual;
		this->p_of_norm = given_p_of_norm;
		this->stop_after_residual_increase = given_stop_after_residual_increase;
		this->verbose = given_verbosity;
	};


	void set_verbosity(const bool given)
	{
		this->verbose = given;
	}


	/*!
	Solves the Poisson's equation in given cells.

	The right hand side (rhs) of the equation must have been initialized
	before calling this function.
	The solution (x) of the equation also must have been initialized
	(with a guess, 0 usually works) before calling this function.

	Given grid must have a default neighborhood of size 0? or ...

	Cells in cells_to_skip and missing neighbors in case the grid is not
	periodic are assumed to be of equal size and have rhs of equal value
	to their current neighbor. If any cell is both in given cells and
	cells_to_skip the result is undefined.

	If the structure of the grid has not changed since the last call
	to this function (new cells haven't been created or existing ones removed,
	cells haven't changed processes, cells_to_skip hasn't changed) then
	cache_is_up_to_date can be set to true in which case the
	caching/preparation step will be skipped speeding up this function.

	If the residual increases to 10x its minimum value while solving
	the last solution will be replaced with the one for which
	the residual was smallest and this function will return.
	*/
	template<class Geometry> void solve(
		/* TODO: overlap computation with communication
		const std::vector<uint64_t>& inner_cells,
		const std::vector<uint64_t>& outer_cells,*/
		const std::vector<uint64_t>& cells,
		dccrg::Dccrg<Poisson_Cell, Geometry>& grid,
		const boost::unordered_set<uint64_t>& cells_to_skip = boost::unordered_set<uint64_t>(),
		const bool cache_is_up_to_date = false
	) {
		// TODO: if a neighbor is not in given cells assume it doesn't exist
		// const boost::unordered_set<uint64_t> unique_cells;
		// unique_cells.reserve(inner_cells.size() + outer_cells.size());
		// unique_cells.insert(inner_cells.begin(), inner_cells.end());
		// unique_cells.insert(outer_cells.begin(), outer_cells.end());

		this->comm = grid.get_communicator();
		this->comm_rank = grid.get_rank();

		if (!cache_is_up_to_date) {
			this->cache_system_info(cells, cells_to_skip, grid);
		}

		double
			// minimum residual reached while solving
			residual_min = std::numeric_limits<double>::max(),
			// local and global values of r0 . r1
			dot_r_l = 0,
			dot_r_g = this->initialize_solver(grid);

		if (this->comm_rank == 0 && this->verbose) {
			std::cout << "r0 . r1: " << dot_r_g << std::endl;
		}

		size_t iteration = 0;
		do {
			iteration++;

			// TODO: update only p0 for calculating alpha and only then p1?
			Poisson_Cell::transfer_switch = Poisson_Cell::SOLVING;
			grid.update_copies_of_remote_neighbors();

			/*
			Calculate alpha
			*/

			// A . p0, cache the result
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;

				data->A_dot_p0 = data->scaling_factor * data->p0;

				BOOST_FOREACH(const neighbor_info_t neigh_info, info.second) {
					Poisson_Cell* const neighbor_data = neigh_info.get<0>();

					double multiplier = 0;
					const int direction = neigh_info.get<1>();
					switch(direction) {
					case +1:
						multiplier = data->f_x_pos;
						break;
					case -1:
						multiplier = data->f_x_neg;
						break;
					case +2:
						multiplier = data->f_y_pos;
						break;
					case -2:
						multiplier = data->f_y_neg;
						break;
					case +3:
						multiplier = data->f_z_pos;
						break;
					case -3:
						multiplier = data->f_z_neg;
						break;
					default:
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Invalid direction: " << direction
							<< std::endl;
						abort();
						break;
					}

					const int rel_ref_lvl = neigh_info.get<2>();
					if (rel_ref_lvl > 0) {
						multiplier /= 4.0;
					}

					data->A_dot_p0 += multiplier * neighbor_data->p0;
				}
			}

			// p1 . (A . p0)
			double dot_p_l = 0, dot_p_g = 0;
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;
				dot_p_l += data->p1 * data->A_dot_p0;
			}
			MPI_Allreduce(&dot_p_l, &dot_p_g, 1, MPI_DOUBLE, MPI_SUM, this->comm);
			if (this->comm_rank == 0 && this->verbose) {
				std::cout << "p1 . (A . p0): " << dot_p_g << std::endl;
			}
			// no sense in continuing with dividing by zero
			if (dot_p_g == 0) {
				break;
			}

			const double alpha = dot_r_g / dot_p_g;
			if (this->comm_rank == 0 && this->verbose) {
				std::cout << "alpha: " << alpha << std::endl;
			}


			// update solution
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;
				data->solution += alpha * data->p0;
			}

			// update residual and possibly stop solving
			const double residual = this->get_residual();
			if (this->comm_rank == 0 && this->verbose) {
				std::cout << "residual: " << residual << std::endl;
			}

			// save solution if at minimum residual so far
			if (residual_min > residual) {
				if (this->comm_rank == 0 && this->verbose) {
					std::cout << "saving solution at residual " << residual << std::endl;
				}
				residual_min = residual;
				BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
					Poisson_Cell* const data = info.first;
					data->best_solution = data->solution;
				}
			}

			if (
				residual <= this->stop_residual
				&& iteration >= this->min_iterations
			) {
				break;
			}
			if (
				residual >= this->stop_after_residual_increase * residual_min
				&& iteration >= this->min_iterations
			) {
				break;
			}

			// update r0
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;
				data->r0 -= alpha * data->A_dot_p0;
			}

			// update r1
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;

				// A . p1
				double A_dot_p1 = data->scaling_factor * data->p1;

				BOOST_FOREACH(const neighbor_info_t neigh_info, info.second) {
					Poisson_Cell* const neighbor_data = neigh_info.get<0>();

					/*
					In transpose(A) . p1 use that multiplier which was calculated
					for updating the current neighbor's value with this cell's,
					i.e. just "reverse" the direction of the update.
					*/
					double multiplier = 0;
					const int direction = neigh_info.get<1>();
					switch(direction) {
					case +1:
						multiplier = neighbor_data->f_x_neg;
						break;
					case -1:
						multiplier = neighbor_data->f_x_pos;
						break;
					case +2:
						multiplier = neighbor_data->f_y_neg;
						break;
					case -2:
						multiplier = neighbor_data->f_y_pos;
						break;
					case +3:
						multiplier = neighbor_data->f_z_neg;
						break;
					case -3:
						multiplier = neighbor_data->f_z_pos;
						break;
					default:
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Invalid direction: " << direction
							<< std::endl;
						abort();
						break;
					}

					const int rel_ref_lvl = neigh_info.get<2>();
					if (rel_ref_lvl > 0) {
						multiplier /= 4.0;
					}

					A_dot_p1 += multiplier * neighbor_data->p1;
				}

				data->r1 -= alpha * A_dot_p1;
			}

			// no sense in continuing with dividing by zero
			if (dot_r_g == 0) {
				break;
			}

			// calculate beta
			const double old_dot_r_g = dot_r_g;
			dot_r_l = dot_r_g = 0;
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;
				dot_r_l += data->r0 * data->r1;
			}
			MPI_Allreduce(&dot_r_l, &dot_r_g, 1, MPI_DOUBLE, MPI_SUM, this->comm);
			if (this->comm_rank == 0 && this->verbose) {
				std::cout << "new r0 . r1: " << dot_r_g << std::endl;
			}

			const double beta = dot_r_g / old_dot_r_g;
			if (this->comm_rank == 0 && this->verbose) {
				std::cout << "beta: " << beta << std::endl;
			}


			// update p0, p1
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;
				data->p0 = data->r0 + beta * data->p0;
				data->p1 = data->r1 + beta * data->p1;
			}

		} while (iteration < this->max_iterations);

		if (this->comm_rank == 0 && this->verbose) {
			std::cout << "iterations: " << iteration << ", residual: " << residual_min << std::endl;
		}

		BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
			Poisson_Cell* const data = info.first;
			data->solution = data->best_solution;
		}

		MPI_Comm_free(&(this->comm));
	}


	/*!
	A more robust but slower and less accurate version of solve().

	Based on http://www.rsmas.miami.edu/personal/miskandarani/Courses/MSC321/Projects/prjpoisson.pdf
	*/
	template<class Geometry> void solve_failsafe(
		const std::vector<uint64_t>& cells,
		dccrg::Dccrg<Poisson_Cell, Geometry>& grid
	) {
		this->comm = grid.get_communicator();
		this->comm_rank = grid.get_rank();

		this->cache_system_info(cells, boost::unordered_set<uint64_t>(), grid);

		// only solution is needed from other processes
		Poisson_Cell::transfer_switch = Poisson_Cell::INIT;

		size_t iteration = 0;
		double norm = std::numeric_limits<double>::max();
		while (iteration++ < this->max_iterations and norm > this->stop_residual) {

			grid.update_copies_of_remote_neighbors();

			double norm_local = 0;
			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;

				const double inv_scaling_factor = -1.0 / data->scaling_factor;

				// use best solution for storing next solution
				data->best_solution = -inv_scaling_factor * data->rhs;

				BOOST_FOREACH(const neighbor_info_t& neigh_info, info.second) {
					Poisson_Cell* const neighbor_data = neigh_info.get<0>();

					double multiplier = 0;
					const int direction = neigh_info.get<1>();
					switch(direction) {
					case +1:
						multiplier = data->f_x_pos;
						break;
					case -1:
						multiplier = data->f_x_neg;
						break;
					case +2:
						multiplier = data->f_y_pos;
						break;
					case -2:
						multiplier = data->f_y_neg;
						break;
					case +3:
						multiplier = data->f_z_pos;
						break;
					case -3:
						multiplier = data->f_z_neg;
						break;
					default:
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Invalid direction: " << direction
							<< std::endl;
						abort();
						break;
					}

					const int rel_ref_lvl = neigh_info.get<2>();
					if (rel_ref_lvl > 0) {
						multiplier /= 4.0;
					}

					data->best_solution
						+= inv_scaling_factor * multiplier * neighbor_data->solution;
				}

				norm_local += std::fabs(data->solution - data->best_solution);
			}

			norm = 0;
			MPI_Allreduce(&norm_local, &norm, 1, MPI_DOUBLE, MPI_SUM, this->comm);
			if (this->comm_rank == 0 && this->verbose) {
				std::cout << "Norm: " << norm << std::endl;
			}

			BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
				Poisson_Cell* const data = info.first;
				data->solution = data->best_solution;
			}

		}

		if (this->comm_rank == 0 && this->verbose) {
			std::cout << "iterations: " << iteration << ", norm: " << norm << std::endl;
		}

		MPI_Comm_free(&(this->comm));
	}



private:

	/*
	Data of cell's neighbor, neighbor's direction from cell
	and neighbor's relative refinement level (if > 0 neighbor is smaller)
	*/
	typedef typename boost::tuple<Poisson_Cell*, int, int> neighbor_info_t;

	// data of a local cell and info of its neighbors
	typedef typename std::pair<Poisson_Cell*, std::vector<neighbor_info_t> > cell_info_t;


	// maximum and minumum number of iterations to do
	unsigned int
		max_iterations,
		min_iterations;

	// stop solving when residual <= stop_residual
	double stop_residual;

	// p to use when calculating the residual as a p-norm
	double p_of_norm;

	// stop solving when residual has increased by this
	// factor from its minimum value encountered so far
	double stop_after_residual_increase;

	// cached and grid agnostic form of the system to solve
	std::vector<cell_info_t> cell_info;

	MPI_Comm comm;
	int comm_rank;

	bool verbose;

	/*!
	Returns global residual of the solution.
	*/
	double get_residual() const
	{
		double local = 0, global = 0;
		BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
			const Poisson_Cell* const data = info.first;
			local += std::pow(std::fabs(data->r0), this->p_of_norm);
		}
		MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, this->comm);
		global = std::pow(global, 1.0 / p_of_norm);
		return global;
	}


	/*!
	Returns the geometrical scaling factor for given cell.

	Neighbors missing from given list are considered to be
	of equal size to given cell.
	*/
	template<class Geometry> void set_scaling_factor(
		const uint64_t cell,
		const std::vector<std::pair<uint64_t, int> >& neighbors,
		const dccrg::Dccrg<Poisson_Cell, Geometry>& grid
	) const {

		Poisson_Cell* const cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< std::endl;
			abort();
		}

		const boost::array<double, 3> cell_length = grid.geometry.get_length(cell);
		const double
			cell_x_half_size = cell_length[0] / 2.0,
			cell_y_half_size = cell_length[1] / 2.0,
			cell_z_half_size = cell_length[2] / 2.0;

		// non-existing or skipped neighbors are equal in size to current cell
		double
			neigh_pos_x_offset = +2 * cell_x_half_size,
			neigh_neg_x_offset = -2 * cell_x_half_size,
			neigh_pos_y_offset = +2 * cell_y_half_size,
			neigh_neg_y_offset = -2 * cell_y_half_size,
			neigh_pos_z_offset = +2 * cell_z_half_size,
			neigh_neg_z_offset = -2 * cell_z_half_size;

		for (size_t i = 0; i < neighbors.size(); i++) {
			const uint64_t neighbor = neighbors[i].first;
			const int direction = neighbors[i].second;
			const boost::array<double, 3> neighbor_length = grid.geometry.get_length(neighbor);
			const double
				neigh_x_half_size = neighbor_length[0] / 2.0,
				neigh_y_half_size = neighbor_length[1] / 2.0,
				neigh_z_half_size = neighbor_length[2] / 2.0;

			// assume rhs and solution are cell-centered
			switch(direction) {
			case +1:
				neigh_pos_x_offset = cell_x_half_size + neigh_x_half_size;
				break;
			case -1:
				neigh_neg_x_offset = -1.0 * (cell_x_half_size + neigh_x_half_size);
				break;
			case +2:
				neigh_pos_y_offset = cell_y_half_size + neigh_y_half_size;
				break;
			case -2:
				neigh_neg_y_offset = -1.0 * (cell_y_half_size + neigh_y_half_size);
				break;
			case +3:
				neigh_pos_z_offset = cell_z_half_size + neigh_z_half_size;
				break;
			case -3:
				neigh_neg_z_offset = -1.0 * (cell_z_half_size + neigh_z_half_size);
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid direction: " << direction
					<< std::endl;
				abort();
				break;
			}
		}

		const double
			total_offset_x = neigh_pos_x_offset - neigh_neg_x_offset,
			total_offset_y = neigh_pos_y_offset - neigh_neg_y_offset,
			total_offset_z = neigh_pos_z_offset - neigh_neg_z_offset;

		// geometry factors are 0 in directions without neighbors
		cell_data->f_x_pos =
		cell_data->f_x_neg =
		cell_data->f_y_pos =
		cell_data->f_y_neg =
		cell_data->f_z_pos =
		cell_data->f_z_neg = 0;

		for (size_t i = 0; i < neighbors.size(); i++) {
			const int direction = neighbors[i].second;

			// don't mind extra work due to 4 smaller face neighbors
			switch(direction) {
			case +1:
				cell_data->f_x_pos = +2.0 / (neigh_pos_x_offset * total_offset_x);
				break;
			case -1:
				cell_data->f_x_neg = -2.0 / (neigh_neg_x_offset * total_offset_x);
				break;
			case +2:
				cell_data->f_y_pos = +2.0 / (neigh_pos_y_offset * total_offset_y);
				break;
			case -2:
				cell_data->f_y_neg = -2.0 / (neigh_neg_y_offset * total_offset_y);
				break;
			case +3:
				cell_data->f_z_pos = +2.0 / (neigh_pos_z_offset * total_offset_z);
				break;
			case -3:
				cell_data->f_z_neg = -2.0 / (neigh_neg_z_offset * total_offset_z);
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid direction: " << direction
					<< std::endl;
				abort();
				break;
			}
		}

		cell_data->scaling_factor
			=
			- cell_data->f_x_pos
			- cell_data->f_x_neg
			- cell_data->f_y_pos
			- cell_data->f_y_neg
			- cell_data->f_z_pos
			- cell_data->f_z_neg;

		return;
	}


	/*!
	Caches a grid agnostic form of the system to solve.

	Prepares geometric factors, scales the initial solution, etc.
	*/
	template<class Geometry> void cache_system_info(
		const std::vector<uint64_t>& cells,
		const boost::unordered_set<uint64_t>& cells_to_skip,
		dccrg::Dccrg<Poisson_Cell, Geometry>& grid
	) {
		/*
		Calculate scaling factors in given cells and
		neighbors of given cells not in cells_to_skip
		*/

		// mark cells_to_skip internally with a scaling factor of 1
		BOOST_FOREACH(const uint64_t cell_to_skip, cells_to_skip) {
			if (grid.is_local(cell_to_skip)) {
				Poisson_Cell* const cell_data = grid[cell_to_skip];
				if (cell_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for cell " << cell_to_skip
						<< std::endl;
					abort();
				}

				cell_data->scaling_factor = 1;
			}
		}
		Poisson_Cell::transfer_switch = Poisson_Cell::GEOMETRY;
		grid.update_copies_of_remote_neighbors();

		boost::unordered_set<uint64_t> scaling_factor_cells;

		// calculate scaling factors in given cells
		BOOST_FOREACH(const uint64_t cell, cells) {

			Poisson_Cell* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			// possibly skip some face neighbors
			const std::vector<std::pair<uint64_t, int> > all_face_neighbors
				= grid.get_face_neighbors_of(cell);

			std::vector<std::pair<uint64_t, int> > face_neighbors;
			for (size_t i = 0; i < all_face_neighbors.size(); i++) {
				const uint64_t neighbor = all_face_neighbors[i].first;
				Poisson_Cell* const neighbor_data = grid[neighbor];
				if (neighbor_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for neighbor " << neighbor
						<< " of cell " << cell
						<< std::endl;
					abort();
				}

				if (neighbor_data->scaling_factor != 1) {
					face_neighbors.push_back(all_face_neighbors[i]);
					if (grid.is_local(neighbor)) {
						scaling_factor_cells.insert(neighbor);
					}
				}
			}

			this->set_scaling_factor(cell, face_neighbors, grid);
		}

		// calculate scaling factor in rest of required cells
		BOOST_FOREACH(const uint64_t cell, cells) {
			scaling_factor_cells.erase(cell);
		}
		BOOST_FOREACH(const uint64_t cell, scaling_factor_cells) {
			const std::vector<std::pair<uint64_t, int> > all_face_neighbors
				= grid.get_face_neighbors_of(cell);

			std::vector<std::pair<uint64_t, int> > face_neighbors;
			for (size_t i = 0; i < all_face_neighbors.size(); i++) {
				const uint64_t neighbor = all_face_neighbors[i].first;
				Poisson_Cell* const neighbor_data = grid[neighbor];
				if (neighbor_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for neighbor " << neighbor
						<< " of cell " << cell
						<< std::endl;
					abort();
				}

				if (neighbor_data->scaling_factor != 1) {
					face_neighbors.push_back(all_face_neighbors[i]);
					if (grid.is_local(neighbor)) {
						scaling_factor_cells.insert(neighbor);
					}
				}
			}

			this->set_scaling_factor(cell, face_neighbors, grid);
		}

		Poisson_Cell::transfer_switch = Poisson_Cell::GEOMETRY;
		grid.update_copies_of_remote_neighbors();

		// cache information about cells' neighbors
		this->cell_info.clear();
		this->cell_info.reserve(cells.size());
		BOOST_FOREACH(const uint64_t cell, cells) {

			cell_info_t temp_cell_info;

			Poisson_Cell* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			temp_cell_info.first = cell_data;

			const std::vector<std::pair<uint64_t, int> > all_face_neighbors
				= grid.get_face_neighbors_of(cell);

			std::vector<std::pair<uint64_t, int> > face_neighbors;
			for (size_t i = 0; i < all_face_neighbors.size(); i++) {
				const uint64_t neighbor = all_face_neighbors[i].first;
				Poisson_Cell* const neighbor_data = grid[neighbor];
				if (neighbor_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for neighbor " << neighbor
						<< " of cell " << cell
						<< std::endl;
					abort();
				}

				if (neighbor_data->scaling_factor != 1) {
					face_neighbors.push_back(all_face_neighbors[i]);
					if (grid.is_local(neighbor)) {
						scaling_factor_cells.insert(neighbor);
					}
				}
			}

			const int cell_ref_lvl = grid.get_refinement_level(cell);
			for (size_t i = 0; i < face_neighbors.size(); i++) {

				neighbor_info_t temp_neigh_info;

				const int direction = face_neighbors[i].second;
				temp_neigh_info.get<1>() = direction;

				const uint64_t neighbor = face_neighbors[i].first;
				Poisson_Cell* const neighbor_data = grid[neighbor];
				if (neighbor_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for neighbor " << neighbor << " of cell " << cell
						<< std::endl;
					abort();
				}
				temp_neigh_info.get<0>() = neighbor_data;

				int relative_ref_lvl = 0;
				const int neigh_ref_lvl = grid.get_refinement_level(neighbor);
				if (neigh_ref_lvl > cell_ref_lvl) {
					relative_ref_lvl = 1;
				} else if (neigh_ref_lvl == cell_ref_lvl) {
					relative_ref_lvl = 0;
				} else {
					relative_ref_lvl = -1;
				}
				temp_neigh_info.get<2>() = relative_ref_lvl;

				temp_cell_info.second.push_back(temp_neigh_info);
			}

			this->cell_info.push_back(temp_cell_info);
		}
	}


	/*!
	Initializes the solution for solve(...).

	Returns the initial value of r0 . r1.
	*/
	template<class Geometry> double initialize_solver(
		dccrg::Dccrg<Poisson_Cell, Geometry>& grid
	) {
		// transfer user's guess for the solution to calculate residual
		Poisson_Cell::transfer_switch = Poisson_Cell::INIT;
		grid.update_copies_of_remote_neighbors();

		// residual == r0 = rhs - A . solution
		BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
			Poisson_Cell* const data = info.first;

			data->r0 = data->rhs - data->scaling_factor * data->solution;

			BOOST_FOREACH(const neighbor_info_t& neigh_info, info.second) {
				Poisson_Cell* const neighbor_data = neigh_info.get<0>();

				// final multiplier to use for current neighbor's data
				double multiplier = 0;
				const int direction = neigh_info.get<1>();
				switch(direction) {
				case +1:
					multiplier = data->f_x_pos;
					break;
				case -1:
					multiplier = data->f_x_neg;
					break;
				case +2:
					multiplier = data->f_y_pos;
					break;
				case -2:
					multiplier = data->f_y_neg;
					break;
				case +3:
					multiplier = data->f_z_pos;
					break;
				case -3:
					multiplier = data->f_z_neg;
					break;
				default:
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Invalid direction: " << direction
						<< std::endl;
					abort();
					break;
				}

				const int rel_ref_lvl = neigh_info.get<2>();
				if (rel_ref_lvl > 0) {
					// average over 4 smaller face neighbors
					multiplier /= 4.0;
				}

				data->r0 -= multiplier * neighbor_data->solution;
			}

			// initially all variables equal to residual
			data->r1 = data->r0;
			data->p0 = data->p1 = data->r0;
		}

		double dot_r_l = 0, dot_r_g = 0;
		BOOST_FOREACH(const cell_info_t& info, this->cell_info) {
			Poisson_Cell* const data = info.first;
			dot_r_l += data->r0 * data->r1;
		}
		MPI_Allreduce(&dot_r_l, &dot_r_g, 1, MPI_DOUBLE, MPI_SUM, this->comm);

		return dot_r_g;
	}

}; // class Poisson_Solve

#endif

