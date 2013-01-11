/*
Parallel Poisson solver build on top of dccrg.

Copyright 2012, 2013 Finnish Meteorological Institute

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
#include "boost/mpi.hpp"
#include "boost/program_options.hpp"
#include "boost/static_assert.hpp"
#include "cmath"
#include "cstdlib"
#include "iostream"
#include "mpi.h"
#include "stdint.h"
#include "vector"

#include "dccrg.hpp"

#include "grid_support.hpp"
#include "reference_poisson_solve.hpp"


/*!
Dccrg cell data class used in the parallel Poisson solver Poisson_Solve(...).
*/
class Poisson_Cell
{
public:

	Poisson_Cell()
	{
		this->rhs = this->solution = 0;
	}

	// rhs is b and solution is x in the equation A . x = b
	double rhs, solution;

	/*
	Same variables are used as in 2.7.6 of numerical recipes
	with the modification that bold r is r0 and bold r with
	dash on top is r1; same for p.
	*/
	double p0, p1, r0, r1;

	/*
	The recurrence in numerical recipes does A . p0 two times,
	this variable caches the result from first to save CPU time.
	*/
	double A_dot_p;

	/*
	Factor by which solution variables in this cell are scaled
	so that diagonal of matrix A has only ones
	*/
	double scaling_factor;

	/*
	Factors that data in neighboring cells is multiplied with
	when calculating A . i for each local cell where i = p0, p1.
	The same factors are used when calculating traspose(A) . i.
	One factor per direction.
	*/
	double f_x_pos, f_x_neg, f_y_pos, f_y_neg, f_z_pos, f_z_neg;

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
	void mpi_datatype(
		void*& address,
		int& count,
		MPI_Datatype& datatype,
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/
	) {
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
	}
};

int Poisson_Cell::transfer_switch = Poisson_Cell::INIT;


/*!
Solves the Poisson's equation.

The right hand side (rhs) of the equation must have been initialized
before calling this function.
The solution (x) of the equation also must have been initialized
(with a guess) before calling this function.

Given grid must have a default neighborhood of size 0? or ...

If grid non-periodic missing neighbors are assumed to be of
equal size and have equal value (of whatever comes after the dot in A . i)

WARNING: doesn't work yet
*/
void Poisson_Solve(
	/* TODO: overlap computation with communication
	const std::vector<uint64_t>& inner_cells,
	const std::vector<uint64_t>& outer_cells,*/
	const std::vector<uint64_t>& cells,
	dccrg::Dccrg<Poisson_Cell>& grid
) {
	MPI_Comm comm = grid.get_comm();

	// make sure copies of remote neighbors exist before caching and
	// transfer user's guess for the solution to calculate residual
	Poisson_Cell::transfer_switch = Poisson_Cell::INIT;
	grid.update_remote_neighbor_data();

	/*
	Cache pointers to local cell data, their neighborhood and sizes
	*/

	/*
	Data of cell's neighbor, neighbor's direction from cell
	and neighbor's relative refinement level (if > 0 neighbor is smaller)
	*/
	typedef typename boost::tuple<Poisson_Cell*, direction_t, int> neighbor_info_t;

	// data of a local cell and info of its neighbors
	typedef typename std::pair<Poisson_Cell*, std::vector<neighbor_info_t> > cell_info_t;

	std::vector<cell_info_t> cell_info;
	cell_info.reserve(cells.size());

	BOOST_FOREACH(const uint64_t cell, cells) {

		cell_info_t temp_cell_info;

		Poisson_Cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< std::endl;
			abort();
		}

		temp_cell_info.first = cell_data;

		const int cell_ref_lvl = grid.get_refinement_level(cell);
		const double
			cell_x_half_size = grid.get_cell_x_size(cell) / 2.0,
			cell_y_half_size = grid.get_cell_y_size(cell) / 2.0,
			cell_z_half_size = grid.get_cell_z_size(cell) / 2.0;

		// get face neighbors of current cell
		std::vector<uint64_t> face_neighbors;
		std::vector<direction_t> directions;
		get_face_neighbors<Poisson_Cell>(cell, grid, face_neighbors, directions);


		/*
		Get cell centers of neighbors relative to current cell
		*/

		// non-existing neighbors are equal in size to current cell
		double
			neigh_pos_x_offset = +2 * cell_x_half_size,
			neigh_neg_x_offset = -2 * cell_x_half_size,
			neigh_pos_y_offset = +2 * cell_y_half_size,
			neigh_neg_y_offset = -2 * cell_y_half_size,
			neigh_pos_z_offset = +2 * cell_z_half_size,
			neigh_neg_z_offset = -2 * cell_z_half_size;

		for (size_t i = 0; i < face_neighbors.size(); i++) {
			const direction_t direction = directions[i];
			const uint64_t neighbor = face_neighbors[i];
			const double
				neigh_x_half_size = grid.get_cell_x_size(neighbor) / 2.0,
				neigh_y_half_size = grid.get_cell_y_size(neighbor) / 2.0,
				neigh_z_half_size = grid.get_cell_z_size(neighbor) / 2.0;

			// assume rhs and solution are cell-centered
			switch(direction) {
			case POS_X:
				neigh_pos_x_offset = cell_x_half_size + neigh_x_half_size;
				break;
			case NEG_X:
				neigh_neg_x_offset = -1.0 * (cell_x_half_size + neigh_x_half_size);
				break;
			case POS_Y:
				neigh_pos_y_offset = cell_y_half_size + neigh_y_half_size;
				break;
			case NEG_Y:
				neigh_neg_y_offset = -1.0 * (cell_y_half_size + neigh_y_half_size);
				break;
			case POS_Z:
				neigh_pos_z_offset = cell_z_half_size + neigh_z_half_size;
				break;
			case NEG_Z:
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

		/*
		Set geometry factors of current cell
		*/

		const double
			total_offset_x = neigh_pos_x_offset - neigh_neg_x_offset,
			total_offset_y = neigh_pos_y_offset - neigh_neg_y_offset,
			total_offset_z = neigh_pos_z_offset - neigh_neg_z_offset;

		// geometry factorsare 0 in directions without neighbors
		cell_data->f_x_pos =
		cell_data->f_x_neg =
		cell_data->f_y_pos =
		cell_data->f_y_neg =
		cell_data->f_z_pos =
		cell_data->f_z_neg = 0;

		for (size_t i = 0; i < face_neighbors.size(); i++) {
			const direction_t direction = directions[i];

			// don't mind extra work due to 4 smaller face neighbors
			switch(direction) {
			case POS_X:
				cell_data->f_x_pos = +2.0 / (neigh_pos_x_offset * total_offset_x);
				break;
			case NEG_X:
				cell_data->f_x_neg = -2.0 / (neigh_neg_x_offset * total_offset_x);
				break;
			case POS_Y:
				cell_data->f_y_pos = +2.0 / (neigh_pos_y_offset * total_offset_y);
				break;
			case NEG_Y:
				cell_data->f_y_neg = -2.0 / (neigh_neg_y_offset * total_offset_y);
				break;
			case POS_Z:
				cell_data->f_z_pos = +2.0 / (neigh_pos_z_offset * total_offset_z);
				break;
			case NEG_Z:
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
			= 0
			- cell_data->f_x_pos
			- cell_data->f_x_neg
			- cell_data->f_y_pos
			- cell_data->f_y_neg
			- cell_data->f_z_pos
			- cell_data->f_z_neg;

		// cache neighbor info
		for (size_t i = 0; i < face_neighbors.size(); i++) {

			neighbor_info_t temp_neigh_info;

			const direction_t direction = directions[i];
			temp_neigh_info.get<1>() = direction;

			const uint64_t neighbor = face_neighbors[i];
			Poisson_Cell* neighbor_data = grid[neighbor];
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

		cell_info.push_back(temp_cell_info);
	}

	/*
	Scale values of A in A . i and due to that also neighbor factors
	*/

	// at this point only correct scaling_factor is needed from remote copies
	Poisson_Cell::transfer_switch = Poisson_Cell::GEOMETRY;
	grid.update_remote_neighbor_data();

	BOOST_FOREACH(cell_info_t& info, cell_info) {
		Poisson_Cell* data = info.first;
		data->solution *= data->scaling_factor;

		// neighbors' values are scaled using their own factor
		bool
			pos_x_done = false,
			neg_x_done = false,
			pos_y_done = false,
			neg_y_done = false,
			pos_z_done = false,
			neg_z_done = false;

		BOOST_FOREACH(neighbor_info_t neigh_info, info.second) {
			Poisson_Cell* neighbor_data = neigh_info.get<0>();
			const direction_t direction = neigh_info.get<1>();

			switch(direction) {
			case POS_X:
				if (!pos_x_done) {
					pos_x_done = true;
					data->f_x_pos /= neighbor_data->scaling_factor;
				}
				break;
			case NEG_X:
				if (!neg_x_done) {
					neg_x_done = true;
					data->f_x_neg /= neighbor_data->scaling_factor;
				}
				break;
			case POS_Y:
				if (!pos_y_done) {
					pos_y_done = true;
					data->f_y_pos /= neighbor_data->scaling_factor;
				}
				break;
			case NEG_Y:
				if (!neg_y_done) {
					neg_y_done = true;
					data->f_y_neg /= neighbor_data->scaling_factor;
				}
				break;
			case POS_Z:
				if (!pos_z_done) {
					pos_z_done = true;
					data->f_z_pos /= neighbor_data->scaling_factor;
				}
				break;
			case NEG_Z:
				if (!neg_z_done) {
					neg_z_done = true;
					data->f_z_neg /= neighbor_data->scaling_factor;
				}
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid direction: " << direction
					<< std::endl;
				abort();
				break;
			}
		}
	}

	// update (rest of) corrected scaling factors
	grid.update_remote_neighbor_data();

	// update scaled solution
	Poisson_Cell::transfer_switch = Poisson_Cell::INIT;
	grid.update_remote_neighbor_data();


	/*
	Initialize solver
	*/

	// residual == r0 = rhs - A . solution
	BOOST_FOREACH(cell_info_t& info, cell_info) {
		Poisson_Cell* data = info.first;

		data->r0 = data->rhs - data->solution;
		std::cout << "r0: " << data->r0 << " ";

		BOOST_FOREACH(neighbor_info_t neigh_info, info.second) {
			Poisson_Cell* neighbor_data = neigh_info.get<0>();

			// final multiplier to use for current neighbor's data
			double multiplier = 0;
			const direction_t direction = neigh_info.get<1>();
			switch(direction) {
			case POS_X:
				multiplier = data->f_x_pos;
				break;
			case NEG_X:
				multiplier = data->f_x_neg;
				break;
			case POS_Y:
				multiplier = data->f_y_pos;
				break;
			case NEG_Y:
				multiplier = data->f_y_neg;
				break;
			case POS_Z:
				multiplier = data->f_z_pos;
				break;
			case NEG_Z:
				multiplier = data->f_z_neg;
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid direction: " << direction
					<< std::endl;
				abort();
				break;
			}

			/* TODO: put into f_*_* ?
			const int rel_ref_lvl = neigh_info.get<2>();
			if (rel_ref_lvl > 0) {
				// average over 4 smaller face neighbors
				multiplier /= 4.0;
			}*/

			data->r0 -= multiplier * neighbor_data->solution;
			std::cout << data->r0 << " ";
		}
		std::cout << std::endl;

		// initially all variables equal to residual
		data->r1 = data->p0 = data->p1 = data->r0;
	}

	// print total residual
	double res_l = 0, res_g = 0;
	BOOST_FOREACH(cell_info_t& info, cell_info) {
		Poisson_Cell* data = info.first;
		res_l += fabs(data->r0);
	}
	MPI_Reduce(&res_l, &res_g, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	const int rank = grid.get_rank();
	if (rank == 0) {
		std::cout << "Global residual after init: " << res_g << std::endl;
	}


	/*
	Solve
	*/

	// local and global versions of scalars used by solver
	double
		alpha = 0, beta = 0,
		// r0 . r1
		dot_r_l = 0, dot_r_g = 0,
		// old r0 . r1, new is needed in the middle of each iteration
		old_dot_r_g = 0;

	// for debugging
	std::vector<double> alphas, betas, p_dots, r_dots, old_r_dots;


	// calculate initial value of r0 . r1
	BOOST_FOREACH(cell_info_t& info, cell_info) {
		Poisson_Cell* data = info.first;
		dot_r_l += data->r0 * data->r1;
	}
	MPI_Allreduce(&dot_r_l, &dot_r_g, 1, MPI_DOUBLE, MPI_SUM, comm);
	//std::cout << "r0 . r1: " << dot_r_g << std::endl;


	size_t step = 0;
	const size_t steps = 10;
	do {
		step++;

		// TODO: update only p0 for calculating alpha and only then p1?
		Poisson_Cell::transfer_switch = Poisson_Cell::SOLVING;
		grid.update_remote_neighbor_data();

		/*
		Calculate alpha
		*/

		// A . p0, cache the result
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;

			data->A_dot_p = data->p0;

			BOOST_FOREACH(neighbor_info_t neigh_info, info.second) {
				Poisson_Cell* neighbor_data = neigh_info.get<0>();

				double multiplier = 0;
				const direction_t direction = neigh_info.get<1>();
				switch(direction) {
				case POS_X:
					multiplier = data->f_x_pos;
					break;
				case NEG_X:
					multiplier = data->f_x_neg;
					break;
				case POS_Y:
					multiplier = data->f_y_pos;
					break;
				case NEG_Y:
					multiplier = data->f_y_neg;
					break;
				case POS_Z:
					multiplier = data->f_z_pos;
					break;
				case NEG_Z:
					multiplier = data->f_z_neg;
					break;
				default:
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Invalid direction: " << direction
						<< std::endl;
					abort();
					break;
				}

				/*const int rel_ref_lvl = neigh_info.get<2>();
				if (rel_ref_lvl > 0) {
					multiplier /= 4.0;
				}*/

				data->A_dot_p += multiplier * neighbor_data->p0;
			}
		}

		// p1 . (A . p0)
		double dot_p_l = 0, dot_p_g = 0;
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;
			// only values used in A . i are scaled
			dot_p_l += data->p1 / data->scaling_factor * data->A_dot_p;
		}
		MPI_Allreduce(&dot_p_l, &dot_p_g, 1, MPI_DOUBLE, MPI_SUM, comm);
		//std::cout << "p1 . (A . p0): " << dot_p_g << std::endl;
		p_dots.push_back(dot_p_g);

		alpha = dot_r_g / dot_p_g;
		//std::cout << "alpha: " << alpha << std::endl;
		alphas.push_back(alpha);


		// update solution
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;
			data->solution += alpha * data->p0;
		}

		// update r0
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;
			data->r0 -= alpha * data->A_dot_p;
		}

		// update r1
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;

			data->r1 -= alpha * data->p1;

			BOOST_FOREACH(neighbor_info_t neigh_info, info.second) {
				Poisson_Cell* neighbor_data = neigh_info.get<0>();

				/*
				In transpose(A) . p1 use that multiplier which was calculated
				for updating the current neighbor's value with this cell's,
				i.e. just "reverse" the direction of the update.
				*/
				double multiplier = 0;
				const direction_t direction = neigh_info.get<1>();
				switch(direction) {
				case POS_X:
					multiplier = neighbor_data->f_x_neg;
					break;
				case NEG_X:
					multiplier = neighbor_data->f_x_pos;
					break;
				case POS_Y:
					multiplier = neighbor_data->f_y_neg;
					break;
				case NEG_Y:
					multiplier = neighbor_data->f_y_pos;
					break;
				case POS_Z:
					multiplier = neighbor_data->f_z_neg;
					break;
				case NEG_Z:
					multiplier = neighbor_data->f_z_pos;
					break;
				default:
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Invalid direction: " << direction
						<< std::endl;
					abort();
					break;
				}

				/*const int rel_ref_lvl = neigh_info.get<2>();
				if (rel_ref_lvl < 0) {
					multiplier /= 4.0;
				}*/

				data->r1 -= alpha * multiplier * neighbor_data->p1;
			}
		}

		// calculate beta
		old_dot_r_g = dot_r_g;
		dot_r_l = dot_r_g = 0;
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;
			dot_r_l += data->r0 * data->r1;
		}
		MPI_Allreduce(&dot_r_l, &dot_r_g, 1, MPI_DOUBLE, MPI_SUM, comm);
		//std::cout << "new r0 . r1: " << dot_r_g << std::endl;
		r_dots.push_back(dot_r_g);
		old_r_dots.push_back(old_dot_r_g);

		beta = dot_r_g / old_dot_r_g;
		//std::cout << "beta: " << beta << std::endl;
		betas.push_back(beta);


		// update p0, p1
		BOOST_FOREACH(cell_info_t& info, cell_info) {
			Poisson_Cell* data = info.first;
			data->p0 = data->r0 + beta * data->p0;
			data->p1 = data->r1 + beta * data->p1;
		}

	} while (step < steps);

	std::cout << "    alpha        beta  old_r_dot  r_dot  p1_dot_A_dot_p0\n";
	for (size_t i = 0; i < alphas.size(); i++) {
		std::cout << alphas[i] << "\t" << betas[i] << "\t" <<  old_r_dots[i] << "\t" << r_dots[i] << "\t" << p_dots[i] << "\n";
	}
}

#endif

