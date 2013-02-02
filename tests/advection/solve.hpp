/*
A solver for the advection tests of dccrg.

Copyright 2012, 2013 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DCCRG_ADVECTION_SOLVE_HPP
#define DCCRG_ADVECTION_SOLVE_HPP


#include "boost/foreach.hpp"
#include "cmath"
#include "iostream"
#include "limits"
#include "string"
#include "vector"

#include "dccrg.hpp"


/*!
Calculates fluxes into and out of given local cells.

The total flux to copies of remote neighbors will be incorrect.
*/
class Solver
{
public:

	template<
		class CellData,
		class Geometry
	> static void calculate_fluxes(
		const double dt,
		const std::vector<uint64_t>& cells,
		dccrg::Dccrg<CellData, Geometry>& grid
	) {

		BOOST_FOREACH(const uint64_t& cell, cells) {

			if (!grid.is_local(cell)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell
					<< " isn't local"
					<< std::endl;
				abort();
			}

			CellData* cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell
					<< std::endl;
				abort();
			}

			const double cell_density = cell_data->density(),
				cell_length_x = grid.get_cell_length_x(cell),
				cell_length_y = grid.get_cell_length_y(cell),
				cell_length_z = grid.get_cell_length_z(cell),
				cell_volume = cell_length_x * cell_length_y * cell_length_z;

			const std::vector<std::pair<uint64_t, int> > neighbors_to_solve
				= grid.get_face_neighbors_of(cell);

			for (size_t i = 0; i < neighbors_to_solve.size(); i++) {

				const uint64_t neighbor = neighbors_to_solve[i].first;
				const int direction = neighbors_to_solve[i].second;

				// solve flux between two local cells only in positive direction
				if (grid.is_local(neighbor) && direction < 0) {
					continue;
				}

				CellData* neighbor_data = grid[neighbor];
				if (neighbor_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No data for cell " << neighbor
						<< std::endl;
					abort();
				}

				const double neighbor_density = neighbor_data->density(),
					neighbor_length_x = grid.get_cell_length_x(neighbor),
					neighbor_length_y = grid.get_cell_length_y(neighbor),
					neighbor_length_z = grid.get_cell_length_z(neighbor),
					neighbor_volume = neighbor_length_x * neighbor_length_y * neighbor_length_z;

				// get area shared between cell and current neighbor
				double min_area = -1;
				switch (direction) {
				case -1:
				case +1:
					min_area = std::min(
						cell_length_y * cell_length_z,
						neighbor_length_y * neighbor_length_z
					);
					break;
				case -2:
				case +2:
					min_area = std::min(
						cell_length_x * cell_length_z,
						neighbor_length_x * neighbor_length_z
					);
					break;
				case -3:
				case +3:
					min_area = std::min(
						cell_length_x * cell_length_y,
						neighbor_length_x * neighbor_length_y
					);
					break;
				}

				/*
				Solve flux
				*/

				// positive flux through a face goes into positive direction
				double flux = 0;

				// velocity interpolated to shared face
				const double
					vx = (cell_length_x * neighbor_data->vx() + neighbor_length_x * cell_data->vx())
						/ (cell_length_x + neighbor_length_x),
					vy = (cell_length_y * neighbor_data->vy() + neighbor_length_y * cell_data->vy())
						/ (cell_length_y + neighbor_length_y),
					vz = (cell_length_z * neighbor_data->vz() + neighbor_length_z * cell_data->vz())
						/ (cell_length_z + neighbor_length_z);

				switch (direction) {
				case +1:
					if (vx >= 0) {
						flux = cell_density * dt * vx * min_area;
					} else {
						flux = neighbor_density * dt * vx * min_area;
					}
					break;

				case +2:
					if (vy >= 0) {
						flux = cell_density * dt * vy * min_area;
					} else {
						flux = neighbor_density * dt * vy * min_area;
					}
					break;

				case +3:
					if (vz >= 0) {
						flux = cell_density * dt * vz * min_area;
					} else {
						flux = neighbor_density * dt * vz * min_area;
					}
					break;

				case -1:
					if (vx >= 0) {
						flux = neighbor_density * dt * vx * min_area;
					} else {
						flux = cell_density * dt * vx * min_area;
					}
					break;

				case -2:
					if (vy >= 0) {
						flux = neighbor_density * dt * vy * min_area;
					} else {
						flux = cell_density * dt * vy * min_area;
					}
					break;

				case -3:
					if (vz >= 0) {
						flux = neighbor_density * dt * vz * min_area;
					} else {
						flux = cell_density * dt * vz * min_area;
					}
					break;
				}

				// save flux
				if (direction > 0) {
					cell_data->flux() -= flux / cell_volume;
					neighbor_data->flux() += flux / neighbor_volume;
				} else {
					cell_data->flux() += flux / cell_volume;
					neighbor_data->flux() -= flux / neighbor_volume;
				}
			}
		}
	}


	/*!
	Applies fluxes to local cells and zeroes the fluxes afterwards.
	*/
	template<
		class CellData,
		class Geometry
	> static void apply_fluxes(dccrg::Dccrg<CellData, Geometry>& grid)
	{
		const std::vector<uint64_t> cells = grid.get_cells();

		BOOST_FOREACH(const uint64_t& cell, cells) {
			CellData* cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "No data for cell " << cell
					<< std::endl;
				abort();
			}

			cell_data->density() += cell_data->flux();
			cell_data->flux() = 0;
		}
	}


	/*!
	Returns the largest allowed global time step.

	Must be called simultaneously on all processes.
	Assumes that cells of same refinement level have the same size
	per dimension.
	*/
	template<
		class CellData,
		class Geometry
	> static double max_time_step(
		MPI_Comm& comm,
		const dccrg::Dccrg<CellData, Geometry>& grid
	) {
		const std::vector<uint64_t> cells = grid.get_cells();

		double min_step = std::numeric_limits<double>::max();

		BOOST_FOREACH(const uint64_t& cell_id, cells) {

			CellData* cell = grid[cell_id];
			if (cell == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell " << cell_id
					<< std::endl;
				abort();
			}

			const double min_step_x = grid.get_cell_length_x(cell_id) / fabs(cell->vx()),
				min_step_y = grid.get_cell_length_y(cell_id) / fabs(cell->vy()),
				min_step_z = grid.get_cell_length_z(cell_id) / fabs(cell->vz()),
				current_min_step = std::min(min_step_x, std::min(min_step_y, min_step_z));

			if (min_step > current_min_step) {
				min_step = current_min_step;
			}
		}

		double ret_val = 0;
		if (
			MPI_Allreduce(
				&min_step,
				&ret_val,
				1,
				MPI_DOUBLE,
				MPI_MIN,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "MPI_Allreduce failed."
				<< std::endl;
			abort();
		}

		return ret_val;
	}

}; // class Solver


class Velocity
{
public:

	static double vx(const double y)
	{
		return -y + 0.5;
	}

	static double vy(const double x)
	{
		return +x - 0.5;
	}

	static double vz(const double /*a*/)
	{
		return 0;
	}
};


class File_Namer
{
public:

	std::string operator()(const double time_step, const std::string& basename)
	{
		std::ostringstream step_string;
		step_string << std::setw(7) << std::setfill('0') << int(time_step * 1000) << "_ms";
		return basename + step_string.str();
	}
};

#endif

