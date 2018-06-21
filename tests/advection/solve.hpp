/*
A solver for the advection tests of dccrg.

Copyright 2012, 2013, 2014, 2015, 2016 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DCCRG_ADVECTION_SOLVE_HPP
#define DCCRG_ADVECTION_SOLVE_HPP


#include "cmath"
#include "iostream"
#include "limits"
#include "string"
#include "vector"

#include "mpi.h"

#include "dccrg.hpp"


/*!
Calculates fluxes into and out of given local cells.

The total flux to copies of remote neighbors will be incorrect.
*/
template<class Cell_Data, class Geometry> void calculate_fluxes(
	const double dt,
	const bool solve_inner,
	const dccrg::Dccrg<Cell_Data, Geometry>& grid
) {
	const auto& cells = [&](){
		if (solve_inner) {
			return grid.inner_cells;
		} else {
			return grid.outer_cells;
		}
	}();

	for (const auto& cell: cells) {
		const std::array<double, 3> cell_length = grid.geometry.get_length(cell.id);

		const double
			cell_density = cell.data->density(),
			cell_volume = cell_length[0] * cell_length[1] * cell_length[2];

		for (const auto& neighbor: cell.neighbors_of) {
			int direction = 0;
			if (neighbor.x == 1 and neighbor.y == 0 and neighbor.z == 0) {
				direction = 1;
			}
			if (neighbor.x == -1 and neighbor.y == 0 and neighbor.z == 0) {
				direction = -1;
			}
			if (neighbor.x == 0 and neighbor.y == -1 and neighbor.z == 0) {
				direction = 2;
			}
			if (neighbor.x == 0 and neighbor.y == -2 and neighbor.z == 0) {
				direction = -2;
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == 1) {
				direction = 3;
			}
			if (neighbor.x == 0 and neighbor.y == 0 and neighbor.z == -1) {
				direction = -3;
			}
			// skip diagonal and other neighbors
			if (direction == 0) {
				continue;
			}

			// solve flux between two local cells only in positive direction
			if (grid.is_local(neighbor.id) && direction < 0) {
				continue;
			}

			const std::array<double, 3> neighbor_length = grid.geometry.get_length(neighbor.id);

			const double
				neighbor_density = neighbor.data->density(),
				neighbor_volume = neighbor_length[0] * neighbor_length[1] * neighbor_length[2];

			// get area shared between cell and current neighbor
			double min_area = -1;
			switch (direction) {
			case -1:
			case +1:
				min_area = std::min(
					cell_length[1] * cell_length[2],
					neighbor_length[1] * neighbor_length[2]
				);
				break;
			case -2:
			case +2:
				min_area = std::min(
					cell_length[0] * cell_length[2],
					neighbor_length[0] * neighbor_length[2]
				);
				break;
			case -3:
			case +3:
				min_area = std::min(
					cell_length[0] * cell_length[1],
					neighbor_length[0] * neighbor_length[1]
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
				vx = (cell_length[0] * neighbor.data->vx() + neighbor_length[0] * cell.data->vx())
					/ (cell_length[0] + neighbor_length[0]),
				vy = (cell_length[1] * neighbor.data->vy() + neighbor_length[1] * cell.data->vy())
					/ (cell_length[1] + neighbor_length[1]),
				vz = (cell_length[2] * neighbor.data->vz() + neighbor_length[2] * cell.data->vz())
					/ (cell_length[2] + neighbor_length[2]);

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
				cell.data->flux() -= flux / cell_volume;
				neighbor.data->flux() += flux / neighbor_volume;
			} else {
				cell.data->flux() += flux / cell_volume;
				neighbor.data->flux() -= flux / neighbor_volume;
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
> void apply_fluxes(dccrg::Dccrg<CellData, Geometry>& grid) {
	for (const auto& cell: grid.local_cells) {
		cell.data->density() += cell.data->flux();
		cell.data->flux() = 0;
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
> double max_time_step(
	MPI_Comm& comm,
	const dccrg::Dccrg<CellData, Geometry>& grid
) {
	double min_step = std::numeric_limits<double>::max();

	for (const auto& cell: grid.local_cells) {
		const std::array<double, 3>
			cell_length = grid.geometry.get_length(cell.id),
			current_steps{{
				cell_length[0] / fabs(cell.data->vx()),
				cell_length[1] / fabs(cell.data->vy()),
				cell_length[2] / fabs(cell.data->vz())
			}};

		const double current_min_step =
			std::min(current_steps[0],
			std::min(current_steps[1], current_steps[2]));

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


double get_vx(const double y) {
	return -y + 0.5;
}

double get_vy(const double x) {
	return +x - 0.5;
}

double get_vz(const double /*a*/) {
	return 0;
}


std::string get_file_name(const double time_step, const std::string& basename) {
	std::ostringstream step_string;
	step_string << std::setw(7) << std::setfill('0') << int(time_step * 1000) << "_ms";
	return basename + step_string.str();
}

#endif
