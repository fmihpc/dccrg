/*
Optimized version of a solver for the advection tests of dccrg.

Copyright 2012, 2013, 2014, 2015, 2016 Finnish Meteorological Institute
Copyright 2015, 2016 Ilja Honkonen

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

#ifndef DCCRG_ADVECTION_SOLVE_OPTIMIZED_HPP
#define DCCRG_ADVECTION_SOLVE_OPTIMIZED_HPP


#include "cmath"
#include "exception"
#include "iostream"
#include "limits"
#include "string"
#include "tuple"
#include "vector"

#include "mpi.h"

#include "dccrg.hpp"


/*!
Calculates fluxes into and out of given local cells.

The total flux to copies of remote neighbors will be incorrect.
*/
class Solver
{
public:

	/*!
	Starts solving at solve_index in dccrg's cell data pointer cache.

	Stops at end of cache or when encouters dccrg::error_cell as cell id.

	Returns last processed index + 1.
	*/
	template<
		class CellData,
		class Geometry
	> static size_t calculate_fluxes(
		const double dt,
		const size_t solve_index,
		dccrg::Dccrg<CellData, Geometry>& grid
	) {
		using std::get;

		const auto& cell_data_pointers = grid.get_cell_data_pointers();
		size_t i = solve_index;
		for ( ; i < cell_data_pointers.size(); i++) {
			const auto& cell_id = get<0>(cell_data_pointers[i]);

			if (cell_id == dccrg::error_cell) {
				break;
			}

			const auto& offset = get<2>(cell_data_pointers[i]);

			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				throw std::runtime_error("Unexpected neighbor cell in cell data pointer cache.");
			}

			#ifdef DEBUG
			if (!grid.is_local(cell_id)) {
				throw std::runtime_error("Unexpected remote cell.");
			}
			#endif

			auto* const cell_data = get<1>(cell_data_pointers[i]);

			const std::array<double, 3> cell_length = grid.geometry.get_length(cell_id);

			const double
				cell_density = cell_data->density(),
				cell_volume = cell_length[0] * cell_length[1] * cell_length[2];

			size_t processed_siblings = 0; // to be able to skip non-face neighbors
			const auto cell_ref_lvl = grid.get_refinement_level(cell_id);

			i++;
			while (i < cell_data_pointers.size()) {

				const auto& neighbor = get<0>(cell_data_pointers[i]);
				if (neighbor == dccrg::error_cell) {
					i--;
					break;
				}

				const auto& neigh_offset = get<2>(cell_data_pointers[i]);

				if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
					i--;
					break;
				}

				int direction = 0;
				bool skip_sibling = true;
				if (neigh_offset[0] == 1 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
					direction = 1;
					/*
					Skip smaller non-face neighbors, they are cached in order
					-x, -y, -z; +x, -y, -z; -x, +y, -z; ...; +x, +y, +z
					*/
					if (processed_siblings % 2 <= 2 / 2 - 1) {
						skip_sibling = false;
					}
				}
				if (neigh_offset[0] == -1 and neigh_offset[1] == 0 and neigh_offset[2] == 0) {
					direction = -1;
					if (processed_siblings % 2 >= 2 / 2) {
						skip_sibling = false;
					}
				}
				if (neigh_offset[0] == 0 and neigh_offset[1] == 1 and neigh_offset[2] == 0) {
					direction = 2;
					if (processed_siblings % 4 <= 4 / 2 - 1) {
						skip_sibling = false;
					}
				}
				if (neigh_offset[0] == 0 and neigh_offset[1] == -1 and neigh_offset[2] == 0) {
					direction = -2;
					if (processed_siblings % 4 >= 4 / 2) {
						skip_sibling = false;
					}
				}
				if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == 1) {
					direction = 3;
					if (processed_siblings % 8 <= 8 / 2 - 1) {
						skip_sibling = false;
					}
				}
				if (neigh_offset[0] == 0 and neigh_offset[1] == 0 and neigh_offset[2] == -1) {
					direction = -3;
					if (processed_siblings % 8 >= 8 / 2) {
						skip_sibling = false;
					}
				}

				if (direction == 0) {
					i++;
					continue;
				}

				// solve flux between two local cells only in positive direction
				if (grid.is_local(neighbor) && direction < 0) {
					i++;
					continue;
				}

				const auto neigh_ref_lvl = grid.get_refinement_level(neighbor);
				if (neigh_ref_lvl <= cell_ref_lvl) {
					processed_siblings = 0;
				} else {
					processed_siblings++;
					if (processed_siblings == 8) {
						processed_siblings = 0;
					}

					if (skip_sibling) {
						i++;
						continue;
					}
				}

				auto* const neighbor_data = get<1>(cell_data_pointers[i]);

				const std::array<double, 3> neighbor_length = grid.geometry.get_length(neighbor);

				const double
					neighbor_density = neighbor_data->density(),
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
					vx = (cell_length[0] * neighbor_data->vx() + neighbor_length[0] * cell_data->vx())
						/ (cell_length[0] + neighbor_length[0]),
					vy = (cell_length[1] * neighbor_data->vy() + neighbor_length[1] * cell_data->vy())
						/ (cell_length[1] + neighbor_length[1]),
					vz = (cell_length[2] * neighbor_data->vz() + neighbor_length[2] * cell_data->vz())
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
					cell_data->flux() -= flux / cell_volume;
					neighbor_data->flux() += flux / neighbor_volume;
				} else {
					cell_data->flux() += flux / cell_volume;
					neighbor_data->flux() -= flux / neighbor_volume;
				}

				i++;
			}
		}

		return i;
	}


	/*!
	Applies fluxes to local cells and zeros fluxes afterwards.
	*/
	template<
		class CellData,
		class Geometry
	> static void apply_fluxes(dccrg::Dccrg<CellData, Geometry>& grid)
	{
		using std::get;
		for (auto& item: grid.get_cell_data_pointers()) {
			const auto& cell_id = get<0>(item);

			if (cell_id == dccrg::error_cell) {
				continue;
			}

			const auto& offset = get<2>(item);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			auto* const cell_data = get<1>(item);
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
		double min_step = std::numeric_limits<double>::max();

		using std::get;

		for (auto& item: grid.get_cell_data_pointers()) {
			const auto& cell_id = get<0>(item);

			if (cell_id == dccrg::error_cell) {
				continue;
			}

			const auto& offset = get<2>(item);
			if (offset[0] != 0 or offset[1] != 0 or offset[2] != 0) {
				continue;
			}

			const auto* const cell_data = get<1>(item);

			const std::array<double, 3>
				cell_length = grid.geometry.get_length(cell_id),
				current_steps{{
					cell_length[0] / fabs(cell_data->vx()),
					cell_length[1] / fabs(cell_data->vy()),
					cell_length[2] / fabs(cell_data->vz())
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

