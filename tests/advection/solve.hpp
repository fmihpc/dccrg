/*
A solver for the advection tests of dccrg.

Copyright 2012, 2013, 2014, 2015, 2016,
2018, 2019 Finnish Meteorological Institute

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
Calculates fluxes into and out of local cells.

Calculates fluxes of inner cells without remote neighbors
if @solve_inner is true, otherwise calculates fluxes of
outer cells which have neighbors on other processes.

The total flux to copies of remote neighbors will be incorrect.
*/
template<class Grid> void calculate_fluxes(
	const double dt,
	const bool solve_inner,
	const Grid& grid,
	const int neighborhood_id = dccrg::default_neighborhood_id
) {
	const auto& cells = [&](){
		if (solve_inner) {
			return grid.inner_cells(neighborhood_id);
		} else {
			return grid.outer_cells(neighborhood_id);
		}
	}();

	for (const auto& cell: cells) {
		const double
			cell_density = cell.data->density(),
			cell_volume
				= cell.data->length_x()
				* cell.data->length_y()
				* cell.data->length_z();
		const int cell_length = grid.mapping.get_cell_length_in_indices(cell.id);

		#ifdef DEBUG
		std::vector<uint64_t> face_neighbors;
		#endif

		for (const auto& neighbor: cell.neighbors_of) {
			const int neighbor_length = grid.mapping.get_cell_length_in_indices(neighbor.id);
			// skip non-face neighbors that overlap in < 2 dimensions
			int overlaps = 0, direction = 0;

			if (neighbor.x < cell_length and neighbor.x > -neighbor_length) {
				overlaps++;
			} else if (neighbor.x == cell_length) {
				direction = 1;
			} else if (neighbor.x == -neighbor_length) {
				direction = -1;
			}

			if (neighbor.y < cell_length and neighbor.y > -neighbor_length) {
				overlaps++;
			} else if (neighbor.y == cell_length) {
				direction = 2;
			} else if (neighbor.y == -neighbor_length) {
				direction = -2;
			}

			if (neighbor.z < cell_length and neighbor.z > -neighbor_length) {
				overlaps++;
			} else if (neighbor.z == cell_length) {
				direction = 3;
			} else if (neighbor.z == -neighbor_length) {
				direction = -3;
			}

			if (overlaps < 2) {
				continue;
			}
			if (overlaps > 2) {
				const auto
					ci = grid.mapping.get_indices(cell.id),
					ni = grid.mapping.get_indices(neighbor.id);
				const auto
					cl = grid.mapping.get_cell_length_in_indices(cell.id),
					nl = grid.mapping.get_cell_length_in_indices(neighbor.id);
				std::cerr << __FILE__ "(" << __LINE__ << "): "
					<< cell.id << " <" << cl << "> (" << ci[0] << "," << ci[1] << "," << ci[2]
					<< ") <=> " << neighbor.id << " <" << nl << "> ("
					<< ni[0] << "," << ni[1] << "," << ni[2] << ")"
					<< std::endl;
				abort();
			}

			if (direction == 0) {
				continue;
			}

			#ifdef DEBUG
			face_neighbors.push_back(neighbor.id);
			#endif

			// solve flux between two local cells only in positive direction
			if (neighbor.is_local && direction < 0) {
				continue;
			}

			const double
				neighbor_density = neighbor.data->density(),
				neighbor_volume = neighbor.data->length_x() * neighbor.data->length_y() * neighbor.data->length_z();

			// get area shared between cell and current neighbor
			double min_area = -1;
			switch (direction) {
			case -1:
			case +1:
				min_area = std::min(
					cell.data->length_y() * cell.data->length_z(),
					neighbor.data->length_y() * neighbor.data->length_z()
				);
				break;
			case -2:
			case +2:
				min_area = std::min(
					cell.data->length_x() * cell.data->length_z(),
					neighbor.data->length_x() * neighbor.data->length_z()
				);
				break;
			case -3:
			case +3:
				min_area = std::min(
					cell.data->length_x() * cell.data->length_y(),
					neighbor.data->length_x() * neighbor.data->length_y()
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
				vx = (cell.data->length_x() * neighbor.data->vx() + neighbor.data->length_x() * cell.data->vx())
					/ (cell.data->length_x() + neighbor.data->length_x()),
				vy = (cell.data->length_y() * neighbor.data->vy() + neighbor.data->length_y() * cell.data->vy())
					/ (cell.data->length_y() + neighbor.data->length_y()),
				vz = (cell.data->length_z() * neighbor.data->vz() + neighbor.data->length_z() * cell.data->vz())
					/ (cell.data->length_z() + neighbor.data->length_z());

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

		#ifdef DEBUG
		const auto ref_face_neighbors = [&](){
			const auto all_neighbors = grid.get_face_neighbors_of(cell.id);
			std::vector<uint64_t> ref_face_neighbors;
			for (const auto& neigh: all_neighbors) {
				if (neigh.first != dccrg::error_cell) {
					ref_face_neighbors.push_back(neigh.first);
				}
			}
			return ref_face_neighbors;
		}();

		if (face_neighbors.size() != ref_face_neighbors.size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Unexpected number of face neighbors for cell " << cell.id
				<< " (ref lvl " << grid.mapping.get_refinement_level(cell.id)
				<< ", child of " << grid.mapping.get_parent(cell.id) << "): " << face_neighbors.size() << "  ";
			for (const auto& n: face_neighbors) {
				std::cerr << n << " (" << grid.mapping.get_refinement_level(n) << "," << grid.mapping.get_parent(n) << "),";
			}
			std::cerr << "   instead of " << ref_face_neighbors.size() << "  ";
			for (const auto& n: ref_face_neighbors) {
				std::cerr << n << " (" << grid.mapping.get_refinement_level(n) << ", " << grid.mapping.get_parent(n) << "), ";
			}
			std::cerr << std::endl;
			abort();
		}
		#endif
	}
}


/*!
Applies fluxes to local cells and zeroes the fluxes afterwards.
*/
template<class Grid> void apply_fluxes(
	Grid& grid, const int neighborhood_id = dccrg::default_neighborhood_id
) {
	for (const auto& cell: grid.local_cells(neighborhood_id)) {
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
template<class Grid> double max_time_step(
	MPI_Comm& comm,
	const Grid& grid,
	const int neighborhood_id = dccrg::default_neighborhood_id
) {
	double min_step = std::numeric_limits<double>::max();

	for (const auto& cell: grid.local_cells(neighborhood_id)) {
		const std::array<double, 3>
			current_steps{{
				cell.data->length_x() / fabs(cell.data->vx()),
				cell.data->length_y() / fabs(cell.data->vy()),
				cell.data->length_z() / fabs(cell.data->vz())
			}};

		if (std::isnormal(current_steps[0])) {
			min_step = std::min(current_steps[0], min_step);
		}
		if (std::isnormal(current_steps[1])) {
			min_step = std::min(current_steps[1], min_step);
		}
		if (std::isnormal(current_steps[2])) {
			min_step = std::min(current_steps[2], min_step);
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
