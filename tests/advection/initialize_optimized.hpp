/*
Optimized version of initializer for advection tests of dccrg.

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

#ifndef DCCRG_ADVECTION_INITIALIZE_OPTIMIZED_HPP
#define DCCRG_ADVECTION_INITIALIZE_OPTIMIZED_HPP

#include "cmath"
#include "iostream"
#include "vector"

#include "dccrg.hpp"

#include "solve_optimized.hpp"

/*!
Initializes the simulation's local cells and the copies of remote neighbors.

Only z direction supported at the moment.
*/
class Initialize
{
public:

	template<
		class CellData,
		class Geometry
	> void operator()(dccrg::Dccrg<CellData, Geometry>& grid)
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

			for (auto& data_item: cell_data->data) {
				data_item = 0;
			}

			const std::array<double, 3> cell_center = grid.geometry.get_center(cell_id);
			const double radius = 0.15;

			// velocities
			cell_data->vx() = Velocity().vx(cell_center[1]);
			cell_data->vy() = Velocity().vy(cell_center[0]);
			cell_data->vz() = Velocity().vz(0);

			/*
			Densities
			*/

			// smooth hump
			const double
				hump_x0 = 0.25,
				hump_y0 = 0.5,
				hump_r
					= std::min(
						std::sqrt(
							std::pow(cell_center[0] - hump_x0, 2.0)
							+ std::pow(cell_center[1] - hump_y0, 2.0)
						),
						radius
					) / radius,
				hump_density = 0.25 * (1 + std::cos(M_PI * hump_r));

			// TODO: slotted disk
			//const double disk_x0 = 0.5, disk_y0 = 0.75;
			// TODO: rotating cone
			//const double cone_x0 = 0.5, cone_y0 = 0.25;

			cell_data->density() = hump_density;
		}

		grid.update_copies_of_remote_neighbors();
	}
};

#endif

