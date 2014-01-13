/*
Particle propagator for dccrg.

Copyright 2012, 2013, 2014 Finnish Meteorological Institute

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

#include "algorithm"
#include "boost/array.hpp"
#include "boost/foreach.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "iomanip"
#include "iostream"
#include "string"
#include "zoltan.h"

#include "../../dccrg.hpp"

#include "cell.hpp"


bool Cell::transfer_particles = false;


using namespace std;
using namespace boost;
using namespace boost::mpi;
using namespace dccrg;


/*!
Propagates particles between local cells and to/from remote cells into local cells.

The velocity field is constant in space and time and only along x direction.
*/
void propagate_particles(Dccrg<Cell>& grid) {

	const double vx = 0.1;

	// propagate particles in local cells and copies of remote neighbors
	vector<uint64_t> cells = grid.get_cells();
	const std::unordered_set<uint64_t>& remote_neighbors
		= grid.get_remote_cells_on_process_boundary_internal();
	cells.insert(cells.begin(), remote_neighbors.begin(), remote_neighbors.end());

	BOOST_FOREACH(const uint64_t& cell, cells) {

		Cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << cell << endl;
			abort();
		}

		for (vector<std::array<double, 3> >::iterator
			particle = cell_data->particles.begin();
			particle != cell_data->particles.end();
			particle++
		) {
			(*particle)[0] += vx;
		}
	}

	// move particles to the particle list of the cell the particles are currently inside of
	BOOST_FOREACH(const uint64_t& previous_cell, cells) {

		Cell* previous_data = grid[previous_cell];
		if (previous_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__ << " No data for cell " << previous_cell << endl;
			abort();
		}

		vector<std::array<double, 3> >::size_type i = 0;
		while (i < previous_data->particles.size()) {

			// hande grid wrap around
			previous_data->particles[i] = grid.geometry.get_real_coordinate(previous_data->particles[i]);

			const uint64_t current_cell = grid.get_existing_cell(previous_data->particles[i]);

			// do nothing if particle hasn't changed cell
			if (current_cell == previous_cell) {
				i++;
				continue;
			}

			// change only local cell data because remote data will be overwritten anyway

			// add particle to the current cell's list
			if (grid.is_local(current_cell)) {
				Cell* current_data = grid[current_cell];
				current_data->particles.push_back(previous_data->particles[i]);
				current_data->number_of_particles = current_data->particles.size();
			}

			// remove particle from its previous cell
			if (grid.is_local(previous_cell)) {
				previous_data->particles.erase(previous_data->particles.begin() + i);
				previous_data->number_of_particles = previous_data->particles.size();
			} else {
				i++;
			}
		}
	}
}


/*!
Saves the local grid and its particles to separate vtk files.

Each process saves its own files.
*/
void save(communicator comm, const Dccrg<Cell>& grid, unsigned int step)
{
	// write the grid
	const string grid_file_name(
		lexical_cast<string>("simple_")
		+ lexical_cast<string>(comm.rank()) + "_"
		+ lexical_cast<string>(step) + "_grid.vtk"
	);
	grid.write_vtk_file(grid_file_name.c_str());

	// write the particles
	const string outname(
		lexical_cast<string>("simple_")
		+ lexical_cast<string>(comm.rank()) + "_"
		+ lexical_cast<string>(step) + ".vtk"
	);
	ofstream outfile;
	outfile.open(outname.c_str());

	const vector<uint64_t> cells = grid.get_cells();

	// write header
	outfile << "# vtk DataFile Version 2.0\n"
		"Particle test\nASCII\nDATASET UNSTRUCTURED_GRID\n";

	// calculate the total number of particles
	uint64_t total_particles = 0;
	BOOST_FOREACH(uint64_t cell, cells) {
		Cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			cerr << __FILE__ << ":" << __LINE__ << ": no data for cell " << cell << endl;
			abort();
		}

		total_particles += cell_data->number_of_particles;
	}
	outfile << "POINTS " << total_particles << " float\n";

	// write out coordinates
	uint64_t written_particles = 0;
	BOOST_FOREACH(uint64_t cell, cells) {

		Cell* cell_data = grid[cell];
		for (vector<std::array<double, 3> >::const_iterator
			coordinates = cell_data->particles.begin();
			coordinates != cell_data->particles.end();
			coordinates++, written_particles++
		) {
			outfile << (*coordinates)[0] << "\t"
				<< (*coordinates)[1] << "\t"
				<< (*coordinates)[2] << "\n";
		}
	}
	if (written_particles != total_particles) {
		cerr << __FILE__ << ":" << __LINE__
			<< ": Incorrect number of particles written: " << written_particles
			<< ", should be " << total_particles
			<< endl;
		abort();
	}

	outfile << "CELLS " << total_particles << " " << 2 * total_particles << "\n";
	for (uint64_t i = 0; i < total_particles; i++) {
		outfile << "1 " << i << "\n";
	}

	outfile << "CELL_TYPES " << total_particles << "\n";
	for (uint64_t i = 0; i < total_particles; i++) {
		outfile << "1 ";
	}

	// process numbers of particles
	outfile << "\nCELL_DATA " << total_particles
		<< "\nSCALARS process int 1\nLOOKUP_TABLE default\n";
	for (uint64_t i = 0; i < total_particles; i++) {
		outfile << comm.rank() << " ";
	}
	outfile << "\n";

	// cell numbers of particles
	outfile << "SCALARS cell int 1\nLOOKUP_TABLE default\n";
	BOOST_FOREACH(uint64_t cell, cells) {
		Cell* cell_data = grid[cell];
		for (vector<std::array<double, 3> >::const_iterator
			coordinates = cell_data->particles.begin();
			coordinates != cell_data->particles.end();
			coordinates++, written_particles++
		) {
			outfile << cell << " ";
		}
	}
	outfile << "\n";

	outfile.close();
}


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}

	// initialize grid
	Dccrg<Cell> grid;
	const std::array<uint64_t, 3> grid_length = {{3, 1, 1}};
	grid.initialize(grid_length, comm, "RANDOM", 1, 0, true, true, true);

	const vector<uint64_t> cells = grid.get_cells();

	// initial condition
	const unsigned int max_particles_per_cell = 5;
	BOOST_FOREACH(uint64_t cell, cells) {

		Cell* cell_data = grid[cell];

		const std::array<double, 3>
			cell_min = grid.geometry.get_min(cell),
			cell_max = grid.geometry.get_max(cell);

		const unsigned int number_of_particles
			= (unsigned int)ceil(max_particles_per_cell * double(rand()) / RAND_MAX);
		for (unsigned int i = 0; i < number_of_particles; i++) {
			std::array<double, 3> coordinates = {{
				cell_min[0] + (cell_max[0] - cell_min[0]) * double(rand()) / RAND_MAX,
				cell_min[1] + (cell_max[1] - cell_min[1]) * double(rand()) / RAND_MAX,
				cell_min[2] + (cell_max[2] - cell_min[2]) * double(rand()) / RAND_MAX
			}};

			cell_data->particles.push_back(coordinates);
			cell_data->number_of_particles = cell_data->particles.size();
		}
	}

	/*
	Visualize the results for example with visit -o simple_particles.visit
	or visit -o simple_grid.visit or overlay them both from the user interface.
	*/
	const string visit_particles_name("simple_particles.visit"),
		visit_grid_name("simple_grid.visit");

	ofstream visit_particles, visit_grid;

	if (comm.rank() == 0) {
		visit_particles.open(visit_particles_name.c_str());
		visit_particles << "!NBLOCKS " << comm.size() << "\n";
		visit_grid.open(visit_grid_name.c_str());
		visit_grid << "!NBLOCKS " << comm.size() << "\n";
	}

	const unsigned int max_steps = 50;
	for (unsigned int step = 0; step < max_steps; step++) {

		// append current output file names to the visit files
		if (comm.rank() == 0) {
			for (int i = 0; i < comm.size(); i++) {
				visit_particles << "simple_"
					<< i << "_"
					<< step << ".vtk\n";
				visit_grid << "simple_"
					<< i << "_"
					<< step << "_grid.vtk\n";
			}
		}

		save(comm, grid, step);

		// reserve space for incoming particles in copies of remote neighbors
		Cell::transfer_particles = false;
		grid.update_copies_of_remote_neighbors();

		const std::unordered_set<uint64_t>& remote_neighbors
			= grid.get_remote_cells_on_process_boundary_internal();

		BOOST_FOREACH(uint64_t remote_neighbor, remote_neighbors) {
			Cell* data = grid[remote_neighbor];
			data->resize();
		}

		// update particle data between neighboring cells on different processes
		Cell::transfer_particles = true;
		grid.update_copies_of_remote_neighbors();

		propagate_particles(grid);
	}

	// append final output file names to the visit files
	if (comm.rank() == 0) {
		for (int i = 0; i < comm.size(); i++) {
				visit_particles << "simple_"
					<< comm.rank() << "_"
					<< max_steps << ".vtk\n";
				visit_grid << "simple_"
					<< comm.rank() << "_"
					<< max_steps << "_grid.vtk\n";
		}
		visit_particles.close();
		visit_grid.close();
	}

	save(comm, grid, max_steps);

	return EXIT_SUCCESS;
}

