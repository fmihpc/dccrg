/*
Program for testing dccrg restart with varying amounts of cell data.

Copyright 2010, 2011, 2012, 2013, 2014,
2015, 2016 Finnish Meteorological Institute
Copyright 2014, 2015, 2016 Ilja Honkonen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "array"
#include "algorithm"
#include "boost/program_options.hpp"
#include "cstdlib"
#include "iostream"
#include "string"
#include "tuple"
#include "vector"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg.hpp"
#include "../../dccrg_cartesian_geometry.hpp"


using namespace std;


struct Cell
{
public:
	uint64_t data_size;
	std::vector<int> data;

	static bool transfer_all, transfer_data;

	std::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype() const
	{
		if (Cell::transfer_all) {
			std::array<int, 2> counts = {{1, int(this->data.size())}};
			std::array<MPI_Aint, 2> displacements = {{
				0,
				(uint8_t*) this->data.data() - (uint8_t*) &(this->data_size)
			}};
			if (this->data.size() == 0) {
				displacements[1] = 0;
			}
			std::array<MPI_Datatype, 2> datatypes = {{MPI_UINT64_T, MPI_INT}};

			MPI_Datatype final_datatype;
			const int ret_val = MPI_Type_create_struct(
				2,
				&counts[0],
				&displacements[0],
				&datatypes[0],
				&final_datatype
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " MPI_Type_create_struct failed: " << dccrg::Error_String()(ret_val)
					<< std::endl;
				abort();
			}

			return std::make_tuple((void*) &(this->data_size), 1, final_datatype);

		} else if (Cell::transfer_data) {
			return std::make_tuple((void*) this->data.data(), int(this->data.size()), MPI_INT);
		} else {
			return std::make_tuple((void*) &(this->data_size), 1, MPI_UINT64_T);
		}
	}
};

bool
	Cell::transfer_all = true,
	Cell::transfer_data = false;


/*
Migrates cells off every third process.
*/
template <class Grid_T> void migrate_cells(Grid_T& grid)
{
	if (grid.get_rank() % 3 == 0) {
		if (grid.get_comm_size() > grid.get_rank() + 1) {
			const std::vector<uint64_t> cells = grid.get_cells();
			for (const auto& cell: cells) {
				grid.pin(cell, grid.get_rank() + 1);
			}
		}
	}

	grid.balance_load();
	grid.unpin_local_cells();
}


int main(int argc, char* argv[])
{
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &comm_size);

	/*
	Options
	*/
	bool restart = false;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("restart", "restart from a file");

	// read options from command line
	boost::program_options::variables_map option_variables;
	boost::program_options::store(
		boost::program_options::parse_command_line(argc, argv, options),
		option_variables
	);
	boost::program_options::notify(option_variables);

	// print a help message if asked
	if (option_variables.count("help") > 0) {
		if (rank == 0) {
			cout << options << endl;
		}
		MPI_Barrier(comm);
		return EXIT_SUCCESS;
	}

	if (option_variables.count("restart") > 0) {
		restart = true;
	}

	// intialize Zoltan
	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
		cerr << "Zoltan_Initialize failed" << endl;
		abort();
	}


	dccrg::Dccrg<Cell, dccrg::Cartesian_Geometry> grid;


	/*
	Common grid parameters
	*/
	const std::array<uint64_t, 3> grid_length = {{20, 1, 1}};
	const unsigned int neighborhood_size = 1;

	dccrg::Cartesian_Geometry::Parameters geom_params;
	geom_params.start[0] =
	geom_params.start[1] =
	geom_params.start[2] = 0;
	geom_params.level_0_cell_length[0] =
	geom_params.level_0_cell_length[1] =
	geom_params.level_0_cell_length[2] = 1.0;

	const std::string restart_name("variable_cell_data.dc");
	// write a restart file...
	if (not restart) {

		grid.initialize(grid_length, comm, "RANDOM", neighborhood_size);

		if (!grid.set_geometry(geom_params)) {
			cerr << "Couldn't set grid geometry" << endl;
			exit(EXIT_FAILURE);
		}

		migrate_cells(grid);

		// initialize
		const std::vector<uint64_t> cells = grid.get_cells();
		for(const auto& cell: cells) {
			Cell* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell: " << cell
					<< std::endl;
				abort();
			}

			if (cell % 4 != 0) {
				for (uint64_t i = 0; i < cell; i++) {
					cell_data->data.push_back(int(i));
				}
			}
			cell_data->data_size = cell_data->data.size();
		}

		std::tuple<void*, int, MPI_Datatype> header;
		std::get<0>(header) = &header;
		std::get<1>(header) = 0;
		std::get<2>(header) = MPI_INT;

		if (!grid.save_grid_data(restart_name, 0, header)) {
			std::cerr << "Process " << rank
				<< " Writing grid to file " << restart_name << " failed"
				<< std::endl;
			abort();
		}

	// ...or restart
	} else {

		std::tuple<void*, int, MPI_Datatype> header;
		std::get<0>(header) = &header;
		std::get<1>(header) = 0;
		std::get<2>(header) = MPI_INT;

		if (!grid.start_loading_grid_data(
			restart_name,
			0,
			header,
			comm,
			"RANDOM"
		)) {
			std::cerr << "Couldn't load initial data" << std::endl;
			abort();
		}

		// 1st load size of data in each cell
		Cell::transfer_all = false;
		Cell::transfer_data = false;

		if (!grid.continue_loading_grid_data()) {
			std::cerr << "Couldn't load initial data" << std::endl;
			abort();
		}

		// verify loaded data
		const std::vector<uint64_t> cells = grid.get_cells();
		for(const auto& cell: cells) {
			const Cell* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell: " << cell
					<< std::endl;
				abort();
			}

			if (cell % 4 != 0) {
				if (cell_data->data_size != cell) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< cell << ": " << cell_data->data_size
						<< std::endl;
					abort();
				}
			} else {
				if (cell_data->data_size != 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< cell << ": " << cell_data->data_size
						<< std::endl;
					abort();
				}
			}
		}

		// resize cell vectors to fit incoming data
		for(const auto& cell: cells) {
			Cell* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell: " << cell
					<< std::endl;
				abort();
			}

			cell_data->data.resize(cell_data->data_size);
		}

		// load rest of cell data
		Cell::transfer_data = true;

		if (!grid.continue_loading_grid_data()) {
			std::cerr << "Couldn't load initial data" << std::endl;
			abort();
		}

		// verify loaded data
		for(const auto& cell: cells) {
			const Cell* const cell_data = grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No data for cell: " << cell
					<< std::endl;
				abort();
			}

			if (cell_data->data.size() != cell_data->data_size) {
				std::cerr << __FILE__ << ":" << __LINE__ << ": "
					<< cell << ", " << cell_data->data_size << ", " << cell_data->data.size()
					<< " but should be " << ((cell % 4 != 0) ? cell : 0)
					<< std::endl;
				abort();
			}
		}

		grid.finish_loading_grid_data();
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

