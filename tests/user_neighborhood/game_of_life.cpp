/*
Tests the grid with some simple game of life patters in 2d using a general neighborhood.
*/

#include "algorithm"
#include "boost/assign/list_of.hpp"
#include "boost/foreach.hpp"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "zoltan.h"

#include "../../dccrg.hpp"

struct game_of_life_cell {
	unsigned int data[2];
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int)
	{
		ar & data;
	}
	#else
	void* at()
	{
		return this;
	}
	const void* at() const
	{
		return this;
	}
	#ifdef DCCRG_USER_MPI_DATA_TYPE
	static MPI_Datatype mpi_datatype()
	{
		MPI_Datatype type;
		MPI_Type_contiguous(1, MPI_UNSIGNED, &type);
		return type;
	}
	#else
	static size_t size()
	{
		return sizeof(unsigned int);
	}
	#endif
	#endif	// ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
};


using namespace boost::mpi;
using namespace dccrg;
using namespace std;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    abort();
	}

	Dccrg<game_of_life_cell> game_grid;

	const uint64_t grid_size = 18;
	const double cell_size = 1.0 / grid_size;
	game_grid.set_geometry(
		grid_size, grid_size, 1,
		0, 0, 0,
		cell_size, cell_size, cell_size
	);

	game_grid.initialize(
		comm,
		"RANDOM",
		2,
		0,
		true, true, true
	);

	/*
	Use a neighborhood like this in the z plane:
	O.O.O
	.....
	O.X.O
	.....
	O.O.O
	and play 4 interlaced games simultaneously.
	*/
	typedef dccrg::Types<3>::neighborhood_item_t neigh_t;
	const std::vector<neigh_t> neighborhood
		= boost::assign::list_of<neigh_t>
			(boost::assign::list_of(-2)(-2)(0))
			(boost::assign::list_of(-2)( 0)(0))
			(boost::assign::list_of(-2)( 2)(0))
			(boost::assign::list_of( 0)(-2)(0))
			(boost::assign::list_of( 0)( 2)(0))
			(boost::assign::list_of( 2)(-2)(0))
			(boost::assign::list_of( 2)( 0)(0))
			(boost::assign::list_of( 2)( 2)(0));

	const int neighborhood_id = 1;
	if (!game_grid.add_remote_update_neighborhood(neighborhood_id, neighborhood)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set neighborhood"
			<< std::endl;
		abort();
	}

	/*
	Patterns at even x and even y cells starting from origin
	*/
	// create a Toad (http://www.conwaylife.com/wiki/Toad)
	const uint64_t toad_start = 12 * grid_size + 6;
	const vector<uint64_t> toad_cells
		= boost::assign::list_of
			(toad_start)
			(toad_start + 2)
			(toad_start + 4)
			(toad_start + 2 + 2 * grid_size)
			(toad_start + 4 + 2 * grid_size)
			(toad_start + 6 + 2 * grid_size);
	BOOST_FOREACH(const uint64_t cell, toad_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->data[0] = 1;
	}

	/*
	Patterns at even x and odd y cells
	*/
	// http://www.conwaylife.com/wiki/Blinker
	const uint64_t blinker_start = 14 * grid_size + 1;
	const vector<uint64_t> blinker_cells
		= boost::assign::list_of
			(blinker_start)
			(blinker_start + 2)
			(blinker_start + 4);
	BOOST_FOREACH(const uint64_t cell, blinker_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->data[0] = 1;
	}

	// http://www.conwaylife.com/wiki/Beehive
	const uint64_t beehive_start = 10 * grid_size + 7;
	const vector<uint64_t> beehive_cells
		= boost::assign::list_of
			(beehive_start + 2)
			(beehive_start + 4)
			(beehive_start - 2 * grid_size)
			(beehive_start - 2 * grid_size + 6)
			(beehive_start - 4 * grid_size + 2)
			(beehive_start - 4 * grid_size + 4);
	BOOST_FOREACH(const uint64_t cell, beehive_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->data[0] = 1;
	}

	/*
	Odd x, even y
	*/
	// http://www.conwaylife.com/wiki/Beacon
	const uint64_t beacon_start = 15 * grid_size + 6;
	const vector<uint64_t> beacon_cells
		= boost::assign::list_of
			(beacon_start)
			(beacon_start + 2)
			(beacon_start - 2 * grid_size)
			(beacon_start + 6 - 4 * grid_size)
			(beacon_start + 4 - 6 * grid_size)
			(beacon_start + 6 - 6 * grid_size);
	BOOST_FOREACH(const uint64_t cell, beacon_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->data[0] = 1;
	}

	// http://www.conwaylife.com/wiki/Block
	const uint64_t block_start = 5 * grid_size + 2;
	const vector<uint64_t> block_cells
		= boost::assign::list_of
			(block_start)
			(block_start + 2)
			(block_start - 2 * grid_size)
			(block_start + 2 - 2 * grid_size);
	BOOST_FOREACH(const uint64_t cell, block_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->data[0] = 1;
	}

	/*
	Odd x, odd y
	*/
	// http://www.conwaylife.com/wiki/Glider
	const uint64_t glider_start = 15 * grid_size + 1;
	const vector<uint64_t> glider_cells
		= boost::assign::list_of
			(glider_start + 2)
			(glider_start + 4 - 2 * grid_size)
			(glider_start - 4 * grid_size)
			(glider_start + 2 - 4 * grid_size)
			(glider_start + 4 - 4 * grid_size);
	BOOST_FOREACH(const uint64_t cell, glider_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data == NULL) {
			continue;
		}
		cell_data->data[0] = 1;
	}


	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "general_neighborhood_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("general_neighborhood.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	#define TIME_STEPS 50
	if (comm.rank() == 0) {
		cout << "step: ";
	}
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.update_remote_neighbor_data(neighborhood_id);
		vector<uint64_t> cells = game_grid.get_cells();

		sort(cells.begin(), cells.end());

		if (comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		// write the game state into a file named according to the current time step
		string current_output_name("");
		ostringstream step_string;
		step_string.fill('0');
		step_string.width(5);
		step_string << step;
		current_output_name += basename.str();
		current_output_name += step_string.str();
		current_output_name += suffix.str();

		// visualize the game with visit -o game_of_life_test.visit
		if (comm.rank() == 0) {
			for (int process = 0; process < comm.size(); process++) {
				visit_file << "general_neighborhood_" << process << "_" << step_string.str() << suffix.str() << endl;
			}
		}

		// write the grid into a file
		game_grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// go through the grids cells and write their state into the file
		outfile << "SCALARS is_alive float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		BOOST_FOREACH(const uint64_t cell, cells) {
			game_of_life_cell* cell_data = game_grid[cell];

			if (cell_data->data[0] == 1) {
				// color live cells of interlaced games differently
				const Types<3>::indices_t indices = game_grid.get_indices(cell);
				outfile << 1 + indices[0] % 2 + (indices[1] % 2) * 2;
			} else {
				outfile << 0;
			}
			outfile << endl;
		}

		// write each cells live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		BOOST_FOREACH(const uint64_t cell, cells) {
			game_of_life_cell* cell_data = game_grid[cell];
			outfile << cell_data->data[1] << endl;
		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		BOOST_FOREACH(const uint64_t cell, cells) {
			const vector<uint64_t>* neighbors = game_grid.get_neighbors(cell);
			outfile << neighbors->size() << endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (size_t i = 0; i < cells.size(); i++) {
			outfile << comm.rank() << endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		BOOST_FOREACH(const uint64_t cell, cells) {
			outfile << cell << endl;
		}
		outfile.close();

		BOOST_FOREACH(const uint64_t cell, cells) {
			game_of_life_cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << endl;
				abort();
			}

			cell_data->data[1] = 0;
			const vector<uint64_t>* neighbors = game_grid.get_neighbors(cell, neighborhood_id);
			if (neighbors == NULL) {
				cerr << "Process " << comm.rank() << ": No neighbor list for cell " << cell << endl;
				abort();
			}

			BOOST_FOREACH(const uint64_t neighbor, *neighbors) {
				if (neighbor == 0) {
					continue;
				}

				game_of_life_cell* neighbor_data = game_grid[neighbor];
				if (neighbor_data == NULL) {
					cout << "Process " << comm.rank()
						<< ": neighbor " << neighbor
						<< " data for cell " << cell
						<< " not available"
						<< endl;
					abort();
				}
				if (neighbor_data->data[0] == 1) {
					cell_data->data[1]++;
				}
			}
		}

		// calculate the next turn
		BOOST_FOREACH(const uint64_t cell, cells) {
			game_of_life_cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << endl;
				abort();
			}

			if (cell_data->data[1] == 3) {
				cell_data->data[0] = 1;
			} else if (cell_data->data[1] != 2) {
				cell_data->data[0] = 0;
			}
		}
	}

	if (comm.rank() == 0) {
		cout << endl;
	}

	return EXIT_SUCCESS;
}

