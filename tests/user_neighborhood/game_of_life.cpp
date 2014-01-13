/*
Tests the grid with some simple game of life patters in 2d using a general neighborhood.
*/

#include "algorithm"
#include "array"
#include "boost/mpi.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "unordered_set"
#include "zoltan.h"

#include "../../dccrg.hpp"

struct game_of_life_cell {
	unsigned int data[2];

	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int)
	{
		ar & data;
	}

	#else

	std::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype() const
	{
		return std::make_tuple((void*) &(this->data[0]), 1, MPI_UNSIGNED);
	}

	#endif
};


using namespace boost::mpi;
using namespace dccrg;
using namespace std;


/*!
The x and y coordinates of moving patterns wander out
of the grid, this puts them back to the proper location.
*/
void wrap_coordinates(
	const int grid_size,
	int& x,
	int& y
) {
	x %= grid_size;
	while (y < 0) {
		y += grid_size;
	}
	y %= grid_size;
}


/*!
Returns cells which are alive at given time step.

Given time step must be divisible by 4.
*/
std::unordered_set<uint64_t> get_live_cells(
	const uint64_t grid_size,
	const int time_step
) {
	if (time_step % 4 > 0) {
		cerr << "Time step must be divisible by 4" << endl;
		abort();
	}

	std::unordered_set<uint64_t> result;

	/*
	The game is interlaced so "neighboring" cells have one
	cell between them, hence 2 * n * grid_size, where n = 1, 2, ...
	*/

	/*
	Patterns at even x and even y cells starting from origin
	*/
	// create a Toad (http://www.conwaylife.com/wiki/Toad)
	const uint64_t toad_start = 12 * grid_size + 6;
	result.insert(toad_start);
	result.insert(toad_start + 2);
	result.insert(toad_start + 4);
	result.insert(toad_start + 2 + 2 * grid_size);
	result.insert(toad_start + 4 + 2 * grid_size);
	result.insert(toad_start + 6 + 2 * grid_size);

	/*
	Patterns at even x and odd y cells
	*/
	// http://www.conwaylife.com/wiki/Blinker
	const uint64_t blinker_start = 14 * grid_size + 1;
	result.insert(blinker_start);
	result.insert(blinker_start + 2);
	result.insert(blinker_start + 4);

	// http://www.conwaylife.com/wiki/Beehive
	const uint64_t beehive_start = 10 * grid_size + 7;
	result.insert(beehive_start + 2);
	result.insert(beehive_start + 4);
	result.insert(beehive_start - 2 * grid_size);
	result.insert(beehive_start - 2 * grid_size + 6);
	result.insert(beehive_start - 4 * grid_size + 2);
	result.insert(beehive_start - 4 * grid_size + 4);

	/*
	Odd x, even y
	*/
	// http://www.conwaylife.com/wiki/Beacon
	const uint64_t beacon_start = 15 * grid_size + 6;
	result.insert(beacon_start);
	result.insert(beacon_start + 2);
	result.insert(beacon_start - 2 * grid_size);
	result.insert(beacon_start + 6 - 4 * grid_size);
	result.insert(beacon_start + 4 - 6 * grid_size);
	result.insert(beacon_start + 6 - 6 * grid_size);

	// http://www.conwaylife.com/wiki/Block
	const uint64_t block_start = 5 * grid_size + 2;
	result.insert(block_start);
	result.insert(block_start + 2);
	result.insert(block_start - 2 * grid_size);
	result.insert(block_start + 2 - 2 * grid_size);

	/*
	Odd x, odd y
	*/
	// http://www.conwaylife.com/wiki/Glider
	std::vector<std::pair<int, int> > glider_coordinates{
		{3 + 2 * (time_step / 4), 15 - 2 * (time_step / 4)},
		{5 + 2 * (time_step / 4), 13 - 2 * (time_step / 4)},
		{1 + 2 * (time_step / 4), 11 - 2 * (time_step / 4)},
		{3 + 2 * (time_step / 4), 11 - 2 * (time_step / 4)},
		{5 + 2 * (time_step / 4), 11 - 2 * (time_step / 4)}
	};
	for (size_t i = 0; i < glider_coordinates.size(); i++) {
		int
			&x = glider_coordinates[i].first,
			&y = glider_coordinates[i].second;
		wrap_coordinates((const int) grid_size, x, y);
		result.insert((uint64_t) y * grid_size + (uint64_t) x);
	}

	return result;
}


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

	const std::array<uint64_t, 3> grid_length = {{18, 18, 1}};
	game_grid.initialize(grid_length, comm, "RANDOM", 2, 0, true, true, true);

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
	const std::vector<neigh_t> neighborhood{
		{-2, -2, 0},
		{-2,  0, 0},
		{-2,  2, 0},
		{ 0, -2, 0},
		{ 0,  2, 0},
		{ 2, -2, 0},
		{ 2,  0, 0},
		{ 2,  2, 0}
	};

	const int neighborhood_id = 1;
	if (!game_grid.add_neighborhood(neighborhood_id, neighborhood)) {
		std::cerr << __FILE__ << ":" << __LINE__
			<< " Couldn't set neighborhood"
			<< std::endl;
		abort();
	}

	// initial condition
	const std::unordered_set<uint64_t> live_cells = get_live_cells(grid_length[0], 0);
	for (const uint64_t cell: live_cells) {
		game_of_life_cell* cell_data = game_grid[cell];
		if (cell_data != NULL) {
			cell_data->data[0] = 1;
		}
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

	#define TIME_STEPS 36
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.update_copies_of_remote_neighbors(neighborhood_id);
		vector<uint64_t> cells = game_grid.get_cells();

		// check that the result is correct
		if (step % 4 == 0) {
			const std::unordered_set<uint64_t> live_cells = get_live_cells(grid_length[0], step);
			for (const uint64_t cell: cells) {
				game_of_life_cell* cell_data = game_grid[cell];
				if (cell_data == NULL) {
					std::cerr << __FILE__ << ":" << __LINE__ << endl;
					abort();
				}

				if (cell_data->data[0] == 0) {
					if (live_cells.count(cell) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Cell " << cell
							<< " shouldn't be alive on step " << step
							<< endl;
						abort();
					}
				} else {
					if (live_cells.count(cell) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Cell " << cell
							<< " should be alive on step " << step
							<< endl;
						abort();
					}
				}
			}
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
		sort(cells.begin(), cells.end());
		game_grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// go through the grids cells and write their state into the file
		outfile << "SCALARS is_alive float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const uint64_t cell: cells) {
			game_of_life_cell* cell_data = game_grid[cell];

			if (cell_data->data[0] == 1) {
				// color live cells of interlaced games differently
				const Types<3>::indices_t indices = game_grid.mapping.get_indices(cell);
				outfile << 1 + indices[0] % 2 + (indices[1] % 2) * 2;
			} else {
				outfile << 0;
			}
			outfile << endl;
		}

		// write each cells live neighbor count
		outfile << "SCALARS live_neighbor_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const uint64_t cell: cells) {
			game_of_life_cell* cell_data = game_grid[cell];
			outfile << cell_data->data[1] << endl;
		}

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const uint64_t cell: cells) {
			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(cell);
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
		for (const uint64_t cell: cells) {
			outfile << cell << endl;
		}
		outfile.close();

		for (const uint64_t cell: cells) {
			game_of_life_cell* cell_data = game_grid[cell];
			if (cell_data == NULL) {
				std::cerr << __FILE__ << ":" << __LINE__ << endl;
				abort();
			}

			cell_data->data[1] = 0;
			const vector<uint64_t>* neighbors = game_grid.get_neighbors_of(cell, neighborhood_id);
			if (neighbors == NULL) {
				cerr << "Process " << comm.rank() << ": No neighbor list for cell " << cell << endl;
				abort();
			}

			for (const uint64_t neighbor: *neighbors) {
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
		for (const uint64_t cell: cells) {
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
		cout << "PASSED" << endl;
	}

	return EXIT_SUCCESS;
}

