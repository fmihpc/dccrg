/*
Tests the grid with some simple game of life patters, returns EXIT_SUCCESS if everything went ok.
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "zoltan.h"

#define DCCRG_ARBITRARY_STRETCH
#include "../../dccrg.hpp"


struct game_of_life_cell {

	// use boost::mpi for data transfers over MPI
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	bool is_alive;
	unsigned int live_neighbour_count;

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/)
	{
		ar & is_alive;
	}

	// use MPI directly for data transfers
	#else

	// data[0] == 1 if cell is alive, data[1] holds the number of live neighbours
	unsigned int data[2];

	void* at(void)
	{
		return this;
	}

	#ifdef DCCRG_USER_MPI_DATA_TYPE
	static MPI_Datatype mpi_datatype(void)
	{
		MPI_Datatype type;
		// processes don't need other processes' live neighbour info
		MPI_Type_contiguous(1, MPI_UNSIGNED, &type);
		return type;
	}
	#else
	static size_t size(void)
	{
		// processes don't need other processes' live neighbour info
		return sizeof(unsigned int);
	}
	#endif

	#endif	// ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
};


using namespace std;
using namespace boost::mpi;


/*!
Returns EXIT_SUCCESS if the state of the given game at given timestep is correct on this process, returns EXIT_FAILURE otherwise.
timestep == 0 means before any turns have been taken.
*/
int check_game_of_life_state(int timestep, dccrg<game_of_life_cell>* grid)
{
	vector<uint64_t> cells = grid->get_cells();
	for (auto cell = cells.cbegin(); cell != cells.cend(); cell++) {

		game_of_life_cell* data = (*grid)[*cell];
		if (data == NULL) {
			cerr << "No data for cell " << *cell << endl;
			return EXIT_FAILURE;
		}

		// check cells that are always supposed to be alive
		switch (*cell) {
		case 22:
		case 23:
		case 32:
		case 33:
		case 36:
		case 39:
		case 47:
		case 48:
		case 52:
		case 53:
		case 94:
		case 95:
		case 110:
		case 122:
		case 137:
		case 138:
		case 188:
		case 199:
		case 206:
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			if (!data->is_alive) {
			#else
			if (!data->data[0]) {
			#endif
				cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		// these are supposed to be alive every other turn
		if (timestep % 2 == 0) {

		switch (*cell) {
		case 109:
		case 123:
		case 189:
		case 190:
		case 198:
		case 200:
		case 204:
		case 205:
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			if (!data->is_alive) {
			#else
			if (!data->data[0]) {
			#endif
				cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		} else {

		switch (*cell) {
		case 174:
		case 184:
		case 214:
		case 220:
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			if (!data->is_alive) {
			#else
			if (!data->data[0]) {
			#endif
				cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
				return EXIT_FAILURE;
			}
			break;
		default:
			break;
		}

		}

		// check that the glider is moving correctly
		switch (timestep) {
		/* can't be bothered manually for cases 1-19, use an automatic method later */
		case 20:
			switch (*cell) {
			case 43:
			case 44:
			case 45:
			case 60:
			case 74:
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (!data->is_alive) {
				#else
				if (!data->data[0]) {
				#endif
					cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 21:
			switch (*cell) {
			case 29:
			case 44:
			case 45:
			case 58:
			case 60:
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (!data->is_alive) {
				#else
				if (!data->data[0]) {
				#endif
					cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 22:
			switch (*cell) {
			case 29:
			case 30:
			case 43:
			case 45:
			case 60:
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (!data->is_alive) {
				#else
				if (!data->data[0]) {
				#endif
					cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 23:
			switch (*cell) {
			case 29:
			case 30:
			case 45:
			case 59:
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (!data->is_alive) {
				#else
				if (!data->data[0]) {
				#endif
					cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;

		case 24:
			switch (*cell) {
			case 29:
			case 30:
			case 45:
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (!data->is_alive) {
				#else
				if (!data->data[0]) {
				#endif
					cerr << "Cell " << *cell << " isn't alive on timestep " << timestep << endl;
					return EXIT_FAILURE;
				}
				break;
			default:
				break;
			}
			break;
		default:
			break;
		}
	}

	return EXIT_SUCCESS;
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
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}


	#define STARTING_CORNER 0.0
	#define GRID_SIZE 15	// in unrefined cells
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
		y_coordinates.push_back(i * CELL_SIZE);
	}
	z_coordinates.push_back(0);
	z_coordinates.push_back(1);
	#define STENCIL_SIZE 1
	dccrg<game_of_life_cell> game_grid(comm, "RANDOM", x_coordinates, y_coordinates, z_coordinates, STENCIL_SIZE);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << game_grid.get_max_refinement_level() << endl;
		cout << "Number of cells: " << (x_coordinates.size() - 1) * (y_coordinates.size() - 1) * (z_coordinates.size() - 1) << endl << endl;
	}


	// create a blinker
	#define BLINKER_START 198
	uint64_t tmp1[] = {BLINKER_START, BLINKER_START + 1, BLINKER_START + 2};
	vector<uint64_t> blinker_cells(tmp1, tmp1 + sizeof(tmp1) / sizeof(uint64_t));
	for (auto cell = blinker_cells.cbegin(); cell != blinker_cells.cend(); cell++) {
		game_of_life_cell* cell_data = game_grid[*cell];
		if (cell_data == NULL) {
			continue;
		}
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->is_alive = true;
		#else
		cell_data->data[0] = 1;
		#endif
	}

	// create a toad
	#define TOAD_START 188
	uint64_t tmp2[] = {TOAD_START, TOAD_START + 1, TOAD_START + 2, TOAD_START + 1 + GRID_SIZE, TOAD_START + 2 + GRID_SIZE, TOAD_START + 3 + GRID_SIZE};
	vector<uint64_t> toad_cells(tmp2, tmp2 + sizeof(tmp2) / sizeof(uint64_t));
	for (auto cell = toad_cells.cbegin(); cell != toad_cells.cend(); cell++) {
		game_of_life_cell* cell_data = game_grid[*cell];
		if (cell_data == NULL) {
			continue;
		}
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->is_alive = true;
		#else
		cell_data->data[0] = 1;
		#endif
	}

	// create a beacon
	#define BEACON_START 137
	uint64_t tmp3[] = {BEACON_START, BEACON_START + 1, BEACON_START - GRID_SIZE, BEACON_START + 1 - GRID_SIZE, BEACON_START + 2 - 2 * GRID_SIZE, BEACON_START + 3 - 2 * GRID_SIZE, BEACON_START + 2 - 3 * GRID_SIZE, BEACON_START + 3 - 3 * GRID_SIZE};
	vector<uint64_t> beacon_cells(tmp3, tmp3 + sizeof(tmp3) / sizeof(uint64_t));
	for (auto cell = beacon_cells.cbegin(); cell != beacon_cells.cend(); cell++) {
		game_of_life_cell* cell_data = game_grid[*cell];
		if (cell_data == NULL) {
			continue;
		}
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->is_alive = true;
		#else
		cell_data->data[0] = 1;
		#endif
	}

	// create a glider
	#define GLIDER_START 143
	uint64_t tmp4[] = {GLIDER_START + 1, GLIDER_START + 2 - GRID_SIZE, GLIDER_START - 2 * GRID_SIZE, GLIDER_START + 1 - 2 * GRID_SIZE, GLIDER_START + 2 - 2 * GRID_SIZE};
	vector<uint64_t> glider_cells(tmp4, tmp4 + sizeof(tmp4) / sizeof(uint64_t));
	for (auto cell = glider_cells.cbegin(); cell != glider_cells.cend(); cell++) {
		game_of_life_cell* cell_data = game_grid[*cell];
		if (cell_data == NULL) {
			continue;
		}
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->is_alive = true;
		#else
		cell_data->data[0] = 1;
		#endif
	}

	// create a block
	#define BLOCK_START 47
	uint64_t tmp5[] = {BLOCK_START, BLOCK_START + 1, BLOCK_START - GRID_SIZE, BLOCK_START + 1 - GRID_SIZE};
	vector<uint64_t> block_cells(tmp5, tmp5 + sizeof(tmp5) / sizeof(uint64_t));
	for (auto cell = block_cells.cbegin(); cell != block_cells.cend(); cell++) {
		game_of_life_cell* cell_data = game_grid[*cell];
		if (cell_data == NULL) {
			continue;
		}
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->is_alive = true;
		#else
		cell_data->data[0] = 1;
		#endif
	}

	// create a beehive
	#define BEEHIVE_START 51
	uint64_t tmp6[] = {BEEHIVE_START - GRID_SIZE, BEEHIVE_START + 1, BEEHIVE_START + 2, BEEHIVE_START + 1 - 2 * GRID_SIZE, BEEHIVE_START + 2 - 2 * GRID_SIZE, BEEHIVE_START + 3 - GRID_SIZE};
	vector<uint64_t> beehive_cells(tmp6, tmp6 + sizeof(tmp6) / sizeof(uint64_t));
	for (auto cell = beehive_cells.cbegin(); cell != beehive_cells.cend(); cell++) {
		game_of_life_cell* cell_data = game_grid[*cell];
		if (cell_data == NULL) {
			continue;
		}
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->is_alive = true;
		#else
		cell_data->data[0] = 1;
		#endif
	}


	// every process outputs the game state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "game_of_life_test_" << comm.rank() << "_";
	ofstream outfile, visit_file;

	// visualize the game with visit -o game_of_life_test.visit
	if (comm.rank() == 0) {
		visit_file.open("game_of_life_test.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}

	#define TIME_STEPS 25
	if (comm.rank() == 0) {
		cout << "step: ";
	}
	for (int step = 0; step < TIME_STEPS; step++) {

		game_grid.balance_load();
		game_grid.start_remote_neighbour_data_update();
		game_grid.wait_neighbour_data_update();
		vector<uint64_t> cells = game_grid.get_cells();

		int result = check_game_of_life_state(step, &game_grid);
		if (GRID_SIZE != 15 || result != EXIT_SUCCESS) {
			cout << "Process " << comm.rank() << ": Game of Life test failed on timestep: " << step << endl;
			abort();
		}

		// the library writes the grid into a file in ascending cell order, do the same for the grid data at every time step
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
				visit_file << "game_of_life_test_" << process << "_" << step_string.str() << suffix.str() << endl;
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
		for (auto cell = cells.cbegin(); cell != cells.cend(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			if (cell_data->is_alive == true) {
			#else
			if (cell_data->data[0] == 1) {
			#endif
				outfile << "1";
			} else {
				outfile << "0";
			}
			outfile << endl;
		}

		// write each cells live neighbour count
		outfile << "SCALARS live_neighbour_count float 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (auto cell = cells.cbegin(); cell != cells.cend(); cell++) {
			game_of_life_cell* cell_data = game_grid[*cell];
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			outfile << cell_data->live_neighbour_count << endl;
			#else
			outfile << cell_data->data[1] << endl;
			#endif
		}

		// write each cells neighbour count
		outfile << "SCALARS neighbours int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			outfile << neighbours->size() << endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << comm.rank() << endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			outfile << *cell << endl;
		}
		outfile.close();


		// get the neighbour counts of every cell, starting with the cells whose neighbour data doesn't come from other processes
		vector<uint64_t> cells_with_local_neighbours = game_grid.get_cells_with_local_neighbours();
		for (auto cell = cells_with_local_neighbours.cbegin(); cell != cells_with_local_neighbours.cend(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->live_neighbour_count = 0;
			#else
			cell_data->data[1] = 0;
			#endif
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			if (neighbours == NULL) {
				cout << "Process " << comm.rank() << ": neighbour list for cell " << *cell << " not available" << endl;
				abort();
			}

			for (auto neighbour = neighbours->cbegin(); neighbour != neighbours->cend(); neighbour++) {

				if (*neighbour == 0) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data == NULL) {
					cout << "Process " << comm.rank() << ": neighbour " << *neighbour << " data for cell " << *cell << " not available" << endl;
					abort();
				}
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (neighbour_data->is_alive) {
					cell_data->live_neighbour_count++;
				}
				#else
				if (neighbour_data->data[0] == 1) {
					cell_data->data[1]++;
				}
				#endif
			}
		}

		vector<uint64_t> cells_with_remote_neighbour = game_grid.get_cells_with_remote_neighbour();
		for (auto cell = cells_with_remote_neighbour.cbegin(); cell != cells_with_remote_neighbour.cend(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->live_neighbour_count = 0;
			#else
			cell_data->data[1] = 0;
			#endif
			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			if (neighbours == NULL) {
				cout << "Process " << comm.rank() << ": neighbour list for cell " << *cell << " not available" << endl;
				abort();
			}

			for (auto neighbour = neighbours->cbegin(); neighbour != neighbours->cend(); neighbour++) {

				if (*neighbour == 0) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
				if (neighbour_data == NULL) {
					cout << "Process " << comm.rank() << ": neighbour " << *neighbour << " data for cell " << *cell << " not available" << endl;
					abort();
				}
				#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
				if (neighbour_data->is_alive) {
					cell_data->live_neighbour_count++;
				}
				#else
				if (neighbour_data->data[0] == 1) {
					cell_data->data[1]++;
				}
				#endif
			}
		}

		// calculate the next turn
		for (auto cell = cells_with_local_neighbours.cbegin(); cell != cells_with_local_neighbours.cend(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			if (cell_data->live_neighbour_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbour_count != 2) {
				cell_data->is_alive = false;
			}
			#else
			if (cell_data->data[1] == 3) {
				cell_data->data[0] = 1;
			} else if (cell_data->data[1] != 2) {
				cell_data->data[0] = 0;
			}
			#endif
		}
		for (auto cell = cells_with_remote_neighbour.cbegin(); cell != cells_with_remote_neighbour.cend(); cell++) {

			game_of_life_cell* cell_data = game_grid[*cell];

			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			if (cell_data->live_neighbour_count == 3) {
				cell_data->is_alive = true;
			} else if (cell_data->live_neighbour_count != 2) {
				cell_data->is_alive = false;
			}
			#else
			if (cell_data->data[1] == 3) {
				cell_data->data[0] = 1;
			} else if (cell_data->data[1] != 2) {
				cell_data->data[0] = 0;
			}
			#endif
		}

	}

	if (comm.rank() == 0) {
		cout << "\nPassed" << endl;
		visit_file.close();
	}

	return EXIT_SUCCESS;
}

