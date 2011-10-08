/*
Tests the scalability of the grid in 2 D
*/

#include "algorithm"
#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "fstream"
#include "iostream"
#include "unistd.h"
#include "zoltan.h"

#include "../../dccrg_arbitrary_geometry.hpp"
#include "../../dccrg.hpp"

// TODO: move this to a separate file
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
	MPI_Datatype mpi_datatype(void)
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
using namespace dccrg;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	time_t before, after, total = 0;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	Dccrg<game_of_life_cell, ArbitraryGeometry> game_grid;

	#define GRID_SIZE 1000	// in unrefined cells
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
		y_coordinates.push_back(i * CELL_SIZE);
	}
	z_coordinates.push_back(0);
	z_coordinates.push_back(1);
	game_grid.set_geometry(x_coordinates, y_coordinates, z_coordinates);

	#define NEIGHBORHOOD_SIZE 1
	#define MAX_REFINEMENT_LEVEL 0
	game_grid.initialize(comm, "RCB", NEIGHBORHOOD_SIZE, MAX_REFINEMENT_LEVEL);
	if (comm.rank() == 0) {
		cout << "Maximum refinement level of the grid: " << game_grid.get_maximum_refinement_level() << endl;
		cout << "Number of cells: " << (x_coordinates.size() - 1) * (y_coordinates.size() - 1) * (z_coordinates.size() - 1) << endl << endl;
	}

	game_grid.balance_load();
	comm.barrier();

	vector<uint64_t> cells_with_local_neighbours = game_grid.get_cells_with_local_neighbours();
	vector<uint64_t> cells_with_remote_neighbour = game_grid.get_cells_with_remote_neighbour();
	cout << "Process " << comm.rank()
		<< ": number of cells with local neighbours: " << cells_with_local_neighbours.size()
		<< ", number of cells with a remote neighbour: " << cells_with_remote_neighbour.size()
		<< endl;

	// initialize the game with a line of living cells in the x direction in the middle
	for (vector<uint64_t>::const_iterator
		cell = cells_with_local_neighbours.begin();
		cell != cells_with_local_neighbours.end();
		cell++
	) {
		game_of_life_cell* cell_data = game_grid[*cell];
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->live_neighbour_count = 0;
		#else
		cell_data->data[1] = 0;
		#endif

		double y = game_grid.get_cell_y(*cell);
		if (fabs(0.5 + 0.1 * game_grid.get_cell_y_size(*cell) - y) < 0.5 * game_grid.get_cell_y_size(*cell)) {
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->is_alive = true;
			#else
			cell_data->data[0] = 1;
			#endif
		} else {
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->is_alive = false;
			#else
			cell_data->data[0] = 0;
			#endif
		}
	}
	for (vector<uint64_t>::const_iterator
		cell = cells_with_remote_neighbour.begin();
		cell != cells_with_remote_neighbour.end();
		cell++
	) {
		game_of_life_cell* cell_data = game_grid[*cell];
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		cell_data->live_neighbour_count = 0;
		#else
		cell_data->data[1] = 0;
		#endif

		double y = game_grid.get_cell_y(*cell);
		if (fabs(0.5 + 0.1 * game_grid.get_cell_y_size(*cell) - y) < 0.5 * game_grid.get_cell_y_size(*cell)) {
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->is_alive = true;
			#else
			cell_data->data[0] = 1;
			#endif
		} else {
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->is_alive = false;
			#else
			cell_data->data[0] = 0;
			#endif
		}
	}

	if (comm.rank() == 0) {
		cout << "step: ";
	}

	#define TIME_STEPS 100
	before = time(NULL);
	for (int step = 0; step < TIME_STEPS; step++) {

		if (comm.rank() == 0) {
			cout << step << " ";
			cout.flush();
		}

		game_grid.start_remote_neighbour_data_update();
		/*
		Get the neighbour counts of every cell, starting with the cells whose neighbour data
		doesn't come from other processes
		*/
		for (vector<uint64_t>::const_iterator
			cell = cells_with_local_neighbours.begin();
			cell != cells_with_local_neighbours.end();
			cell++
		) {
			game_of_life_cell* cell_data = game_grid[*cell];
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->live_neighbour_count = 0;
			#else
			cell_data->data[1] = 0;
			#endif

			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			for (vector<uint64_t>::const_iterator
				neighbour = neighbours->begin();
				neighbour != neighbours->end();
				neighbour++
			) {
				if (*neighbour == 0) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
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

		// wait for neighbour data updates to this process to finish and go through the rest of the cells
		game_grid.wait_neighbour_data_update_receives();
		for (vector<uint64_t>::const_iterator
			cell = cells_with_remote_neighbour.begin();
			cell != cells_with_remote_neighbour.end();
			cell++
		) {
			game_of_life_cell* cell_data = game_grid[*cell];
			#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
			cell_data->live_neighbour_count = 0;
			#else
			cell_data->data[1] = 0;
			#endif

			const vector<uint64_t>* neighbours = game_grid.get_neighbours(*cell);
			for (vector<uint64_t>::const_iterator
				neighbour = neighbours->begin();
				neighbour != neighbours->end();
				neighbour++
			) {
				if (*neighbour == 0) {
					continue;
				}

				game_of_life_cell* neighbour_data = game_grid[*neighbour];
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
		/*
		Wait for neighbour data updates from this process to finish until
		updating live status of own cells
		*/
		game_grid.wait_neighbour_data_update_sends();

		// calculate the next turn
		for (vector<uint64_t>::const_iterator
			cell = cells_with_local_neighbours.begin();
			cell != cells_with_local_neighbours.end();
			cell++
		) {
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
		for (vector<uint64_t>::const_iterator
			cell = cells_with_remote_neighbour.begin();
			cell != cells_with_remote_neighbour.end();
			cell++
		) {
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
	after = time(NULL);
	total += after - before;
	if (comm.rank() == 0) {
		cout << endl;
	}
	comm.barrier();

	int number_of_cells = cells_with_local_neighbours.size() + cells_with_remote_neighbour.size();
	cout << "Process " << comm.rank()
		<< ": " << number_of_cells * TIME_STEPS << " cells processed at the speed of "
		<< double(number_of_cells * TIME_STEPS) / total << " cells / second"
		<< endl;

	return EXIT_SUCCESS;
}
