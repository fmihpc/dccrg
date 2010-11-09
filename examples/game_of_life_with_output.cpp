/*
A simple 2 D game of life program to demonstrate the efficient usage of dccrg and shows an example of how to output dccrg grid data into a file
*/

#include "boost/mpi.hpp"
#include "cstdio"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "iostream"
#include "mpi.h"
#include "stdint.h"
#include "zoltan.h"

#include "../dccrg.hpp"

using namespace std;
using namespace boost::mpi;


// store in every cell of the grid whether the cell is alive and the number of live neighbours it has
struct game_of_life_cell {

	// boost requires this from user data
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
		/* live_neighbour_count from neighbouring cells is not used
		ar & live_neighbour_count;*/
	}

	bool is_alive;
	unsigned int live_neighbour_count;
};


/*!
Initializes the given cells, all of which must be local
*/
void initialize_game(const vector<uint64_t>* cells, dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];
		cell_data->live_neighbour_count = 0;

		if (double(rand()) / RAND_MAX < 0.2) {
			cell_data->is_alive = true;
		} else {
			cell_data->is_alive = false;
		}
	}
}


/*!
Calculates the number of live neihgbours for every cell given, all of which must be local
*/
void get_live_neighbour_counts(const vector<uint64_t>* cells, dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];

		cell_data->live_neighbour_count = 0;
		const vector<uint64_t>* neighbours = game_grid->get_neighbours(*cell);

		for (vector<uint64_t>::const_iterator neighbour = neighbours->begin(); neighbour != neighbours->end(); neighbour++) {
			game_of_life_cell* neighbour_data = (*game_grid)[*neighbour];
			if (neighbour_data->is_alive) {
				cell_data->live_neighbour_count++;
			}
		}
	}
}


/*!
Applies the game of life rules to every given cell, all of which must be local
*/
void apply_rules(const vector<uint64_t>* cells, dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];

		if (cell_data->live_neighbour_count == 3) {
			cell_data->is_alive = true;
		} else if (cell_data->live_neighbour_count != 2) {
			cell_data->is_alive = false;
		}
	}
}


/*!
Writes the game state into a file named game_of_life_, postfixed with the timestep and .dc
Fileformat:
double x_start
double y_start
double z_start
double cell_size
uint64_t x_length
uint64_t y_length
uint64_t z_length
uint64_t cell1
uint64_t is_alive1
uint64_t cell2
uint64_t is_alive2
uint64_t cell3
uint64_t is_alive3
...
*/
bool write_game_data(const int step, communicator comm, dccrg<game_of_life_cell>* game_grid)
{
	int result;

	// get the output filename
	ostringstream basename("game_of_life_"), step_string, suffix(".dc");
	step_string.width(3);
	step_string.fill('0');
	step_string << step;

	string output_name("");
	output_name += basename.str();
	output_name += step_string.str();
	output_name += suffix.str();

	// MPI_File_open wants a non-constant string
	char* output_name_c_string = new char [output_name.size() + 1];
	strncpy(output_name_c_string, output_name.c_str(), output_name.size() + 1);

	/*
	Contrary to what http://www.open-mpi.org/doc/v1.4/man3/MPI_File_open.3.php writes, MPI_File_open doesn't truncate the file with OpenMPI 1.4.1 on Ubuntu, so use a fopen call first (http://www.opengroup.org/onlinepubs/009695399/functions/fopen.html)
	*/
	if (comm.rank() == 0) {
		FILE* i = fopen(output_name_c_string, "w");
		fflush(i);
		fclose(i);
	}
	comm.barrier();

	MPI_File outfile;
	result = MPI_File_open(comm, output_name_c_string, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfile);
	if (result != MPI_SUCCESS) {
		char mpi_error_string[MPI_MAX_ERROR_STRING + 1];
		int mpi_error_string_length;
		MPI_Error_string(result, mpi_error_string, &mpi_error_string_length);
		mpi_error_string[mpi_error_string_length + 1] = '\0';
		cerr << "Couldn't open file " << output_name_c_string << ": " << mpi_error_string << endl;
		return false;
	}

	MPI_Status status;
	MPI_Offset offset;

	// first process writes the header
	MPI_File_set_view(outfile, 0, MPI_DOUBLE, MPI_DOUBLE, (char*)"native", MPI_INFO_NULL);
	if (comm.rank() == 0) {
		// for file writes the offset is in units of sizeof(MPI_DOUBLE)
		offset = 0;
		double x_start = game_grid->get_x_start();
		MPI_File_write_at(outfile, offset, &x_start, 1, MPI_DOUBLE, &status);

		offset = 1;
		double y_start = game_grid->get_y_start();
		MPI_File_write_at(outfile, offset, &y_start, 1, MPI_DOUBLE, &status);

		offset = 2;
		double z_start = game_grid->get_z_start();
		MPI_File_write_at(outfile, offset, &z_start, 1, MPI_DOUBLE, &status);

		offset = 3;
		double cell_size = game_grid->get_cell_x_size(1);
		MPI_File_write_at(outfile, offset, &cell_size, 1, MPI_DOUBLE, &status);
	}

	// for file views the offset is in bytes
	MPI_File_set_view(outfile, 4 * sizeof(double), MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, (char*)"native", MPI_INFO_NULL);
	if (comm.rank() == 0) {
		// for file writes the offset is in units of sizeof(uint64_t), starting from the offset given in set_view
		offset = 0;
		uint64_t x_length = game_grid->get_x_length();
		MPI_File_write_at(outfile, offset, &x_length, 1, MPI_UNSIGNED_LONG_LONG, &status);

		offset = 1;
		uint64_t y_length = game_grid->get_y_length();
		MPI_File_write_at(outfile, offset, &y_length, 1, MPI_UNSIGNED_LONG_LONG, &status);

		offset = 2;
		uint64_t z_length = game_grid->get_z_length();
		MPI_File_write_at(outfile, offset, &z_length, 1, MPI_UNSIGNED_LONG_LONG, &status);
	}

	// figure out how many bytes every process will write
	vector<uint64_t> cells = game_grid->get_cells();
	vector<uint64_t> all_bytes;
	all_gather(comm, 2 * sizeof(uint64_t) * cells.size(), all_bytes);

	// add the header to every process' offset
	offset = 4 * sizeof(double) + 3 * sizeof(uint64_t);
	for (int i = 0; i < comm.rank(); i++) {
		offset += all_bytes[i];
	}

	// every process writes its cells after the previous processes
	MPI_File_set_view(outfile, offset, MPI_UNSIGNED_LONG_LONG, MPI_UNSIGNED_LONG_LONG, (char*)"native", MPI_INFO_NULL);
	offset = 0;
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		uint64_t current_cell = *cell;
		MPI_File_write_at(outfile, offset, &current_cell, 1, MPI_UNSIGNED_LONG_LONG, &status);
		offset++;

		uint64_t current_is_alive = uint64_t((*game_grid)[*cell]->is_alive);
		MPI_File_write_at(outfile, offset, &current_is_alive, 1, MPI_UNSIGNED_LONG_LONG, &status);
		offset++;
	}

	MPI_File_close(&outfile);
	return true;
}


int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (comm.rank() == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}


	// create the grid
	#define GRID_X_SIZE 10	// in unrefined cells
	#define GRID_Y_SIZE 10
	#define GRID_Z_SIZE 1
	#define CELL_SIZE 1.0
	#define STENCIL_SIZE 1	// the cells that share a vertex are considered neighbours
	#define MAX_REFINEMENT_LEVEL 0
	dccrg<game_of_life_cell> game_grid(comm, "RCB", 0, 0, 0, CELL_SIZE, CELL_SIZE, CELL_SIZE, GRID_X_SIZE, GRID_Y_SIZE, GRID_Z_SIZE, STENCIL_SIZE, MAX_REFINEMENT_LEVEL);	// use the recursive coordinate bisection method for load balancing (http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg_rcb.html)

	// since the grid doesn't change (isn't refined / unrefined) during the game, workload can be balanced just once in the beginning
	game_grid.balance_load();

	/*
	Get the cells on this process just once, since the grid doesn't change during the game
	To make the game scale better, separate local cells into those without even one neighbour on another process and those that do.
	While updating cell data between processes, start calculating the next turn for cells which don't have neighbours on other processes
	*/
	vector<uint64_t> cells_with_local_neighbours = game_grid.get_cells_with_local_neighbours();
	vector<uint64_t> cells_with_remote_neighbour = game_grid.get_cells_with_remote_neighbour();

	initialize_game(&cells_with_local_neighbours, &game_grid);
	initialize_game(&cells_with_remote_neighbour, &game_grid);

	#define TURNS 10
	for (int turn = 0; turn < TURNS; turn++) {

		write_game_data(turn, comm, &game_grid);

		// start updating cell data from other processes and calculate the next turn for cells without neighbours on other processes in the meantime
		game_grid.start_remote_neighbour_data_update();
		get_live_neighbour_counts(&cells_with_local_neighbours, &game_grid);

		// wait for neighbour data updates to finish and the calculate the next turn for rest of the cells on this process
		game_grid.wait_neighbour_data_update();
		get_live_neighbour_counts(&cells_with_remote_neighbour, &game_grid);

		// update the state of life for all local cells
		apply_rules(&cells_with_local_neighbours, &game_grid);
		apply_rules(&cells_with_remote_neighbour, &game_grid);
	}
	write_game_data(TURNS, comm, &game_grid);

	return EXIT_SUCCESS;
}
