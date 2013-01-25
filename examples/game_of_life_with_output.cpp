/*
A simple 2 D game of life program to demonstrate the efficient usage of dccrg and shows an example of how to output dccrg grid data into a file
*/

#include "boost/mpi.hpp"
#include "cstddef"
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
using namespace dccrg;

// store in every cell of the grid whether the cell is alive and the number of live neighbors it has
struct game_of_life_cell {

	// boost requires this from user data
	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & is_alive;
		/* live_neighbor_count from neighboring cells is not used
		ar & live_neighbor_count;*/
	}

	bool is_alive;
	unsigned int live_neighbor_count;
};


/*!
Initializes the given cells, all of which must be local
*/
void initialize_game(const vector<uint64_t>* cells, Dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];
		cell_data->live_neighbor_count = 0;

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
void get_live_neighbor_counts(const vector<uint64_t>* cells, Dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];

		cell_data->live_neighbor_count = 0;
		const vector<uint64_t>* neighbors = game_grid->get_neighbors(*cell);

		for (vector<uint64_t>::const_iterator neighbor = neighbors->begin(); neighbor != neighbors->end(); neighbor++) {
			if (*neighbor == 0) {
				continue;
			}

			game_of_life_cell* neighbor_data = (*game_grid)[*neighbor];
			if (neighbor_data->is_alive) {
				cell_data->live_neighbor_count++;
			}
		}
	}
}


/*!
Applies the game of life rules to every given cell, all of which must be local
*/
void apply_rules(const vector<uint64_t>* cells, Dccrg<game_of_life_cell>* game_grid)
{
	for (vector<uint64_t>::const_iterator cell = cells->begin(); cell != cells->end(); cell++) {

		game_of_life_cell* cell_data = (*game_grid)[*cell];

		if (cell_data->live_neighbor_count == 3) {
			cell_data->is_alive = true;
		} else if (cell_data->live_neighbor_count != 2) {
			cell_data->is_alive = false;
		}
	}
}


/*!
Writes the game state into a file named game_of_life_, postfixed with the timestep and .dc
Fileformat:
uint64_t time step
double   start_x
double   start_y
double   start_z
double   cell_x_size
double   cell_y_size
double   cell_z_size
uint64_t x_length
uint64_t y_length
uint64_t z_length
int8_t   maximum_refinement_level
uint64_t 1st cell
uint64_t  1st cell is_alive
uint64_t 2nd cell
uint64_t  2nd cell is_alive
...
*/
bool write_game_data(const uint64_t step, communicator comm, Dccrg<game_of_life_cell>* game_grid)
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

	delete [] output_name_c_string;

	// figure out how many bytes each process will write
	size_t bytes = 0;

	// header
	if (comm.rank() == 0) {
		bytes += sizeof(int) + 4 * sizeof(uint64_t) + 6 * sizeof(double);
	}
	vector<uint64_t> cells = game_grid->get_cells();
	bytes += cells.size() * (sizeof(uint64_t) + sizeof(uint64_t));

	// collect data from this process into one buffer
	uint8_t* buffer = new uint8_t [bytes];

	size_t offset = 0;

	// header
	if (comm.rank() == 0) {

		memcpy(buffer + offset, &step, sizeof(uint64_t));
		offset += sizeof(uint64_t);

		{
		double value = game_grid->get_start_x();
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);
		value = game_grid->get_start_y();
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);
		value = game_grid->get_start_z();
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);
		value = game_grid->get_cell_x_size(1);
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);
		value = game_grid->get_cell_y_size(1);
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);
		value = game_grid->get_cell_z_size(1);
		memcpy(buffer + offset, &value, sizeof(double));
		offset += sizeof(double);
		}
		{
		uint64_t value = game_grid->get_x_length();
		memcpy(buffer + offset, &value, sizeof(uint64_t));
		offset += sizeof(uint64_t);
		value = game_grid->get_y_length();
		memcpy(buffer + offset, &value, sizeof(uint64_t));
		offset += sizeof(uint64_t);
		value = game_grid->get_z_length();
		memcpy(buffer + offset, &value, sizeof(uint64_t));
		offset += sizeof(uint64_t);
		}
		{
		int value = game_grid->get_maximum_refinement_level();
		memcpy(buffer + offset, &value, sizeof(int));
		offset += sizeof(int);
		}
	}

	for (uint64_t i = 0; i < cells.size(); i++) {
		const uint64_t cell = cells[i];
		memcpy(buffer + offset, &cell, sizeof(uint64_t));
		offset += sizeof(uint64_t);

		game_of_life_cell* data = (*game_grid)[cells[i]];
		const uint64_t alive = data->is_alive ? 1 : 0;
		memcpy(buffer + offset, &alive, sizeof(uint64_t));
		offset += sizeof(uint64_t);
	}

	vector<size_t> all_bytes;
	all_gather(comm, bytes, all_bytes);

	// calculate offset of this process in the file
	MPI_Offset mpi_offset = 0;
	for (int i = 0; i < comm.rank(); i++) {
		mpi_offset += all_bytes[i];
	}
	MPI_File_set_view(outfile, mpi_offset, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL);

	MPI_Status status;
	MPI_File_write_at_all(outfile, 0, (void*)buffer, bytes, MPI_BYTE, &status);
	//if (status...

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

	Dccrg<game_of_life_cell> game_grid;

	#define X_LENGTH 10	// in unrefined cells
	#define Y_LENGTH 10
	#define Z_LENGTH 1
	#define CELL_SIZE 1.0
	game_grid.set_geometry(X_LENGTH, Y_LENGTH, Z_LENGTH, 0, 0, 0, CELL_SIZE, CELL_SIZE, CELL_SIZE);

	// the cells that share a vertex are considered neighbors
	#define NEIGHBORHOOD_SIZE 1
	#define MAX_REFINEMENT_LEVEL 0
	// use the recursive coordinate bisection method for load balancing (http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg_rcb.html)
	game_grid.initialize(comm, "RCB", NEIGHBORHOOD_SIZE, MAX_REFINEMENT_LEVEL);

	// since the grid doesn't change (isn't refined / unrefined) during the game, workload can be balanced just once in the beginning
	game_grid.balance_load();

	/*
	Get the cells on this process just once, since the grid doesn't change during the game
	To make the game scale better, separate local cells into those without even one neighbor on another process and those that do.
	While updating cell data between processes, start calculating the next turn for cells which don't have neighbors on other processes
	*/
	vector<uint64_t> cells_with_local_neighbors = game_grid.get_cells_with_local_neighbors();
	vector<uint64_t> cells_with_remote_neighbor = game_grid.get_cells_with_remote_neighbor();

	initialize_game(&cells_with_local_neighbors, &game_grid);
	initialize_game(&cells_with_remote_neighbor, &game_grid);

	#define TURNS 10
	for (unsigned int turn = 0; turn < TURNS; turn++) {

		write_game_data(turn, comm, &game_grid);

		// start updating cell data from other processes and calculate the next turn for cells without neighbors on other processes in the meantime
		game_grid.start_remote_neighbor_data_update();
		get_live_neighbor_counts(&cells_with_local_neighbors, &game_grid);

		// wait for neighbor data updates to finish and the calculate the next turn for rest of the cells on this process
		game_grid.wait_neighbor_data_update();
		get_live_neighbor_counts(&cells_with_remote_neighbor, &game_grid);

		// update the state of life for all local cells
		apply_rules(&cells_with_local_neighbors, &game_grid);
		apply_rules(&cells_with_remote_neighbor, &game_grid);
	}
	write_game_data(TURNS, comm, &game_grid);

	return EXIT_SUCCESS;
}
