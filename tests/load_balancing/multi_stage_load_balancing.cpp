/*
Tests load balancing when cell data has to be sent more than once.

Required for example when the amount of data to be received is not
known in advance.
*/

#include "algorithm"
#include "boost/program_options.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "fstream"
#include "iostream"
#include "mpi.h"
#include "vector"
#include "stdint.h"
#include "zoltan.h"

#include "../../dccrg.hpp"

using namespace std;
using namespace dccrg;

/*!
Cell data class for testing multi stage load balancing.
*/
class Cell
{
public:
	// has as many items as there will
	// be load balancing stages
	std::vector<int> data;

	// load balancing stage currently in progress
	static size_t stage;

	boost::tuple<
		void*,
		int,
		MPI_Datatype
	> get_mpi_datatype(
		const uint64_t /*cell_id*/,
		const int /*sender*/,
		const int /*receiver*/,
		const bool /*receiving*/
	) {
		return boost::make_tuple(&(this->data[Cell::stage]), 1, MPI_INT);
	}
};

size_t Cell::stage = 0;

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
	size_t lb_stages = 0;
	boost::program_options::options_description options("Usage: program_name [options], where options are:");
	options.add_options()
		("help", "print this help message")
		("lb-stages",
			boost::program_options::value<size_t>(&lb_stages)->default_value(5),
			"Send cell data in arg stages");

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


	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	Dccrg<Cell> grid;
	const boost::array<uint64_t, 3> grid_length = {{3, 3, 3}};
	grid.initialize(grid_length, comm, "RANDOM", 3);

	// local cells must be obtained before starting load balancing
	const vector<uint64_t> cells = grid.get_cells();
	BOOST_FOREACH(const uint64_t cell, cells) {
		Cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< endl;
			abort();
		}

		cell_data->data.resize(lb_stages);
		BOOST_FOREACH(int& data_item, cell_data->data) {
			data_item = -1;
		}
	}

	grid.initialize_balance_load(true);

	/*
	Prepare cell data for sending.

	Resize local cells' data to the number of stages.
	All items are initially == -1.
	*/

	// cells that will be moved to this process
	const boost::unordered_set<uint64_t>& added_cells
		= grid.get_cells_added_by_balance_load();
	BOOST_FOREACH(const uint64_t cell, added_cells) {
		Cell* cell_data = grid[cell];
		if (cell_data == NULL) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No data for cell " << cell
				<< endl;
			abort();
		}

		cell_data->data.resize(lb_stages);
	}

	// cells that will be moved from this process
	const boost::unordered_set<uint64_t>& removed_cells
		= grid.get_cells_removed_by_balance_load();

	// allow going through all cells with one loop
	vector<uint64_t> all_cells;
	all_cells.insert(all_cells.end(), cells.begin(), cells.end());
	all_cells.insert(all_cells.end(), added_cells.begin(), added_cells.end());

	/*
	Test by creating and transferring the following data:
	 0, -1, -1, -1, ... if current stage is 0
	-1,  1, -1, -1, ... if current stage is 1
	-1, -1,  2, -1, ... if current stage is 2
	...
	*/
	for (size_t current_stage = 0; current_stage < lb_stages; current_stage++) {
		Cell::stage = current_stage;

		// set all data of all cells to -1
		BOOST_FOREACH(const uint64_t cell, all_cells) {
			Cell* cell_data = grid[cell];
			BOOST_FOREACH(int& data_item, cell_data->data) {
				data_item = -1;
			}
		}

		// set one item of outgoing cells to current_stage
		BOOST_FOREACH(const uint64_t cell, removed_cells) {
			Cell* cell_data = grid[cell];
			cell_data->data[current_stage] = (int) current_stage;
		}

		// transfer data
		grid.continue_balance_load();

		// check that data arrived correctly
		BOOST_FOREACH(const uint64_t cell, added_cells) {
			Cell* cell_data = grid[cell];
			if (cell_data->data.size() != lb_stages) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Rank " << rank
					<< ": Incorrect data size in cell " << cell
					<< std::endl;
				abort();
			}
			for (size_t i = 0; i < cell_data->data.size(); i++) {
				if (i == current_stage) {
					if (cell_data->data[i] != (int) current_stage) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Rank " << rank
							<< ", in stage " << current_stage
							<< ": incorrect data item in cell " << cell
							<< " at index " << i
							<< ": " << cell_data->data[i]
							<< ", should be " << current_stage
							<< std::endl;
						abort();
					}
				} else {
					if (cell_data->data[i] != -1) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Rank " << rank
							<< ", in stage " << current_stage
							<< ": incorrect data item in cell " << cell
							<< " at index " << i
							<< ": " << cell_data->data[i]
							<< ", should be -1"
							<< std::endl;
						abort();
					}
				}
			}
		}
	}

	// check that correct cells will be moved
	boost::unordered_set<uint64_t> new_cells_reference;
	BOOST_FOREACH(const uint64_t cell, cells) {
		new_cells_reference.insert(cell);
	}
	BOOST_FOREACH(const uint64_t cell, added_cells) {
		new_cells_reference.insert(cell);
	}
	BOOST_FOREACH(const uint64_t cell, removed_cells) {
		new_cells_reference.erase(cell);
	}

	grid.finish_balance_load();

	const vector<uint64_t> new_local = grid.get_cells();
	BOOST_FOREACH(const uint64_t cell, new_local) {
		if (new_cells_reference.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Rank " << rank
				<< ": Cell " << cell
				<< " should not be on process " << rank
				<< std::endl;
			abort();
		}
	}

	if (rank == 0) {
		cout << "PASSED" << endl;
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

