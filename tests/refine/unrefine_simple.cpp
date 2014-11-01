/*
Tests the grid using simple refining and unrefining
*/

#include "cstdlib"
#include "fstream"
#include "iostream"
#include "sstream"

#include "mpi.h"
#include "zoltan.h"

#include "../../dccrg_stretched_cartesian_geometry.hpp"
#include "../../dccrg.hpp"


using namespace std;
using namespace dccrg;

struct Cell {
	std::tuple<void*, int, MPI_Datatype> get_mpi_datatype()
	{
		return std::make_tuple(this, 0, MPI_BYTE);
	}
};

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

	clock_t before, after;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}
	if (rank == 0) {
		cout << "Using Zoltan version " << zoltan_version << endl;
	}

	Dccrg<Cell, Stretched_Cartesian_Geometry> grid;

	#define GRID_SIZE 2
	const std::array<uint64_t, 3> grid_length{{GRID_SIZE, 1, 1}};
	#define NEIGHBORHOOD_SIZE 1
	grid.initialize(grid_length, comm, "RANDOM", NEIGHBORHOOD_SIZE, 5);

	#define CELL_SIZE (1.0 / GRID_SIZE)
	Stretched_Cartesian_Geometry::Parameters geom_params;
	for (size_t i = 0; i <= GRID_SIZE; i++) {
		geom_params.coordinates[0].push_back(i * CELL_SIZE);
	}
	geom_params.coordinates[1].push_back(0);
	geom_params.coordinates[1].push_back(1);
	geom_params.coordinates[2].push_back(0);
	geom_params.coordinates[2].push_back(1);
	grid.set_geometry(geom_params);

	if (rank == 0) {
		cout << "Maximum refinement level of the grid: "
			<< grid.get_maximum_refinement_level()
			<< "\nNumber of cells: "
			<< (geom_params.coordinates.size() - 1)
				* (geom_params.coordinates.size() - 1)
				* (geom_params.coordinates.size() - 1)
			<< endl << endl;
	}

	// every process outputs state into its own file
	ostringstream basename, suffix(".vtk");
	basename << "unrefine_simple_" << rank << "_";
	ofstream outfile, visit_file;

	// visualize with visit -o game_of_life_test.visit
	if (rank == 0) {
		visit_file.open("unrefine_simple.visit");
		visit_file << "!NBLOCKS " << comm_size << endl;
	}

	#define TIME_STEPS 8
	for (int step = 0; step < TIME_STEPS; step++) {

		if (rank == 0) {
			cout << "step " << step << endl;
		}

		grid.balance_load();
		auto cells = grid.get_cells();
		sort(cells.begin(), cells.end());

		// write state into a file named according to the current time step
		string current_output_name("");
		ostringstream step_string;
		step_string.fill('0');
		step_string.width(5);
		step_string << step;
		current_output_name += basename.str();
		current_output_name += step_string.str();
		current_output_name += suffix.str();

		// visualize with visit -o game_of_life_test.visit
		if (rank == 0) {
			for (int process = 0; process < comm_size; process++) {
				visit_file << "unrefine_simple_" << process
					<< "_" << step_string.str() << suffix.str()
					<< endl;
			}
		}

		// write the grid into a file
		grid.write_vtk_file(current_output_name.c_str());
		// prepare to write the game data into the same file
		outfile.open(current_output_name.c_str(), ofstream::app);
		outfile << "CELL_DATA " << cells.size() << endl;

		// write each cells neighbor count
		outfile << "SCALARS neighbors int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			const auto* const neighbors = grid.get_neighbors_of(cell);
			outfile << neighbors->size() << endl;
		}

		// write each cells process
		outfile << "SCALARS process int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (size_t i = 0; i < cells.size(); i++) {
			outfile << rank << endl;
		}

		// write each cells id
		outfile << "SCALARS id int 1" << endl;
		outfile << "LOOKUP_TABLE default" << endl;
		for (const auto& cell: cells) {
			outfile << cell << endl;
		}
		outfile.close();

		before = clock();

		// refine / unrefine the smallest cell that is closest to the grid starting corner
		const std::array<double, 3> adapt_coord{{
			0.0001 * CELL_SIZE,
			0.0001 * CELL_SIZE,
			0.0001 * CELL_SIZE
		}};
		if (step < 4) {
			grid.refine_completely_at(adapt_coord);
		} else {
			grid.unrefine_completely_at(adapt_coord);
		}

		auto new_cells = grid.stop_refining();

		after = clock();
		cout << "Process " << rank
			<<": Refining / unrefining took " << double(after - before) / CLOCKS_PER_SEC
			<< " seconds, " << new_cells.size()
			<< " new cells created"
			<< endl;
	}

	if (rank == 0) {
		visit_file.close();
	}

	return EXIT_SUCCESS;
}

