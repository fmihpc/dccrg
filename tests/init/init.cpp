/*
Tests the scalability of initializing the grid
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "zoltan.h"

#define DCCRG_ARBITRARY_STRETCH
#include "../../dccrg.hpp"


struct CellData {

	double variables[3];

	template<typename Archiver> void serialize(Archiver& ar, const unsigned int /*version*/) {
		ar & variables[0] & variables[1] & variables[2];
	}
};


using namespace std;
using namespace boost::mpi;

int main(int argc, char* argv[])
{
	environment env(argc, argv);
	communicator comm;

	float zoltan_version;
	if (Zoltan_Initialize(argc, argv, &zoltan_version) != ZOLTAN_OK) {
	    cout << "Zoltan_Initialize failed" << endl;
	    exit(EXIT_FAILURE);
	}

	#define GRID_SIZE 100	// in unrefined cells
	#define CELL_SIZE (1.0 / GRID_SIZE)
	vector<double> x_coordinates, y_coordinates, z_coordinates;
	for (int i = 0; i <= GRID_SIZE; i++) {
		x_coordinates.push_back(i * CELL_SIZE);
		y_coordinates.push_back(i * CELL_SIZE);
		z_coordinates.push_back(i * CELL_SIZE);
	}
	clock_t before = clock();
	dccrg<CellData> game_grid(comm, "RCB", x_coordinates, y_coordinates, z_coordinates, 1, 0);
	clock_t after = clock();
	cout << "Process " << comm.rank() << ": grid initialization took " << double(after - before) / CLOCKS_PER_SEC << " seconds (total grid size " << (x_coordinates.size() - 1) * (y_coordinates.size() - 1) * (z_coordinates.size() - 1) << ")" << endl;

	// don't write too large a grid to disk
	if ((x_coordinates.size() - 1) * (y_coordinates.size() - 1) * (z_coordinates.size() - 1) > 1000) {
		return EXIT_SUCCESS;
	}

	// write the grid onto disk
	vector<uint64_t> cells = game_grid.get_cells();
	sort(cells.begin(), cells.end());
	ostringstream basename, suffix(".vtk");
	basename << "init_" << comm.rank();
	ofstream outfile, visit_file;
	if (comm.rank() == 0) {
		visit_file.open("init.visit");
		visit_file << "!NBLOCKS " << comm.size() << endl;
	}
	string current_output_name("");
	current_output_name += basename.str();
	current_output_name += suffix.str();
	if (comm.rank() == 0) {
		for (int process = 0; process < comm.size(); process++) {
			visit_file << "init_" << process << suffix.str() << endl;
		}
	}
	game_grid.write_vtk_file(current_output_name.c_str());
	outfile.open(current_output_name.c_str(), ofstream::app);
	outfile << "CELL_DATA " << cells.size() << endl;
	outfile << "SCALARS process int 1" << endl;
	outfile << "LOOKUP_TABLE default" << endl;
	for (vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
		outfile << comm.rank() << endl;
	}
	outfile.close();
	if (comm.rank() == 0) {
		visit_file.close();
	}

	return EXIT_SUCCESS;
}
