/*
Tests the grid with some simple game of life patters
*/

#include "boost/mpi.hpp"
#include "boost/unordered_set.hpp"
#include "cstdlib"
#include "ctime"
#include "../../dccrg.hpp"
#include "zoltan.h"


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

	clock_t before = clock();
	#define GRID_SIZE 100	// in unrefined cells
	dccrg<Data> game_grid(comm, "RCB", 0.0, 0.0, 0.0, 1, GRID_SIZE, GRID_SIZE, GRID_SIZE, 1, 0);
	clock_t after = clock();
	cout << "Process " << comm.rank() << ": grid initialization took " << double(after - before) / CLOCKS_PER_SEC << " seconds" << endl;

	if (GRID_SIZE * GRID_SIZE * GRID_SIZE > 1000) {
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
