/*
Turns the files saved from the GoL restart test into .vtk files,
visualize for example with VisIt (https://wci.llnl.gov/codes/visit/)
*/

#include "algorithm"
#include "boost/unordered_map.hpp"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "iostream"
#include "mpi.h"
#include "stdint.h"

#include "../../dccrg_cartesian_geometry.hpp"

using namespace std;
using namespace dccrg;

int main(int argc, char* argv[])
{
	int ret_val;

	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		cerr << "Coudln't initialize MPI." << endl;
		abort();
	}

	MPI_Comm comm = MPI_COMM_WORLD;

	int rank = 0, comm_size = 0;
	MPI_Comm_rank(comm, &rank);
	if (rank < 0) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}
	MPI_Comm_size(comm, &comm_size);
	if (comm_size < 0) {
		std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
		abort();
	}

	// each process converts its own files
	for (int arg = 1; arg < argc; arg++) {

		if ((arg - 1) % comm_size != rank) {
			continue;
		}

		const std::string name(argv[arg]);

		MPI_File file;
		MPI_Offset offset = 0;

		ret_val = MPI_File_open(
			MPI_COMM_SELF,
			argv[arg],
			MPI_MODE_RDONLY,
			MPI_INFO_NULL,
			&file
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't open file " << name
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			continue;
		}

		// read time step
		uint64_t time_step;
		ret_val = MPI_File_read_at(
			file,
			offset,
			&time_step,
			1,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't read time step from file " << name
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			continue;
		}
		offset += sizeof(uint64_t);
		//cout << "time step: " << time_step << endl;

		// check endianness
		uint64_t
			endianness_original = 0x1234567890abcdef,
			endianness_read = 0;

		ret_val = MPI_File_read_at(
			file,
			offset,
			&endianness_read,
			1,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't read endianness check from file " << name
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			continue;
		}
		offset += sizeof(uint64_t);

		if (endianness_original != endianness_read) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " File " << name
				<< " has wrong endianness, value from file is "
				<< endianness_read
				<< " but should be " << endianness_original
				<< std::endl;
			continue;
		}

		// initialize mapping
		Mapping mapping;
		if (!mapping.read(file, offset)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't read mapping from file " << name
				<< std::endl;
			continue;
		}
		offset += mapping.data_size();

		//cout << "x_length: " << mapping.length.get()[0] << "\n"
		//	<< "y_length: " << mapping.length.get()[1] << "\n"
		//	<< "z_length: " << mapping.length.get()[2] << "\n"
		//	<< "max_ref_level: " << mapping.get_maximum_refinement_level()
		//	<< endl;

		// read neighborhood length
		unsigned int neighborhood_length = 0;
		if (MPI_File_read_at(
			file,
			offset,
			(void*) &neighborhood_length,
			1,
			MPI_UNSIGNED,
			MPI_STATUS_IGNORE
		) != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't read length of cells' neighborhood from file " << name
				<< std::endl;
			continue;
		}
		offset += sizeof(unsigned int);

		// initialize topology
		Grid_Topology topology;
		if (!topology.read(file, offset)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't read grid topology from file " << name
				<< std::endl;
			continue;
		}
		offset += topology.data_size();

		//cout << "periodic in x: " << topology.is_periodic(0) << "\n"
		//	<< "periodic in y: " << topology.is_periodic(1) << "\n"
		//	<< "periodic in z: " << topology.is_periodic(2) << "\n"
		//	<< endl;

		// initialize geometry
		Cartesian_Geometry geometry(mapping.length, mapping, topology);
		if (!geometry.read(file, offset)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << rank
				<< " Couldn't read geometry from file"
				<< std::endl;
			continue;
		}
		offset += geometry.data_size();

		//cout << "x_start: " << geometry.get_start()[0] << "\n"
		//	<< "y_start: " << geometry.get_start()[1] << "\n"
		//	<< "z_start: " << geometry.get_start()[2] << "\n"
		//	<< "level 0 cell x size: " << geometry.get_level_0_cell_length()[0] << "\n"
		//	<< "level 0 cell y size: " << geometry.get_level_0_cell_length()[1] << "\n"
		//	<< "level 0 cell z size: " << geometry.get_level_0_cell_length()[2]
		//	<< endl;

		uint64_t number_of_cells;
		if (MPI_File_read_at(
			file,
			offset,
			&number_of_cells,
			1,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		) != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		offset += sizeof(uint64_t);
		//cout << "number of cells: " << number_of_cells << endl;

		// read in game data (1st cell, is_alive 1st), (2nd cell, is_alive 2nd), ...
		std::vector<std::pair<uint64_t, uint64_t> > cell_data(number_of_cells);

		// read in cell list
		for (uint64_t i = 0; i < number_of_cells; i++) {
			if (MPI_File_read_at(
				file,
				offset,
				&(cell_data[i].first),
				1,
				MPI_UINT64_T,
				MPI_STATUS_IGNORE
			) != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			offset += sizeof(uint64_t);
			// skip data offset
			offset += sizeof(uint64_t);
		}

		// read in cell data
		for (uint64_t i = 0; i < number_of_cells; i++) {
			if (MPI_File_read_at(
				file,
				offset,
				&(cell_data[i].second),
				1,
				MPI_UINT64_T,
				MPI_STATUS_IGNORE
			) != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}
			offset += sizeof(uint64_t);
		}
		MPI_File_close(&file);

		// write the game data to a .vtk file
		const string input_name(argv[arg]),
			current_output_name(input_name.substr(0, input_name.size() - 2) + "vtk");

		std::ofstream outfile(current_output_name.c_str());
		if (!outfile.is_open()) {
			std::cerr << "Couldn't open file " << current_output_name << std::endl;
			exit(EXIT_FAILURE);
		}

		outfile
			<< "# vtk DataFile Version 2.0\n"
			<< "Game of Life data\n"
			<< "ASCII\n"
			<< "DATASET UNSTRUCTURED_GRID"
			<< std::endl;

		// write separate points for every cells' corners
		outfile << "POINTS " << cell_data.size() * 8 << " float" << std::endl;
		for (uint64_t i = 0; i < cell_data.size(); i++) {
			const uint64_t cell = cell_data[i].first;
			const boost::array<double, 3>
				cell_min = geometry.get_min(cell),
				cell_max = geometry.get_max(cell);

			outfile
				<< cell_min[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
				<< cell_max[0] << " " << cell_min[1] << " " << cell_min[2] << "\n"
				<< cell_min[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
				<< cell_max[0] << " " << cell_max[1] << " " << cell_min[2] << "\n"
				<< cell_min[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
				<< cell_max[0] << " " << cell_min[1] << " " << cell_max[2] << "\n"
				<< cell_min[0] << " " << cell_max[1] << " " << cell_max[2] << "\n"
				<< cell_max[0] << " " << cell_max[1] << " " << cell_max[2] << "\n";
		}

		// map cells to written points
		outfile << "CELLS " << cell_data.size() << " " << cell_data.size() * 9 << "\n";
		for (unsigned int j = 0; j < cell_data.size(); j++) {
			outfile << "8 ";
			for (int i = 0; i < 8; i++) {
				 outfile << j * 8 + i << " ";
			}
			outfile << "\n";
		}

		// cell types
		outfile << "CELL_TYPES " << cell_data.size() << "\n";
		for (unsigned int i = 0; i < cell_data.size(); i++) {
			outfile << 11 << "\n";
		}

		outfile << "CELL_DATA " << cell_data.size() << "\n";
		outfile << "SCALARS is_alive int 1" << "\n";
		outfile << "LOOKUP_TABLE default" << "\n";
		for (uint64_t i = 0; i < cell_data.size(); i++) {
			outfile << cell_data[i].second << "\n";
		}
	}

	MPI_Finalize();

	return EXIT_SUCCESS;
}

