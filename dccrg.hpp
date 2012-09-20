/*
A distributed cartesian cell-refinable grid

Copyright 2009, 2010, 2011, 2012 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_HPP
#define DCCRG_HPP


/*
If compilation fails with a complaint about
MPI_UNSIGNED_LONG_LONG try replacing it with
MPI_UNSIGNED_LONG in the following:
*/
#ifndef MPI_UINT64_T
#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
#endif


/*!
\mainpage Distributed Cartesian Cell-Refinable Grid.

\section intro_sec Introduction
dccrg is a grid library for simulations using the finite volume method.
See the examples directory for some simple examples and the tests directory
for more advanced usage of dccrg.
*/



/*
If the size of the data in every cell is known in advance by the user, neighbor data updates can be optimized by defining DCCRG_CELL_DATA_SIZE_FROM_USER, in which case:
	-UserData class must have a static function size() which returns the size of data in bytes of all cells.
	-UserData instances must have a function at() which returns the starting address of their data.

Additionally if DCCRG_USER_MPI_DATA_TYPE is defined:
	-UserData class must have a static function mpi_datatype() which returns the MPI_Datatype of all cells.
	-UserData function size() is not needed (size() == 1 is assumed).
	-Cells can have non-contiguous data.

If DCCRG_SEND_SINGLE_CELLS is defined then cell data is sent one cell at a time.

*/
#ifdef DCCRG_USER_MPI_DATA_TYPE
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		#error DCCRG_CELL_DATA_SIZE_FROM_USER must defined when using DCCRG_USER_MPI_DATA_TYPE
	#endif
#endif

#include "algorithm"
#include "boost/array.hpp"
#include "boost/foreach.hpp"
#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
#include "boost/mpi.hpp"
#endif
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cassert"
#include "cstdio"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "functional"
#include "limits"
#include "stdint.h"
#include "utility"
#include "vector"
#include "zoltan.h"

#ifdef USE_SFC
#include "sfc++.hpp"
#endif

#include "dccrg_index.hpp"
#include "dccrg_mpi_support.hpp"
#include "dccrg_types.hpp"
#include "dccrg_constant_geometry.hpp"


namespace dccrg
{

template <class UserData, class UserGeometry = ConstantGeometry> class Dccrg : public UserGeometry
{

public:

	// helper type for iterating over this->cells using BOOST_FOREACH
	typedef typename std::pair<const uint64_t&, const UserData&> cell_and_data_pair_t;


	/*!
	Creates an uninitialized instance of the grid.

	The instance's set_geometry and initialize functions must be called
	before doing anything else, otherwise the results will be undefined.
	*/
	Dccrg()
	{
		this->initialized = false;
	}


	/*!
	Initializes the instance of the grid with given parameters.

	The geometry of the grid instance must have been set using set_geometry
	before calling this function.

	Zoltan_Initialize must have been called before calling this function.

	comm: the grid will span all the processes in the communicator comm

	load_balancing_method:
		- The method that Zoltan will use for load balancing, given as a string.
		- All methods except REFTREE are supported, see this page for a list of available methods:
		- http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg.html#LB_METHOD

	neighborhood_size:
		- Determines which cells are considered neighbors.
		- When calculating the neighbors of a given cell a cube of length
		  2 * neighborhood_size + 1 in every direction is considered, centered
		  at the cell for which neighbors are being calculated.
		- The unit lenght of the cube is the cell for which neighbors are being calculated.
		- If neighborhood_size == 0, only cells (or children within the volume of
		  cells of the same size as the current cell) that share a face are considered.

	maximum_refinement_level:
		- The maximum number of times an unrefined cell can be refined
		  (replacingit with 8 smaller cells).
		- If not given the maximum refinement level is maximized based on the grids initial size.

	periodic_in_x, y and z:
		- The grid neighborhoods wrap around in periodic directions, e.g. if periodic in some
		  direction cells on the opposite sides of the grid in that direction can be neighbors.
	*/
	void initialize(
		const MPI_Comm& given_comm,
		const char* load_balancing_method,
		const unsigned int neighborhood_size,
		const int maximum_refinement_level = -1,
		const bool periodic_in_x = false,
		const bool periodic_in_y = false,
		const bool periodic_in_z = false,
		const uint64_t sfc_caching_batches = 1
	) {
		if (this->initialized) {
			std::cerr << "Initialize function called for an already initialized dccrg" << std::endl;
			// TODO: throw an exception instead
			abort();
		}

		if (sfc_caching_batches == 0) {
			std::cerr << "sfc_caching_batches must be > 0" << std::endl;
			abort();
		}

		if (MPI_Comm_dup(given_comm, &this->comm) != MPI_SUCCESS) {
			std::cerr << "Couldn't duplicate given communicator" << std::endl;
			abort();
		}

		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->boost_comm = boost::mpi::communicator(this->comm, boost::mpi::comm_attach);
		#endif

		int temp_size = 0;
		if (MPI_Comm_size(this->comm, &temp_size) != MPI_SUCCESS) {
			std::cerr << "Couldn't get size of communicator" << std::endl;
			abort();
		}
		if (temp_size < 0) {
			std::cerr << "Negative MPI comm size not supported" << std::endl;
			abort();
		}
		this->comm_size = (uint64_t) temp_size;

		int temp_rank = 0;
		if (MPI_Comm_rank(this->comm, &temp_rank) != MPI_SUCCESS) {
			std::cerr << "Couldn't get rank for communicator" << std::endl;
			abort();
		}
		if (temp_rank < 0) {
			std::cerr << "Negative MPI rank not supported" << std::endl;
			abort();
		}
		this->rank = (uint64_t) temp_rank;

		this->max_ref_lvl_diff = 1;

		/*
		Setup Zoltan
		*/
		this->zoltan = Zoltan_Create(this->comm);
		if (this->zoltan == NULL) {
			std::cerr << "Zoltan_Create failed"  << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}

		// check whether Zoltan_LB_Partition is expected to fail
		if (strncmp(load_balancing_method, "NONE", sizeof("NONE")) == 0) {
			this->no_load_balancing = true;
		} else {
			this->no_load_balancing = false;
		}

		// reserved options that the user cannot change
		this->reserved_options.insert("EDGE_WEIGHT_DIM");
		this->reserved_options.insert("NUM_GID_ENTRIES");
		this->reserved_options.insert("NUM_LID_ENTRIES");
		this->reserved_options.insert("OBJ_WEIGHT_DIM");
		this->reserved_options.insert("RETURN_LISTS");
		this->reserved_options.insert("NUM_GLOBAL_PARTS");
		this->reserved_options.insert("NUM_LOCAL_PARTS");
		this->reserved_options.insert("AUTO_MIGRATE");

		// set reserved options
		Zoltan_Set_Param(this->zoltan, "EDGE_WEIGHT_DIM", "0");	// 0 because Zoltan crashes in hierarchial with larger values
		Zoltan_Set_Param(this->zoltan, "NUM_GID_ENTRIES", "1");
		Zoltan_Set_Param(this->zoltan, "NUM_LID_ENTRIES", "0");
		Zoltan_Set_Param(this->zoltan, "OBJ_WEIGHT_DIM", "1");
		Zoltan_Set_Param(this->zoltan, "RETURN_LISTS", "ALL");

		// set other options
		Zoltan_Set_Param(this->zoltan, "DEBUG_LEVEL", "0");
		Zoltan_Set_Param(this->zoltan, "HIER_DEBUG_LEVEL", "0");
		Zoltan_Set_Param(this->zoltan, "HIER_CHECKS", "0");
		Zoltan_Set_Param(this->zoltan, "LB_METHOD", load_balancing_method);
		Zoltan_Set_Param(this->zoltan, "REMAP", "1");


		// size of cells id in unsigned ints, but has to be 1 even when global id is uint64_t, for some reason
		/*char global_id_length_string[10];
		snprintf(global_id_length_string, 10, "%0i", int(sizeof(uint64_t) / sizeof(unsigned int)));*/

		// set the grids callback functions in Zoltan
		Zoltan_Set_Num_Obj_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::get_number_of_cells, this);
		Zoltan_Set_Obj_List_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_cell_list, this);
		Zoltan_Set_Num_Geom_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::get_grid_dimensionality, NULL);
		Zoltan_Set_Geom_Multi_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_with_cell_coordinates, this);
		Zoltan_Set_Num_Edges_Multi_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_number_of_neighbors_for_cells, this);
		Zoltan_Set_Edge_List_Multi_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_neighbor_lists, this);
		Zoltan_Set_HG_Size_CS_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_number_of_hyperedges, this);
		Zoltan_Set_HG_CS_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_hyperedge_lists, this);
		Zoltan_Set_HG_Size_Edge_Wts_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_number_of_edge_weights, this);
		Zoltan_Set_HG_Edge_Wts_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::fill_edge_weights, this);
		Zoltan_Set_Hier_Num_Levels_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::get_number_of_load_balancing_hierarchies, this);
		Zoltan_Set_Hier_Part_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::get_part_number, this);
		Zoltan_Set_Hier_Method_Fn(this->zoltan, &Dccrg<UserData, UserGeometry>::set_partitioning_options, this);


		/*
		Set grid parameters
		*/

		this->set_periodicity(0, periodic_in_x);
		this->set_periodicity(1, periodic_in_y);
		this->set_periodicity(2, periodic_in_z);

		// set / check neighborhood_of
		this->neighborhood_size = neighborhood_size;
		if (this->neighborhood_size == 0) {
			{
			Types<3>::neighborhood_item_t item = {{0, 0, -1}};
			this->neighborhood_of.push_back(item);
			}
			{
			Types<3>::neighborhood_item_t item = {{0, -1, 0}};
			this->neighborhood_of.push_back(item);
			}
			{
			Types<3>::neighborhood_item_t item = {{-1, 0, 0}};
			this->neighborhood_of.push_back(item);
			}
			{
			Types<3>::neighborhood_item_t item = {{1, 0, 0}};
			this->neighborhood_of.push_back(item);
			}
			{
			Types<3>::neighborhood_item_t item = {{0, 1, 0}};
			this->neighborhood_of.push_back(item);
			}
			{
			Types<3>::neighborhood_item_t item = {{0, 0, 1}};
			this->neighborhood_of.push_back(item);
			}
		} else {
			for (int z = -neighborhood_size; (unsigned int) abs(z) < neighborhood_size + 1; z++)
			for (int y = -neighborhood_size; (unsigned int) abs(y) < neighborhood_size + 1; y++)
			for (int x = -neighborhood_size; (unsigned int) abs(x) < neighborhood_size + 1; x++) {
				if (x == 0 && y == 0 && z == 0) {
					continue;
				}
				const Types<3>::neighborhood_item_t item = {{x, y, z}};
				this->neighborhood_of.push_back(item);
			}
		}

		// set neighborhood_to
		BOOST_FOREACH(const Types<3>::neighborhood_item_t& offset, this->neighborhood_of) {
			Types<3>::neighborhood_item_t item = {{-offset[0], -offset[1], -offset[2]}};
			this->neighborhood_to.push_back(item);
		}


		if (maximum_refinement_level < 0) {
			this->set_maximum_refinement_level(this->get_maximum_possible_refinement_level());
		} else if (!this->set_maximum_refinement_level(maximum_refinement_level)) {
			std::cerr << "Couldn't set maximum refinement level to " << maximum_refinement_level << std::endl;
			abort();
		}

		// TODO: check that the last index in the grid in every direction is less than error_index


		// create unrefined cells
		uint64_t cells_per_process = 0;
		if (this->grid_length < this->comm_size) {
			cells_per_process = 1;
		} else if (this->grid_length % this->comm_size > 0) {
			cells_per_process = this->grid_length / this->comm_size + 1;
		} else {
			cells_per_process = this->grid_length / this->comm_size;
		}

		// some processes get fewer cells if grid size not divisible by this->comm_size
		uint64_t procs_with_fewer = cells_per_process * this->comm_size - this->grid_length;

		#ifndef USE_SFC

		uint64_t cell_to_create = 1;
		for (unsigned int process = 0; process < this->comm_size; process++) {

			uint64_t cells_to_create;
			if (process < procs_with_fewer) {
				cells_to_create = cells_per_process - 1;
			} else {
				cells_to_create = cells_per_process;
			}

			for (uint64_t i = 0; i < cells_to_create; i++) {
				this->cell_process[cell_to_create] = process;
				if (process == this->rank) {
					this->cells[cell_to_create];
				}
				cell_to_create++;
			}
		}
		assert(cell_to_create == this->grid_length + 1);

		#else

		const dccrg::Types<3>::indices_t length = {{
			this->x_length,
			this->y_length,
			this->z_length
		}};
		sfc::Sfc<3, uint64_t> mapping(length);

		/*
		Cache only batch_size number of sfc indices at a time.
		Saves memory and can even be faster than caching everything at once
		*/
		uint64_t batch_size;
		if (mapping.size() % sfc_caching_batches > 0) {
			batch_size = 1 + mapping.size() / sfc_caching_batches;
		} else {
			batch_size = mapping.size() / sfc_caching_batches;
		}

		uint64_t cache_start = 0, cache_end = batch_size - 1;
		mapping.cache_sfc_index_range(cache_start, cache_end);

		uint64_t sfc_index = 0;
		for (uint64_t process = 0; process < this->comm_size; process++) {

			uint64_t cells_to_create;
			if (process < procs_with_fewer) {
				cells_to_create = cells_per_process - 1;
			} else {
				cells_to_create = cells_per_process;
			}

			for (uint64_t i = 0; i < cells_to_create; i++) {

				// cache new sfc index batch
				if (sfc_index > cache_end) {
					cache_start = cache_end;
					cache_end = cache_start + batch_size;

					if (cache_end >= mapping.size()) {
						cache_end = mapping.size() - 1;
					}

					mapping.clear();
					mapping.cache_sfc_index_range(cache_start, cache_end);
				}

				dccrg::Types<3>::indices_t indices = mapping.get_indices(sfc_index);
				// transform indices to those of refinement level 0 cells
				indices[0] *= uint64_t(1) << this->max_refinement_level;
				indices[1] *= uint64_t(1) << this->max_refinement_level;
				indices[2] *= uint64_t(1) << this->max_refinement_level;
				const uint64_t cell_to_create = this->get_cell_from_indices(indices, 0);

				this->cell_process[cell_to_create] = process;
				if (process == this->rank) {
					this->cells[cell_to_create];
				}

				sfc_index++;
			}
		}
		mapping.clear();

		if (sfc_index != this->grid_length) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Process " << this->rank
				<< ": Incorrect number of cells created: " << sfc_index
				<< ", should be " << this->grid_length
				<< std::endl;
			abort();
		}

		#endif

		// update neighbor lists of created cells
		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
			this->neighbors[item.first]
				= this->find_neighbors_of(item.first, this->neighborhood_of, this->max_ref_lvl_diff);
			this->neighbors_to[item.first]
				= this->find_neighbors_to(item.first, this->neighborhood_to);
		}
		#ifdef DEBUG
		if (!this->verify_neighbors()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbor lists are inconsistent" << std::endl;
			abort();
		}
		#endif

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
			this->update_remote_neighbor_info(item.first);
		}
		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info is not consistent"
				<< std::endl;
			abort();
		}
		#endif

		this->recalculate_neighbor_update_send_receive_lists();

		this->initialized = true;
	}


	/*!
	Returns all cells on this process that don't have children (e.g. leaf cells)
	*/
	std::vector<uint64_t> get_cells() const
	{
		std::vector<uint64_t> all_cells;
		all_cells.reserve(this->cells.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			#ifdef DEBUG
			if (this->cell_process.count(item.first) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << item.first
					<< " shouldn't exist"
					<< std::endl;
				abort();
			}

			if (this->cell_process.at(item.first) != this->rank) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Cell " << item.first
					<< " should be on process " << this->cell_process.at(item.first)
					<< std::endl;
				abort();
			}

			const uint64_t child = this->get_child(item.first);
			if (child == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Child == 0 for cell " << item.first
					<< std::endl;
				abort();
			}

			if (child != item.first) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Cell " << item.first
					<< " has a child"
					<< std::endl;
				abort();
			}
			#endif

			all_cells.push_back(item.first);
		}

		return all_cells;
	}


	/*!
	Writes the current grid data into a file with given name.

	Returns true on success, false otherwise (one one or more processes)
	Data stored in local cells is also written.
	The file is written in parallel using MPI_IO.
	Must be called by all processes.
	Data is written starting at start_offset given in bytes
	(e.g. write global simulation data yourself into the
	beginning of the file).
	Requires at least as much additional memory as local cell data.

	Does nothing unless DCCRG_CELL_DATA_SIZE_FROM_USER and
	DCCRG_USER_MPI_DATA_TYPE are defined
	*/
	bool write_grid(const std::string& name, const MPI_Offset& start_offset)
	{
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		#ifdef DCCRG_USER_MPI_DATA_TYPE

		/*
		File format:

		uint8_t * start_offset, data skipped by this function
		uint64_t number of cells
		uint64_t 1st cell id
		uint64_t offset in bytes to data of 1st cell in file
		uint64_t 2nd cell id
		uint64_t offset...
		...
		uint64_t offset to data of last cell in file
		uint8_t * N, data of 1st cell
		uint8_t * M, data of 2nd cell
		...
		*/

		// MPI_File_open wants a non-constant string
		char* name_c_string = new char [name.size() + 1];
		strncpy(name_c_string, name.c_str(), name.size() + 1);

		MPI_File outfile;

		int result = MPI_File_open(
			this->comm,
			name_c_string,
			MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL,
			&outfile
		);

		delete [] name_c_string;

		if (result != MPI_SUCCESS) {
			std::cerr << "Couldn't open file " << name_c_string
				<< ": " << Error_String()(result)
				<< std::endl;
			return false;
		}

		// datatypes of local cells
		std::vector<MPI_Datatype> datatypes;
		datatypes.reserve(this->cells.size());

		// 1 datatype per local cell
		std::vector<int> block_lengths(this->cells.size(), 1);

		// address of local cell data relative to data of first cell
		std::vector<MPI_Aint> displacements(this->cells.size(), 0);

		// offsets at which local cells' data starts in the file
		std::vector<std::pair<uint64_t, uint64_t> > local_cell_offsets(this->cells.size());

		// initialize offset info with dummy offsets, update when gathering datatypes
		size_t i = 0;
		BOOST_FOREACH(const cell_and_data_pair_t& cell_and_data, this->cells) {
			const uint64_t cell = cell_and_data.first;
			local_cell_offsets[i] = std::make_pair(cell, 0);
			i++;
		}

		// gather datatypes, etc. of local cells
		uint64_t local_number_of_bytes = 0;
		for (size_t i = 0; i < local_cell_offsets.size(); i++) {
			const uint64_t cell = local_cell_offsets[i].first;
			datatypes.push_back(this->cells.at(cell).mpi_datatype());

			// cell data displacement relative to data of first cell
			displacements[i]
				= (uint8_t*) this->cells.at(cell).at()
				- (uint8_t*) this->cells.at(local_cell_offsets[0].first).at();

			// set cell data offsets in file to be relative to local cells
			local_cell_offsets[i].second += local_number_of_bytes;
			int current_bytes;
			MPI_Type_size(datatypes[i], &current_bytes);
			local_number_of_bytes += (uint64_t) current_bytes;
		}

		// get number of cells each process will write
		std::vector<uint64_t> number_of_cells(this->comm_size, 0);
		uint64_t local_cells = this->cells.size();
		if (
			MPI_Allgather(
				&local_cells,
				1,
				MPI_UINT64_T,
				&(number_of_cells[0]),
				1,
				MPI_UINT64_T,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__ << "MPI_Allgather failed." << std::endl;
			abort();
		}

		MPI_Comm non_boost_comm = this->comm;
		const uint64_t total_number_of_cells =
			All_Reduce()(this->cells.size(), non_boost_comm);

		// process 0 writes the total number of cells
		if (this->rank == 0) {
			result = MPI_File_write_at_all(
				outfile,
				start_offset,
				(void*) &total_number_of_cells,
				sizeof(uint64_t),
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);

			if (result != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't write cell list to file " << name
					<< ": " << Error_String()(result)
					<< std::endl;
				return false;
			}
		} else {
			result = MPI_File_write_at_all(
				outfile,
				start_offset,
				(void*) &total_number_of_cells,
				0,
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);
		}

		// get offset of local cell list in output file
		size_t cell_list_offset = (size_t) start_offset + sizeof(uint64_t);
		for (size_t i = 0; i < (size_t) this->rank; i++) {
			cell_list_offset += number_of_cells[i] * sizeof(std::pair<uint64_t, uint64_t>);
		}

		// get number of bytes each process will write
		std::vector<uint64_t> number_of_bytes(this->comm_size, 0);
		if (
			MPI_Allgather(
				&local_number_of_bytes,
				1,
				MPI_UINT64_T,
				&(number_of_bytes[0]),
				1,
				MPI_UINT64_T,
				comm
			) != MPI_SUCCESS
		) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		// get offset of local cell data in output file
		size_t cell_data_offset =
			(size_t) start_offset
			+ sizeof(uint64_t)
			+ total_number_of_cells * sizeof(std::pair<uint64_t, uint64_t>);

		for (size_t i = 0; i < (size_t) this->rank; i++) {
			cell_data_offset += number_of_bytes[i];
		}

		for (size_t i = 0; i < local_cell_offsets.size(); i++) {
			local_cell_offsets[i].second += cell_data_offset;
		}

		if (this->cells.size() > 0) {

			// write cell list
			result = MPI_File_write_at_all(
				outfile,
				(MPI_Aint) cell_list_offset,
				(void*) &local_cell_offsets[0],
				local_cell_offsets.size() * sizeof(std::pair<uint64_t, uint64_t>),
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);

			if (result != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't write cell list to file " << name
					<< ": " << Error_String()(result)
					<< std::endl;
				return false;
			}

			// represent local cell data with one type
			MPI_Datatype final_type;

			if (MPI_Type_create_struct(
				datatypes.size(),
				&block_lengths[0],
				&displacements[0],
				&datatypes[0],
				&final_type) != MPI_SUCCESS
			) {
				std::cerr << "Process " << this->rank
					<< " Couldn't create final datatype"
					<< std::endl;
				abort();
			}

			if (MPI_Type_commit(&final_type) != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't commit final datatype"
					<< std::endl;
				abort();
			}

			// copy cell data to a contiguous buffer
			uint8_t* buffer = new uint8_t [local_number_of_bytes];
			if (buffer == NULL) {
				std::cerr << "Couldn't reserve memory for temporary buffer" << std::endl;
			}

			int buffer_position = 0;
			assert(MPI_Pack(
				this->cells.at(local_cell_offsets[0].first).at(),
				1,
				final_type,
				buffer,
				local_number_of_bytes,
				&buffer_position,
				this->comm
			) == MPI_SUCCESS);

			// deallocate datatypes
			if (MPI_Type_free(&final_type) != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "MPI_Type_free failed"
					<< std::endl;
				abort();
			}

			BOOST_FOREACH(MPI_Datatype& type, datatypes) {
				if (MPI_Type_free(&type) != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< "MPI_Type_free failed"
						<< std::endl;
					abort();
				}
			}

			// write cell data
			if (MPI_File_write_at_all(
				outfile,
				(MPI_Aint) cell_data_offset,
				buffer,
				local_number_of_bytes,
				MPI_BYTE,
				MPI_STATUS_IGNORE) != MPI_SUCCESS
			) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Couldn't write cell data"
					<< std::endl;
				abort();
			}

			delete [] buffer;

		} else {
			uint8_t temp;

			// write empty cell list
			MPI_File_write_at_all(
				outfile,
				0,
				&temp,
				0,
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);

			// write no cell data
			MPI_File_write_at_all(
				outfile,
				0,
				&temp,
				0,
				MPI_BYTE,
				MPI_STATUS_IGNORE
			);
		}

		if (MPI_File_close(&outfile) != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't close file " << name
				<< std::endl;
			return false;
		}
		#endif
		#endif

		return true;
	}


	/*!
	Reads data of local cells from given file.

	Returns true on success, false otherwise (one one or more processes)
	The file is read in parallel using MPI_IO.
	Must be called by all processes and all cells in the grid must be of
	refinement level 0 prior to calling this function.
	Data is read starting at start_offset given in bytes
	(e.g. read global simulation data yourself from the
	beginning of the file).
	TODO: ?Requires at least as much additional memory as local cell data?

	Does nothing unless DCCRG_CELL_DATA_SIZE_FROM_USER and
	DCCRG_USER_MPI_DATA_TYPE are defined

	TODO: Reads at most number_of_cells number of cell data offsets at a time,
	give a smaller number if all cell ids won't fit	into memory at once.
	*/
	bool read_grid(
		const std::string& name,
		const MPI_Offset start_offset,
		const uint64_t /*number_of_cells*/ = ~uint64_t(0)
	) {
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		#ifdef DCCRG_USER_MPI_DATA_TYPE

		// MPI_File_open wants a non-constant string
		char* name_c_string = new char [name.size() + 1];
		strncpy(name_c_string, name.c_str(), name.size() + 1);

		MPI_File infile;

		int result = MPI_File_open(
			this->comm,
			name_c_string,
			MPI_MODE_RDONLY,
			MPI_INFO_NULL,
			&infile
		);

		delete [] name_c_string;

		if (result != MPI_SUCCESS) {
			std::cerr << "Couldn't open file " << name_c_string
				<< ": " << Error_String()(result)
				<< std::endl;
			return false;
		}

		// read the total number of cells in the file
		uint64_t total_number_of_cells;
		if (MPI_File_read_at_all(
			infile,
			start_offset,
			&total_number_of_cells,
			sizeof(uint64_t),
			MPI_BYTE,
			MPI_STATUS_IGNORE) != MPI_SUCCESS
		) {
			std::cerr << "Couldn't read total number of cells" << std::endl;
			return false;
		}

		// read cell data offsets
		// TODO: read in batches of size number_of_cells
		std::vector<std::pair<uint64_t, uint64_t> > all_data_offsets(total_number_of_cells);
		if (MPI_File_read_at_all(
			infile,
			start_offset + sizeof(uint64_t),
			&all_data_offsets[0],
			total_number_of_cells * sizeof(std::pair<uint64_t, uint64_t>),
			MPI_BYTE,
			MPI_STATUS_IGNORE) != MPI_SUCCESS
		) {
			std::cerr << "Couldn't read number of cells" << std::endl;
			return false;
		}

		// remove all but local cell data offsets
		std::vector<std::pair<uint64_t, uint64_t> > data_offsets;
		data_offsets.reserve(this->cells.size());

		for (uint64_t i = 0; i < all_data_offsets.size(); i++) {
			const uint64_t cell = all_data_offsets[i].first;
			if (this->cell_overlaps_local(cell)) {
				data_offsets.push_back(all_data_offsets[i]);
			}
		}
		all_data_offsets.clear();

		// refine the grid
		std::vector<uint64_t> local_cells;
		local_cells.reserve(data_offsets.size());
		for (std::vector<std::pair<uint64_t, uint64_t> >::const_iterator
			item = data_offsets.begin();
			item != data_offsets.end();
			item++
		) {
			local_cells.push_back(item->first);
		}

		if (!this->load(local_cells)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Couldn't load grid"
				<< std::endl;
			abort();
		}
		local_cells.clear();

		#ifdef DEBUG
		if (this->cells.size() != data_offsets.size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of cell data offsets: " << data_offsets.size()
				<< ", should be " << this->cells.size()
				<< std::endl;
			abort();
		}
		#endif

		// create a file view representing local cell data
		MPI_Datatype file_type;

		if (this->cells.size() > 0) {

			// datatypes of local cells
			std::vector<MPI_Datatype> datatypes(this->cells.size());

			// 1 datatype per local cell
			std::vector<int> block_lengths(this->cells.size(), 1);

			// address of local cell data relative to start of file
			std::vector<MPI_Aint> file_displacements(this->cells.size(), 0);

			for (uint64_t i = 0; i < data_offsets.size(); i++) {
				const uint64_t cell = data_offsets[i].first;
				datatypes[i] = this->cells.at(cell).mpi_datatype();
				file_displacements[i] = data_offsets[i].second;
			}

			int result
				= MPI_Type_create_struct(
					datatypes.size(),
					&block_lengths[0],
					&file_displacements[0],
					&datatypes[0],
					&file_type
				);
			if (result != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't create datatype for file view: " << Error_String()(result)
					<< std::endl;
				abort();
			}

			result = MPI_Type_commit(&file_type);
			if (result != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't commit datatype for file view: " << Error_String()(result)
					<< std::endl;
				abort();
			}

		// create an empty file type
		} else {
			MPI_Type_contiguous(0, MPI_BYTE, &file_type);
			if (MPI_Type_commit(&file_type) != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't commit datatype for file view"
					<< std::endl;
				abort();
			}
		}

		// set the file view corresponding to local cell data
		result = MPI_File_set_view(
			infile,
			0,
			MPI_BYTE,
			file_type,
			const_cast<char*>("native"),
			MPI_INFO_NULL
		);

		if (result != MPI_SUCCESS) {
			std::cerr << "Couldn't set file view: "
				<< Error_String()(result)
				<< std::endl;
			abort();
		}

		MPI_Type_free(&file_type);


		// create a datatype representing local cell data
		MPI_Datatype local_type;

		if (this->cells.size() > 0) {

			// datatypes of local cells
			std::vector<MPI_Datatype> datatypes(this->cells.size());

			// 1 datatype per local cell
			std::vector<int> block_lengths(this->cells.size(), 1);

			// address of local cell data relative to data of 1st cell
			std::vector<MPI_Aint> local_displacements(this->cells.size(), 0);

			for (uint64_t i = 0; i < data_offsets.size(); i++) {
				const uint64_t
					cell = data_offsets[i].first,
					first_cell = data_offsets[0].first;

				datatypes[i] = this->cells.at(cell).mpi_datatype();
				local_displacements[i]
					= (uint8_t*) this->cells.at(cell).at()
					- (uint8_t*) this->cells.at(first_cell).at();
			}

			if (MPI_Type_create_struct(
				datatypes.size(),
				&block_lengths[0],
				&local_displacements[0],
				&datatypes[0],
				&local_type) != MPI_SUCCESS
			) {
				std::cerr << "Process " << this->rank
					<< " Couldn't create datatype for local cells"
					<< std::endl;
				abort();
			}

			if (MPI_Type_commit(&local_type) != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't commit datatype for file view"
					<< std::endl;
				abort();
			}

			BOOST_FOREACH(MPI_Datatype& type, datatypes) {
				MPI_Type_free(&type);
			}

		// create an empty type for local data
		} else {
			MPI_Type_contiguous(0, MPI_BYTE, &local_type);
			if (MPI_Type_commit(&local_type) != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't commit datatype for file view"
					<< std::endl;
				abort();
			}
		}

		// read local cell data from file
		if (this->cells.size() > 0) {
			result = MPI_File_read_at_all(
				infile,
				0,
				this->cells.at(data_offsets[0].first).at(),
				1,
				local_type,
				MPI_STATUS_IGNORE
			);
		} else {
			uint64_t temp;
			result = MPI_File_read_at_all(
				infile,
				0,
				&temp,
				1,
				local_type,
				MPI_STATUS_IGNORE
			);
		}

		if (result != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " couldn't read local cell data: "
				<< Error_String()(result)
				<< std::endl;
			abort();
		}

		MPI_Type_free(&local_type);

		MPI_File_close(&infile);

		#endif
		#endif

		return true;
	}


	/*!
	Returns a begin const_iterator to the internal storage of local cells and their data.
	*/
	typename boost::unordered_map<uint64_t, UserData>::const_iterator begin() const
	{
		return this->cells.begin();
	}

	/*!
	Returns an end const_iterator to the internal storage of local cells and their data.
	*/
	typename boost::unordered_map<uint64_t, UserData>::const_iterator end() const
	{
		return this->cells.end();
	}


	/*!
	Returns the number of local cells without children, e.g. leaf cells.
	*/
	size_t size() const
	{
		return this->cells.size();
	}


	/*!
	Returns local cells without remote neighbors.

	Returns cell in this process that don't have children (e.g. leaf cells)
	and don't have neighbors on other processes.
	*/
	std::vector<uint64_t> get_cells_with_local_neighbors() const
	{
		std::vector<uint64_t> return_cells;
		return_cells.reserve(this->cells.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			uint64_t child = this->get_child(item.first);
			assert(child > 0);

			if (child != item.first) {
				continue;
			}

			bool has_remote_neighbor = false;

			assert(this->neighbors.count(item.first) > 0);

			BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(item.first)) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					has_remote_neighbor = true;
					break;
				}
			}

			if (!has_remote_neighbor) {
				return_cells.push_back(item.first);
			}
		}

		return return_cells;
	}

	/*!
	Returns local cells without remote neighbors for given cell neighborhood.

	Returns nothing if given neighborhood doesn't exist.
	*/
	std::vector<uint64_t> get_cells_with_local_neighbors(const int neighborhood_id) const
	{
		std::vector<uint64_t> return_cells;

		if (this->user_neigh_of.count(neighborhood_id) == 0) {
			return return_cells;
		}

		return_cells.reserve(this->cells.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			const uint64_t
				cell = item.first,
				child = this->get_child(cell);

			if (child != cell) {
				continue;
			}

			#ifdef DEBUG
			if (child == error_cell) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid child for cell " << cell
					<< std::endl;
				abort();
			}

			if (this->user_neigh_of.at(neighborhood_id).count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No neighbors for cell " << cell
					<< " in neighborhood " << neighborhood_id
					<< std::endl;
				abort();
			}
			#endif

			bool has_remote_neighbor = false;

			BOOST_FOREACH(
				const uint64_t neighbor,
				this->user_neigh_of.at(neighborhood_id).at(cell)
			) {

				if (neighbor == error_cell) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					has_remote_neighbor = true;
					break;
				}
			}

			if (!has_remote_neighbor) {
				return_cells.push_back(cell);
			}
		}

		return return_cells;
	}


	/*!
	Returns local cells with at least one remote neighbor.

	Returns all cells on this process that don't have children (e.g. leaf cells)
	and have at least one neighbor on another processes.
	*/
	std::vector<uint64_t> get_cells_with_remote_neighbor() const
	{
		std::vector<uint64_t> return_cells;
		return_cells.reserve(this->cells.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			uint64_t child = this->get_child(item.first);
			assert(child > 0);

			if (child != item.first) {
				continue;
			}

			bool has_remote_neighbor = false;

			assert(this->neighbors.count(item.first) > 0);

			BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(item.first)) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					has_remote_neighbor = true;
					break;
				}
			}

			if (has_remote_neighbor) {
				return_cells.push_back(item.first);
			}
		}

		return return_cells;
	}

	/*!
	Returns local cells with at least one remote neighbor in given cell neighborhood

	Returns nothing if given neighborhood doesn't exist.
	*/
	std::vector<uint64_t> get_cells_with_remote_neighbor(const int neighborhood_id) const
	{
		std::vector<uint64_t> return_cells;

		if (this->user_neigh_of.count(neighborhood_id) == 0) {
			return return_cells;
		}

		return_cells.reserve(this->cells.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			const uint64_t
				cell = item.first,
				child = this->get_child(cell);

			if (child != cell) {
				continue;
			}

			#ifdef DEBUG
			if (child == error_cell) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid child for cell " << cell
					<< std::endl;
				abort();
			}

			if (this->user_neigh_of.at(neighborhood_id).count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No neighbors for cell " << cell
					<< " in neighborhood " << neighborhood_id
					<< std::endl;
				abort();
			}
			#endif

			bool has_remote_neighbor = false;

			BOOST_FOREACH(
				const uint64_t& neighbor,
				this->user_neigh_of.at(neighborhood_id).at(cell)
			) {
				if (neighbor == error_cell) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					has_remote_neighbor = true;
					break;
				}
			}

			if (has_remote_neighbor) {
				return_cells.push_back(cell);
			}
		}

		return return_cells;
	}


	/*!
	Returns all cells in the grid that don't have children (e.g. leaf cells)
	*/
	std::vector<uint64_t> get_all_cells() const
	{
		std::vector<uint64_t> all_cells;
		all_cells.reserve(this->cell_process.size());

		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			item = this->cell_process.begin();
			item != this->cell_process.end();
			item++
		) {

			const uint64_t child = this->get_child(item->first);

			if (child == item->first) {
				all_cells.push_back(item->first);
			}
		}

		return all_cells;
	}


	/*!
	Returns a pointer to the user supplied data of given cell.

	The data of local cells is always available, including refined cells
	before the next call to stop_refining.
	The data of cells which are on other processes can also be available if:
		- the cells are neighbors to a local cell and remote neighbor data has been updated
		- the cells were unrefined and their parent is now a local cell
	*/
	UserData* operator [] (const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			return (UserData*) &(this->cells.at(cell));
		} else if (this->remote_neighbors.count(cell) > 0) {
			return (UserData*) &(this->remote_neighbors.at(cell));
		} else if (this->refined_cell_data.count(cell) > 0) {
			return (UserData*) &(this->refined_cell_data.at(cell));
		} else if (this->unrefined_cell_data.count(cell) > 0) {
			return (UserData*) &(this->unrefined_cell_data.at(cell));
		} else {
			return NULL;
		}
	}


	/*!
	Returns true if given cell overlaps a local cell.

	Returns true if given cell is either a local cell of
	refinement level 0 or is a (grandgrand...)child of a
	local cell with refinement level 0.
	Returns false otherwise.
	*/
	bool cell_overlaps_local(const uint64_t cell)
	{
		const uint64_t parent = this->get_level_0_parent(cell);

		if (this->cell_process.count(parent) > 0
		&& this->cell_process.at(parent) == this->rank) {
			return true;
		} else {
			return false;
		}
	}


	/*!
	Returns those cells from given that are either local or a child of local.

	Assumes local cells are of refinement level 0.
	*/
	boost::unordered_set<uint64_t> get_cells_overlapping_local(const std::vector<uint64_t>& given_cells)
	{
		boost::unordered_set<uint64_t> result;

		BOOST_FOREACH(const uint64_t cell, given_cells) {

			if (this->cell_overlaps_local(cell)) {
				result.insert(cell);
			}
		}

		return result;
	}


	/*!
	Refines the current grid so that the given cells exist.

	Must be called by all processes and only cells of
	refinement level 0 must exist in the grid at that time.
	Ignores cells in given list that aren't a (grand...)
	child of a local cell.
	Returns true on success and false otherwise.
	*/
	bool load(const std::vector<uint64_t>& given_cells)
	{
		boost::unordered_set<uint64_t> overlapping
			= this->get_cells_overlapping_local(given_cells);

		int max_ref_lvl_of_overlapping = 0;
		BOOST_FOREACH(const uint64_t cell, overlapping) {
			const int refinement_level = this->get_refinement_level(cell);
			max_ref_lvl_of_overlapping
				= std::max(refinement_level, max_ref_lvl_of_overlapping);
		}

		/*
		Starting from refinement level 0 refine each local cell
		that has a child in given_cells
		*/
		std::vector<boost::unordered_set<uint64_t> > cells_and_parents(max_ref_lvl_of_overlapping);

		BOOST_FOREACH(const uint64_t cell, overlapping) {
			uint64_t current_child = cell;
			const int refinement_level = this->get_refinement_level(current_child);
			for (int i = refinement_level - 1; i >= 0; i--) {
				const uint64_t parent = this->get_parent_for_removed(current_child);
				cells_and_parents[i].insert(parent);
				current_child = parent;
			}
		}

		for (int current_ref_lvl = 0; current_ref_lvl < max_ref_lvl_of_overlapping; current_ref_lvl++) {
			BOOST_FOREACH(const uint64_t cell, cells_and_parents[current_ref_lvl]) {
				this->refine_completely(cell);
			}
			this->stop_refining();
			this->clear_refined_unrefined_data();
		}

		return true;
	}


	/*!
	Load balances the grid's cells among processes.

	Must be called simultaneously on all processes.
	Cells which haven't been pinned are moved as suggested by Zoltan,
	pinned cells are moved as requested by the user.
	Does not update remote neighbor data between processes afterward.
	Discards refines / unrefines.

	If prepare_to_balance_load was called before this then has_been_prepared must be true
	otherwise it must be false.
	*/
	void balance_load(const bool has_been_prepared = false)
	{
		if (!has_been_prepared) {
			this->make_new_partition(true);
		}
		this->move_cells();
		this->added_cells.clear();
		this->removed_cells.clear();
	}

	/*!
	Moves pinned grid cells as requested by the user.

	Must be called simultaneously on all processes.
	Cells which haven't been pinned are not moved.
	Does not update remote neighbor data between processes afterward.
	Discards refines / unrefines.

	If prepare_to_migrate_cells was called before this then
	has_been_prepared must be true otherwise it must be false.
	*/
	void migrate_cells(const bool has_been_prepared = false)
	{
		if (!has_been_prepared) {
			this->make_new_partition(false);
		}
		this->move_cells();
		this->added_cells.clear();
		this->removed_cells.clear();
	}


	/*!
	Same as balance_load but only prepares to move cells with balance_load.

	Must be used when cells contain variable mpi_datatypes so that when
	cells are moved receiving processes can construct the receiving
	datatypes in balance_load based on cells data transferred by this function.
	The next dccrg function to be called after this one must be balance_load.

	When calling balance_load after this function has_been_prepared must be true.
	*/
	void prepare_to_balance_load()
	{
		this->make_new_partition(true);
		this->prepare_to_move_cells();
	}

	/*!
	Same as migrate_cells but only prepares to move cells with migrate_cells.

	Must be used when cells contain variable mpi_datatypes so that when
	cells are moved receiving processes can construct the receiving
	datatypes in migrate_cells based on cells data transferred by this function.
	The next dccrg function to be called after this one must be migrate_cells.

	When calling migrate_cells after this function has_been_prepared must be true.
	*/
	void prepare_to_migrate_cells()
	{
		this->make_new_partition(true);
		this->prepare_to_move_cells();
	}


	/*!
	Updates the user data of neighboring cells between processes.

	User data of any local cell which is considered as a neighbor
	to a cell on another process is sent to that process.
	User data of any cell that is considered as a neighbor of a
	cell on this process is received from the other process.
	Cells' user data is only sent to / received from a process once.
	Must be called simultaneously on all processes
	*/
	void update_remote_neighbor_data()
	{
		this->start_remote_neighbor_data_update();
		this->wait_neighbor_data_update();
	}

	/*!
	Same as the version without id but uses the given neighborhood.

	Must be called with the same id by all processes.
	Can be used to transfer the cell data of less number of cells
	between processes than the full neighborhood would.
	add_remote_update_neighborhood must be used once prior to
	calling this with the same id in order to create the reduced
	neighborhood.
	Does nothing in case a neighborhood with given id doesn't exist.
	*/
	void update_remote_neighbor_data(const int id)
	{
		this->start_remote_neighbor_data_update(id);
		this->wait_neighbor_data_update(id);
	}


	/*!
	Starts the update of neighbor data between processes, returns immediately.

	Must be called simultaneously on all processes
	*/
	void start_remote_neighbor_data_update()
	{
		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->remote_neighbors, this->cells_to_receive, this->cells_to_send
		#elif defined (DCCRG_CELL_DATA_SIZE_FROM_USER)
		this->remote_neighbors, this->cells_to_receive, this->cells_to_send
		#else
		this->cells_to_receive, this->cells_to_send
		#endif
		);
	}

	/*!
	Same as the version without id but uses the given neighborhood.

	add_remote_update_neighborhood must been used once with the
	same id prior to calling this.
	Does nothing in case a neighborhood with given id doesn't exist.
	*/
	void start_remote_neighbor_data_update(const int id)
	{
		if (this->user_hood_of.count(id) == 0) {

			#ifdef DEBUG
			if (this->user_hood_to.count(id) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should not have id " << id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_of.count(id) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should not have id " << id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_to.count(id) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should not have id " << id
					<< std::endl;
				abort();
			}
			#endif

			return;
		}

		#ifdef DEBUG
		if (this->user_hood_to.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_of.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_to.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_send.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_receive.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << id
				<< std::endl;
			abort();
		}
		#endif

		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->remote_neighbors,
		this->user_neigh_cells_to_receive.at(id),
		this->user_neigh_cells_to_send.at(id)
		#elif defined (DCCRG_CELL_DATA_SIZE_FROM_USER)
		this->remote_neighbors,
		this->user_neigh_cells_to_receive.at(id),
		this->user_neigh_cells_to_send.at(id)
		#else
		this->user_neigh_cells_to_receive.at(id),
		this->user_neigh_cells_to_send.at(id)
		#endif
		);
	}


	/*!
	Waits for remote neighbor data transfers to/from this process to complete.

	Must be called simultaneously on all processes.
	*/
	void wait_neighbor_data_update()
	{
		this->wait_neighbor_data_update_receives();
		this->wait_neighbor_data_update_sends();
	}

	/*!
	Same as the version without id but uses the given neighborhood.

	Does nothing if neighborhood with given id doesn't exist.
	*/
	void wait_neighbor_data_update(const int id)
	{
		this->wait_neighbor_data_update_receives(id);
		this->wait_neighbor_data_update_sends();
	}


	/*!
	Waits for remote neighbor data transfers to other processes to complete.

	Waits until all sends associated with neighbor data update transfers
	between processes have completed.
	Must be called simultaneously on all processes and probably must be
	called after wait...update_receives().
	*/
	void wait_neighbor_data_update_sends()
	{
		this->wait_user_data_transfer_sends();
	}


	/*!
	Waits for remote neighbor data transfers from other processes to complete.

	Waits until all receives associated with neighbor data update transfers
	between processes have completed and incorporates that data.
	Must be called simultaneously on all processes and probably must be
	called before wait...update_sends().
	*/
	void wait_neighbor_data_update_receives()
	{
		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->remote_neighbors, this->cells_to_receive
		#endif
		#endif
		);
	}

	/*!
	Same as the version without id but uses the given neighborhood.

	Does nothing if neighborhood with given id doesn't exist.
	*/
	void wait_neighbor_data_update_receives(const int id)
	{
		if (this->user_hood_of.count(id) == 0) {
			return;
		}

		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->remote_neighbors, this->user_neigh_cells_to_receive.at(id)
		#endif
		#endif
		);
	}


	/*!
	Returns the number of cells whose data this process has to send during a neighbor data update.

	The total amount of cells to be sent is returned so if a cell's data will be sent to
	N processes it is counted N times.
	*/
	uint64_t get_number_of_update_send_cells() const
	{
		uint64_t result = 0;
		for (
			#ifdef DCCRG_SEND_SINGLE_CELLS
			boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			#else
			boost::unordered_map<int, std::vector<uint64_t> >::const_iterator
			#endif
			receiver = cells_to_send.begin();
			receiver != cells_to_send.end();
			receiver++
		) {
			result += receiver->second.size();
		}
		return result;
	}

	/*!
	Returns the number of cells whose data this process has to receive during a neighbor data update.
	*/
	uint64_t get_number_of_update_receive_cells() const
	{
		uint64_t result = 0;
		for (
			#ifdef DCCRG_SEND_SINGLE_CELLS
			boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			#else
			boost::unordered_map<int, std::vector<uint64_t> >::const_iterator
			#endif
			sender = cells_to_receive.begin();
			sender != cells_to_receive.end();
			sender++
		) {
			result += sender->second.size();
		}
		return result;
	}


	/*!
	Returns a pointer to the neighbors of given cell.

	In case the grid is not periodic in one or more directions,
	neighbors that would be outside of the grid are error_cell.
	Some neighbors might be on another process, but have a copy of their data on this process.
	The local copy of remote neighbors' data is updated, for example, by calling
	update_remote_neighbor_data().
	Returns NULL if given cell doesn't exist or is on another process.

	The neighbors are always in the following order:
		- if all neighbors are of the same size then they are in z order, e.g.
		  with a neighborhood size of 2 the first neighbor is at offset (-2, -2, -2)
		  from the given cell, the second one is at (-1, -2, -2), etc, in size units
		  of the given cell.
		- if one or more of the cells in 1) is refined then instead of one cell
		  there are 8 which are again in z order.
	For example with maximum refinement level 1 and neighborhood size of 1
	the neighbors of a cell of refinement level 0 at indices (2, 2, 2) could
	be in the following order: (0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0),
	rest of refined cells..., (2, 0, 0), (3, 0, 0), (0, 2, 0), (2, 2, 0), ...

	Offset (0, 0, 0) is skipped in all neighbor lists, so with a neighborhood
	size of 2 the minimum length of neighbors lists is 124 and not 5^3 = 125.
	*/
	const std::vector<uint64_t>* get_neighbors(const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			#ifdef DEBUG
			if (this->neighbors.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Neighbor list for cell " << cell
					<< " doesn't exist"
					<< std::endl;
				abort();
			}
			#endif
			return &(this->neighbors.at(cell));
		} else {
			return NULL;
		}
	}

	/*!
	Same as the function without id but returns neighbors for the given neighborhood.

	Neighbors are in the same order as the offsets in the neighborhood with given id.
	*/
	const std::vector<uint64_t>* get_neighbors(const uint64_t cell, const int id) const
	{
		if (this->cells.count(cell) > 0) {
			#ifdef DEBUG
			if (this->user_neigh_of.at(id).count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Neighbor list for cell " << cell
					<< " doesn't exist for neighborhood id " << id
					<< std::endl;
				abort();
			}
			#endif
			return &(this->user_neigh_of.at(id).at(cell));
		} else {
			return NULL;
		}
	}


	/*!
	Returns a pointer to the cells that consider given cell as a neighbor.

	This list doesn't include 0s even if the grid isn't periodic in some direction.
	Returns NULL if given cell doesn't exist or is on another process.
	Neighbors returned by this function are in no particular order.
	*/
	const std::vector<uint64_t>* get_neighbors2(const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			#ifdef DEBUG
			if (this->neighbors_to.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbors_to list for cell " << cell
					<< " doesn't exist"
					<< std::endl;
				abort();
			}
			#endif
			return &(this->neighbors_to.at(cell));
		} else {
			return NULL;
		}
	}

	/*!
	Same as the functions without id but with respect to the given neighborhood.
	*/
	const std::vector<uint64_t>* get_neighbors2(const uint64_t cell, const int id) const
	{
		if (this->cells.count(cell) > 0) {
			#ifdef DEBUG
			if (this->user_neigh_to.at(id).count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbors_to list for cell " << cell
					<< " doesn't exist for neighborhood id " << id
					<< std::endl;
				abort();
			}
			#endif
			return &(this->user_neigh_to.at(id).at(cell));
		} else {
			return NULL;
		}
	}


	/*!
	Returns the size of cells' neihgbourhood in every direction.
	*/
	unsigned int get_neighborhood_size() const
	{
		return this->neighborhood_size;
	}


	/*!
	Returns all neighbors of given cell that are at given offsets from it.

	Returns nothing in the following cases:
		- given cell doesn't exist
		- given cell is on another process
		- any of given offsets is larger in absolute value than the neighborhood
		  size or larger than 1 if neihgborhood size == 0
		- i == 0 && j == 0 && k == 0
	*/
	std::vector<uint64_t> get_neighbors_of(const uint64_t cell, const int i, const int j, const int k) const
	{
		std::vector<uint64_t> return_neighbors;
		if (this->cell_process.count(cell) == 0
		|| this->cell_process.at(cell) != this->rank
		|| (i == 0 && j == 0 && k == 0)) {
			return return_neighbors;
		}

		const int refinement_level = this->get_refinement_level(cell);

		// find cell(s) at given indices in the stored neighbor list
		const int last_offset = (this->neighborhood_size > 0) ? int(this->neighborhood_size) : 1;
		int index = 0;
		for (int
			current_k = (this->neighborhood_size > 0) ? -int(this->neighborhood_size) : -1;
			current_k <= last_offset;
			current_k++
		)
		for (int
			current_j = (this->neighborhood_size > 0) ? -int(this->neighborhood_size) : -1;
			current_j <= last_offset;
			current_j++
		)
		for (int
			current_i = (this->neighborhood_size > 0) ? -int(this->neighborhood_size) : -1;
			current_i <= last_offset;
			current_i++
		) {
			if (current_i == 0 && current_j == 0 && current_k == 0) {
				continue;
			}

			if (this->neighborhood_size == 0) {
				// skip diagonal offsets
				const int zero_offsets_in_current =
					  ((current_i == 0) ? 1 : 0)
					+ ((current_j == 0) ? 1 : 0)
					+ ((current_k == 0) ? 1 : 0);
				if (zero_offsets_in_current != 2) {
					continue;
				}
			}

			const int current_refinement_level = this->get_refinement_level(this->neighbors.at(cell)[index]);
			if (i == current_i && j == current_j && k == current_k) {

				// TODO check for 0 neighbor instead of error from get_refinement_level
				if (current_refinement_level == -1) {
					return_neighbors.push_back(0);
				} else {
					return_neighbors.push_back(this->neighbors.at(cell)[index]);

					if (current_refinement_level > refinement_level) {
						return_neighbors.reserve(8);
						for (int i = 1; i < 8; i++) {
							index++;
							return_neighbors.push_back(this->neighbors.at(cell)[index]);
						}
					}
				}

				current_i = current_j = current_k = last_offset + 1;

			} else {
				if (current_refinement_level > refinement_level) {
					index += 7;
				}
			}

			index++;
		}

		return return_neighbors;
	}


	/*!
	Returns the given cell's neighbors that are on another process.

	Returns nothing if given cell doesn't exist or is on another process
	or doesn't have remote neighbors.
	*/
	std::vector<uint64_t> get_remote_neighbors(const uint64_t cell) const
	{
		std::vector<uint64_t> result;

		if (this->cells.count(cell) == 0
		|| this->neighbors.count(cell) == 0) {
			return result;
		}

		BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(cell)) {

			if (neighbor == 0) {
				continue;
			}

			if (this->cell_process.at(neighbor) != this->rank) {
				result.push_back(neighbor);
			}
		}

		return result;
	}


	/*!
	Returns remote neighbors of local cells.

	Returns cells on other processes which at least one local cell
	considers as neighbors.
	Use update_remote_neighbor_data to make sure that a local copy
	of remote neighbors' data exists.
	*/
	std::vector<uint64_t> get_remote_neighbors() const
	{
		return this->get_list_of_remote_cells_with_local_neighbors();
	}


	/*!
	Returns true if given cell is on this process and false otherwise.
	*/
	bool is_local(const uint64_t cell) const
	{
		if (this->cell_process.count(cell) > 0
		&& this->cell_process.at(cell) == this->rank) {
			return true;
		} else {
			return false;
		}
	}


	/*!
	Writes the cells on this process into a vtk file with given name in ASCII format.

	The cells are written in ascending order.
	Must be called simultaneously on all processes.
	*/
	void write_vtk_file(const char* file_name) const
	{
		std::ofstream outfile(file_name);
		if (!outfile.is_open()) {
			std::cerr << "Couldn't open file " << file_name << std::endl;
			// TODO: throw an exception instead
			exit(1);
		}

		std::vector<uint64_t> leaf_cells = this->get_cells();
		std::sort(leaf_cells.begin(), leaf_cells.end());
		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << "Cartesian cell refinable grid" << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;

		// write separate points for every cells corners
		outfile << "POINTS " << leaf_cells.size() * 8 << " float" << std::endl;
		for (unsigned int i = 0; i < leaf_cells.size(); i++) {
			outfile << this->get_cell_x_min(leaf_cells[i]) << " "
				<< this->get_cell_y_min(leaf_cells[i]) << " "
				<< this->get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_max(leaf_cells[i]) << " "
				<< this->get_cell_y_min(leaf_cells[i]) << " "
				<< this->get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_min(leaf_cells[i]) << " "
				<< this->get_cell_y_max(leaf_cells[i]) << " "
				<< this->get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_max(leaf_cells[i]) << " "
				<< this->get_cell_y_max(leaf_cells[i]) << " "
				<< this->get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_min(leaf_cells[i]) << " "
				<< this->get_cell_y_min(leaf_cells[i]) << " "
				<< this->get_cell_z_max(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_max(leaf_cells[i]) << " "
				<< this->get_cell_y_min(leaf_cells[i]) << " "
				<< this->get_cell_z_max(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_min(leaf_cells[i]) << " "
				<< this->get_cell_y_max(leaf_cells[i]) << " "
				<< this->get_cell_z_max(leaf_cells[i]) << std::endl;
			outfile << this->get_cell_x_max(leaf_cells[i]) << " "
				<< this->get_cell_y_max(leaf_cells[i]) << " "
				<< this->get_cell_z_max(leaf_cells[i]) << std::endl;
		}

		// map cells to written points
		outfile << "CELLS " << leaf_cells.size() << " " << leaf_cells.size() * 9 << std::endl;
		for (unsigned int j = 0; j < leaf_cells.size(); j++) {
			outfile << "8 ";
			for (int i = 0; i < 8; i++) {
				 outfile << j * 8 + i << " ";
			}
			outfile << std::endl;
		}

		// cell types
		outfile << "CELL_TYPES " << leaf_cells.size() << std::endl;
		for (unsigned int i = 0; i < leaf_cells.size(); i++) {
			outfile << 11 << std::endl;
		}

		if (!outfile.good()) {
			std::cerr << "Writing of vtk file probably failed" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}

		outfile.close();
	}


	/*!
	Creates all children of given cell (and possibly of other cells due to induced refinement).

	Takes priority over unrefining.
	Refines / unrefines take effect only after a call to stop_refining() and are lost
	after a call to balance_load().
	Does nothing in any of the following cases:
		- given cell has already been refined (including induced refinement)
		  and stop_refining() has not been called afterwards
		- given cell doesn't exist on this process
		- given cell's children already exist
	Children are created on their parent's process.

	If given cell is at maximum refinement level dont_unrefine will be invoked instead.
	 */
	void refine_completely(const uint64_t cell)
	{
		if (cell == error_cell) {
			return;
		}

		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cells.count(cell) == 0) {
			return;
		}

		const int refinement_level = this->get_refinement_level(cell);

		if (refinement_level > this->max_refinement_level) {
			return;
		}

		// not if cell has children
		if (cell != this->get_child(cell)) {
			return;
		}

		if (refinement_level == this->max_refinement_level) {
			this->dont_unrefine(cell);
			return;
		}

		this->cells_to_refine.insert(cell);

		// override local unrefines
		const std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(cell));
		BOOST_FOREACH(const uint64_t& sibling, siblings) {
			this->cells_to_unrefine.erase(sibling);
		}

		BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(cell)) {
			if (this->get_refinement_level(neighbor) <= refinement_level) {
				const std::vector<uint64_t> neighbor_siblings = this->get_all_children(this->get_parent(neighbor));
				BOOST_FOREACH(const uint64_t& sibling, neighbor_siblings) {
					this->cells_to_unrefine.erase(sibling);
				}
			}
		}

		BOOST_FOREACH(const uint64_t& neighbor, this->neighbors_to.at(cell)) {
			if (this->get_refinement_level(neighbor) <= refinement_level) {
				const std::vector<uint64_t> neighbor_siblings = this->get_all_children(this->get_parent(neighbor));
				BOOST_FOREACH(const uint64_t& sibling, neighbor_siblings) {
					this->cells_to_unrefine.erase(sibling);
				}
			}
		}
	}

	/*!
	As refine_completely, but uses the smallest existing cell at given coordinates.

	Does nothing in the same cases as refine_completely and additionally
	if the coordinate is outside of the grid.
	*/
	void refine_completely_at(const double x, const double y, const double z)
	{
		const uint64_t cell = this->get_existing_cell(x, y, z);
		if (cell == 0) {
			return;
		}

		this->refine_completely(cell);
	}


	/*!
	Removes the given cell and its siblings from the grid.

	Refining (including induced refining) takes priority over unrefining.
	Refines / unrefines take effect only after a call to stop_refining()
	and are lost after a call to balance_load().
	Does nothing in any of the following cases:
		- dont_unrefine was called previously for given cell or its siblings
		- given cell or one of its siblings has already been unrefined
		  and stop_refining() has not been called
		- given cell doesn't exist on this process
		- given cell has children
		- given cells refinement level is 0

	After a cell and its siblings have been unrefined, their data has been moved
	to their parent's process.
	When no longer needed that data can be freed using clear_refined_unrefined_data.
	*/
	void unrefine_completely(const uint64_t cell)
	{
		if (cell == error_cell) {
			return;
		}

		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cells.count(cell) == 0) {
			return;
		}

		if (this->get_refinement_level(cell) == 0) {
			return;
		}

		const std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(cell));

		// don't unrefine if any sibling...
		BOOST_FOREACH(const uint64_t& sibling, siblings) {

			// ...has children
			if (sibling != this->get_child(sibling)) {
				return;
			}

			// ...cannot be unrefined
			if (this->cells_to_refine.count(sibling) > 0
			|| this->cells_not_to_unrefine.count(sibling) > 0) {
				return;
			}
		}

		// unrefinement succeeds if parent of unrefined will fulfill requirements
		const uint64_t parent = this->get_parent(cell);
		const int refinement_level = this->get_refinement_level(parent);

		#ifdef DEBUG
		if (parent == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent" << std::endl;
			abort();
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid refinement level for parent" << std::endl;
			abort();
		}
		#endif

		const std::vector<uint64_t> neighbors
			= this->find_neighbors_of(
				parent,
				this->neighborhood_of,
				2 * this->max_ref_lvl_diff,
				true
			);

		BOOST_FOREACH(const uint64_t& neighbor, neighbors) {

			const int neighbor_ref_lvl = this->get_refinement_level(neighbor);

			if (neighbor_ref_lvl > refinement_level + this->max_ref_lvl_diff) {
				return;
			}

			if (neighbor_ref_lvl == refinement_level + this->max_ref_lvl_diff
			&& this->cells_to_refine.count(neighbor) > 0) {
				return;
			}
		}

		// record only one sibling to unrefine / process
		BOOST_FOREACH(const uint64_t& sibling, siblings) {
			if (this->cells_to_unrefine.count(sibling) > 0) {
				return;
			}
		}

		this->cells_to_unrefine.insert(cell);
	}


	/*!
	As unrefine_completely, but uses the smallest existing cell at given coordinates.

	Does nothing in the same cases as unrefine_completely and additionally
	if the coordinate is outside of the grid.
	*/
	void unrefine_completely_at(const double x, const double y, const double z)
	{
		const uint64_t cell = this->get_existing_cell(x, y, z);
		if (cell == 0) {
			return;
		}

		this->unrefine_completely(cell);
	}


	/*!
	Prevents the given cell or its siblings from being unrefined.

	Has an effect only during the next call to stop_refining().
	Has no effect if balance_load() is called before stop_refining().
	Does nothing in any of the following cases:
		- given cell doesn't exist on this process
		- given cell has children
		- given cell's refinement level is 0
	*/
	void dont_unrefine(const uint64_t cell)
	{
		if (cell == error_cell) {
			return;
		}

		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cells.count(cell) == 0) {
			return;
		}

		if (this->get_refinement_level(cell) == 0) {
			return;
		}

		if (cell != this->get_child(cell)) {
			// cell already has children
			return;
		}

		// record only one sibling / process
		const std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(cell));
		BOOST_FOREACH(const uint64_t& sibling, siblings) {
			if (this->cells_not_to_unrefine.count(sibling) > 0) {
				return;
			}
		}

		// override local unrefines
		BOOST_FOREACH(const uint64_t& sibling, siblings) {
			this->cells_to_unrefine.erase(sibling);
		}

		this->cells_not_to_unrefine.insert(cell);
	}


	/*!
	As dont_unrefine but uses the smallest existing cell at given coordinates.

	Does nothing in the same cases as dont_unrefine and additionally if the
	coordinate is outside of the grid.
	*/
	void dont_unrefine_at(const double x, const double y, const double z)
	{
		const uint64_t cell = this->get_existing_cell(x, y, z);
		if (cell == error_cell) {
			return;
		}

		this->dont_unrefine(cell);
	}


	/*!
	Executes refines / unrefines that have been requested so far.

	Must be called simultaneously on all processes.
	Returns cells that were created by refinement on this process.
	Moves user data of unrefined cells to the process of their parent.
	*/
	std::vector<uint64_t> stop_refining()
	{
		this->induce_refines();

		// update dont_refines between processes
		this->all_to_all_set(this->cells_not_to_unrefine);

		this->override_unrefines();
		this->cells_not_to_unrefine.clear();

		return this->execute_refines();
	}


	/*!
	Returns cells that were removed by unrefinement whose parent is on this process
	Removed cells data is also on this process, but only until balance_load() is called
	*/
	std::vector<uint64_t> get_removed_cells() const
	{
		std::vector<uint64_t> unref_removed_cells;
		unref_removed_cells.reserve(this->unrefined_cell_data.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->unrefined_cell_data) {
			unref_removed_cells.push_back(item.first);
		}

		return unref_removed_cells;
	}


	/*!
	Given a cell that exists and has a parent returns the parent cell
	Returns the given cell if it doesn't have a parent or 0 if given cell doesn't exist
	*/
	uint64_t get_parent(const uint64_t cell) const
	{
		if (this->cell_process.count(cell) == 0) {
			return 0;
		}

		// given cell cannot have a parent
		if (this->get_refinement_level(cell) == 0) {
			return cell;
		}

		const uint64_t parent = this->get_cell_from_indices(
			this->get_indices(cell),
			this->get_refinement_level(cell) - 1
		);

		if (this->cell_process.count(parent) > 0) {
			return parent;
		} else {
			return cell;
		}
	}


	/*!
	Returns the indices corresponding to the given neighborhood at given indices.

	Neighborhood is returned in units of size_in_indices.
	If grid is not periodic in some direction then indices which would fall outside
	of the grid in that direction are returned as error_index.
	*/
	std::vector<Types<3>::indices_t> indices_from_neighborhood(
		const Types<3>::indices_t indices,
		const uint64_t size_in_indices,
		const std::vector<Types<3>::neighborhood_item_t>& neighborhood
	) const {
		// TODO: make neighborhood a const reference
		std::vector<Types<3>::indices_t> return_indices;
		return_indices.reserve(neighborhood.size());

		// grid length in indices
		const uint64_t grid_length[3] = {
			this->get_x_length() * (uint64_t(1) << this->max_refinement_level),
			this->get_y_length() * (uint64_t(1) << this->max_refinement_level),
			this->get_z_length() * (uint64_t(1) << this->max_refinement_level)
		};

		#ifdef DEBUG
		for (unsigned int dimension = 0; dimension < 3; dimension++) {
			if (indices[dimension] >= grid_length[dimension]) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Given indices outside of the grid in dimension " << dimension
					<< std::endl;
				abort();
			}
		}
		#endif

		BOOST_FOREACH(const Types<3>::neighborhood_item_t& offsets, neighborhood) {

			Types<3>::indices_t temp_indices = {{indices[0], indices[1], indices[2]}};

			for (unsigned int dimension = 0; dimension < 3; dimension++) {
				if (offsets[dimension] < 0) {

					if (this->is_periodic(dimension)) {

						// neighborhood might wrap around the grid several times
						for (int i = 0; i > offsets[dimension]; i--) {

							#ifdef DEBUG
							if (temp_indices[dimension] < size_in_indices - 1
							&& temp_indices[dimension] > 0) {
								std::cerr << __FILE__ << ":" << __LINE__
									<< " Cells aren't supposed to wrap around the grid."
									<< std::endl;
								abort();
							}
							#endif

							if (temp_indices[dimension] >= size_in_indices) {
								temp_indices[dimension] -= size_in_indices;
							} else {
								temp_indices[dimension] = grid_length[dimension] - size_in_indices;
							}
						}
					// use error_indices to signal that this neighborhood item is outside of the grid
					} else {
						if (indices[dimension] < abs(offsets[dimension]) * size_in_indices) {
							temp_indices[0] = error_index;
							temp_indices[1] = error_index;
							temp_indices[2] = error_index;
							break;
						}

						temp_indices[dimension] += offsets[dimension] * size_in_indices;
					}

				} else {

					if (this->is_periodic(dimension)) {
						for (int i = 0; i < offsets[dimension]; i++) {

							#ifdef DEBUG
							if (temp_indices[dimension] > grid_length[dimension] - size_in_indices) {
								std::cerr << __FILE__ << ":" << __LINE__ << " Cells aren't supposed to wrap around the grid." << std::endl;
								abort();
							}
							#endif

							if (temp_indices[dimension] < grid_length[dimension] - size_in_indices) {
								temp_indices[dimension] += size_in_indices;
							} else {
								temp_indices[dimension] = 0;
							}
						}
					} else {
						if (indices[dimension] + offsets[dimension] * size_in_indices >= grid_length[dimension]) {
							temp_indices[0] = error_index;
							temp_indices[1] = error_index;
							temp_indices[2] = error_index;
							break;
						}

						temp_indices[dimension] += offsets[dimension] * size_in_indices;
					}
				}
			}

			return_indices.push_back(temp_indices);
		}

		return return_indices;
	}


	/*!
	Returns the existing neighbors (that don't have children) of given cell.

	Uses given neighborhood when searching for neighbors.
	max_diff is the distance to search in refinement level from given cell inclusive.
	Returns nothing if the following cases:
		- given cell has children
		- given doesn't exist
	*/
	// TODO: make private?
	std::vector<uint64_t> find_neighbors_of(
		const uint64_t cell,
		const std::vector<Types<3>::neighborhood_item_t>& neighborhood,
		const int max_diff,
		const bool has_children = false
	) const {
		std::vector<uint64_t> return_neighbors;

		const int refinement_level = this->get_refinement_level(cell);

		#ifdef DEBUG
		if (max_diff < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " max_diff must not be negative" << std::endl;
			abort();
		}

		if (cell == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given: " << cell << std::endl;
			abort();
		}

		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of given cell (" << cell << ") is too large: " << refinement_level << std::endl;
			abort();
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid refinement level for cell " << cell << ": " << refinement_level << std::endl;
			abort();
		}
		#endif

		if (this->cell_process.count(cell) == 0) {
			return return_neighbors;
		}

		if (!has_children && cell != this->get_child(cell)) {
			return return_neighbors;
		}

		const uint64_t cell_size = this->get_cell_size_in_indices(cell);

		const std::vector<Types<3>::indices_t> indices_of = this->indices_from_neighborhood(
			this->get_indices(cell),
			cell_size,
			neighborhood
		);

		BOOST_FOREACH(const Types<3>::indices_t& index_of, indices_of) {

			if (index_of[0] == error_index) {
				return_neighbors.push_back(0);
				continue;
			}

			const uint64_t neighbor = this->get_existing_cell(
				index_of,
				(refinement_level < max_diff)
					? 0 : refinement_level - max_diff,
				(refinement_level <= this->max_refinement_level - max_diff)
					? refinement_level + max_diff : this->max_refinement_level
			);

			#ifdef DEBUG
			if (neighbor == 0) {
				const Types<3>::indices_t indices = this->get_indices(cell);
				const uint64_t smallest = this->get_existing_cell(index_of, 0, this->max_refinement_level);
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor not found for cell " << cell
					<< " (at indices " << indices[0]
					<< "," << indices[1]
					<< "," << indices[2]
					<< "; ref. lvl. " << refinement_level
					<< ", child of " << this->get_parent(cell)
					<< ") within refinement levels [" << refinement_level - max_diff
					<< ", " << refinement_level + max_diff
					<< "], smallest cell found at indices " << index_of[0]
					<< "," << index_of[1]
					<< "," << index_of[2]
					<< " was " << smallest
					<< " with refinement level " << this->get_refinement_level(smallest)
					<< std::endl;
				abort();
			}

			if (this->cell_process.count(neighbor) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor " << neighbor
					<< " doesn't exist"
					<< std::endl;
				abort();
			}
			#endif

			const int neighbor_ref_lvl = this->get_refinement_level(neighbor);

			#ifdef DEBUG
			if (neighbor_ref_lvl < 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid refinement level for neighbor " << neighbor
					<< " of cell " << cell
					<< std::endl;
				abort();
			}
			#endif

			if (neighbor_ref_lvl <= refinement_level) {
				return_neighbors.push_back(neighbor);
			// add all cells at current search indices within size of given cell
			} else {

				const Types<3>::indices_t index_max = {{
					index_of[0] + cell_size - 1,
					index_of[1] + cell_size - 1,
					index_of[2] + cell_size - 1
				}};

				const std::vector<uint64_t> current_neighbors = this->find_cells(
					index_of,
					index_max,
					std::max(0, refinement_level - max_diff),
					std::min(this->max_refinement_level, refinement_level + max_diff)
				);

				#ifdef DEBUG
				if (current_neighbors.size() == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " No neighbors for cell " << cell
						<< " starting at indices " << index_of[0]
						<< ", " << index_of[1]
						<< ", " << index_of[2]
						<< " between refinement levels " << refinement_level - max_diff
						<< ", " << refinement_level + max_diff
						<< std::endl;
					abort();
				}

				if (current_neighbors.size() < 8) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Too few neighbors for cell " << cell
						<< " of size " << cell_size
						<< " with max_diff " << max_diff
						<< std::endl;

					std::cerr << "Found: ";
					BOOST_FOREACH(const uint64_t& found, current_neighbors) {
						std::cerr << found << " ";
					}

					std::cerr << "\nShould be: ";
					const std::vector<uint64_t> real_neighbors = this->find_cells(
						index_of,
						index_max,
						0,
						this->max_refinement_level
					);
					BOOST_FOREACH(const uint64_t& real, real_neighbors) {
						std::cerr << real << " ";
					}
					std::cerr << std::endl;

					abort();
				}

				BOOST_FOREACH(const uint64_t& current_neighbor, current_neighbors) {
					if (this->cell_process.count(current_neighbor) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Neighbor " << current_neighbor
							<< " doesn't exist between refinement levels " << refinement_level - max_diff
							<< ", " << refinement_level + max_diff
							<< std::endl;
						abort();
					}
				}
				#endif

				return_neighbors.insert(
					return_neighbors.end(),
					current_neighbors.begin(),
					current_neighbors.end()
				);
			}
		}

		return return_neighbors;
	}



	/*!
	Returns cells (which don't have children) that consider given cell as a neighbor.

	Returns nothing if the following cases:
		- given cell has children
		- given cell doesn't exist on any process
	Returned cells are not in any particular order.
	Assumes a maximum refinement level difference of one between neighbors
	(both cases: neighbors_of, neighbors_to).
	*/
	std::vector<uint64_t> find_neighbors_to(
		const uint64_t cell,
		const std::vector<Types<3>::neighborhood_item_t>& neighborhood
	) const {
		std::vector<uint64_t> return_neighbors;

		if (cell == 0
		|| cell > this->last_cell
		|| this->cell_process.count(cell) == 0
		|| cell != this->get_child(cell)) {
			return return_neighbors;
		}

		const int refinement_level = this->get_refinement_level(cell);

		#ifdef DEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << cell << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			abort();
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << cell << " is less than 0: " << refinement_level << std::endl;
			abort();
		}
		#endif

		/*
		FIXME: neighbors should be unique only in within
		each offset in user neighborhood?
		Would allow users to have the same offset
		multiple times in their neighborhoods.
		*/
		boost::unordered_set<uint64_t> unique_neighbors;

		// neighbors_to larger than given cell
		if (refinement_level > 0) {
			const uint64_t parent = this->get_parent(cell);
			const Types<3>::indices_t indices = this->get_indices(parent);
			const uint64_t size_in_indices = this->get_cell_size_in_indices(parent);

			std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
				indices,
				size_in_indices,
				neighborhood
			);

			BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

				if (search_index[0] == error_index) {
					continue;
				}

				const uint64_t found = this->get_cell_from_indices(search_index, refinement_level - 1);
				// only add if found cell doesn't have children
				if (found == this->get_child(found)) {
					unique_neighbors.insert(found);
				}
			}
		}

		// neighbors_to smaller than given cell
		if (refinement_level < this->max_refinement_level) {

			const std::vector<uint64_t> children = this->get_all_children(cell);
			#ifdef DEBUG
			if (children.size() == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Got no children for cell " << cell << std::endl;
				abort();
			}
			#endif

			const uint64_t size_in_indices = this->get_cell_size_in_indices(children[0]);

			BOOST_FOREACH(const uint64_t& child, children) {

				const Types<3>::indices_t indices = this->get_indices(child);

				std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
					indices,
					size_in_indices,
					neighborhood
				);

				BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

					if (search_index[0] == error_index) {
						continue;
					}

					const uint64_t found = this->get_cell_from_indices(search_index, refinement_level + 1);

					if (found == this->get_child(found)) {
						unique_neighbors.insert(found);
					}
				}
			}
		}

		// neighbors_to of the same size as given cell
		const Types<3>::indices_t indices = this->get_indices(cell);
		const uint64_t size_in_indices = this->get_cell_size_in_indices(cell);

		std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
			indices,
			size_in_indices,
			neighborhood
		);

		BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

			if (search_index[0] == error_index) {
				continue;
			}

			const uint64_t found = this->get_cell_from_indices(search_index, refinement_level);
			if (found == this->get_child(found)) {
				unique_neighbors.insert(found);
			}
		}

		return_neighbors.reserve(unique_neighbors.size());

		BOOST_FOREACH(const uint64_t& neighbor, unique_neighbors) {
			return_neighbors.push_back(neighbor);
		}

		return return_neighbors;
	}


	/*!
	As find_neighbors_to(cell) but uses the given neighbors_of list.

	Given list is assumed to have all neighbors_to of given cell that are
	as large or smaller than given cell.
	*/
	std::vector<uint64_t> find_neighbors_to(
		const uint64_t cell,
		const std::vector<uint64_t>& found_neighbors_of
	) const {
		std::vector<uint64_t> return_neighbors;

		if (cell == 0
		|| cell > this->last_cell
		|| this->cell_process.count(cell) == 0
		|| cell != this->get_child(cell)) {
			return return_neighbors;
		}

		// get neighbors_to of given cell, first from its neighbors_of
		boost::unordered_set<uint64_t> unique_neighbors_to;

		BOOST_FOREACH(const uint64_t& neighbor_of, found_neighbors_of) {
			// neighbors_to doesn't store cells that would be outside of the grid
			if (neighbor_of == 0) {
				continue;
			}

			if (this->is_neighbor(neighbor_of, cell)) {
				unique_neighbors_to.insert(neighbor_of);
			}
		}

		const int refinement_level = this->get_refinement_level(cell);
		#ifdef DEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Refinement level (" << refinement_level
				<< ") of cell " << cell
				<< " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")"
				<< std::endl;
			abort();
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Refinement level of cell " << cell
				<< " is less than 0: " << refinement_level
				<< std::endl;
			abort();
		}
		#endif

		// find cells larger than given cell for neighbors_to list
		if (refinement_level > 0) {

			const uint64_t parent = this->get_parent(cell);
			#ifdef DEBUG
			if (parent == cell) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent for cell " << cell << std::endl;
				abort();
			}
			#endif

			const Types<3>::indices_t indices = this->get_indices(parent);
			const uint64_t size_in_indices = this->get_cell_size_in_indices(parent);
			const std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
				indices,
				size_in_indices,
				this->neighborhood_to
			);

			BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

				if (search_index[0] == error_index) {
					continue;
				}

				const uint64_t found = this->get_cell_from_indices(search_index, refinement_level - 1);
				// only add if found cell doesn't have children
				if (found == this->get_child(found)) {
					unique_neighbors_to.insert(found);
				}
			}
		}

		return_neighbors.reserve(unique_neighbors_to.size());
		return_neighbors.insert(
			return_neighbors.begin(),
			unique_neighbors_to.begin(),
			unique_neighbors_to.end()
		);

		return return_neighbors;
	}


	/*!
	Returns unique cells within given rectangular box and refinement levels (both inclusive).

	Cells within given volume are always returned in the following order:
	Starting from the corner closest to the starting corner of the grid cells are returned
	first in the positive x direction then y direction and finally z direction.
	*/
	// TODO: make private?
	std::vector<uint64_t> find_cells
	(
		const Types<3>::indices_t indices_min,
		const Types<3>::indices_t indices_max,
		const int minimum_refinement_level,
		const int maximum_refinement_level
	) const {
		// size of cells in indices of given maximum_refinement_level
		const uint64_t index_increase = uint64_t(1) << (this->max_refinement_level - maximum_refinement_level);

		#ifdef DEBUG
		if (minimum_refinement_level > maximum_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid refinement levels given" << std::endl;
			abort();
		}

		if (maximum_refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid maximum refinement level given" << std::endl;
			abort();
		}

		// check that outer shell makes sense
		if (indices_min[0] > indices_max[0]) {
			std::cerr << __FILE__ << ":" << __LINE__ << " minimum x index > maximum x index" << std::endl;
			abort();
		}

		if (indices_min[1] > indices_max[1]) {
			std::cerr << __FILE__ << ":" << __LINE__ << " minimum y index > maximum y index" << std::endl;
			abort();
		}

		if (indices_min[2] > indices_max[2]) {
			std::cerr << __FILE__ << ":" << __LINE__ << " minimum z index > maximum z index" << std::endl;
			abort();
		}
		#endif

		std::vector<uint64_t> result;
		boost::unordered_set<uint64_t> uniques;

		Types<3>::indices_t indices = {{0, 0, 0}};
		for (indices[2] = indices_min[2]; indices[2] <= indices_max[2]; indices[2] += index_increase)
		for (indices[1] = indices_min[1]; indices[1] <= indices_max[1]; indices[1] += index_increase)
		for (indices[0] = indices_min[0]; indices[0] <= indices_max[0]; indices[0] += index_increase) {

			const uint64_t cell = this->get_existing_cell(indices, minimum_refinement_level, maximum_refinement_level);

			#ifdef DEBUG
			if (cell == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No cell found between refinement levels [" << minimum_refinement_level
					<< ", " << maximum_refinement_level
					<< "] at indices " << indices[0]
					<< " " << indices[1]
					<< " " << indices[2]
					<< std::endl;

				const uint64_t smallest = this->get_existing_cell(indices, 0, this->max_refinement_level);
				std::cerr << __FILE__ << ":" << __LINE__
					<< " smallest cell there is " << smallest
					<< " with refinement level " << this->get_refinement_level(smallest)
					<< std::endl;
				abort();
			}

			if (this->cell_process.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell
					<< " doesn't exist"
					<< std::endl;
				abort();
			}

			if (cell > this->last_cell) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell can't exist" << std::endl;
				abort();
			}
			#endif

			/*
			When searching for neighbors_to, cells may exist with larger refinement
			level than given in find_neighbors_to and shouldn't be considered.
			*/
			if (cell != this->get_child(cell)) {
				continue;
			}

			// don't add the same cell twice
			if (uniques.count(cell) == 0) {
				uniques.insert(cell);
				result.push_back(cell);
			}
		}

		return result;
	}


	/*!
	Removes user data of refined and unrefined cells from this process.
	*/
	void clear_refined_unrefined_data()
	{
		this->refined_cell_data.clear();
		this->unrefined_cell_data.clear();
	}


	/*!
	Sets the given option for non-hierarchial partitioning.

	Does nothing if option name is one of:
		- RETURN_LISTS
		- EDGE_WEIGHT_DIM
		- NUM_GID_ENTRIES
		- OBJ_WEIGHT_DIM
	Call this with name = LB_METHOD and value = HIER to use hierarchial
	partitioning and set those options using the other function with this name.
	*/
	void set_partitioning_option(const std::string name, const std::string value)
	{
		if (this->reserved_options.count(name) > 0) {
			#ifdef DEBUG
			std::cerr << "User tried to set an option reserved for dccrg (" << name << ": " << value << ")" << std::endl;
			#endif
			return;
		}

		Zoltan_Set_Param(this->zoltan, name.c_str(), value.c_str());
	}


	/*!
	Adds a new level for hierarchial partitioning with each part having given number of processes.

	Assigns default partitioning options for the added level.
	Does nothing if processes_per_part < 1.
	*/
	void add_partitioning_level(const int processes)
	{
		if (processes < 1) {
			#ifdef DEBUG
			std::cerr << "User tried to assign " << processes << " processes per part for a new hierarchial partitioning level" << std::endl;
			#endif
			return;
		}

		this->processes_per_part.push_back(processes);

		// create default partitioning options for the level
		boost::unordered_map<std::string, std::string> default_load_balance_options;
		default_load_balance_options["LB_METHOD"] = "HYPERGRAPH";
		default_load_balance_options["PHG_CUT_OBJECTIVE"] = "CONNECTIVITY";
		this->partitioning_options.push_back(default_load_balance_options);
	}


	/*!
	Rremoves the given hierarhchial partitioning level.

	Level numbering starts from 0.
	Does nothing if given level doesn't exist.
	*/
	void remove_partitioning_level(const int hierarchial_partitioning_level)
	{
		if (hierarchial_partitioning_level < 0
		|| hierarchial_partitioning_level >= int(this->processes_per_part.size())) {
			return;
		}

		this->processes_per_part.erase(this->processes_per_part.begin() + hierarchial_partitioning_level);
		this->partitioning_options.erase(this->partitioning_options.begin() + hierarchial_partitioning_level);
	}


	/*!
	Adds (or overwrites) the given option and its value for hierarchial partitioning of given level.

	Level numbering starts from 0.
	Does nothing in the following cases:
		- option name is one of: RETURN_LISTS, ...
		- given level doesn't exist
	*/
	void add_partitioning_option(const int hierarchial_partitioning_level, const std::string name, const std::string value)
	{
		if (hierarchial_partitioning_level < 0
		|| hierarchial_partitioning_level >= int(this->processes_per_part.size())) {
			return;
		}

		if (this->reserved_options.count(name) > 0) {
			#ifdef DEBUG
			std::cerr << "User tried to set an option reserved for dccrg (" << name << ": " << value << ") for level " << hierarchial_partitioning_level << " of hierarchial partitioning" << std::endl;
			#endif
			return;
		}

		this->partitioning_options[hierarchial_partitioning_level][name] = value;
	}


	/*!
	Removes the given option from the given level of hierarchial partitioning.

	Level numbering starts from 0.
	Does nothing if given level doesn't exist.
	*/
	void remove_partitioning_option(const int hierarchial_partitioning_level, const std::string name)
	{
		if (hierarchial_partitioning_level < 0
		|| hierarchial_partitioning_level >= int(this->processes_per_part.size())) {
			return;
		}

		this->partitioning_options[hierarchial_partitioning_level].erase(name);
	}


	/*!
	Returns the names of partitioning options for hierarchial partitioning at given level.

	Returns nothing if given level doesn't exist.
	*/
	std::vector<std::string> get_partitioning_options(const int hierarchial_partitioning_level) const
	{
		std::vector<std::string> partitioning_options;

		if (hierarchial_partitioning_level < 0
		|| hierarchial_partitioning_level >= int(this->processes_per_part.size())) {
			return partitioning_options;
		}

		for (boost::unordered_map<std::string, std::string>::const_iterator
			option = this->partitioning_options.at(hierarchial_partitioning_level).begin();
			option != this->partitioning_options.at(hierarchial_partitioning_level).end();
			option++
		) {
			partitioning_options.push_back(option->first);
		}

		return partitioning_options;
	}


	/*!
	Returns the value of given non-hierarchial partitioning option.

	Returns an empty string if given option or given level doesn't exist.
	*/
	std::string get_partitioning_option_value(const int hierarchial_partitioning_level, const std::string name) const
	{
		std::string value;

		if (hierarchial_partitioning_level < 0
		|| hierarchial_partitioning_level >= int(this->processes_per_part.size())) {
			return value;
		}

		if (this->partitioning_options.count(name) > 0) {
			value = this->partitioning_options.at(name);
		}

		return value;
	}


	/*!
	Returns the process which has the given cell or -1 if the cell doesn't exist.
	*/
	int get_process(const uint64_t cell) const
	{
		if (this->cell_process.count(cell) == 0) {
			return -1;
		} else {
			return this->cell_process.at(cell);
		}
	}


	/*!
	Given cell is kept on this process during subsequent load balancing.

	Does nothing in the same cases as pin(cell, process).
	*/
	void pin(const uint64_t cell)
	{
		this->pin(cell, this->rank);
	}

	/*!
	Given cell is sent to the given process and kept there during subsequent load balancing.

	Does nothing in the following cases:
		- given cell doesn't exist
		- given cell exists on another process
		- given cell has children
		- given process doesn't exist
	*/
	void pin(const uint64_t cell, const int process)
	{
		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		if (process < 0 || process >= (int) this->comm_size) {
			return;
		}

		// do nothing if the request already exists
		if (this->pin_requests.count(cell) > 0
		&& (int) this->pin_requests.at(cell) == process) {
			return;
		}

		this->new_pin_requests[cell] = (uint64_t) process;
	}

	/*!
	Allows the given cell to be moved to another process during subsequent load balancing.

	Does nothing in the following cases:
		- given cell has children
		- given cell doesn't exist
		- given cell exists on another process
	*/
	void unpin(const uint64_t cell)
	{
		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		if (this->pin_requests.count(cell) > 0) {
			// non-existing process means unpin
			this->new_pin_requests[cell] = this->comm_size;
		} else {
			this->new_pin_requests.erase(cell);
		}
	}

	/*!
	Executes unpin(cell) for all cells on this process.
	*/
	void unpin_local_cells()
	{
		#ifdef DEBUG
		// check that all child cells on this process are also in this->cells.
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			i = this->cell_process.begin();
			i != this->cell_process.end();
			i++
		) {
			const uint64_t cell = i->first;

			if (this->cell_process.at(cell) != this->rank) {
				return;
			}

			if (cell == this->get_child(cell)) {
				if (this->cells.count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " should be in this->cells of process " << this->rank << std::endl;
					abort();
				}
			} else {
				if (this->cells.count(cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " shouldn't be in this->cells of process " << this->rank << std::endl;
					abort();
				}
			}
		}
		#endif

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
			this->unpin(item.first);
		}
	}

	/*!
	All cells in the grid are free to move between processes.

	Must be called simultaneously on all processes.
	*/
	void unpin_all_cells()
	{
		this->new_pin_requests.clear();
		this->pin_requests.clear();
	}


	/*!
	Returns a pointer to the send lists of this process.

	These lists record which cells' user data this process will send during neighbor data updates.
	The key is the target process.
	*/
	#ifdef DCCRG_SEND_SINGLE_CELLS
	const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >*
	#else
	const boost::unordered_map<int, std::vector<uint64_t> >*
	#endif
	get_send_lists()
	{
		return &(this->cells_to_send);
	}

	/*!
	Returns a pointer to the receive lists of this process.

	These lists record which cells' user data this process will receive during neighbor data updates.
	The key is the source process.
	*/
	#ifdef DCCRG_SEND_SINGLE_CELLS
	const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >*
	#else
	const boost::unordered_map<int, std::vector<uint64_t> >*
	#endif
	get_receive_lists()
	{
		return &(this->cells_to_receive);
	}


	/*!
	Returns a pointer to the set of local cells which have at least one neighbor
	on another process.
	*/
	const boost::unordered_set<uint64_t>* get_cells_with_remote_neighbors() const
	{
		return &(this->cells_with_remote_neighbors);
	}

	/*!
	Returns a vector of local cells which have at least one neighbor on another process.
	*/
	std::vector<uint64_t> get_list_of_cells_with_remote_neighbors() const
	{
		std::vector<uint64_t> result(
			this->cells_with_remote_neighbors.begin(),
			this->cells_with_remote_neighbors.end()
		);
		return result;
	}

	/*!
	Returns a pointer to the set of remote cells which have at least one local neighbor.
	*/
	const boost::unordered_set<uint64_t>* get_remote_cells_with_local_neighbors() const
	{
		return &(this->remote_cells_with_local_neighbors);
	}

	/*!
	See get_remote_neighbors().
	*/
	std::vector<uint64_t> get_list_of_remote_cells_with_local_neighbors() const
	{
		std::vector<uint64_t> result(
			this->remote_cells_with_local_neighbors.begin(),
			this->remote_cells_with_local_neighbors.end()
		);
		return result;
	}


	/*!
	Sets the weight of given local existing cell without children.

	Does nothing if above conditions are not met.
	Cell weights are given to Zoltan when balancing the load.
	Unset cell weights are assumed to be 1.

	User set cell weights are removed when balance_load is called.
	Children of refined cells inherit their parent's weight.
	Parents of unrefined cells do not inherit the moved cells' weights.
	*/
	void set_cell_weight(const uint64_t cell, const double weight)
	{
		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		this->cell_weights[cell] = weight;
	}

	/*!
	Returns the weight of given local existing cell without children.

	Returns a quiet nan if above conditions are not met.
	Unset cell weights are assumed to be 1.
	*/
	double get_cell_weight(const uint64_t cell) const
	{
		if (this->cell_process.count(cell) == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (this->cell_process.at(cell) != this->rank) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (cell != this->get_child(cell)) {
			return std::numeric_limits<double>::quiet_NaN();
		}

		if (this->cell_weights.count(cell) == 0) {
			return 1;
		} else {
			return this->cell_weights.at(cell);
		}
	}


	/*!
	Returns a pointer to the list of cells that will be added
	to this process after preparing to migrate or load balance cells.
	*/
	const boost::unordered_set<uint64_t>* get_balance_added_cells() const
	{
		return &(this->added_cells);
	}

	/*!
	Returns a pointer to the list of cells that will be removed
	from this process after preparing to migrate or load balance cells.
	*/
	const boost::unordered_set<uint64_t>* get_balance_removed_cells() const
	{
		return &(this->removed_cells);
	}


	/*!
	Returns the smallest existing cell at the given coordinate.

	Returns error_cell if the coordinate is outside of the grid or the cell is on another process.
	*/
	uint64_t get_existing_cell(const double x, const double y, const double z) const
	{
		const Types<3>::indices_t indices = {{
			this->get_x_index_of_coord(x),
			this->get_y_index_of_coord(y),
			this->get_z_index_of_coord(z)
		}};

		if (indices[0] == error_index
		|| indices[1] == error_index
		|| indices[2] == error_index) {
			return error_cell;
		}

		return this->get_existing_cell(indices, 0, this->max_refinement_level);
	}


 	/*!
	Returns the siblings of given cell regardless of whether they exist.

	If given a cell of refinement level 0 returns the given cell.
	Returns nothing if given cell's refinement level exceeds the maximum of this grid.
	*/
	std::vector<uint64_t> get_siblings(const uint64_t cell) const
	{
		std::vector<uint64_t> siblings;

		const int refinement_level = this->get_refinement_level(cell);
		if (refinement_level < 0
		|| refinement_level > this->max_refinement_level) {
			return siblings;
		}

		if (refinement_level == 0) {
			siblings.push_back(cell);
			return siblings;
		}

		return this->get_all_children(this->get_parent(cell));
	}


	/*!
	Returns all children of given cell regardless of whether they exist.

	Returns nothing if childrens' refinement level would exceed max_refinement_level or
	given cell doesn't exist.
	 */
	std::vector<uint64_t> get_all_children(const uint64_t cell) const
	{
		std::vector<uint64_t> children;

		if (cell == error_cell) {
			return children;
		}

		if (this->cell_process.count(cell) == 0) {
			return children;
		}

		// given cell cannot have children
		int refinement_level = this->get_refinement_level(cell);
		if (refinement_level >= this->max_refinement_level) {
			return children;
		}

		children.reserve(8);

		Types<3>::indices_t indices = this->get_indices(cell);

		// get indices of next refinement level within this cell
		refinement_level++;
		const uint64_t index_offset = (uint64_t(1) << (this->max_refinement_level - refinement_level));
		for (uint64_t z_index_offset = 0; z_index_offset < 2 * index_offset; z_index_offset += index_offset)
		for (uint64_t y_index_offset = 0; y_index_offset < 2 * index_offset; y_index_offset += index_offset)
		for (uint64_t x_index_offset = 0; x_index_offset < 2 * index_offset; x_index_offset += index_offset) {
			children.push_back(
				this->get_cell_from_indices(
					indices[0] + x_index_offset,
					indices[1] + y_index_offset,
					indices[2] + z_index_offset,
					refinement_level
				)
			);
		}

		return children;
	}


	boost::unordered_map<uint64_t, double>::size_type get_number_of_cell_weights() const
	{
		return this->cell_weights.size();
	}


	/*!
	Adds a new neighborhood for remote neighbor updates.

	Must be called with identical parameters on all processes.
	No neighborhood_item_t should have all offsets equal to 0.
	Returns true on success and false in any of the following cases:
		- given id already exists, use remove_remote_... before calling this
		- (part of) given neighborhood is outside of initial neighborhood size

	Use this to reduce data transfers between processes when the full
	neighborhood doesn't have to be used.
	*/
	bool add_remote_update_neighborhood(
		const int id,
		const std::vector<Types<3>::neighborhood_item_t>& given_neigh
	) {
		if (this->user_hood_of.count(id) > 0) {

			#ifdef DEBUG
			if (this->user_hood_to.count(id) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should have id " << id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_of.count(id) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should have id " << id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_to.count(id) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should have id " << id
					<< std::endl;
				abort();
			}
			#endif

			return false;
		}

		#ifdef DEBUG
		if (this->user_hood_to.count(id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_of.count(id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_to.count(id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_send.count(id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_receive.count(id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << id
				<< std::endl;
			abort();
		}
		#endif

		if (this->neighborhood_size > 0) {

			BOOST_FOREACH(const Types<3>::neighborhood_item_t& neigh_item, given_neigh) {
				for (size_t i = 0; i < 3; i++) {
					if ((unsigned int) abs(neigh_item[i]) > this->neighborhood_size) {
						return false;
					}
				}

				if (neigh_item[0] == 0 && neigh_item[1] == 0 && neigh_item[2] == 0) {
					return false;
				}
			}

		} else {

			BOOST_FOREACH(const Types<3>::neighborhood_item_t& neigh_item, given_neigh) {
				int zero_offsets = 0;

				for (size_t i = 0; i < 3; i++) {
					if (neigh_item[i] == 0) {
						zero_offsets++;
					}
					if (abs(neigh_item[i]) > 1) {
						return false;
					}
				}

				if (zero_offsets != 2) {
					return false;
				}
			}

		}

		// set user_hood_of and _to
		this->user_hood_of[id] = given_neigh;
		this->user_hood_to[id];
		BOOST_FOREACH(const Types<3>::neighborhood_item_t& neigh_item, given_neigh) {
			const Types<3>::neighborhood_item_t neigh_item_to = {{
				-neigh_item[0],
				-neigh_item[1],
				-neigh_item[2]
			}};
			this->user_hood_to[id].push_back(neigh_item_to);
		}

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
			this->update_user_neighbors(item.first, id);
		}

		this->recalculate_neighbor_update_send_receive_lists(id);

		return true;
	}


	/*!
	Removes the given neighborhood from remote neighbor updates.

	Must be called with identical id on all processes.
	Frees local neighbor lists and other resources associated
	with given neighborhood.
	*/
	void remove_remote_update_neighborhood(const int id)
	{
		this->user_hood_of.erase(id);
		this->user_hood_to.erase(id);
		this->user_neigh_of.erase(id);
		this->user_neigh_to.erase(id);
		this->user_neigh_cells_to_send.erase(id);
		this->user_neigh_cells_to_receive.erase(id);
	}


	/*!
	Returns dccrg's communicator.

	Returns a duplicate of the MPI communicator of dccrg
	which itself is a duplicate of the communicator that
	was given to dccrg's initialize function.
	*/
	MPI_Comm get_communicator() const
	{
		MPI_Comm result;
		if (MPI_Comm_dup(this->comm, &result) != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__ << " MPI_Comm_dup failed." << std::endl;
			abort();
		}
		return result;
	}



private:

	bool initialized;

	// size of the neighbor stencil of a cells in cells (of the same size as the cell itself)
	unsigned int neighborhood_size;

	// maximum difference in refinement level between neighbors
	int max_ref_lvl_diff;

	// the grid is distributed between these processes
	MPI_Comm comm;
	uint64_t rank, comm_size;
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::mpi::communicator boost_comm;
	#endif

	// cells and their data on this process
	boost::unordered_map<uint64_t, UserData> cells;

	// cell on this process and its neighbors
	boost::unordered_map<uint64_t, std::vector<uint64_t> > neighbors;

	/*
	Offsets of cells that are considered as neighbors of a cell and
	offsets of cells that consider a cell as a neighbor
	*/
	std::vector<Types<3>::neighborhood_item_t> neighborhood_of, neighborhood_to;

	/*
	User defined versions of neighborhood_of and _to.
	*/
	boost::unordered_map<
		int, // user defined id of neighborhood
		std::vector<Types<3>::neighborhood_item_t>
	> user_hood_of, user_hood_to;

	/*!
	Cell on this process and those cells that aren't neighbors of this cell but whose neighbor this cell is.
	For example with a stencil size of 1 in the following grid:
\verbatim
|-----------|
|     |5 |6 |
|  1  |--|--|
|     |9 |10|
|-----------|
\endverbatim
	neighbors_to[6] = 1 because neighbors[6] = 5, 9, 10 while
	neighbors_to[5] is empty because neighbors[5] = 1, 6, 9, 10
	*/
	boost::unordered_map<uint64_t, std::vector<uint64_t> > neighbors_to;

	/*
	User defined versions of neighbors_of and _to
	*/
	boost::unordered_map<
		int, // user defined id of neighbor lists
		boost::unordered_map<uint64_t, std::vector<uint64_t> >
	> user_neigh_of, user_neigh_to;

	// on which process every cell in the grid is
	boost::unordered_map<uint64_t, uint64_t> cell_process;

	// cells on this process that have a neighbor on another process or are considered as a neighbor of a cell on another process
	boost::unordered_set<uint64_t> cells_with_remote_neighbors;

	// cells on other processes that have a neighbor on this process or are considered as a neighbor of a cell on this process
	boost::unordered_set<uint64_t> remote_cells_with_local_neighbors;

	// remote neighbors and their data, of cells on this process
	boost::unordered_map<uint64_t, UserData> remote_neighbors;

	#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::unordered_map<int, std::vector<MPI_Request> > send_requests, receive_requests;
	#else
	boost::unordered_map<int, std::vector<boost::mpi::request> > send_requests, receive_requests;
	#endif

	// cells whose data has to be received / sent by this process from the process as the key
	#ifdef DCCRG_SEND_SINGLE_CELLS
	// store cell, tag pairs so users can also send the data themselves easily
	boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > > cells_to_send, cells_to_receive;
	#else
	boost::unordered_map<int, std::vector<uint64_t> > cells_to_send, cells_to_receive;
	#endif

	/*
	User defined neighborhood versions of cells_to_send and _receive.
	*/
	boost::unordered_map<
		int, // user defined id of neighbor lists
		boost::unordered_map<
			int, // process to send to / receive from
			#ifdef DCCRG_SEND_SINGLE_CELLS
			std::vector<std::pair<uint64_t, int> >
			#else
			std::vector<uint64_t>
			#endif
		>
	> user_neigh_cells_to_send, user_neigh_cells_to_receive;

	// cells added to / removed from this process by load balancing
	boost::unordered_set<uint64_t> added_cells, removed_cells;

	#ifndef DCCRG_SEND_SINGLE_CELLS
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	// storage for cells' user data that awaits transfer to or from this process
	boost::unordered_map<int, std::vector<UserData> > incoming_data, outgoing_data;
	#endif
	#endif

	// cells to be refined / unrefined after a call to stop_refining()
	boost::unordered_set<uint64_t> cells_to_refine, cells_to_unrefine;

	// cells whose siblings shouldn't be unrefined
	boost::unordered_set<uint64_t> cells_not_to_unrefine;

	// stores user data of cells whose children were created while refining
	boost::unordered_map<uint64_t, UserData> refined_cell_data;
	// stores user data of cells that were removed while unrefining
	boost::unordered_map<uint64_t, UserData> unrefined_cell_data;

	// cell that should be kept on a particular process
	boost::unordered_map<uint64_t, uint64_t> pin_requests;
	// pin requests given since that last time load was balanced
	boost::unordered_map<uint64_t, uint64_t> new_pin_requests;

	// variables for load balancing using Zoltan
	Zoltan_Struct* zoltan;
	// number of processes per part in a hierarchy level (numbering starts from 0)
	std::vector<unsigned int> processes_per_part;
	// options for each level of hierarchial load balancing (numbering start from 0)
	std::vector<boost::unordered_map<std::string, std::string> > partitioning_options;
	// record whether Zoltan_LB_Partition is expected to fail (when the user selects NONE as the load balancing algorithm)
	bool no_load_balancing;
	// reserved options that the user cannot change
	boost::unordered_set<std::string> reserved_options;

	// optional user-given weights of cells on this process
	boost::unordered_map<uint64_t, double> cell_weights;

	// processes which have cells close enough from cells of this process
	boost::unordered_set<int> neighbor_processes;


	/*!
	Moves cells between processes due to load balancing or user request.

	Recalculates neighbor lists, etc.
	Must be called simultaneously on all processes.
	Clears user-given weights of all cells.
	*/
	void move_cells() {
		// TODO: get rid of added_cells and removed_cells and use cells_to_send and receive instead?

		this->cell_weights.clear();
		this->cells_with_remote_neighbors.clear();
		this->remote_cells_with_local_neighbors.clear();
		this->remote_neighbors.clear();
		this->cells_to_refine.clear();
		this->refined_cell_data.clear();
		this->cells_to_unrefine.clear();
		this->unrefined_cell_data.clear();
		this->cells_not_to_unrefine.clear();

		/*
		Calculate where cells have migrated to update internal data structures
		Any cell can end up on any process and any neighbor of any cell can end up on yet another process
		*/

		// removed cells on all processes
		std::vector<uint64_t> temp_removed_cells(
			this->removed_cells.begin(),
			this->removed_cells.end()
		);
		std::sort(temp_removed_cells.begin(), temp_removed_cells.end());

		std::vector<std::vector<uint64_t> > all_removed_cells;
		All_Gather()(temp_removed_cells, all_removed_cells, this->comm);

		// created cells on all processes
		std::vector<uint64_t> temp_added_cells(
			this->added_cells.begin(),
			this->added_cells.end()
		);
		std::sort(temp_added_cells.begin(), temp_added_cells.end());

		std::vector<std::vector<uint64_t> > all_added_cells;
		All_Gather()(temp_added_cells, all_added_cells, this->comm);

		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->cells, this->cells_to_receive, this->cells_to_send
		#elif defined (DCCRG_CELL_DATA_SIZE_FROM_USER)
		this->cells, this->cells_to_receive, this->cells_to_send
		#else
		this->cells_to_receive, this->cells_to_send
		#endif
		);

		#ifdef DEBUG
		// check that there are no duplicate adds / removes
		boost::unordered_set<uint64_t> all_adds, all_removes;

		BOOST_FOREACH(const std::vector<uint64_t>& item, all_removed_cells) {
			BOOST_FOREACH(const uint64_t& removed_cell, item) {

				if (all_removes.count(removed_cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << removed_cell
						<< " was already removed"
						<< std::endl;
					abort();
				}
				all_removes.insert(removed_cell);
			}
		}

		BOOST_FOREACH(const std::vector<uint64_t>& item, all_added_cells) {
			BOOST_FOREACH(const uint64_t& added_cell, item) {

				if (all_adds.count(added_cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << added_cell
						<< " was already removed"
						<< std::endl;
					abort();
				}
				all_adds.insert(added_cell);
			}
		}

		// check that cells were removed by their process
		for (uint64_t cell_remover = 0; cell_remover < all_removed_cells.size(); cell_remover++) {

			BOOST_FOREACH(const uint64_t& removed_cell, all_removed_cells.at(cell_remover)) {

				if (this->cell_process.at(removed_cell) != cell_remover) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << removed_cell
						<< " doesn't belong to process " << cell_remover
						<< std::endl;
					abort();
				}
			}
		}
		#endif

		// update cell to process mappings
		for (uint64_t cell_creator = 0; cell_creator < all_added_cells.size(); cell_creator++) {

			BOOST_FOREACH(const uint64_t& created_cell, all_added_cells.at(cell_creator)) {
				this->cell_process.at(created_cell) = cell_creator;
			}
		}

		#ifdef DEBUG
		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Grid is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->pin_requests_succeeded()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Pin requests didn't succeed" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// create neighbor lists for cells without children that came to this process
		BOOST_FOREACH(const uint64_t& added_cell, this->added_cells) {

			if (added_cell != this->get_child(added_cell)) {
				continue;
			}

			// TODO: use this->update_neighbors(added_cell)
			this->neighbors[added_cell]
				= this->find_neighbors_of(added_cell, this->neighborhood_of, this->max_ref_lvl_diff);
			this->neighbors_to[added_cell]
				= this->find_neighbors_to(added_cell, this->neighborhood_to);

			// also update user neighbor lists
			for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
				item = this->user_hood_of.begin();
				item != this->user_hood_of.end();
				item++
			) {
				this->update_user_neighbors(added_cell, item->first);
			}
		}

		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->cells, this->cells_to_receive
		#endif
		#endif
		);
		this->wait_user_data_transfer_sends();
		this->cells_to_send.clear();
		this->cells_to_receive.clear();

		// free user data and neighbor lists of cells removed from this process
		BOOST_FOREACH(const uint64_t& removed_cell, this->removed_cells) {
			this->cells.erase(removed_cell);
			this->neighbors.erase(removed_cell);
			this->neighbors_to.erase(removed_cell);

			// also user neighbor lists
			for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
				item = this->user_hood_of.begin();
				item != this->user_hood_of.end();
				item++
			) {
				this->user_neigh_of.at(item->first).erase(removed_cell);
				this->user_neigh_to.at(item->first).erase(removed_cell);
			}
		}

		this->update_remote_neighbor_info();

		this->recalculate_neighbor_update_send_receive_lists();

		#ifdef DEBUG
		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " The grid is inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_neighbors()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbor lists are incorrect" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbor info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " virhe" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}

	/*!
	Prepares to move cells between processes with move_cells.

	Sends user data of cells between processes once before move_cells so that
	the correct mpi datatype can be constructed when actually moving cells.

	Must be called simultaneously on all processes.
	move_cells must be the next dccrg function to be called after this one.
	*/
	void prepare_to_move_cells() {

		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbor info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " virhe" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->cells_with_remote_neighbors.clear();
		this->remote_cells_with_local_neighbors.clear();
		this->remote_neighbors.clear();
		this->cells_to_refine.clear();
		this->refined_cell_data.clear();
		this->cells_to_unrefine.clear();
		this->unrefined_cell_data.clear();

		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->cells, this->cells_to_receive, this->cells_to_send
		#elif defined (DCCRG_CELL_DATA_SIZE_FROM_USER)
		this->cells, this->cells_to_receive, this->cells_to_send
		#else
		this->cells_to_receive, this->cells_to_send
		#endif
		);

		#ifdef DEBUG
		// check that there are no duplicate adds / removes
		// removed cells on all processes
		std::vector<uint64_t> temp_removed_cells(
			this->removed_cells.begin(),
			this->removed_cells.end()
		);
		std::sort(temp_removed_cells.begin(), temp_removed_cells.end());

		std::vector<std::vector<uint64_t> > all_removed_cells;
		All_Gather()(temp_removed_cells, all_removed_cells, this->comm);

		// created cells on all processes
		std::vector<uint64_t> temp_added_cells(
			this->added_cells.begin(),
			this->added_cells.end()
		);
		std::sort(temp_added_cells.begin(), temp_added_cells.end());

		std::vector<std::vector<uint64_t> > all_added_cells;
		All_Gather()(temp_added_cells, all_added_cells, this->comm);

		boost::unordered_set<uint64_t> all_adds, all_removes;

		BOOST_FOREACH(const std::vector<uint64_t>& item, all_removed_cells) {
			BOOST_FOREACH(const uint64_t& removed_cell, item) {

				if (all_removes.count(removed_cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << removed_cell
						<< " was already removed"
						<< std::endl;
					abort();
				}
				all_removes.insert(removed_cell);
			}
		}

		BOOST_FOREACH(const std::vector<uint64_t>& item, all_added_cells) {
			BOOST_FOREACH(const uint64_t& added_cell, item) {

				if (all_adds.count(added_cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << added_cell
						<< " was already removed"
						<< std::endl;
					abort();
				}
				all_adds.insert(added_cell);
			}
		}

		// check that cells were removed by their process
		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {

			BOOST_FOREACH(const uint64_t& removed_cell, all_removed_cells.at(cell_remover)) {

				if (this->cell_process.at(removed_cell) != cell_remover) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << removed_cell
						<< " doesn't belong to process " << cell_remover
						<< std::endl;
					abort();
				}
			}
		}
		#endif

		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->cells, this->cells_to_receive
		#endif
		#endif
		);
		this->wait_user_data_transfer_sends();
	}


	/*!
	Updates user pin requests globally based on new_pin_requests.

	Must be called simultaneously on all processes.
	*/
	void update_pin_requests()
	{
		std::vector<uint64_t> new_pinned_cells, new_pinned_processes;

		new_pinned_cells.reserve(this->new_pin_requests.size());
		new_pinned_processes.reserve(this->new_pin_requests.size());
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			item = this->new_pin_requests.begin();
			item != this->new_pin_requests.end();
			item++
		) {
			new_pinned_cells.push_back(item->first);
			new_pinned_processes.push_back(item->second);
		}

		std::vector<std::vector<uint64_t> > all_new_pinned_cells, all_new_pinned_processes;
		All_Gather()(new_pinned_cells, all_new_pinned_cells, this->comm);
		All_Gather()(new_pinned_processes, all_new_pinned_processes, this->comm);

		for (uint64_t process = 0; process < all_new_pinned_cells.size(); process++) {
			for (uint64_t i = 0; i < all_new_pinned_cells.at(process).size(); i++) {

				const uint64_t requested_process = all_new_pinned_processes[process][i];

				if (requested_process >= this->comm_size) {
					this->pin_requests.erase(all_new_pinned_cells[process][i]);
				} else {
					this->pin_requests[all_new_pinned_cells[process][i]] = requested_process;
				}

				#ifdef DEBUG
				if (this->cell_process.at(all_new_pinned_cells[process][i]) != process) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Process " << process << " tried pin cell " << all_new_pinned_cells[process][i] << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}
		}

		this->new_pin_requests.clear();
	}


	/*!
	Repartitions cells across processes based on user requests and
	Zoltan if use_zoltan is true.

	Updates send & receive lists.
	*/
	void make_new_partition(const bool use_zoltan)
	{
		this->update_pin_requests();

		int partition_changed, global_id_size, local_id_size, number_to_receive, number_to_send;
		ZOLTAN_ID_PTR global_ids_to_receive, local_ids_to_receive, global_ids_to_send, local_ids_to_send;
		int *sender_processes, *receiver_processes;

		if (use_zoltan && Zoltan_LB_Balance(
			this->zoltan,
			&partition_changed,
			&global_id_size,
			&local_id_size,
			&number_to_receive,
			&global_ids_to_receive,
			&local_ids_to_receive,
			&sender_processes,
			&number_to_send,
			&global_ids_to_send,
			&local_ids_to_send,
			&receiver_processes
			) != ZOLTAN_OK
		) {
			if (!this->no_load_balancing) {
				if (this->rank == 0) {
					std::cerr << "Zoltan_LB_Partition failed" << std::endl;
				}
				Zoltan_Destroy(&this->zoltan);
				// TODO: throw an exception instead
				abort();
			}

			#ifdef DEBUG
			// check that processes have the cells they're supposed to send
			for (int i = 0; i < number_to_send; i++) {
				if (this->cells.count(global_ids_to_send[i]) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cannot send cell " << global_ids_to_send[i] << " to process " << receiver_processes[i] << std::endl;
					abort();
				}
			}

			// check that cells to be received are on the sending process
			for (int i = 0; i < number_to_receive; i++) {
				if (this->cell_process.at(global_ids_to_receive[i]) != (uint64_t)sender_processes[i]) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cannot receive cell " << global_ids_to_receive[i]
						<< " from process " << sender_processes[i]
						<< std::endl;
					abort();
				}
			}
			#endif
		}

		this->cells_to_receive.clear();
		this->cells_to_send.clear();

		/*
		Processes and the cells for which data has to be received by this process
		*/

		// migration from user
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			pin_request = this->pin_requests.begin();
			pin_request != this->pin_requests.end();
			pin_request++
		) {
			const uint64_t current_process_of_cell = this->cell_process.at(pin_request->first);

			if (pin_request->second == this->rank
			&& current_process_of_cell != this->rank) {
				this->cells_to_receive[current_process_of_cell].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(pin_request->first, -1)
					#else
					pin_request->first
					#endif
				);
				this->added_cells.insert(pin_request->first);
			}
		}

		// migration from Zoltan
		if (use_zoltan) {
			for (int i = 0; i < number_to_receive; i++) {

				// don't send / receive from self
				if ((uint64_t)sender_processes[i] == this->rank) {
					continue;
				}

				// skip user-migrated cells
				if (this->pin_requests.count(global_ids_to_receive[i]) > 0) {
					continue;
				}

				this->cells_to_receive[sender_processes[i]].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(global_ids_to_receive[i], -1)
					#else
					global_ids_to_receive[i]
					#endif
				);

				#ifdef DEBUG
				if (added_cells.count(global_ids_to_receive[i]) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << global_ids_to_receive[i]
						<< " has already been received from process " << this->rank
						<< std::endl;
					abort();
				}
				#endif

				this->added_cells.insert(global_ids_to_receive[i]);
			}
		}

		// receive cells in known order and add message tags
		for (
			#ifdef DCCRG_SEND_SINGLE_CELLS
			boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			#else
			boost::unordered_map<int, std::vector<uint64_t> >::iterator
			#endif
			sender = this->cells_to_receive.begin();
			sender != this->cells_to_receive.end();
			sender++
		) {
			std::sort(sender->second.begin(), sender->second.end());
			#ifdef DCCRG_SEND_SINGLE_CELLS
			// TODO: check that message tags don't overflow
			for (unsigned int i = 0; i < sender->second.size(); i++) {
				sender->second[i].second = i + 1;
			}
			#endif
		}


		/*
		Processes and the cells for which data has to be sent by this process
		*/

		// migration from user
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			pin_request = this->pin_requests.begin();
			pin_request != this->pin_requests.end();
			pin_request++
		) {
			const uint64_t current_process_of_cell = this->cell_process.at(pin_request->first);
			const uint64_t destination_process = pin_request->second;

			if (destination_process != this->rank
			&& current_process_of_cell == this->rank) {
				this->cells_to_send[destination_process].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(pin_request->first, -1)
					#else
					pin_request->first
					#endif
				);
				this->removed_cells.insert(pin_request->first);
			}
		}

		// migration from Zoltan
		if (use_zoltan) {
			for (int i = 0; i < number_to_send; i++) {

				// don't send / receive from self
				if ((uint64_t) receiver_processes[i] == this->rank) {
					continue;
				}

				// skip user-migrated cells
				if (this->pin_requests.count(global_ids_to_send[i]) > 0) {
					continue;
				}

				this->cells_to_send[receiver_processes[i]].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(global_ids_to_send[i], -1)
					#else
					global_ids_to_send[i]
					#endif
				);

				#ifdef DEBUG
				if (removed_cells.count(global_ids_to_send[i]) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << global_ids_to_send[i]
						<< " has already been sent from process " << this->rank
						<< std::endl;
					abort();
				}
				#endif

				this->removed_cells.insert(global_ids_to_send[i]);
			}

			Zoltan_LB_Free_Data(
				&global_ids_to_receive,
				&local_ids_to_receive,
				&sender_processes,
				&global_ids_to_send,
				&local_ids_to_send,
				&receiver_processes
			);
		}

		// send cells in known order and add message tags
		for (
			#ifdef DCCRG_SEND_SINGLE_CELLS
			boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			#else
			boost::unordered_map<int, std::vector<uint64_t> >::iterator
			#endif
			receiver = this->cells_to_send.begin();
			receiver != this->cells_to_send.end();
			receiver++
		) {
			std::sort(receiver->second.begin(), receiver->second.end());
			#ifdef DCCRG_SEND_SINGLE_CELLS
			// TODO: check that message tags don't overflow
			for (unsigned int i = 0; i < receiver->second.size(); i++) {
				receiver->second[i].second = i + 1;
			}
			#endif
		}
	}


	/*!
	Calculates what to send and where during a remote neighbor data update.

	Assumes up-to-date internal and user neighbor lists,
	clears previous send / receive lists.
	*/
	void recalculate_neighbor_update_send_receive_lists()
	{
		// clear previous lists
		this->cells_to_send.clear();
		this->cells_to_receive.clear();

		// only send a cell to a process once
		boost::unordered_map<
			int, // process to send to / receive from
			boost::unordered_set<uint64_t>
		> unique_cells_to_send, unique_cells_to_receive;

		// calculate new lists for neighbor data updates
		BOOST_FOREACH(const uint64_t& cell, this->cells_with_remote_neighbors) {

			#ifdef DEBUG
			if (cell != this->get_child(cell)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell << " has children"
					<< std::endl;
				abort();
			}
			#endif

			uint64_t current_process = this->rank;

			// data must be received from neighbors_of
			BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(cell)) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) != current_process) {
					unique_cells_to_receive[this->cell_process.at(neighbor)].insert(neighbor);
				}
			}

			// data must be sent to neighbors_to
			BOOST_FOREACH(const uint64_t& neighbor, this->neighbors_to.at(cell)) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) != current_process) {
					unique_cells_to_send[this->cell_process.at(neighbor)].insert(cell);
				}
			}
		}

		// populate final send list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			receiver = unique_cells_to_send.begin();
			receiver != unique_cells_to_send.end();
			receiver++
		) {
			this->cells_to_send[receiver->first].reserve(receiver->second.size());

			BOOST_FOREACH(const uint64_t& cell, receiver->second) {
				this->cells_to_send[receiver->first].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(cell, -1)
					#else
					cell
					#endif
				);
			}

			std::sort(
				this->cells_to_send[receiver->first].begin(),
				this->cells_to_send[receiver->first].end()
			);

			#ifdef DCCRG_SEND_SINGLE_CELLS
			// sequential tags for messages: 1, 2, ...
			for (unsigned int i = 0; i < this->cells_to_send[receiver->first].size(); i++) {
				this->cells_to_send[receiver->first][i].second = i + 1;
			}
			#endif
		}

		// populate final receive list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			sender = unique_cells_to_receive.begin();
			sender != unique_cells_to_receive.end();
			sender++
		) {
			this->cells_to_receive[sender->first].reserve(sender->second.size());

			BOOST_FOREACH(const uint64_t& cell, sender->second) {

				this->cells_to_receive[sender->first].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(cell, -1)
					#else
					cell
					#endif
				);
			}

			std::sort(
				this->cells_to_receive[sender->first].begin(),
				this->cells_to_receive[sender->first].end()
			);

			#ifdef DCCRG_SEND_SINGLE_CELLS
			// sequential tags for messages: 1, 2, ...
			for (unsigned int i = 0; i < this->cells_to_receive[sender->first].size(); i++) {
				this->cells_to_receive[sender->first][i].second = i + 1;
			}
			#endif
		}

		for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
			item = this->user_hood_of.begin();
			item != this->user_hood_of.end();
			item++
		) {
			this->recalculate_neighbor_update_send_receive_lists(item->first);
		}
	}

	/*!
	Same as the version without an id but for user defined neighborhood.

	Updates send/receive lists of cells using only the given
	neighborhood id.
	*/
	void recalculate_neighbor_update_send_receive_lists(const int id)
	{
		// clear previous lists
		this->user_neigh_cells_to_send[id].clear();
		this->user_neigh_cells_to_receive[id].clear();

		boost::unordered_map<int, boost::unordered_set<uint64_t> >
			user_neigh_unique_sends,
			user_neigh_unique_receives;

		// calculate new lists for neighbor data updates
		BOOST_FOREACH(const uint64_t cell, this->cells_with_remote_neighbors) {

			#ifdef DEBUG
			if (cell != this->get_child(cell)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell << " has children"
					<< std::endl;
				abort();
			}
			#endif

			// data must be received from neighbors_of
			BOOST_FOREACH(const uint64_t& neighbor, this->user_neigh_of.at(id).at(cell)) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					user_neigh_unique_receives[this->cell_process.at(neighbor)].insert(neighbor);
				}
			}

			// data must be sent to neighbors_to
			BOOST_FOREACH(const uint64_t& neighbor, this->user_neigh_to.at(id).at(cell)) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					user_neigh_unique_sends[this->cell_process.at(neighbor)].insert(cell);
				}
			}
		}

		// populate final send list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			receiver = user_neigh_unique_sends.begin();
			receiver != user_neigh_unique_sends.end();
			receiver++
		) {
			this->user_neigh_cells_to_send[id][receiver->first].reserve(receiver->second.size());

			BOOST_FOREACH(const uint64_t& cell, receiver->second) {
				this->user_neigh_cells_to_send.at(id)[receiver->first].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(cell, -1)
					#else
					cell
					#endif
				);
			}

			std::sort(
				this->user_neigh_cells_to_send.at(id)[receiver->first].begin(),
				this->user_neigh_cells_to_send.at(id)[receiver->first].end()
			);

			#ifdef DCCRG_SEND_SINGLE_CELLS
			// sequential tags for messages: 1, 2, ...
			for (unsigned int i = 0; i < this->cells_to_send[receiver->first].size(); i++) {
				this->user_neigh_cells_to_send.at(id)[receiver->first][i].second = i + 1;
			}
			#endif
		}

		// populate final receive list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			sender = user_neigh_unique_receives.begin();
			sender != user_neigh_unique_receives.end();
			sender++
		) {
			this->user_neigh_cells_to_receive[id][sender->first].reserve(sender->second.size());

			BOOST_FOREACH(const uint64_t& cell, sender->second) {

				this->user_neigh_cells_to_receive.at(id)[sender->first].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(cell, -1)
					#else
					cell
					#endif
				);
			}

			std::sort(
				this->user_neigh_cells_to_receive.at(id)[sender->first].begin(),
				this->user_neigh_cells_to_receive.at(id)[sender->first].end()
			);

			#ifdef DCCRG_SEND_SINGLE_CELLS
			// sequential tags for messages: 1, 2, ...
			for (unsigned int i = 0; i < this->cells_to_receive[sender->first].size(); i++) {
				this->user_neigh_cells_to_receive.at(id)[sender->first][i].second = i + 1;
			}
			#endif
		}
	}


	/*!
	Updates neighbor and neighbor_to lists around given cell's neighborhood.

	Does nothing in the following cases:
		- given cell doesn't exist in the grid
		- given cell has children
	Assumes that the refinement level difference between given cell and its neighborhood is no larger than 1.
	*/
	void update_neighbors(const uint64_t cell)
	{
		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		this->neighbors.at(cell) = this->find_neighbors_of(cell, this->neighborhood_of, this->max_ref_lvl_diff);
		this->neighbors_to.at(cell) = this->find_neighbors_to(cell, this->neighbors.at(cell));

		#ifdef DEBUG
		if (
			!this->verify_neighbors(
				cell,
				this->neighborhood_of,
				this->neighborhood_to,
				this->neighbors,
				this->neighbors_to
			)
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Neighbor update failed for cell " << cell
				<< " (child of " << this->get_parent(cell) << ")"
				<< std::endl;
			abort();
		}
		#endif
	}


	/*!
	Updates the neighbors and _to of given cell based on given neighborhood.

	Does nothing in the following cases:
		- given cell doesn't exist in the grid
		- given cell has children
	Assumes that update_neighbors(cell) has been called prior to this.
	*/
	void update_user_neighbors(const uint64_t cell, const int id)
	{
		if (this->user_hood_of.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No user neighborhood with id " << id
				<< std::endl;
			abort();
		}

		#ifdef DEBUG
		if (this->user_hood_to.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No user neighborhood to with id " << id
				<< std::endl;
			abort();
		}
		#endif

		// find neighbors_of, should be in order given by user
		this->user_neigh_of[id][cell].clear();
		BOOST_FOREACH(const Types<3>::neighborhood_item_t& item, this->user_hood_of[id]) {
			std::vector<uint64_t> cells_at_offset
				= this->get_neighbors_of(cell, item[0], item[1], item[2]);
			this->user_neigh_of[id][cell].insert(
				this->user_neigh_of[id][cell].end(),
				cells_at_offset.begin(),
				cells_at_offset.end()
			);
		}

		// find neighbors_to
		this->user_neigh_to[id][cell]
			= this->find_neighbors_to(cell, this->user_hood_to[id]);
	}


	/*!
	Updates the remote neighbor info of given cell on this process without children.

	Uses current neighbor lists.
	Does nothing if given cell doesn't exist on this process or has children
	*/
	void update_remote_neighbor_info(const uint64_t cell)
	{
		if (this->cells.count(cell) == 0) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		// TODO: also update remote_cells_with_local_neighbors
		this->cells_with_remote_neighbors.erase(cell);

		#ifdef DEBUG
		if (this->neighbors.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Neighbor list for cell " << cell
				<< " doesn't exist"
				<< std::endl;
			abort();
		}

		if (this->neighbors_to.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Neighbors_to list for cell " << cell
				<< " doesn't exist"
				<< std::endl;
			abort();
		}
		#endif

		// neighbors of given cell
		BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(cell)) {

			if (neighbor == 0) {
				continue;
			}

			#ifdef DEBUG
			if (this->cell_process.count(neighbor) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor " << neighbor
					<< " doesn't exist in cell_process"
					<< std::endl;
				abort();
			}
			#endif

			if (this->cell_process.at(neighbor) != this->rank) {
				this->cells_with_remote_neighbors.insert(cell);
				this->remote_cells_with_local_neighbors.insert(neighbor);
			}
		}
		// cells with given cell as neighbor
		BOOST_FOREACH(const uint64_t& neighbor_to, this->neighbors_to.at(cell)) {

			#ifdef DEBUG
			if (this->cell_process.count(neighbor_to) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor_to " << neighbor_to
					<< " doesn't exist in cell_process"
					<< std::endl;
				abort();
			}
			#endif

			if (this->cell_process.at(neighbor_to) != this->rank) {
				this->cells_with_remote_neighbors.insert(cell);
				this->remote_cells_with_local_neighbors.insert(neighbor_to);
			}
		}

		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info(cell)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbor info for cell " << cell << " is not consistent" << std::endl;
			abort();
		}
		#endif
	}


	/*!
	Updates the remote neighbor info of all cells on this process without children.

	Uses current neighbor lists.
	*/
	void update_remote_neighbor_info()
	{
		// TODO this probably can't be optimized without storing neighbor lists also for remote neighbors
		this->cells_with_remote_neighbors.clear();
		this->remote_cells_with_local_neighbors.clear();

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			if (item.first != this->get_child(item.first)) {
				continue;
			}

			this->update_remote_neighbor_info(item.first);

			#ifdef DEBUG
			if (!this->verify_remote_neighbor_info(item.first)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Remote neighbor info for cell " << item.first
					<< " is not consistent"
					<< std::endl;
				abort();
			}
			#endif
		}

		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info is not consistent"
				<< std::endl;
			abort();
		}
		#endif
	}


	/*!
	Returns true if cell1 considers cell2 as a neighbor, even if neither of them exists
	*/
	bool is_neighbor(const uint64_t cell1, const uint64_t cell2) const
	{
		#ifdef DEBUG
		if (cell1 == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell1 given." << std::endl;
			abort();
		}

		if (cell1 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Impossible cell1 given." << std::endl;
			abort();
		}

		if (cell2 == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell2 given." << std::endl;
			abort();
		}

		if (cell2 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Impossible cell2 given." << std::endl;
			abort();
		}
		#endif

		const Types<3>::indices_t indices1 = this->get_indices(cell1);
		const Types<3>::indices_t indices2 = this->get_indices(cell2);
		const uint64_t cell1_size = this->get_cell_size_in_indices(cell1);
		const uint64_t cell2_size = this->get_cell_size_in_indices(cell2);

		// distance in indices between given cells
		Types<3>::indices_t distance = {{0, 0, 0}};

		const uint64_t grid_length[3] = {
			this->get_x_length() * (uint64_t(1) << this->max_refinement_level),
			this->get_y_length() * (uint64_t(1) << this->max_refinement_level),
			this->get_z_length() * (uint64_t(1) << this->max_refinement_level)
		};

		uint64_t max_distance = 0;

		for (unsigned int i = 0; i < 3; i++) {
			if (indices1[i] <= indices2[i]) {
				if (indices2[i] <= indices1[i] + cell1_size) {
					distance[i] = 0;
				} else {
					distance[i] = indices2[i] - (indices1[i] + cell1_size);
				}

				if (this->is_periodic(i)) {
					const uint64_t distance_to_end = grid_length[i] - (indices2[i] + cell2_size);
					distance[i] = std::min(distance[i], indices1[i] + distance_to_end);
				}
			} else {
				if (indices1[i] <= indices2[i] + cell2_size) {
					distance[i] = 0;
				} else {
					distance[i] = indices1[i] - (indices2[i] + cell2_size);
				}

				if (this->is_periodic(i)) {
					const uint64_t distance_to_end = grid_length[i] - (indices1[i] + cell1_size);
					distance[i] = std::min(distance[i], indices2[i] + distance_to_end);
				}
			}

			max_distance = std::max(max_distance, distance[i]);
		}

		if (this->neighborhood_size == 0) {
			if (max_distance < cell1_size
			&& this->overlapping_indices(cell1, cell2) >= 2) {
				return true;
			// diagonal cell isn't a neighbor
			} else {
				return false;
			}
		}

		if (max_distance < this->neighborhood_size * cell1_size) {
			return true;
		} else {
			return false;
		}
	}


	/*!
	Given a cell that exists and has children returns one of the children.

	Returns the given cell if it doesn't have children or error_cell if the cell doesn't exist.
	*/
	uint64_t get_child(const uint64_t cell) const
	{
		if (this->cell_process.count(cell) == 0) {
			return error_cell;
		}

		const int refinement_level = this->get_refinement_level(cell);

		// given cell cannot have children
		if (refinement_level == this->max_refinement_level) {
			return cell;
		}

		const uint64_t child = this->get_cell_from_indices(
			this->get_indices(cell),
			refinement_level + 1
		);

		if (this->cell_process.count(child) > 0) {
			return child;
		} else {
			return cell;
		}
	}


	/*!
	Enforces maximum refinement level difference between neighbors.

	Adds new cells to cells_to_refine in order to enforce maximum refinement
	level difference of max_ref_lvl_diff between neighbors (also across processes).
	After this function cells_to_refine will contain the refines of all processes.
	*/
	void induce_refines()
	{
		std::vector<uint64_t> new_refines(this->cells_to_refine.begin(), this->cells_to_refine.end());
		MPI_Comm non_boost_comm = this->comm;
		while (All_Reduce()(new_refines.size(), non_boost_comm) > 0) {

			std::vector<std::vector<uint64_t> > all_new_refines;
			All_Gather()(new_refines, all_new_refines, this->comm);
			new_refines.clear();

			boost::unordered_set<uint64_t> unique_induced_refines;

			// induced refines on this process
			BOOST_FOREACH(const uint64_t& refined, all_new_refines.at(this->rank)) {

				// refine local neighbors that are too large
				BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(refined)) {

					if (neighbor == 0) {
						continue;
					}

					#ifdef DEBUG
					if (this->cell_process.count(neighbor) == 0) {
						std::cerr << "Process " << this->rank
							<< ": Cell " << refined
							<< " had a non-existing neighbor in neighbor list: " << neighbor
							<< std::endl;
					}
					#endif

					if (this->cell_process.at(neighbor) != this->rank) {
						continue;
					}

					if (this->get_refinement_level(neighbor) < this->get_refinement_level(refined)) {
						if (this->cells_to_refine.count(neighbor) == 0) {
							unique_induced_refines.insert(neighbor);
						}
					}
				}

				BOOST_FOREACH(const uint64_t& neighbor_to, this->neighbors_to.at(refined)) {

					if (neighbor_to == 0) {
						continue;
					}

					#ifdef DEBUG
					if (this->cell_process.count(neighbor_to) == 0) {
						std::cerr << "Process " << this->rank
							<< ": Cell " << refined
							<< " had a non-existing neighbor in neighbor list: " << neighbor_to
							<< std::endl;
					}
					#endif

					if (this->cell_process.at(neighbor_to) != this->rank) {
						continue;
					}

					if (this->get_refinement_level(neighbor_to) < this->get_refinement_level(refined)) {
						if (this->cells_to_refine.count(neighbor_to) == 0) {
							unique_induced_refines.insert(neighbor_to);
						}
					}
				}
			}

			// refines induced here by other processes
			for (unsigned int process = 0; process < this->comm_size; process++) {

				if (process == this->rank) {
					continue;
				}

				BOOST_FOREACH(const uint64_t& refined, all_new_refines.at(process)) {

					if (this->remote_cells_with_local_neighbors.count(refined) == 0) {
						continue;
					}

					// refine all local cells that are too large and neighboring the refined cell
					/*
					TODO: probably faster to search for local neighbors of refined
					cell, even faster would be to also store neighbors lists of
					remote cells with local neighbors
					*/
					BOOST_FOREACH(const uint64_t& local, this->cells_with_remote_neighbors) {

						if (this->is_neighbor(local, refined)
						&& this->get_refinement_level(local) < this->get_refinement_level(refined)
						&& this->cells_to_refine.count(local) == 0) {
							unique_induced_refines.insert(local);
						}
					}
				}
			}
			all_new_refines.clear();

			new_refines.insert(new_refines.end(), unique_induced_refines.begin(), unique_induced_refines.end());
			this->cells_to_refine.insert(unique_induced_refines.begin(), unique_induced_refines.end());
			unique_induced_refines.clear();
		}

		// add refines from all processes to cells_to_refine
		std::vector<uint64_t> refines(this->cells_to_refine.begin(), this->cells_to_refine.end());
		std::vector<std::vector<uint64_t> > all_refines;
		All_Gather()(refines, all_refines, this->comm);

		for (unsigned int process = 0; process < this->comm_size; process++) {
			this->cells_to_refine.insert(all_refines[process].begin(), all_refines[process].end());
		}

		#ifdef DEBUG
		// check that all required refines have been induced
		BOOST_FOREACH(const uint64_t& refined, this->cells_to_refine) {

			// neighbors_of
			std::vector<uint64_t> neighbors_of
				= this->find_neighbors_of(refined, this->neighborhood_of, this->max_ref_lvl_diff);

			BOOST_FOREACH(const uint64_t& neighbor_of, neighbors_of) {

				if (neighbor_of == 0) {
					continue;
				}

				if (this->get_refinement_level(neighbor_of) < this->get_refinement_level(refined)
				&& this->cells_to_refine.count(neighbor_of) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbor (" << neighbor_of
						<< ") of cell that will be refined (" << refined
						<< ", ref lvl " << this->get_refinement_level(refined)
						<< ") has too small refinement level: " << this->get_refinement_level(neighbor_of)
						<< std::endl;
					abort();
				}
			}

			// neighbors_to
			std::vector<uint64_t> neighbors_to
				= this->find_neighbors_to(refined, this->neighborhood_to);
			BOOST_FOREACH(const uint64_t& neighbor_to, neighbors_to) {

				if (neighbor_to == 0) {
					continue;
				}

				if (this->get_refinement_level(neighbor_to) < this->get_refinement_level(refined)
				&& this->cells_to_refine.count(neighbor_to) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbor (" << neighbor_to
						<< ") of cell that will be refined (" << refined
						<< ", ref lvl " << this->get_refinement_level(refined)
						<< ") has too small refinement level: " << this->get_refinement_level(neighbor_to)
						<< std::endl;
					abort();
				}
			}
		}

		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Grid isn't consistent" << std::endl;
			abort();
		}
		#endif
	}


	/*!
	Sends the numbers in s to all other processes and adds the numbers sent by all others to s.
	*/
	void all_to_all_set(boost::unordered_set<uint64_t>& s)
	{
		std::vector<uint64_t> local_s(s.begin(), s.end());

		std::vector<std::vector<uint64_t> > all_s;
		All_Gather()(local_s, all_s, this->comm);

		BOOST_FOREACH(const std::vector<uint64_t>& i, all_s) {
			BOOST_FOREACH(const uint64_t& cell, i) {
				s.insert(cell);
			}
		}
	}


	/*!
	Overrides local unrefines based on global refines.

	Removes cells from cells_to_unrefine in order to enforce maximum refinement level
	difference of one between neighbors.
	cells_to_refine and cells_not_to_unrefine must be identical between processes.
	After this function cells_to_unrefine will contain the unrefines of all processes.
	*/
	void override_unrefines()
	{
		// unrefines that were not overridden
		boost::unordered_set<uint64_t> final_unrefines;

		// don't unrefine if...
		BOOST_FOREACH(const uint64_t& unrefined, this->cells_to_unrefine) {

			bool can_unrefine = true;

			// ...any sibling cannot be
			const uint64_t parent = this->get_parent(unrefined);
			const std::vector<uint64_t> siblings = this->get_all_children(parent);

			BOOST_FOREACH(const uint64_t& sibling, siblings) {
				if (this->cells_to_refine.count(sibling) > 0
				|| this->cells_not_to_unrefine.count(sibling) > 0) {
					can_unrefine = false;
					break;
				}
			}

			if (!can_unrefine) {
				continue;
			}

			// ...parent of unrefined wouldn't fulfill requirements
			const int refinement_level = this->get_refinement_level(parent);

			#ifdef DEBUG
			if (parent == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent" << std::endl;
				abort();
			}

			if (refinement_level < 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid refinement level for parent" << std::endl;
				abort();
			}
			#endif

			const std::vector<uint64_t> neighbors
				= this->find_neighbors_of(
					parent,
					this->neighborhood_of,
					2 * this->max_ref_lvl_diff,
					true
				);

			BOOST_FOREACH(const uint64_t& neighbor, neighbors) {

				const int neighbor_ref_lvl = this->get_refinement_level(neighbor);

				if (neighbor_ref_lvl == refinement_level + this->max_ref_lvl_diff
				&& this->cells_to_refine.count(neighbor) > 0) {
					can_unrefine = false;
					break;
				}
			}

			if (can_unrefine) {
				final_unrefines.insert(unrefined);
			}
		}
		this->cells_to_unrefine.clear();

		// add unrefines from all processes to cells_to_unrefine
		std::vector<uint64_t> unrefines(final_unrefines.begin(), final_unrefines.end());
		std::vector<std::vector<uint64_t> > all_unrefines;
		All_Gather()(unrefines, all_unrefines, this->comm);

		for (unsigned int process = 0; process < this->comm_size; process++) {
			this->cells_to_unrefine.insert(all_unrefines[process].begin(), all_unrefines[process].end());
		}

		#ifdef DEBUG
		// check that maximum refinement level difference between future neighbors <= 1
		BOOST_FOREACH(const uint64_t& unrefined, this->cells_to_unrefine) {

			if (unrefined != this->get_child(unrefined)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << unrefined
					<< " has children"
					<< std::endl;
				abort();
			}

			if (this->cell_process.count(unrefined) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << unrefined
					<< " to be unrefined doesn't exist"
					<< std::endl;
				abort();
			}

			if (this->cell_process.at(unrefined) == this->rank
			&& this->cells.count(unrefined) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << unrefined
					<< " to be unrefined has no data"
					<< std::endl;
				abort();
			}

			const int ref_lvl = this->get_refinement_level(unrefined);

			// neighbors_of
			const std::vector<uint64_t> neighbors
				= this->find_neighbors_of(
					this->get_parent(unrefined),
					this->neighborhood_of,
					2 * this->max_ref_lvl_diff,
					true
				);

			BOOST_FOREACH(const uint64_t& neighbor, neighbors) {

				if (neighbor == 0) {
					continue;
				}

				const int neighbor_ref_lvl = this->get_refinement_level(neighbor);

				if (neighbor_ref_lvl > ref_lvl) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbor " << neighbor
						<< " of cell that will be unrefined (" << unrefined
						<< ", ref lvl " << ref_lvl
						<< ") has too large refinement level: " << neighbor_ref_lvl
						<< std::endl;
					abort();
				}

				if (neighbor_ref_lvl == ref_lvl
				&& this->cells_to_refine.count(neighbor) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbor " << neighbor
						<< " of cell that will be unrefined (" << unrefined
						<< ", ref lvl " << ref_lvl
						<< ") is identical in size and will be refined"
						<< std::endl;
					abort();
				}
			}
		}

		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Grid isn't consistent" << std::endl;
			abort();
		}
		#endif
	}


	/*!
	Adds refined cells to the grid, removes unrefined cells from the grid.

	cells_to_refine and cells_to_unrefine must contain the cells to refine/unrefine of all processes.
	Returns new cells created on this process by refinement.
	Moves unrefined cell data to the process of their parent.
	*/
	std::vector<uint64_t> execute_refines()
	{
		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbor info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " User data is inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		std::vector<uint64_t> new_cells;

		this->remote_neighbors.clear();
		this->cells_to_send.clear();
		this->cells_to_receive.clear();
		this->refined_cell_data.clear();
		this->unrefined_cell_data.clear();
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->incoming_data.clear();
		this->outgoing_data.clear();
		#endif
		#endif

		#ifdef DEBUG
		// check that cells_to_refine is identical between processes
		std::vector<uint64_t> ordered_cells_to_refine(this->cells_to_refine.begin(), this->cells_to_refine.end());
		std::sort(ordered_cells_to_refine.begin(), ordered_cells_to_refine.end());

		std::vector<std::vector<uint64_t> > all_ordered_cells_to_refine;
		All_Gather()(ordered_cells_to_refine, all_ordered_cells_to_refine, this->comm);

		for (unsigned int process = 0; process < this->comm_size; process++) {
			if (!std::equal(all_ordered_cells_to_refine[process].begin(), all_ordered_cells_to_refine[process].end(), all_ordered_cells_to_refine[0].begin())) {
				std::cerr << __FILE__ << ":" << __LINE__ << " cells_to_refine differ between processes 0 and " << process << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// check that cells_to_unrefine is identical between processes
		std::vector<uint64_t> ordered_cells_to_unrefine(this->cells_to_unrefine.begin(), this->cells_to_unrefine.end());
		std::sort(ordered_cells_to_unrefine.begin(), ordered_cells_to_unrefine.end());

		std::vector<std::vector<uint64_t> > all_ordered_cells_to_unrefine;
		All_Gather()(ordered_cells_to_unrefine, all_ordered_cells_to_unrefine, this->comm);

		for (unsigned int process = 0; process < this->comm_size; process++) {
			if (!std::equal(all_ordered_cells_to_unrefine[process].begin(), all_ordered_cells_to_unrefine[process].end(), all_ordered_cells_to_unrefine[0].begin())) {
				std::cerr << __FILE__ << ":" << __LINE__ << " cells_to_unrefine differ between processes 0 and " << process << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		#endif

		// cells whose neighbor lists have to be updated afterwards
		boost::unordered_set<uint64_t> update_neighbors;

		// a separate neighborhood update function has to be used for cells whose children were removed by unrefining
		boost::unordered_set<uint64_t> update_neighbors_unrefined;

		// refines
		BOOST_FOREACH(const uint64_t& refined, this->cells_to_refine) {

			#ifdef DEBUG
			if (this->cell_process.count(refined) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << refined
					<< " doesn't exist"
					<< std::endl;
				abort();
			}

			if (this->rank == this->cell_process.at(refined)
			&& this->cells.count(refined) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Data for cell " << refined
					<< " doesn't exist"
					<< std::endl;
				abort();
			}


			if (this->cell_process.at(refined) == this->rank
			&& this->neighbors.count(refined) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor list for cell " << refined
					<< " doesn't exist"
					<< std::endl;
				abort();
			}

			if (this->cell_process.at(refined) == this->rank
			&& this->neighbors_to.count(refined) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor_to list for cell " << refined
					<< " doesn't exist"
					<< std::endl;
				abort();
			}
			#endif

			const uint64_t process_of_refined = this->cell_process.at(refined);

			// move user data of refined cells into refined_cell_data
			if (this->rank == process_of_refined) {
				// TODO: move data instead of copying, using boost::move or c++0x move?
				this->refined_cell_data[refined] = this->cells.at(refined);
				this->cells.erase(refined);
			}

			// add children of refined cells into the grid
			const std::vector<uint64_t> children = this->get_all_children(refined);
			BOOST_FOREACH(const uint64_t& child, children) {
				this->cell_process[child] = process_of_refined;

				if (this->rank == process_of_refined) {
					this->cells[child];
					this->neighbors[child];
					this->neighbors_to[child];
					new_cells.push_back(child);
				}
			}

			// children of refined cells inherit their pin request status
			if (this->pin_requests.count(refined) > 0) {
				BOOST_FOREACH(const uint64_t& child, children) {
					this->pin_requests[child] = this->pin_requests.at(refined);
				}
				this->pin_requests.erase(refined);
			}
			if (this->new_pin_requests.count(refined) > 0) {
				BOOST_FOREACH(const uint64_t& child, children) {
					this->new_pin_requests[child] = this->new_pin_requests.at(refined);
				}
				this->new_pin_requests.erase(refined);
			}

			// children of refined cells inherit their weight
			if (this->rank == process_of_refined
			&& this->cell_weights.count(refined) > 0) {
				BOOST_FOREACH(const uint64_t& child, children) {
					this->cell_weights[child] = this->cell_weights.at(refined);
				}
				this->cell_weights.erase(refined);
			}

			// use local neighbor lists to find cells whose neighbor lists have to updated
			if (this->rank == process_of_refined) {
				// update the neighbor lists of created local cells
				BOOST_FOREACH(const uint64_t& child, children) {
					update_neighbors.insert(child);
				}

				// update neighbor lists of all the parent's neighbors
				BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(refined)) {
					if (neighbor == 0) {
						continue;
					}

					if (this->cell_process.at(neighbor) == this->rank) {
						update_neighbors.insert(neighbor);
					}
				}

				BOOST_FOREACH(const uint64_t& neighbor_to, this->neighbors_to.at(refined)) {
					if (this->cell_process.at(neighbor_to) == this->rank) {
						update_neighbors.insert(neighbor_to);
					}
				}
			}

			// without using local neighbor lists figure out rest of the neighbor lists that need updating
			if (this->remote_cells_with_local_neighbors.count(refined) > 0) {

				/*
				No need to update local neighbors_to of refined cell, if they are larger
				they will also be refined and updated.
				*/
				const std::vector<uint64_t> neighbors
					= this->find_neighbors_of(
						refined,
						this->neighborhood_of,
						2 * this->max_ref_lvl_diff,
						true
					);

				BOOST_FOREACH(const uint64_t& neighbor, neighbors) {
					if (neighbor == 0) {
						continue;
					}

					if (this->is_local(neighbor)) {
						update_neighbors.insert(neighbor);
					}
				}
			}
		}

		// needed for checking which neighborhoods to update due to unrefining
		boost::unordered_set<uint64_t> parents_of_unrefined;

		// initially only one sibling is recorded per process when unrefining, insert the rest of them now
		boost::unordered_set<uint64_t> all_to_unrefine;
		BOOST_FOREACH(const uint64_t& unrefined, this->cells_to_unrefine) {

			const uint64_t parent_of_unrefined = this->get_parent(unrefined);
			#ifdef DEBUG
			if (unrefined != this->get_child(unrefined)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << unrefined
					<< " has children"
					<< std::endl;
				abort();
			}

			if (parent_of_unrefined == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent cell" << std::endl;
				abort();
			}

			if (parent_of_unrefined == unrefined) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << unrefined
					<< " has no parent"
					<< std::endl;
				abort();
			}
			#endif

			parents_of_unrefined.insert(parent_of_unrefined);

			const std::vector<uint64_t> siblings = this->get_all_children(parent_of_unrefined);

			#ifdef DEBUG
			bool unrefined_in_siblings = false;
			BOOST_FOREACH(const uint64_t& sibling, siblings) {

				if (this->cell_process.count(sibling) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << sibling
						<< " doesn't exist"
						<< std::endl;
					abort();
				}

				if (sibling != this->get_child(sibling)) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << sibling
						<< " has has children"
						<< std::endl;
					abort();
				}

				if (this->cell_process.at(sibling) == this->rank
				&& this->cells.count(sibling) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << sibling
						<< " has no data"
						<< std::endl;
					abort();
				}

				if (unrefined == sibling) {
					unrefined_in_siblings = true;
				}
			}

			if (!unrefined_in_siblings) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell to unrefine isn't its parent's child" << std::endl;
				abort();
			}
			#endif

			all_to_unrefine.insert(siblings.begin(), siblings.end());
		}

		// unrefines
		BOOST_FOREACH(const uint64_t& unrefined, all_to_unrefine) {

			const uint64_t parent_of_unrefined = this->get_parent(unrefined);
			#ifdef DEBUG
			if (parent_of_unrefined == unrefined) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << unrefined
					<< " has no parent"
					<< std::endl;
				abort();
			}
			#endif

			const uint64_t process_of_parent = this->cell_process.at(parent_of_unrefined);
			const uint64_t process_of_unrefined = this->cell_process.at(unrefined);

			// remove unrefined cells and their siblings from the grid, but don't remove user data yet
			this->cell_process.erase(unrefined);
			update_neighbors.erase(unrefined);
			this->pin_requests.erase(unrefined);
			this->new_pin_requests.erase(unrefined);
			this->cell_weights.erase(unrefined);

			// don't send unrefined cells' user data to self
			if (this->rank == process_of_unrefined
			&& this->rank == process_of_parent) {

				#ifdef DEBUG
				if (this->cells.count(unrefined) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << unrefined
						<< " to be unrefined has no data"
						<< std::endl;
					abort();
				}
				#endif

				// TODO move data instead of copying
				this->unrefined_cell_data[unrefined] = this->cells.at(unrefined);
				this->cells.erase(unrefined);

			// send user data of removed cell to the parent's process
			} else if (this->rank == process_of_unrefined) {

				this->cells_to_send[process_of_parent].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(unrefined, -1)
					#else
					unrefined
					#endif
				);

			// receive user data of removed cell from its process
			} else if (this->rank == process_of_parent) {
				this->cells_to_receive[process_of_unrefined].push_back(
					#ifdef DCCRG_SEND_SINGLE_CELLS
					std::make_pair(unrefined, -1)
					#else
					unrefined
					#endif
				);
			}
		}

		// receive cells in known order and add message tags
		for (
			#ifdef DCCRG_SEND_SINGLE_CELLS
			boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			#else
			boost::unordered_map<int, std::vector<uint64_t> >::iterator
			#endif
			sender = this->cells_to_receive.begin();
			sender != this->cells_to_receive.end();
			sender++
		) {
			std::sort(sender->second.begin(), sender->second.end());
			#ifdef DCCRG_SEND_SINGLE_CELLS
			// TODO: merge with identical code in make_new_partition
			for (unsigned int i = 0; i < sender->second.size(); i++) {
				 sender->second[i].second = i + 1;
			}
			#endif
		}

		// send cells in known order and add message tags
		for (
			#ifdef DCCRG_SEND_SINGLE_CELLS
			boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			#else
			boost::unordered_map<int, std::vector<uint64_t> >::iterator
			#endif
			receiver = this->cells_to_send.begin();
			receiver != this->cells_to_send.end();
			receiver++
		) {
			std::sort(receiver->second.begin(), receiver->second.end());
			#ifdef DCCRG_SEND_SINGLE_CELLS
			// TODO: check that message tags don't overflow
			for (unsigned int i = 0; i < receiver->second.size(); i++) {
				receiver->second[i].second = i + 1;
			}
			#endif
		}


		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->unrefined_cell_data, this->cells_to_receive, this->cells_to_send
		#elif defined (DCCRG_CELL_DATA_SIZE_FROM_USER)
		this->unrefined_cell_data, this->cells_to_receive, this->cells_to_send
		#else
		this->cells_to_receive, this->cells_to_send
		#endif
		);

		// update data for parents (and their neighborhood) of unrefined cells
		BOOST_FOREACH(const uint64_t& parent, parents_of_unrefined) {

			/* TODO: skip unrefined cells far enough away
			std::vector<uint64_t> children = this->get_all_children(*parent);
			*/

			#ifdef DEBUG
			if (this->cell_process.count(parent) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Parent " << parent
					<< " doesn't exist"
					<< std::endl;
				abort();
			}

			if (parent != this->get_child(parent)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Parent " << parent
					<< " still has children"
					<< std::endl;
				abort();
			}
			#endif

			const std::vector<uint64_t> new_neighbors_of
				= this->find_neighbors_of(parent, this->neighborhood_of, this->max_ref_lvl_diff);

			BOOST_FOREACH(const uint64_t& neighbor, new_neighbors_of) {

				if (neighbor == 0) {
					continue;
				}

				if (this->cell_process.at(neighbor) == this->rank) {
					update_neighbors.insert(neighbor);
				}
			}

			const std::vector<uint64_t> new_neighbors_to
				= this->find_neighbors_to(parent, this->neighborhood_to);
			BOOST_FOREACH(const uint64_t& neighbor, new_neighbors_to) {
				if (this->cell_process.at(neighbor) == this->rank) {
					update_neighbors.insert(neighbor);
				}
			}

			// add user data and neighbor lists of local parents of unrefined cells
			if (this->cell_process.at(parent) == this->rank) {
				this->cells[parent];
				this->neighbors[parent] = new_neighbors_of;
				this->neighbors_to[parent] = new_neighbors_to;

				// add user neighbor lists
				for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
					item = this->user_hood_of.begin();
					item != this->user_hood_of.end();
					item++
				) {
					this->update_user_neighbors(parent, item->first);
				}
			}
		}

		// update neighbor lists of cells affected by refining / unrefining
		BOOST_FOREACH(const uint64_t& cell, update_neighbors) {
			this->update_neighbors(cell);
			//update also neighbor lists of user neighborhoods
			for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
				item = this->user_hood_of.begin();
				item != this->user_hood_of.end();
				item++
			) {
				this->update_user_neighbors(cell, item->first);
			}
		}

		// remove neighbor lists of added cells' parents
		BOOST_FOREACH(const uint64_t& refined, this->cells_to_refine) {

			if (this->cell_process.at(refined) == this->rank) {

				#ifdef DEBUG
				if (this->neighbors.count(refined) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbor list for cell " << refined
						<< " doesn't exist"
						<< std::endl;
					abort();
				}

				if (this->neighbors_to.count(refined) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbor_to list for cell " << refined
						<< " doesn't exist"
						<< std::endl;
					abort();
				}
				#endif

				this->neighbors.erase(refined);
				this->neighbors_to.erase(refined);

				// remove also from user's neighborhood
				for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
					item = this->user_hood_of.begin();
					item != this->user_hood_of.end();
					item++
				) {
					this->user_neigh_of.at(item->first).erase(refined);
					this->user_neigh_to.at(item->first).erase(refined);
				}
			}
		}

		// remove neighbor lists of removed cells
		BOOST_FOREACH(const uint64_t& unrefined, all_to_unrefine) {
			this->neighbors.erase(unrefined);
			this->neighbors_to.erase(unrefined);
			// also from user neighborhood
			for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
				item = this->user_hood_of.begin();
				item != this->user_hood_of.end();
				item++
			) {
				this->user_neigh_of.at(item->first).erase(unrefined);
				this->user_neigh_to.at(item->first).erase(unrefined);
			}
		}

		this->update_remote_neighbor_info();

		#ifdef DEBUG
		if (!this->verify_neighbors()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbor lists are inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->unrefined_cell_data, this->cells_to_receive
		#endif
		#endif
		);
		this->wait_user_data_transfer_sends();
		this->cells_to_send.clear();
		this->cells_to_receive.clear();

		#ifdef DEBUG
		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " virhe" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// remove user data of unrefined cells from this->cells
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = all_to_unrefine.begin(); unrefined != all_to_unrefine.end(); unrefined++) {
			this->cells.erase(*unrefined);
		}

		this->cells_to_refine.clear();
		this->cells_to_unrefine.clear();

		this->recalculate_neighbor_update_send_receive_lists();

		return new_cells;
	}


	/*!
	Default constructs local copy of the data of remote neighbor cells.

	Use this function to modify the default constructed incoming cells
	before they start receiving data from other processes.

	For example if the default constructed version of cell data class
	doesn't return a correct MPI_Datatype when updating remote neighbor
	data then call this beforehand and use get_remote_neighbors() to
	get a list of cells which you should change.
	*/
	void create_copies_of_remote_neighbors()
	{
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			sender = this->cells_to_receive.begin();
			sender != this->cells_to_receive.end();
			sender++
		) {
			for (std::vector<std::pair<uint64_t, int> >::const_iterator
				item = sender->second.begin();
				item != sender->second.end();
				item++
			) {
				const uint64_t cell = item->first;
				this->remote_neighbors[cell];
			}
		}
	}


	/*!
	Starts user data transfers between processes based on cells_to_send and cells_to_receive.

	User data arriving to this process is saved in given destination.
	*/
	void start_user_data_transfers(
	#ifdef DCCRG_SEND_SINGLE_CELLS
	boost::unordered_map<uint64_t, UserData>& destination,
	const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >& receive_item,
	const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >& send_item
	#elif defined (DCCRG_CELL_DATA_SIZE_FROM_USER)
	boost::unordered_map<uint64_t, UserData>& destination,
	boost::unordered_map<int, std::vector<uint64_t> >& receive_item,
	boost::unordered_map<int, std::vector<uint64_t> >& send_item
	#else
	boost::unordered_map<int, std::vector<uint64_t> >& receive_item,
	boost::unordered_map<int, std::vector<uint64_t> >& send_item
	#endif
	) {
		#ifdef DCCRG_SEND_SINGLE_CELLS

		// post all receives, messages are unique between different senders so just iterate over processes in random order
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			sender = receive_item.begin();
			sender != receive_item.end();
			sender++
		) {
			const int process = sender->first;

			#ifdef DEBUG
			if (process == (int) this->rank
			&& sender->second.size() > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " trying to transfer to self"
					<< std::endl;
				abort();
			}
			#endif

			for (std::vector<std::pair<uint64_t, int> >::const_iterator
				item = sender->second.begin();
				item != sender->second.end();
				item++
			) {
				const uint64_t cell = item->first;

				if (destination.count(cell) == 0) {
					destination[cell];
				}

				// TODO: preallocate MPI_Requests?
				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->receive_requests[process].push_back(MPI_Request());

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Datatype user_datatype = destination.at(cell).mpi_datatype();
				MPI_Type_commit(&user_datatype);
				#endif

				MPI_Irecv(
					destination.at(cell).at(),
					#ifdef DCCRG_USER_MPI_DATA_TYPE
					1,
					user_datatype,
					#else
					UserData::size(),
					MPI_BYTE,
					#endif
					process,
					item->second,
					this->comm,
					&(this->receive_requests[process].back())
				);

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Type_free(&user_datatype);
				#endif

				#else // ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

				this->receive_requests[process].push_back(
					this->boost_comm.irecv(
						process,
						item->second,
						destination.at(cell)
					)
				);
				#endif
			}
		}

		// post all sends
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			receiver = send_item.begin();
			receiver != send_item.end();
			receiver++
		) {
			const int process = receiver->first;

			#ifdef DEBUG
			if (process == (int) this->rank
			&& receiver->second.size() > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Trying to transfer to self" << std::endl;
				abort();
			}
			#endif

			for (std::vector<std::pair<uint64_t, int> >::const_iterator
				item = receiver->second.begin();
				item != receiver->second.end();
				item++
			) {
				const uint64_t cell = item->first;

				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->send_requests[process].push_back(MPI_Request());

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Datatype user_datatype = this->cells.at(cell).mpi_datatype();
				MPI_Type_commit(&user_datatype);
				#endif

				// FIXME: check the return value
				MPI_Isend(
					this->cells.at(cell).at(),
					#ifdef DCCRG_USER_MPI_DATA_TYPE
					1,
					user_datatype,
					#else
					UserData::size(),
					MPI_BYTE,
					#endif
					process,
					item->second,
					this->comm,
					&(this->send_requests[process].back())
				);

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Type_free(&user_datatype);
				#endif

				#else

				this->send_requests[process].push_back(
					this->boost_comm.isend(
						process,
						item->second,
						this->cells.at(cell)
					)
				);
				#endif
			}
		}

		// all user data is sent using one MPI message / process
		#else	// ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		// receive one MPI datatype per process
		for (boost::unordered_map<int, std::vector<uint64_t> >::iterator
			sender = receive_item.begin();
			sender != receive_item.end();
			sender++
		) {
			// TODO: still needed?
			std::sort(sender->second.begin(), sender->second.end());

			// reserve space for incoming user data at our end
			for (uint64_t i = 0; i < sender->second.size(); i++) {
				if (destination.count(sender->second[i]) == 0) {
					destination[sender->second[i]];
				}
			}

			// get displacements in bytes for incoming user data
			std::vector<MPI_Aint> displacements(sender->second.size(), 0);
			for (uint64_t i = 0; i < sender->second.size(); i++) {
				displacements[i]
					= (uint8_t*) destination.at(sender->second[i]).at()
					- (uint8_t*) destination.at(sender->second[0]).at();
			}

			MPI_Datatype receive_datatype;

			#ifdef DCCRG_USER_MPI_DATA_TYPE
			std::vector<int> block_lengths(displacements.size(), 1);
			std::vector<MPI_Datatype> datatypes(displacements.size());
			for (uint64_t i = 0; i < sender->second.size(); i++) {
				datatypes[i] = destination.at(sender->second[i]).mpi_datatype();
			}
			#else
			std::vector<int> block_lengths(displacements.size(), UserData::size());
			#endif

			#ifdef DCCRG_USER_MPI_DATA_TYPE
			MPI_Type_create_struct(
				displacements.size(),
				&block_lengths[0],
				&displacements[0],
				&datatypes[0],
				&receive_datatype
			);
			#else
			MPI_Type_create_hindexed(
				displacements.size(),
				&block_lengths[0],
				&displacements[0],
				MPI_BYTE,
				&receive_datatype
			);
			#endif

			MPI_Type_commit(&receive_datatype);

			int receive_tag = sender->first * this->comm_size + this->rank;

			this->receive_requests[sender->first].push_back(MPI_Request());

			MPI_Irecv(
				destination.at(sender->second[0]).at(),
				1,
				receive_datatype,
				sender->first,
				receive_tag,
				this->comm,
				&(this->receive_requests[sender->first].back())
			);

			MPI_Type_free(&receive_datatype);
			#ifdef DCCRG_USER_MPI_DATA_TYPE
			BOOST_FOREACH(MPI_Datatype& type, datatypes) {
				if (MPI_Type_free(&type) != MPI_SUCCESS) {
					std::cout << __FILE__ << ":" << __LINE__
						<< "Couldn't free MPI_Datatype"
						<< std::endl;
						abort();
				}
			}
			#endif
		}

		// send one MPI datatype per process
		for (boost::unordered_map<int, std::vector<uint64_t> >::iterator
			receiver = send_item.begin();
			receiver != send_item.end();
			receiver++
		) {
			std::sort(receiver->second.begin(), receiver->second.end());

			// get displacements in bytes for outgoing user data
			std::vector<MPI_Aint> displacements(receiver->second.size(), 0);
			for (uint64_t i = 0; i < receiver->second.size(); i++) {
				displacements[i]
					= (uint8_t*) this->cells.at(receiver->second[i]).at()
					- (uint8_t*) this->cells.at(receiver->second[0]).at();
			}

			MPI_Datatype send_datatype;

			#ifdef DCCRG_USER_MPI_DATA_TYPE
			std::vector<int> block_lengths(displacements.size(), 1);
			std::vector<MPI_Datatype> datatypes(displacements.size());
			for (uint64_t i = 0; i < receiver->second.size(); i++) {
				datatypes[i] = this->cells.at(receiver->second[i]).mpi_datatype();
			}
			#else
			std::vector<int> block_lengths(displacements.size(), UserData::size());
			#endif

			#ifdef DCCRG_USER_MPI_DATA_TYPE
			MPI_Type_create_struct(
				displacements.size(),
				&block_lengths[0],
				&displacements[0],
				&datatypes[0],
				&send_datatype
			);
			#else
			MPI_Type_create_hindexed(
				displacements.size(),
				&block_lengths[0],
				&displacements[0],
				MPI_BYTE,
				&send_datatype
			);
			#endif

			MPI_Type_commit(&send_datatype);

			int send_tag = this->rank * this->comm_size + receiver->first;

			this->send_requests[receiver->first].push_back(MPI_Request());

			MPI_Isend(
				this->cells.at(receiver->second[0]).at(),
				1,
				send_datatype,
				receiver->first,
				send_tag,
				this->comm,
				&(this->send_requests[receiver->first].back())
			);

			MPI_Type_free(&send_datatype);
			#ifdef DCCRG_USER_MPI_DATA_TYPE
			BOOST_FOREACH(MPI_Datatype& type, datatypes) {
				if (MPI_Type_free(&type) != MPI_SUCCESS) {
					std::cout << __FILE__ << ":" << __LINE__
						<< "Couldn't free MPI_Datatype"
						<< std::endl;
						abort();
				}
			}
			#endif
		}

		#else	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		// post all receives
		for (int sender = 0; sender < (int) this->comm_size; sender++) {

			if ((uint64_t) sender == this->rank) {
				continue;
			}

			if (receive_item.count(sender) == 0) {
				// no data to send / receive
				continue;
			}

			int receive_tag = sender * this->comm_size + this->rank;

			this->receive_requests[sender].push_back(
				this->boost_comm.irecv(
					sender,
					receive_tag,
					this->incoming_data[sender]
				)
			);
		}

		// gather all data to send
		for (int receiver = 0; receiver < (int) this->comm_size; receiver++) {

			if ((uint64_t) receiver == this->rank) {
				// don't send to self
				continue;
			}

			if (send_item.count(receiver) == 0) {
				// no data to send / receive
				continue;
			}

			std::sort(send_item.at(receiver).begin(), send_item.at(receiver).end());
			// construct the outgoing data vector
			BOOST_FOREACH(const uint64_t& cell, send_item.at(receiver)) {
				UserData* user_data = (*this)[cell];
				assert(user_data != NULL);
				this->outgoing_data[receiver].push_back(*user_data);
			}
		}

		// post all sends
		for (int receiver = 0; receiver < (int) this->comm_size; receiver++) {

			if ((uint64_t) receiver == this->rank) {
				continue;
			}

			if (send_item.count(receiver) == 0) {
				// no data to send / receive
				continue;
			}

			int send_tag = this->rank * this->comm_size + receiver;

			this->send_requests[receiver].push_back(
				this->boost_comm.isend(
					receiver,
					send_tag,
					this->outgoing_data[receiver]
				)
			);
		}
		#endif	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		#endif	// ifdef DCCRG_SEND_SINGLE_CELLS
	}


	/*!
	Waits for the receives of user data transfers between processes to complete.

	User data arriving to this process is saved in given destination.
	*/
	void wait_user_data_transfer_receives(
	#ifndef DCCRG_SEND_SINGLE_CELLS
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::unordered_map<uint64_t, UserData>& destination,
	boost::unordered_map<int, std::vector<uint64_t> >& receive_item
	#endif
	#endif
	) {
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator
			process = this->receive_requests.begin();
			process != this->receive_requests.end();
			process++
		) {
			std::vector<MPI_Status> statuses;
			statuses.resize(process->second.size());

			if (MPI_Waitall(process->second.size(), &(process->second[0]), &(statuses[0])) != MPI_SUCCESS) {
				BOOST_FOREACH(const MPI_Status& status, statuses) {
					if (status.MPI_ERROR != MPI_SUCCESS) {
						std::cerr << "MPI receive failed from process " << status.MPI_SOURCE
							<< " with tag " << status.MPI_TAG
							<< std::endl;
					}
				}
			}
		}

		#else	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator
			process = this->receive_requests.begin();
			process != this->receive_requests.end();
			process++
		) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}

		#ifndef DCCRG_SEND_SINGLE_CELLS

		// incorporate received data
		for (typename boost::unordered_map<int, std::vector<UserData> >::const_iterator
			sender = this->incoming_data.begin();
			sender != this->incoming_data.end();
			sender++
		) {
			std::sort(receive_item.at(sender->first).begin(), receive_item.at(sender->first).end());

			int i = 0;
			BOOST_FOREACH(const uint64_t& cell, receive_item.at(sender->first)) {
				// TODO move data instead of copying
				destination[cell] = this->incoming_data.at(sender->first)[i];
				i++;
			}
		}
		this->incoming_data.clear();

		#endif	// ifndef DCCRG_SEND_SINGLE_CELLS
		#endif	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		this->receive_requests.clear();
	}


	/*!
	Waits for the sends of user data transfers between processes to complete.
	*/
	void wait_user_data_transfer_sends()
	{
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator
			process = this->send_requests.begin();
			process != this->send_requests.end();
			process++
		) {
			std::vector<MPI_Status> statuses;
			statuses.resize(process->second.size());

			if (MPI_Waitall(process->second.size(), &(process->second[0]), &(statuses[0])) != MPI_SUCCESS) {
				BOOST_FOREACH(const MPI_Status& status, statuses) {
					if (status.MPI_ERROR != MPI_SUCCESS) {
						std::cerr << "MPI receive failed from process " << status.MPI_SOURCE
							<< " with tag " << status.MPI_TAG
							<< std::endl;
					}
				}
			}
		}

		#else	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator
			process = this->send_requests.begin();
			process != this->send_requests.end();
			process++
		) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}

		#ifndef DCCRG_SEND_SINGLE_CELLS

		this->outgoing_data.clear();

		#endif	// ifndef DCCRG_SEND_SINGLE_CELLS
		#endif	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		this->send_requests.clear();
	}


	/*!
	Returns true if cells with given index properties overlap.

	Sizes are also given in indices.
	*/
	bool indices_overlap(const uint64_t index1, const uint64_t size1, const uint64_t index2, const uint64_t size2) const
	{
		#ifdef DEBUG
		if (index1 >= this->get_x_length() * (uint64_t(1) << this->max_refinement_level)
		&& index1 >= this->get_y_length() * (uint64_t(1) << this->max_refinement_level)
		&& index1 >= this->get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid index given" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (index2 >= this->get_x_length() * (uint64_t(1) << this->max_refinement_level)
		&& index2 >= this->get_y_length() * (uint64_t(1) << this->max_refinement_level)
		&& index2 >= this->get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid index given" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (size1 > (uint64_t(1) << this->max_refinement_level)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid size given" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (size2 > (uint64_t(1) << this->max_refinement_level)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid size given" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		if (index1 + size1 > index2 && index1 < index2 + size2) {
			return true;
		} else {
			return false;
		}
	}

	/*!
	Same as the uint64_t version but in 3d, returns true only if all indices overlap.
	*/
	bool indices_overlap(const Types<3>::indices_t indices1, const uint64_t size1, const Types<3>::indices_t indices2, const uint64_t size2) const
	{
		for (int i = 0; i < 3; i++) {
			if (indices1[i] + size1 <= indices2[i] || indices1[i] >= indices2[i] + size2) {
				return false;
			}
		}
		return true;
	}

	/*!
	Returns true if x indices of given cells overlap, even if they don't exist
	*/
	bool x_indices_overlap(const uint64_t cell1, const uint64_t cell2) const
	{
		assert(cell1 > 0);
		assert(cell1 <= this->last_cell);
		assert(cell2 > 0);
		assert(cell2 <= this->last_cell);

		const uint64_t index1 = this->get_x_index(cell1);
		const uint64_t index2 = this->get_x_index(cell2);
		const uint64_t size1 = this->get_cell_size_in_indices(cell1);
		const uint64_t size2 = this->get_cell_size_in_indices(cell2);

		return this->indices_overlap(index1, size1, index2, size2);
	}

	/*!
	Returns true if y indices of given cells overlap, even if they don't exist
	*/
	bool y_indices_overlap(const uint64_t cell1, const uint64_t cell2) const
	{
		assert(cell1 > 0);
		assert(cell1 <= this->last_cell);
		assert(cell2 > 0);
		assert(cell2 <= this->last_cell);

		const uint64_t index1 = this->get_y_index(cell1);
		const uint64_t index2 = this->get_y_index(cell2);
		const uint64_t size1 = this->get_cell_size_in_indices(cell1);
		const uint64_t size2 = this->get_cell_size_in_indices(cell2);

		return this->indices_overlap(index1, size1, index2, size2);
	}

	/*!
	Returns true if z indices of given cells overlap, even if they don't exist
	*/
	bool z_indices_overlap(const uint64_t cell1, const uint64_t cell2) const
	{
		assert(cell1 > 0);
		assert(cell1 <= this->last_cell);
		assert(cell2 > 0);
		assert(cell2 <= this->last_cell);

		const uint64_t index1 = this->get_z_index(cell1);
		const uint64_t index2 = this->get_z_index(cell2);
		const uint64_t size1 = this->get_cell_size_in_indices(cell1);
		const uint64_t size2 = this->get_cell_size_in_indices(cell2);

		return this->indices_overlap(index1, size1, index2, size2);
	}


	/*!
	Returns the number of directions in which given cells' indices overlap
	Returns 0 if even one of given cells doesn't exist
	*/
	int overlapping_indices(const uint64_t cell1, const uint64_t cell2) const
	{
		#ifdef DEBUG
		if (cell1 == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			abort();
		}
		if (cell2 == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			abort();
		}

		if (cell1 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			abort();
		}
		if (cell2 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			abort();
		}
		#endif

		if (this->cell_process.count(cell1) == 0 || this->cell_process.count(cell2) == 0) {
			return 0;
		}

		const Types<3>::indices_t indices1 = this->get_indices(cell1);
		const Types<3>::indices_t indices2 = this->get_indices(cell2);

		const uint64_t size1 = this->get_cell_size_in_indices(cell1);
		const uint64_t size2 = this->get_cell_size_in_indices(cell2);

		int ret = 0;
		if (this->indices_overlap(indices1[0], size1, indices2[0], size2)) {
			ret++;
		}
		if (this->indices_overlap(indices1[1], size1, indices2[1], size2)) {
			ret++;
		}
		if (this->indices_overlap(indices1[2], size1, indices2[2], size2)) {
			ret++;
		}

		return ret;
	}


	/*!
	Returns the smallest existing cell at given indices between given refinement levels inclusive.

	Returns error_cell if no cell between given refinement ranges exists or an index is outside of
	the grid or minimum_refinement_level > maximum_refinement_level.
	*/
	uint64_t get_existing_cell(
		const Types<3>::indices_t& indices,
		const int minimum_refinement_level,
		const int maximum_refinement_level
	) const
	{
		if (indices[0] >= this->x_length * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (indices[1] >= this->y_length * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (indices[2] >= this->z_length * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (minimum_refinement_level > maximum_refinement_level) {
			return error_cell;
		}

		int average_refinement_level = (maximum_refinement_level + minimum_refinement_level) / 2;
		const uint64_t average_cell = this->get_cell_from_indices(indices, average_refinement_level);

		// use binary search recursively (assumes that all cells refine to 8 children)
		if (this->cell_process.count(average_cell) == 0) {

			// search for larger cell
			if (average_refinement_level > minimum_refinement_level) {

				uint64_t larger_cell = this->get_existing_cell(indices, minimum_refinement_level, average_refinement_level - 1);

				if (this->cell_process.count(larger_cell) == 0) {
					return 0;
				} else {
					return larger_cell;
				}
			} else {
				return 0;
			}
		} else {
			// search for smaller cell
			if (average_refinement_level < maximum_refinement_level) {
				uint64_t smaller_cell = this->get_existing_cell(indices, average_refinement_level + 1, maximum_refinement_level);

				if (this->cell_process.count(smaller_cell) == 0) {
					return average_cell;
				} else {
					return smaller_cell;
				}
			} else {
				return average_cell;
			}
		}
	}


	/*!
	Returns the number of values needed to represent the coordinate of a cell
	*/
	static int get_grid_dimensionality(void* /*data*/, int* error)
	{
		*error = ZOLTAN_OK;
		return 3;
	}


	/*!
	Fills geom_vec with the coordinates of cells given in global_id
	*/
	static void fill_with_cell_coordinates(void *data, int /*global_id_size*/, int /*local_id_size*/, int number_of_cells, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int /*number_of_dimensions*/, double *geom_vec, int *error)
	{
		Dccrg<UserData, UserGeometry>* dccrg_instance = reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		for (int i = 0; i < number_of_cells; i++) {
			uint64_t cell = uint64_t(global_ids[i]);
			if (dccrg_instance->cells.count(cell) == 0) {
				*error = ZOLTAN_FATAL;
				std::cerr << "Process " << dccrg_instance->rank
					<< ": Zoltan wanted the coordinates of a non-existing cell " << cell
					<< std::endl;
				return;
			}

			geom_vec[3 * i + 0] = dccrg_instance->get_cell_x(cell);
			geom_vec[3 * i + 1] = dccrg_instance->get_cell_y(cell);
			geom_vec[3 * i + 2] = dccrg_instance->get_cell_z(cell);
		}
	}


	/*!
	Returns the number of cells on this process
	*/
	static int get_number_of_cells(void* data, int* error)
	{
		Dccrg<UserData, UserGeometry>* dccrg_instance = reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;
		return dccrg_instance->cells.size();
	}


	/*!
	Writes all cell ids on this process to the global_ids array
	*/
	static void fill_cell_list(void* data, int /*global_id_size*/, int /*local_id_size*/, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int number_of_weights_per_object, float* object_weights, int* error)
	{
		Dccrg<UserData, UserGeometry>* dccrg_instance = reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		int i = 0;
		BOOST_FOREACH(const cell_and_data_pair_t& item, dccrg_instance->cells) {

			#ifdef DEBUG
			if (item.first == 0) {
				std::cerr << "User data exist for an illegal cell" << std::endl;
				abort();
			}
			#endif

			global_ids[i] = item.first;

			if (number_of_weights_per_object > 0) {
				if (dccrg_instance->cell_weights.count(item.first) > 0) {
					object_weights[i] = dccrg_instance->cell_weights.at(item.first);
				} else {
					object_weights[i] = 1;
				}
			}

			i++;
		}
	}


	/*!
	Writes the number of neighbors into number_of_neighbors for all cells given in global_ids.
	*/
	static void fill_number_of_neighbors_for_cells(
		void* data,
		int /*global_id_size*/,
		int /*local_id_size*/,
		int number_of_cells,
		ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR /*local_ids*/,
		int* number_of_neighbors,
		int* error
	) {
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		for (int i = 0; i < number_of_cells; i++) {
			uint64_t cell = uint64_t(global_ids[i]);
			if (dccrg_instance->cells.count(cell) == 0) {
				*error = ZOLTAN_FATAL;
				std::cerr << "Process " << dccrg_instance->rank
					<< ": Zoltan wanted the number of neighbors of a non-existing cell " << cell
					<< std::endl;
				return;
			}

			number_of_neighbors[i] = 0;
			BOOST_FOREACH(const uint64_t& neighbor, dccrg_instance->neighbors.at(cell)) {
				if (neighbor != 0
				/* Zoltan 3.501 crashes in hierarchial
				if a cell is a neighbor to itself */
				&& neighbor != cell) {
					number_of_neighbors[i]++;
				}
			}
		}
	}


	/*!
	Writes neighbor lists of given cells into neighbors, etc.
	*/
	static void fill_neighbor_lists(
		void* data,
		int /*global_id_size*/,
		int /*local_id_size*/,
		int number_of_cells,
		ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR /*local_ids*/,
		int* number_of_neighbors,
		ZOLTAN_ID_PTR neighbors,
		int* processes_of_neighbors,
		int number_of_weights_per_edge,
		float* edge_weights, int* error
	) {
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		int current_neighbor_number = 0;
		for (int i = 0; i < number_of_cells; i++) {
			uint64_t cell = uint64_t(global_ids[i]);
			if (dccrg_instance->cells.count(cell) == 0) {
				*error = ZOLTAN_FATAL;
				std::cerr << "Process " << dccrg_instance->rank
					<< ": Zoltan wanted neighbor list of a non-existing cell " << cell
					<< std::endl;
				return;
			}

			number_of_neighbors[i] = 0;

			BOOST_FOREACH(const uint64_t& neighbor, dccrg_instance->neighbors.at(cell)) {

				if (neighbor == 0
				/* Zoltan 3.501 crashes in hierarchial
				if a cell is a neighbor to itself */
				|| neighbor == cell) {
					continue;
				}

				number_of_neighbors[i]++;

				neighbors[current_neighbor_number] = neighbor;
				processes_of_neighbors[current_neighbor_number]
					= dccrg_instance->cell_process.at(neighbor);

				// weight of edge from cell to *neighbor
				if (number_of_weights_per_edge > 0) {
					edge_weights[current_neighbor_number] = 1.0;
				}

				current_neighbor_number++;
			}
		}
	}


	/*!
	Writes the number of hyperedges (self + one per neighbor cell) in the grid for all cells on this process.
	*/
	static void fill_number_of_hyperedges(
		void* data,
		int* number_of_hyperedges,
		int* number_of_connections,
		int* format, int* error
	) {
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		*number_of_hyperedges = dccrg_instance->cells.size();
		*format = ZOLTAN_COMPRESSED_EDGE;

		*number_of_connections = 0;
		BOOST_FOREACH(const cell_and_data_pair_t& item, dccrg_instance->cells) {

			(*number_of_connections)++;

			BOOST_FOREACH(const uint64_t& neighbor, dccrg_instance->neighbors.at(item.first)) {
				if (neighbor != 0
				/* Zoltan 3.501 crashes in hierarchial
				if a cell is a neighbor to itself */
				&& neighbor != item.first) {
					(*number_of_connections)++;
				}
			}
		}
	}


	/*!
	Writes the hypergraph in compressed edge format
	*/
	static void fill_hyperedge_lists(
		void* data,
		int /*global_id_size*/,
		int number_of_hyperedges,
		int number_of_connections,
		int format,
		ZOLTAN_ID_PTR hyperedges,
		int* hyperedge_connection_offsets,
		ZOLTAN_ID_PTR connections,
		int* error
	) {
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		if (format != ZOLTAN_COMPRESSED_EDGE) {
			std::cerr << "Only compressed edge format supported for hypergraph partitioning"
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		if ((unsigned int) number_of_hyperedges != dccrg_instance->cells.size()) {
			std::cerr << "Zoltan is expecting wrong number of hyperedges: " << number_of_hyperedges
				<< " instead of " << dccrg_instance->cells.size()
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		int i = 0;
		int connection_number = 0;
		BOOST_FOREACH(const cell_and_data_pair_t& item, dccrg_instance->cells) {

			hyperedges[i] = item.first;
			hyperedge_connection_offsets[i] = connection_number;

			// add a connection to the cell itself from its hyperedge
			connections[connection_number++] = item.first;

			BOOST_FOREACH(const uint64_t& neighbor, dccrg_instance->neighbors.at(item.first)) {
				if (neighbor == 0
				/* Zoltan 3.501 crashes in hierarchial
				if a cell is a neighbor to itself */
				|| neighbor == item.first) {
					continue;
				}

				connections[connection_number++] = neighbor;
			}

			i++;
		}

		if (connection_number != number_of_connections) {
			std::cerr << "Zoltan is expecting wrong number of connections from hyperedges: "
				<< number_of_connections
				<< " instead of " << connection_number
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}
	}


	/*!
	Writes the number of hyperedge weights (one per hyperedge) on this process
	*/
	static void fill_number_of_edge_weights(void* data, int* number_of_edge_weights, int* error)
	{
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		*number_of_edge_weights = dccrg_instance->cells.size();
		return;
	}


	/*!
	Writes hyperedge weights (one per hyperedge) on this process
	*/
	static void fill_edge_weights(
		void* data,
		int /*global_id_size*/,
		int /*local_id_size*/,
		int number_of_hyperedges,
		int number_of_weights_per_hyperedge,
		ZOLTAN_ID_PTR hyperedges,
		ZOLTAN_ID_PTR /*hyperedges_local_ids*/,
		float* hyperedge_weights,
		int* error
	) {
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;

		if ((unsigned int) number_of_hyperedges != dccrg_instance->cells.size()) {
			std::cerr << "Zoltan is expecting wrong number of hyperedges: " << number_of_hyperedges
				<< " instead of " << dccrg_instance->cells.size()
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		int i = 0;
		BOOST_FOREACH(const cell_and_data_pair_t& item, dccrg_instance->cells) {

			hyperedges[i] = item.first;

			if (number_of_weights_per_hyperedge > 0) {
				int number_of_hyperedges = 0;

				BOOST_FOREACH(const uint64_t& neighbor, dccrg_instance->neighbors.at(item.first)) {
					if (neighbor != 0
					/* Zoltan 3.501 crashes in hierarchial
					if a cell is a neighbor to itself (periodic grid) */
					&& neighbor != item.first) {
						number_of_hyperedges++;
					}
				}

				hyperedge_weights[i] = 1.0 * number_of_hyperedges;
			}

			i++;
		}
	}


	/*!
	Returns the number of hierarchies to use for load balancing.
	*/
	static int get_number_of_load_balancing_hierarchies(void* data, int* error)
	{
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);
		*error = ZOLTAN_OK;
		return dccrg_instance->processes_per_part.size();
	}


	/*!
	Returns the part number of this process on given hierarchy level of load balancing.
	*/
	static int get_part_number(void* data, int level, int* error)
	{
		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);

		if (level < 0 || level >= int(dccrg_instance->processes_per_part.size())) {
			std::cerr << "Zoltan wanted a part number for an invalid hierarchy level (should be [0, "
				<< dccrg_instance->processes_per_part.size() - 1
				<< "]): " << level
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return -1;
		} else {
			*error = ZOLTAN_OK;
		}

		int process = dccrg_instance->rank;
		int part;

		for (int i = 0; i <= level; i++) {
			part = process / dccrg_instance->processes_per_part[i];
			process %= dccrg_instance->processes_per_part[i];
		}

		return part;
	}


	/*!
	Sets the partitioning options of given zoltan instance for given level.
	*/
	static void set_partitioning_options(void* data, int level, struct Zoltan_Struct* zz, int* error)
	{
		if (zz == NULL) {
			std::cerr << "Zoltan gave a NULL pointer for zz" << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		Dccrg<UserData, UserGeometry>* dccrg_instance
			= reinterpret_cast<Dccrg<UserData, UserGeometry> *>(data);

		if (level < 0 || level >= int(dccrg_instance->processes_per_part.size())) {
			std::cerr
				<< "Zoltan wanted partitioning options for an invalid hierarchy level (level should be between 0 and "
				<< dccrg_instance->processes_per_part.size() - 1
				<< " inclusive): " << level
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return;
		} else {
			*error = ZOLTAN_OK;
		}

		for (boost::unordered_map<std::string, std::string>::const_iterator
			option = dccrg_instance->partitioning_options.at(level).begin();
			option != dccrg_instance->partitioning_options.at(level).end();
			option++
		) {
			Zoltan_Set_Param(zz, option->first.c_str(), option->second.c_str());
		}
	}



	#ifdef DEBUG
	/*!
	Returns false if the same cells don't exist on the same process for all processes.
	*/
	bool is_consistent()
	{
		// sort existing cells from this process
		std::vector<uint64_t> local_cells;
		local_cells.reserve(this->cell_process.size());
		for (typename boost::unordered_map<uint64_t, uint64_t>::const_iterator
			cell = this->cell_process.begin();
			cell != this->cell_process.end();
			cell++
		) {
			local_cells.push_back(cell->first);
		}
		std::sort(local_cells.begin(), local_cells.end());

		// processes of existing cells from this process
		std::vector<uint64_t> local_processes;
		local_processes.reserve(this->cell_process.size());
		for (std::vector<uint64_t>::const_iterator
			cell = local_cells.begin();
			cell != local_cells.end();
			cell++
		) {
			local_processes.push_back((uint64_t)this->cell_process[*cell]);
		}

		// compare the above between processes
		std::vector<std::vector<uint64_t> > all_cells;
		All_Gather()(local_cells, all_cells, this->comm);
		std::vector<std::vector<uint64_t> > all_processes;
		All_Gather()(local_processes, all_processes, this->comm);

		for (uint64_t process = 0; process < this->comm_size; process++) {
			if (!std::equal(
				all_cells[process].begin(),
				all_cells[process].end(),
				all_cells[0].begin()
			)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Grid has different cells between processes 0 and " << process
					<< std::endl;
				return false;
			}

			if (!std::equal(
				all_processes[process].begin(),
				all_processes[process].end(),
				all_processes[0].begin()
			)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Grid's cells have different processes between processes 0 and " << process
					<< std::endl;
				return false;
			}
		}

		return true;
	}


	/*!
	Return false if neighbors lists of the given cell aren't consistent
	*/
	bool verify_neighbors(
		const uint64_t cell,
		const std::vector<Types<3>::neighborhood_item_t>& hood_of,
		const std::vector<Types<3>::neighborhood_item_t>& hood_to,
		const boost::unordered_map<uint64_t, std::vector<uint64_t> >& neighbor_of_lists,
		boost::unordered_map<uint64_t, std::vector<uint64_t> >& neighbor_to_lists
	) {
		if (cell == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			return false;
		}

		if (cell > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ <<
				" Cell " << cell << " shouldn't exist"
				<< std::endl;
			return false;
		}

		if (this->cell_process.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Cell " << cell << " doesn't exist"
				<< std::endl;
			return false;
		}

		if (cell == this->get_child(cell)) {

			if (neighbor_of_lists.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No neighbor list for cell " << cell
					<< std::endl;
				return false;
			}

			if (neighbor_to_lists.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No neighbor_to list for cell " << cell
					<< std::endl;
				return false;
			}

		} else {

			if (neighbor_of_lists.count(cell) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor list for cell " << cell
					<< " shouldn't exist"
					<< std::endl;
				return false;
			}

			if (neighbor_to_lists.count(cell) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor_to list for cell " << cell
					<< " shouldn't exist"
					<< std::endl;
				return false;
			}

		}

		if (cell != this->get_child(cell)) {
			return true;
		}

		// neighbors
		std::vector<uint64_t> compare_neighbors
			= this->find_neighbors_of(cell, hood_of, this->max_ref_lvl_diff);

		if (
			neighbor_of_lists.at(cell).size() != compare_neighbors.size()
		|| (
			neighbor_of_lists.at(cell).size() > 0
			&& compare_neighbors.size() > 0
			&& !std::equal(
				neighbor_of_lists.at(cell).begin(),
				neighbor_of_lists.at(cell).end(),
				compare_neighbors.begin()
			)
		)
		) {
			std::cerr << "Process " << this->rank
				<< " neighbor counts for cell " << cell
				<< " (child of " << this->get_parent(cell)
				<< ") don't match "
				<< neighbor_of_lists.at(cell).size() << ": ";

			BOOST_FOREACH(const uint64_t& c, neighbor_of_lists.at(cell)) {
				std::cerr << c << " ";
			}
			std::cerr << ", should be (+ child of) " << compare_neighbors.size() << ": ";
			BOOST_FOREACH(const uint64_t& c, compare_neighbors) {
				std::cerr << c << "(" << this->get_parent(c) << ") ";
			}
			std::cerr << std::endl;
			return false;
		}

		// neighbors_to
		if (neighbor_to_lists.at(cell).size() > 0) {
			std::sort(neighbor_to_lists.at(cell).begin(), neighbor_to_lists.at(cell).end());
		}
		std::vector<uint64_t> compare_neighbors_to
			= this->find_neighbors_to(cell, hood_to);
		if (compare_neighbors_to.size() > 0) {
			std::sort(compare_neighbors_to.begin(), compare_neighbors_to.end());
		}

		if (
			!std::equal(
				neighbor_to_lists.at(cell).begin(),
				neighbor_to_lists.at(cell).end(),
				compare_neighbors_to.begin()
			)
		) {
			std::cerr << "Process " << this->rank
				<< " neighbor_to counts for cell " << cell
				<< " (child of " << this->get_parent(cell)
				<< ") don't match: " << neighbor_to_lists.at(cell).size()
				<< " (";

			BOOST_FOREACH(const uint64_t& c, neighbor_to_lists.at(cell)) {
				std::cerr << c;
				if (c != this->get_child(c)) {
					std::cerr << " [has a child " << this->get_child(c) << "], ";
				} else {
					std::cerr << ", ";
				}
			}
			std::cerr << ") should be " << compare_neighbors_to.size() << " (";
			BOOST_FOREACH(const uint64_t& c, compare_neighbors_to) {
				std::cerr << c << ", ";
			}
			std::cerr << ")" << std::endl;
			return false;
		}

		return true;
	}


	/*!
	Returns false if neighbor lists on this process aren't consistent
	*/
	bool verify_neighbors()
	{
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			cell = this->cell_process.begin();
			cell != this->cell_process.end();
			cell++
		) {
			if (cell->second != this->rank) {
				continue;
			}

			// verify default neighbor lists
			if (
				!this->verify_neighbors(
					cell->first,
					this->neighborhood_of,
					this->neighborhood_to,
					this->neighbors,
					this->neighbors_to
				)
			) {
				return false;
			}

			// verify user neighbor lists
			for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
				item = this->user_hood_of.begin();
				item != this->user_hood_of.end();
				item++
			) {
				const int hood_id = item->first;
				if (
					!this->verify_neighbors(
						cell->first,
						this->user_hood_of.at(hood_id),
						this->user_hood_to.at(hood_id),
						this->user_neigh_of.at(hood_id),
						this->user_neigh_to.at(hood_id)
					)
				) {
					return false;
				}
			}
		}

		return true;
	}


	/*!
	Returns false if remote neighbor info for given cell is inconsistent.

	Remote neighbor info consists of cells_with_remote_neighbors and remote_cells_with_local_neighbors.
	*/
	bool verify_remote_neighbor_info(const uint64_t cell)
	{
		if (
			!this->verify_neighbors(
				cell,
				this->neighborhood_of,
				this->neighborhood_to,
				this->neighbors,
				this->neighbors_to
			)
		) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Cell " << cell
				<< " has inconsistent neighbors"
				<< std::endl;
			return false;
		}

		if (cell != this->get_child(cell)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Cell " << cell << " has children"
				<< std::endl;
			return true;
		}

		std::vector<uint64_t> all_neighbors(
			this->neighbors.at(cell).begin(),
			this->neighbors.at(cell).end()
		);
		all_neighbors.insert(
			all_neighbors.end(),
			this->neighbors_to.at(cell).begin(),
			this->neighbors_to.at(cell).end()
		);

		BOOST_FOREACH(const uint64_t& neighbor, all_neighbors) {

			if (neighbor == 0) {
				continue;
			}

			if (this->cell_process.at(neighbor) != this->rank) {

				if (this->cells_with_remote_neighbors.count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Local cell " << cell
						<< " should be in cells_with_remote_neighbors"
						<< std::endl;
					return false;
				}

				if (this->remote_cells_with_local_neighbors.count(neighbor) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Remote cell " << neighbor
						<< " should be in remote_cells_with_local_neighbors"
						<< std::endl;
					return false;
				}
			}
		}

		return true;
	}


	/*!
	Returns false if remote neighbor info on this process is inconsistent.

	Remote neighbor info consists of cells_with_remote_neighbors
	and remote_cells_with_local_neighbors.
	*/
	bool verify_remote_neighbor_info()
	{
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			item = this->cell_process.begin();
			item != this->cell_process.end();
			item++
		) {

			if (item->first != this->get_child(item->first)) {
				continue;
			}

			// check whether this cell should be in remote_cells_with_local_neighbors
			if (item->second != this->rank) {

				bool should_be_in_remote_cells = false;

				BOOST_FOREACH(const cell_and_data_pair_t& cell, this->cells) {

					if (cell.first != this->get_child(cell.first)) {
						continue;
					}

					if (item->first == cell.first) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Same cell." << std::endl;
						abort();
					}

					if (this->is_neighbor(item->first, cell.first)
					|| this->is_neighbor(cell.first, item->first)) {
						should_be_in_remote_cells = true;
					}
				}

				if (should_be_in_remote_cells) {

					if (this->remote_cells_with_local_neighbors.count(item->first) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Remote cell " << item->first
							<< " should be in remote_cells_with_local_neighbors because:"
							<< std::endl;

						BOOST_FOREACH(const cell_and_data_pair_t& cell, this->cells) {
							if (item->first == cell.first) {
								std::cerr << __FILE__ << ":" << __LINE__ << " Same cell." << std::endl;
								abort();
							}

							if (this->is_neighbor(item->first, cell.first)
							|| this->is_neighbor(cell.first, item->first)) {
								std::cerr << "\tremote cell " << item->first
									<< " has a local neighbor " << cell.first
									<< std::endl;
							}
						}
						return false;
					}

				} else {

					if (this->remote_cells_with_local_neighbors.count(item->first) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Remote cell " << item->first
							<< " should not be in remote_cells_with_local_neighbors"
							<< std::endl;
						return false;
					}
				}

			// check whether this cell should be in cells_with_remote_neighbor
			} else {

				bool no_remote_neighbor = true;

				// search in neighbors_of
				const std::vector<uint64_t> neighbors_of
					= this->find_neighbors_of(item->first, this->neighborhood_of, this->max_ref_lvl_diff);

				BOOST_FOREACH(const uint64_t& neighbor, neighbors_of) {

					if (neighbor == 0) {
						continue;
					}

					if (this->cell_process.at(neighbor) != this->rank) {
						no_remote_neighbor = false;
					}

					if (!this->is_neighbor(item->first, neighbor)) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Cell " << neighbor
							<< " should be a neighbor of cell " << item->first
							<< std::endl;
						abort();
					}
				}

				// search in neighbors_to
				std::vector<uint64_t> neighbors_to
					= this->find_neighbors_to(item->first, this->neighborhood_to);
				BOOST_FOREACH(const uint64_t& neighbor, neighbors_to) {

					if (neighbor == 0) {
						continue;
					}

					if (this->cell_process.at(neighbor) != this->rank) {
						no_remote_neighbor = false;
					}

					if (!this->is_neighbor(neighbor, item->first)) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Cell " << item->first
							<< " should be a neighbor of cell " << neighbor
							<< std::endl;
						exit(EXIT_FAILURE);
					}
				}

				if (no_remote_neighbor) {
					if (this->cells_with_remote_neighbors.count(item->first) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Local cell " << item->first
							<< " should not be in cells_with_remote_neighbors"
							<< std::endl;
						return false;
					}
				} else {
					if (this->cells_with_remote_neighbors.count(item->first) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Local cell " << item->first
							<< " should be in cells_with_remote_neighbors"
							<< std::endl;
						return false;
					}
				}
			}
		}

		return true;
	}


	/*!
	Returns true if user data exists for local cells.
	*/
	bool verify_user_data()
	{
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			item = this->cell_process.begin();
			item != this->cell_process.end();
			item++
		) {
			if (item->second == this->rank
			&& item->first == this->get_child(item->first)
			&& this->cells.count(item->first) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " User data for local cell " << item->first
					<< " should exist"
					<< std::endl;
				return false;
			}
			if (item->second != this->rank
			&& this->cells.count(item->first) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " User data for local cell " << item->first
					<< " shouldn't exist"
					<< std::endl;
				return false;
			}
		}

		return true;
	}


	/*!
	Returns true if all cells are where pin reqests should have placed them.
	*/
	bool pin_requests_succeeded()
	{
		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			pin_request = this->pin_requests.begin();
			pin_request != this->pin_requests.end();
			pin_request++
		) {
			if (this->cell_process.at(pin_request->first) != pin_request->second) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << pin_request->first
					<< " not at requested process " << pin_request->second
					<< " but at " << this->cell_process.at(pin_request->first)
					<< std::endl;
				return false;
			}
		}

		return true;
	}
	#endif

};

}	// namespace

#endif

