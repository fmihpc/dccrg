/*
A distributed cartesian cell-refinable grid

Copyright 2009, 2010, 2011 Finnish Meteorological Institute

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
By default dccrg creates unrefined cells of constant size in x, y and z directions.
Cells of arbitrary size in x, y and z directions can be created by defining DCCRG_ARBITRARY_STRETCH.
DCCRG_CONSTANT_STRETCH is not supported at the moment.
*/
#ifdef DCCRG_ARBITRARY_STRETCH
	#ifdef DCCRG_CONSTANT_STRETCH
		#error Only one type of grid stretching can be used at a time
	#endif
#endif


/*
If the size of the data in every cell is known in advance by the user, neighbour data updates can be optimized by defining DCCRG_CELL_DATA_SIZE_FROM_USER, in which case:
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
#include "boost/mpi.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cassert"
#include "cstdio"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "functional"
#include "stdint.h"
#include "utility"
#include "vector"
#include "zoltan.h"

// make CellGeometry serializable
#include "boost/serialization/serialization.hpp"
#include "dccrg_cell_geometry.hpp"


template <class UserData> class dccrg
{

public:

	/*!
	Creates an uninitialized instance of the grid.

	The instance's initialize function must be called before doing anything else, otherwise the results will be undefined.
	*/
	dccrg()
	{
		this->initialized = false;
	}


	/*!
	Initializes the grid with given parameters, see the initialize function for their description.

	Zoltan_Initialize must have been called before calling this constructor.
	The instance's initialize function must not be called after using this constuctor to create an instance of the grid.
	 */
	dccrg(
		boost::mpi::communicator comm,
		const char* load_balancing_method,
		#ifdef DCCRG_ARBITRARY_STRETCH
		const std::vector<double> x_coordinates,
		const std::vector<double> y_coordinates,
		const std::vector<double> z_coordinates,
		#else
		const double x_start,
		const double y_start,
		const double z_start,
		const double cell_x_size,
		const double cell_y_size,
		const double cell_z_size,
		const uint64_t x_length,
		const uint64_t y_length,
		const uint64_t z_length,
		#endif
		const unsigned int neighbourhood_size,
		const int maximum_refinement_level = -1
	)
	{
		this->initialized = false;

		this->initialize(
			comm,
			load_balancing_method,
			#ifdef DCCRG_ARBITRARY_STRETCH
			x_coordinates,
			y_coordinates,
			z_coordinates,
			#else
			x_start,
			y_start,
			z_start,
			cell_x_size,
			cell_y_size,
			cell_z_size,
			x_length,
			y_length,
			z_length,
			#endif
			neighbourhood_size,
			maximum_refinement_level
		);
	}


	/*!
	Initializes this instance of the grid with given parameters.

	Zoltan_Initialize must have been called before calling this function.

	comm: the grid will span all the processes in the communicator comm

	load_balancing_method:
		The method that Zoltan will use for load balancing, given as a string.
		All methods except REFTREE are supported, see this page for a list of available methods: http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg.html#LB_METHOD

	neighbourhood_size:
		Determines which cells are considered neighbours.
		When calculating the neighbours of a given cell a cube of length neighbourhood_size + 1 in every direction is considered, centered at the cell for which neighbours are being calculated.
		The unit lenght of the cube is the cell for which neighbours are being calculated.
		If neighbourhood_size == 0, only cells (or children within the volume of cells of the same size as the current cell) that share a face are considered.

	maximum_refinement_level:
		The maximum number of times an unrefined cell can be refined (replacing it with 8 smaller cells).
		If not given the maximum refinement level is maximized based on the grids initial size.

	Depending on the type of geometry selected when compiling programs using the grid, one of the following sets of parameters is needed:

	#ifdef DCCRG_ARBITRARY_STRETCH

	x, y and z_coordinates:
		The coordinates of unrefined cells in the respective directions.
		First coordinate is the starting point of the grid, the following ith value is the endpoint of the ith unrefined cell.

	#else

	x_start, y_start, z_start:
		The starting corner of the grid.
	cell_size:
		The size of each unrefined cell in every direction.
	x_length, y_length, z_length:
		The number of cells in the grid in x, y and z direction.

	#endif

	*/
	void initialize(
		boost::mpi::communicator comm,
		const char* load_balancing_method,
		#ifdef DCCRG_ARBITRARY_STRETCH
		const std::vector<double> x_coordinates,
		const std::vector<double> y_coordinates,
		const std::vector<double> z_coordinates,
		#else
		const double x_start,
		const double y_start,
		const double z_start,
		const double cell_x_size,
		const double cell_y_size,
		const double cell_z_size,
		const uint64_t x_length,
		const uint64_t y_length,
		const uint64_t z_length,
		#endif
		const unsigned int neighbourhood_size,
		const int maximum_refinement_level = -1
	)
	{
		if (this->initialized) {
			std::cerr << "Initialize function called for an already initialized dccrg" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}

		this->comm = comm;

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
		this->reserved_options.insert("OBJ_WEIGHT_DIM");
		this->reserved_options.insert("RETURN_LISTS");
		this->reserved_options.insert("NUM_GLOBAL_PARTS");
		this->reserved_options.insert("NUM_LOCAL_PARTS");
		this->reserved_options.insert("AUTO_MIGRATE");

		// set reserved options
		Zoltan_Set_Param(this->zoltan, "EDGE_WEIGHT_DIM", "0");	// 0 because Zoltan crashes in hierarchial with larger values
		Zoltan_Set_Param(this->zoltan, "NUM_GID_ENTRIES", "1");
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
		Zoltan_Set_Num_Obj_Fn(this->zoltan, &dccrg<UserData>::get_number_of_cells, this);
		Zoltan_Set_Obj_List_Fn(this->zoltan, &dccrg<UserData>::fill_cell_list, this);
		Zoltan_Set_Num_Geom_Fn(this->zoltan, &dccrg<UserData>::get_grid_dimensionality, NULL);
		Zoltan_Set_Geom_Multi_Fn(this->zoltan, &dccrg<UserData>::fill_with_cell_coordinates, this);
		Zoltan_Set_Num_Edges_Multi_Fn(this->zoltan, &dccrg<UserData>::fill_number_of_neighbours_for_cells, this);
		Zoltan_Set_Edge_List_Multi_Fn(this->zoltan, &dccrg<UserData>::fill_neighbour_lists, this);
		Zoltan_Set_HG_Size_CS_Fn(this->zoltan, &dccrg<UserData>::fill_number_of_hyperedges, this);
		Zoltan_Set_HG_CS_Fn(this->zoltan, &dccrg<UserData>::fill_hyperedge_lists, this);
		Zoltan_Set_HG_Size_Edge_Wts_Fn(this->zoltan, &dccrg<UserData>::fill_number_of_edge_weights, this);
		Zoltan_Set_HG_Edge_Wts_Fn(this->zoltan, &dccrg<UserData>::fill_edge_weights, this);
		Zoltan_Set_Hier_Num_Levels_Fn(this->zoltan, &dccrg<UserData>::get_number_of_load_balancing_hierarchies, this);
		Zoltan_Set_Hier_Part_Fn(this->zoltan, &dccrg<UserData>::get_part_number, this);
		Zoltan_Set_Hier_Method_Fn(this->zoltan, &dccrg<UserData>::set_partitioning_options, this);


		/*
		Set grid parameters
		*/

		#ifdef DCCRG_ARBITRARY_STRETCH
		if (!this->geometry.set_coordinates(x_coordinates, y_coordinates, z_coordinates)) {
			std::cerr << "Failed to set grid geometry" << std::endl;
			exit(EXIT_FAILURE);
		}
		#else
		this->geometry.set_x_start(x_start);
		this->geometry.set_y_start(y_start);
		this->geometry.set_z_start(z_start);
		this->geometry.set_cell_x_size(cell_x_size);
		this->geometry.set_cell_y_size(cell_y_size);
		this->geometry.set_cell_z_size(cell_z_size);
		this->geometry.set_x_length(x_length);
		this->geometry.set_y_length(y_length);
		this->geometry.set_z_length(z_length);

		if (x_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the x direction" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}
		this->geometry.set_x_length(x_length);

		if (y_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the y direction" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}
		this->geometry.set_y_length(y_length);

		if (z_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the z direction" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}
		this->geometry.set_z_length(z_length);

		#endif

		this->neighbourhood_size = neighbourhood_size;

		// get the maximum refinement level based on the size of the grid when using uint64_t for cell ids
		double max_id = uint64_t(~0), last_id = this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length();
		int refinement_level = 0;
		while (last_id / max_id < 1) {
			refinement_level++;
			last_id += double(this->geometry.get_x_length()) * this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << refinement_level * 3);
		}
		refinement_level--;

		// grid is too large even without refinement
		if (refinement_level < 0) {
			std::cerr << "Given grid would contain more than 2^64 - 1 unrefined cells" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);
		}


		if (maximum_refinement_level > refinement_level) {

			std::cerr << "Given max_refinement_level (" << maximum_refinement_level << ") is too large: " << "x_length * this->geometry.get_y_length() * this->geometry.get_z_length() * 8^max_refinement_level / (2^64 - 1) >= " << this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << maximum_refinement_level * 3) / max_id << " but must be < 1" << std::endl;
			// TODO: throw an exception instead
			exit(EXIT_FAILURE);

		} else if (maximum_refinement_level < 0) {
			this->max_refinement_level = refinement_level;
		} else {
			this->max_refinement_level = maximum_refinement_level;
		}

		this->geometry.set_maximum_refinement_level(this->max_refinement_level);

		// the number of the last cell at maximum refinement level
		uint64_t id = 0;
		for (refinement_level = 0; refinement_level <= this->max_refinement_level; refinement_level++) {
			id += this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << refinement_level * 3);
		}
		this->max_cell_number = id;

		// create unrefined cells
		uint64_t cells_per_process = 1 + this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() / uint64_t(comm.size());
		for (uint64_t id = 1; id <= this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length(); id++) {
			if (id / cells_per_process == uint64_t(comm.rank())) {
				this->cells[id];
			}
			this->cell_process[id] = id / cells_per_process;
		}

		// update neighbour lists of created cells
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
			this->neighbours[cell->first] = this->find_neighbours_of(cell->first);
			this->neighbours_to[cell->first] = this->find_neighbours_to(cell->first);
		}
		#ifdef DEBUG
		if (!this->verify_neighbours()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour lists are inconsistent" << std::endl;
			// TODO: throw an exception instead when debugging?
			exit(EXIT_FAILURE);
		}
		#endif

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
			this->update_remote_neighbour_info(cell->first);
		}
		#ifdef DEBUG
		if (!this->verify_remote_neighbour_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->recalculate_neighbour_update_send_receive_lists();

		this->initialized = true;
	}


	/*!
	Returns all cells on this process that don't have children (e.g. leaf cells)
	*/
	std::vector<uint64_t> get_cells(void) const
	{
		std::vector<uint64_t> all_cells;
		all_cells.reserve(this->cells.size());

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

			#ifdef DEBUG
			if (this->cell_process.count(cell->first) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell->first << " shouldn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (this->cell_process.at(cell->first) != this->comm.rank()) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Process " << this->comm.rank() << ": Cell " << cell->first << " should be on process " << this->cell_process.at(cell->first) << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			const uint64_t child = this->get_child(cell->first);
			assert(child > 0);

			if (child == cell->first) {
				all_cells.push_back(cell->first);
			}
		}

		return all_cells;
	}


	/*!
	Returns all cells on this process that don't have children (e.g. leaf cells) and don't have neighbours on other processes
	*/
	std::vector<uint64_t> get_cells_with_local_neighbours(void) const
	{
		std::vector<uint64_t> return_cells;
		return_cells.reserve(this->cells.size());

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

			uint64_t child = this->get_child(cell->first);
			assert(child > 0);

			if (child != cell->first) {
				continue;
			}

			bool has_remote_neighbour = false;

			assert(this->neighbours.count(cell->first) > 0);
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell->first).begin(); neighbour != this->neighbours.at(cell->first).end(); neighbour++) {
				if (this->cell_process.at(*neighbour) != this->comm.rank()) {
					has_remote_neighbour = true;
					break;
				}
			}

			if (!has_remote_neighbour) {
				return_cells.push_back(cell->first);
			}
		}

		return return_cells;
	}


	/*!
	Returns all cells on this process that don't have children (e.g. leaf cells) and have at least one neighbour on another processes
	*/
	std::vector<uint64_t> get_cells_with_remote_neighbour(void) const
	{
		std::vector<uint64_t> return_cells;
		return_cells.reserve(this->cells.size());

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

			uint64_t child = this->get_child(cell->first);
			assert(child > 0);

			if (child != cell->first) {
				continue;
			}

			bool has_remote_neighbour = false;

			assert(this->neighbours.count(cell->first) > 0);
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell->first).begin(); neighbour != this->neighbours.at(cell->first).end(); neighbour++) {
				if (this->cell_process.at(*neighbour) != this->comm.rank()) {
					has_remote_neighbour = true;
					break;
				}
			}

			if (has_remote_neighbour) {
				return_cells.push_back(cell->first);
			}
		}

		return return_cells;
	}


	/*!
	Returns all cells in the grid that don't have children (e.g. leaf cells)
	*/
	std::vector<uint64_t> get_all_cells(void) const
	{
		std::vector<uint64_t> all_cells;
		all_cells.reserve(this->cell_process.size());

		for (auto item = this->cell_process.cbegin(); item != this->cell_process.cend(); item++) {

			const uint64_t child = this->get_child(item->first);

			if (child == item->first) {
				all_cells.push_back(item->first);
			}
		}

		return all_cells;
	}


	/*!
	Returns a pointer to the user supplied data of given cell
	Return NULL if the given cell isn't on this process and if the given cell isn't a neighbour of any cell on this process
	*/
	UserData* operator [] (const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			return (UserData*) &(this->cells.at(cell));
		} else if (this->remote_neighbours.count(cell) > 0) {
			return (UserData*) &(this->remote_neighbours.at(cell));
		} else if (this->refined_cell_data.count(cell) > 0) {
			return (UserData*) &(this->refined_cell_data.at(cell));
		} else if (this->unrefined_cell_data.count(cell) > 0) {
			return (UserData*) &(this->unrefined_cell_data.at(cell));
		} else {
			return NULL;
		}
	}


	/*!
	Refines the grid so that at least the given cells whose parents are on this process will exist in the grid.

	Must be called simultaneously on all processes.
	Does not store the user data of any refined cell.
	Returns true on this process if successful and false if given an invalid cell (0 or a cell with a too large refinement level).
	*/
	bool load(const std::vector<uint64_t>& cells)
	{
		// see for example http://www.informit.com/articles/article.aspx?p=376878&rll=1 for an explanation about template<template...
		this->comm.barrier();

		// calculate which cells must be refined...
		boost::unordered_set<uint64_t> cells_to_refine;

		// ...and check for invalid cells
		for (std::vector<uint64_t>::const_iterator cell = cells.begin(); cell != cells.end(); cell++) {
			if (*cell == 0) {
				return false;
			}

			if (this->get_refinement_level(*cell) < 0) {
				return false;
			}

			// refine all parents of current cell
			uint64_t parent = this->get_parent_for_removed(*cell);
			while (parent != this->get_parent_for_removed(parent)) {
				cells_to_refine.insert(parent);
			}
			cells_to_refine.insert(parent);
		}

		// keep refining until no more refines on any process, TODO

		return true;
	}


	/*!
	Load balances the grid's cells among processes.

	Must be called simultaneously on all processes.
	Cells which haven't been pinned are moved as suggested by Zoltan, pinned cells are moved as requested by the user.
	Does not update remote neighbour data between processes afterward.
	Discards refines / unrefines.
	*/
	void balance_load(void)
	{
		this->comm.barrier();

		this->update_pin_requests();

		int partition_changed, global_id_size, local_id_size, number_to_receive, number_to_send;
		ZOLTAN_ID_PTR global_ids_to_receive, local_ids_to_receive, global_ids_to_send, local_ids_to_send;
		int *sender_processes, *receiver_processes;

		if (Zoltan_LB_Balance(this->zoltan, &partition_changed, &global_id_size, &local_id_size, &number_to_receive, &global_ids_to_receive, &local_ids_to_receive, &sender_processes, &number_to_send, &global_ids_to_send, &local_ids_to_send, &receiver_processes) != ZOLTAN_OK) {
			if (!this->no_load_balancing) {
				std::cerr << "Zoltan_LB_Partition failed" << std::endl;
				Zoltan_Destroy(&this->zoltan);
				// TODO: throw an exception instead
				abort();
			}
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
			if (this->cell_process.at(global_ids_to_receive[i]) != sender_processes[i]) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cannot receive cell " << global_ids_to_receive[i] << " from process " << sender_processes[i] << std::endl;
				abort();
			}
		}
		#endif

		this->cells_to_receive.clear();
		this->cells_to_send.clear();

		/*
		Processes and the cells for which data has to be received by this process
		*/

		// cells added to / removed from this process by load balancing
		boost::unordered_set<uint64_t> added_cells, removed_cells;

		// migration from user
		for (auto pin_request = this->pin_requests.cbegin(); pin_request != this->pin_requests.cend(); pin_request++) {

			const int current_process_of_cell = this->cell_process.at(pin_request->first);

			if (pin_request->second == this->comm.rank()
			&& current_process_of_cell != this->comm.rank()) {
				this->cells_to_receive[current_process_of_cell].push_back(pin_request->first);
				added_cells.insert(pin_request->first);
			}
		}

		// migration from Zoltan
		for (int i = 0; i < number_to_receive; i++) {

			// don't send / receive from self
			if (sender_processes[i] == this->comm.rank()) {
				continue;
			}

			// skip user-migrated cells
			if (this->pin_requests.count(global_ids_to_receive[i]) > 0) {
				continue;
			}

			this->cells_to_receive[sender_processes[i]].push_back(global_ids_to_receive[i]);

			#ifdef DEBUG
			if (added_cells.count(global_ids_to_receive[i]) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << global_ids_to_receive[i] << " has already been received from process " << this->comm.rank() << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			added_cells.insert(global_ids_to_receive[i]);
		}

		/*
		Processes and the cells for which data has to be sent by this process
		*/

		// migration from user
		for (auto pin_request = this->pin_requests.cbegin(); pin_request != this->pin_requests.cend(); pin_request++) {

			const int current_process_of_cell = this->cell_process.at(pin_request->first);
			const int destination_process = pin_request->second;

			if (destination_process != this->comm.rank()
			&& current_process_of_cell == this->comm.rank()) {
				this->cells_to_send[destination_process].push_back(pin_request->first);
				removed_cells.insert(pin_request->first);
			}
		}

		// migration from Zoltan
		for (int i = 0; i < number_to_send; i++) {

			// don't send / receive from self
			if (receiver_processes[i] == this->comm.rank()) {
				continue;
			}

			// skip user-migrated cells
			if (this->pin_requests.count(global_ids_to_send[i]) > 0) {
				continue;
			}

			this->cells_to_send[receiver_processes[i]].push_back(global_ids_to_send[i]);

			#ifdef DEBUG
			if (removed_cells.count(global_ids_to_send[i]) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << global_ids_to_send[i] << " has already been sent from process " << this->comm.rank() << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			removed_cells.insert(global_ids_to_send[i]);
		}

		Zoltan_LB_Free_Data(&global_ids_to_receive, &local_ids_to_receive, &sender_processes, &global_ids_to_send, &local_ids_to_send, &receiver_processes);

		this->move_cells(added_cells, removed_cells);
	}

	/*!
	Moves pinned grid cells as requested by the user.

	Must be called simultaneously on all processes.
	Cells which haven't been pinned are not moved.
	Does not update remote neighbour data between processes afterward.
	Discards refines / unrefines.
	*/
	void migrate_cells(void)
	{
		this->comm.barrier();

		this->update_pin_requests();

		this->cells_to_receive.clear();
		this->cells_to_send.clear();

		// cells added to / removed from this process by load balancing
		boost::unordered_set<uint64_t> added_cells, removed_cells;

		// processes and the cells for which data has to be received by this process by user request
		for (auto pin_request = this->pin_requests.cbegin(); pin_request != this->pin_requests.cend(); pin_request++) {

			const int current_process_of_cell = this->cell_process.at(pin_request->first);

			if (pin_request->second == this->comm.rank()
			&& current_process_of_cell != this->comm.rank()) {
				this->cells_to_receive[current_process_of_cell].push_back(pin_request->first);
				added_cells.insert(pin_request->first);
			}
		}

		// processes and the cells for which data has to be sent by this process by user request
		for (auto pin_request = this->pin_requests.cbegin(); pin_request != this->pin_requests.cend(); pin_request++) {

			const int current_process_of_cell = this->cell_process.at(pin_request->first);
			const int destination_process = pin_request->second;

			if (destination_process != this->comm.rank()
			&& current_process_of_cell == this->comm.rank()) {
				this->cells_to_send[destination_process].push_back(pin_request->first);
				removed_cells.insert(pin_request->first);
			}
		}

		this->move_cells(added_cells, removed_cells);
	}


	/*!
	Updates the user data between processes of those cells that have at least one neighbour or are considered as a neighbour of a cell on another process
	Must be called simultaneously on all processes
	*/
	void update_remote_neighbour_data(void)
	{
		this->start_remote_neighbour_data_update();
		this->wait_neighbour_data_update();
	}


	/*!
	Starts the update of neighbour data between processes and returns before (probably) it has completed
	Must be called simultaneously on all processes
	*/
	void start_remote_neighbour_data_update(void)
	{
		this->comm.barrier();
		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->remote_neighbours
		#else
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->remote_neighbours
		#endif
		#endif
		);
	}


	/*!
	Waits until all neighbour data update transfers between processes have completed and incorporates that data.
	Must be called simultaneously on all processes.
	*/
	void wait_neighbour_data_update(void)
	{
		this->wait_neighbour_data_update_receives();
		this->wait_neighbour_data_update_sends();
	}


	/*!
	Waits until all sends associated with neighbour data update transfers between processes have completed.
	Must be called simultaneously on all processes and probably must be called after wait...update_receives(void).
	*/
	void wait_neighbour_data_update_sends(void)
	{
		this->wait_user_data_transfer_sends();
	}


	/*!
	Waits until all receives associated with neighbour data update transfers between processes have completed and incorporates that data.
	Must be called simultaneously on all processes and probably must be called before wait...update_sends(void).
	*/
	void wait_neighbour_data_update_receives(void)
	{
		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->remote_neighbours
		#endif
		#endif
		);
	}


	/*!
	Returns the number of cells whose data this process has to send during a neighbour data update.

	The total amount of cells to be sent is returned so if a cell's data will be sent to N processes it is counted N times.
	*/
	uint64_t get_number_of_update_send_cells(void) const
	{
		uint64_t result = 0;
		for (boost::unordered_map<int, std::vector<uint64_t> >::const_iterator receiver = cells_to_send.begin(); receiver != cells_to_send.end(); receiver++) {
			result += receiver->second.size();
		}
		return result;
	}

	/*!
	Returns the number of cells whose data this process has to receive during a neighbour data update.
	*/
	uint64_t get_number_of_update_receive_cells(void) const
	{
		uint64_t result = 0;
		for (boost::unordered_map<int, std::vector<uint64_t> >::const_iterator sender = cells_to_receive.begin(); sender != cells_to_receive.end(); sender++) {
			result += sender->second.size();
		}
		return result;
	}


	/*!
	Returns a pointer to the neighbours of given cell
	Some neighbours might be on another process, but have a copy of their data on this process
	The local copy of remote neighbours' data is updated, for example, by calling update_remote_neighbour_data()
	Returns NULL if given cell doesn't exist or is on another process
	*/
	const std::vector<uint64_t>* get_neighbours(const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			#ifdef DEBUG
			if (this->neighbours.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Process " << this->comm.rank() << ": Neighbour list for cell " << cell << " doesn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif
			return &(this->neighbours.at(cell));
		} else {
			return NULL;
		}
	}

	/*!
	Returns a pointer to the cells that consider given cell as a neighbour but that are not neighbours of given cell
	Some cell might be on another process or the structure might be empty
	Returns NULL if given cell doesn't exist or is on another process
	*/
	const std::vector<uint64_t>* get_neighbours2(const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			#ifdef DEBUG
			if (this->neighbours_to.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbours_to list for cell " << cell << " doesn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif
			return &(this->neighbours_to.at(cell));
		} else {
			return NULL;
		}
	}


	/*!
	Returns the neighbour(s) of given cell in the positive or negative x direction, e.g. when viewed from that direction returns all neighbours that overlap the given cell
	Returns neighbours in positive x direction if given direction > 0 and negative direction otherwise
	Returns nothing if given cell doesn't exist or exists on another process
	*/
	std::vector<uint64_t> get_neighbours_x(const uint64_t cell, const double direction) const
	{
		std::vector<uint64_t> return_neighbours;

		if (this->cells.count(cell) == 0) {
			return return_neighbours;
		}

		uint64_t y_index = this->get_y_index(cell), z_index = this->get_z_index(cell);
		uint64_t size_i = this->get_cell_size_in_indices(cell);
		double x = this->get_cell_x(cell);

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {

			if (direction > 0) {
				if (this->get_cell_x(*neighbour) < x) {
					continue;
				}
			} else {
				if (this->get_cell_x(*neighbour) > x) {
					continue;
				}
			}

			uint64_t neigh_y_i = this->get_y_index(*neighbour), neigh_z_i = this->get_z_index(*neighbour);
			uint64_t neigh_size_i = this->get_cell_size_in_indices(*neighbour);

			// return only neighbours whose indices overlap with given cell in y and z directions
			if (neigh_y_i + neigh_size_i > y_index
			 && neigh_y_i < y_index + size_i
			 && neigh_z_i + neigh_size_i > z_index
			 && neigh_z_i < z_index + size_i) {
				return_neighbours.push_back(*neighbour);
			}
		}

		return return_neighbours;
	}
	/*!
	Same as get_neighbours_x but in y direction
	*/
	std::vector<uint64_t> get_neighbours_y(const uint64_t cell, const double direction) const
	{
		std::vector<uint64_t> return_neighbours;

		if (this->cells.count(cell) == 0) {
			return return_neighbours;
		}

		uint64_t x_index = this->get_x_index(cell), z_index = this->get_z_index(cell);
		uint64_t size_i = this->get_cell_size_in_indices(cell);
		double y = this->get_cell_y(cell);

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {

			if (direction > 0) {
				if (this->get_cell_y(*neighbour) < y) {
					continue;
				}
			} else {
				if (this->get_cell_y(*neighbour) > y) {
					continue;
				}
			}

			uint64_t neigh_x_i = this->get_x_index(*neighbour), neigh_z_i = this->get_z_index(*neighbour);
			uint64_t neigh_size_i = this->get_cell_size_in_indices(*neighbour);

			// return only neighbours whose indices overlap with given cell in y and z directions
			if (neigh_x_i + neigh_size_i > x_index
			 && neigh_x_i < x_index + size_i
			 && neigh_z_i + neigh_size_i > z_index
			 && neigh_z_i < z_index + size_i) {
				return_neighbours.push_back(*neighbour);
			}
		}

		return return_neighbours;
	}
	/*!
	Same as get_neighbours_x but in z direction
	*/
	std::vector<uint64_t> get_neighbours_z(const uint64_t cell, const double direction) const
	{
		std::vector<uint64_t> return_neighbours;

		if (this->cells.count(cell) == 0) {
			return return_neighbours;
		}

		uint64_t x_index = this->get_x_index(cell), y_index = this->get_y_index(cell);
		uint64_t size_i = this->get_cell_size_in_indices(cell);
		double z = this->get_cell_z(cell);

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {

			if (direction > 0) {
				if (this->get_cell_z(*neighbour) < z) {
					continue;
				}
			} else {
				if (this->get_cell_z(*neighbour) > z) {
					continue;
				}
			}

			uint64_t neigh_x_i = this->get_x_index(*neighbour), neigh_y_i = this->get_y_index(*neighbour);
			uint64_t neigh_size_i = this->get_cell_size_in_indices(*neighbour);

			// return only neighbours whose indices overlap with given cell in y and z directions
			if (neigh_x_i + neigh_size_i > x_index
			 && neigh_x_i < x_index + size_i
			 && neigh_y_i + neigh_size_i > y_index
			 && neigh_y_i < y_index + size_i) {
				return_neighbours.push_back(*neighbour);
			}
		}

		return return_neighbours;
	}


	/*!
	Returns the given cells neighbours that are on another process
	Returns nothing if given cell doesn't exist or is on another process or doesn't have remote neighbours
	*/
	std::vector<uint64_t> get_remote_neighbours(const uint64_t cell) const
	{
		std::vector<uint64_t> result;

		if (this->cells.count(cell) == 0) {
			return result;
		}

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {
			if (this->cell_process[*neighbour] != this->comm.rank()) {
				result.push_back(*neighbour);
			}
		}

		return result;
	}


	bool is_local(const uint64_t cell) const
	{
		if (this->cell_process.at(cell) == this->comm.rank()) {
			return true;
		} else {
			return false;
		}
	}


	/*!
	Returns the maximum possible refinement level of any cell in the grid (0 means unrefined)
	*/
	int get_max_refinement_level(void) const
	{
		return this->max_refinement_level;
	}


	/*!
	The following return the x, y or z coordinate of the center of given cell regardless of whether it exists or has children
	*/
	double get_cell_x(const uint64_t cell) const
	{
		return this->geometry.get_cell_x(cell);
	}
	double get_cell_y(const uint64_t cell) const
	{
		return this->geometry.get_cell_y(cell);
	}
	double get_cell_z(const uint64_t cell) const
	{
		return this->geometry.get_cell_z(cell);
	}


	double get_cell_x_min(const uint64_t cell) const
	{
		return this->geometry.get_cell_x_min(cell);
	}
	double get_cell_y_min(const uint64_t cell) const
	{
		return this->geometry.get_cell_y_min(cell);
	}
	double get_cell_z_min(const uint64_t cell) const
	{
		return this->geometry.get_cell_z_min(cell);
	}

	double get_cell_x_max(const uint64_t cell) const
	{
		return this->geometry.get_cell_x_max(cell);
	}
	double get_cell_y_max(const uint64_t cell) const
	{
		return this->geometry.get_cell_y_max(cell);
	}
	double get_cell_z_max(const uint64_t cell) const
	{
		return this->geometry.get_cell_z_max(cell);
	}

	/*!
	The following return the length of given cell in x, y or z direction regardless of whether it exists or has children
	*/
	double get_cell_x_size(const uint64_t cell) const
	{
		return this->geometry.get_cell_x_size(cell);
	}
	double get_cell_y_size(const uint64_t cell) const
	{
		return this->geometry.get_cell_y_size(cell);
	}
	double get_cell_z_size(const uint64_t cell) const
	{
		return this->geometry.get_cell_z_size(cell);
	}


	/*!
	Returns the refinement level of given cell, even if it doesn't exist
	Returns -1 if cell == 0 or cell would exceed maximum refinement level
	*/
	int get_refinement_level(const uint64_t cell) const
	{
		return this->geometry.get_refinement_level(cell);
	}


	/*!
	Writes the cells on this process into a vtk file with given name in ASCII format
	The cells are written in ascending order
	Must be called simultaneously on all processes
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
			outfile << geometry.get_cell_x_min(leaf_cells[i]) << " " << geometry.get_cell_y_min(leaf_cells[i]) << " " << geometry.get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_max(leaf_cells[i]) << " " << geometry.get_cell_y_min(leaf_cells[i]) << " " << geometry.get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_min(leaf_cells[i]) << " " << geometry.get_cell_y_max(leaf_cells[i]) << " " << geometry.get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_max(leaf_cells[i]) << " " << geometry.get_cell_y_max(leaf_cells[i]) << " " << geometry.get_cell_z_min(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_min(leaf_cells[i]) << " " << geometry.get_cell_y_min(leaf_cells[i]) << " " << geometry.get_cell_z_max(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_max(leaf_cells[i]) << " " << geometry.get_cell_y_min(leaf_cells[i]) << " " << geometry.get_cell_z_max(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_min(leaf_cells[i]) << " " << geometry.get_cell_y_max(leaf_cells[i]) << " " << geometry.get_cell_z_max(leaf_cells[i]) << std::endl;
			outfile << geometry.get_cell_x_max(leaf_cells[i]) << " " << geometry.get_cell_y_max(leaf_cells[i]) << " " << geometry.get_cell_z_max(leaf_cells[i]) << std::endl;
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

	Takes priority over unrefining. Refines / unrefines take effect only after a call to stop_refining() and are lost after a call to balance_load(). Does nothing in any of the following cases:
		-given cell has already been refined (including induced refinement) and stop_refining() has not been called afterwards
		-given cell doesn't exist on this process
		-given cells children already exist
		-the created childrens' refinement level would exceed max_refinement_level
	Children are created on their parent's process.
	 */
	void refine_completely(const uint64_t cell)
	{
		if (this->cells.count(cell) == 0) {
			return;
		}

		if (this->get_refinement_level(cell) >= this->max_refinement_level) {
			return;
		}

		if (cell != this->get_child(cell)) {
			// cell already has children
			return;
		}

		this->cells_to_refine.insert(cell);
	}

	/*!
	As refine_completely, but uses the smallest existing cell at given coordinates.
	Does nothing in the same cases as refine_completely and additionally if the coordinate is outside of the grid.
	*/
	void refine_completely_at(const double x, const double y, const double z)
	{
		const uint64_t cell = this->get_smallest_cell_from_coordinate(x, y, z);
		if (cell == 0) {
			return;
		}

		this->refine_completely(cell);
	}


	/*!
	Removes the given cell and its siblings from the grid.

	Refining (including induced refining) takes priority over unrefining. Refines / unrefines take effect only after a call to stop_refining() and are lost after a call to balance_load(). Does nothing in any of the following cases:
		-given cell or one of its siblings has already been unrefined and stop_refining() has not been called
		-given cell doesn't exist on this process
		-given cell has children
		-given cells refinement level is 0

	After a cell and its siblings have been unrefined, their data has been moved to their parent's process. When no longer needed that data can be freed using ?.
	*/
	void unrefine_completely(const uint64_t cell)
	{
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

		// record only one sibling to unrefine / process
		const std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(cell));
		for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {
			if (this->cells_to_unrefine.count(*sibling) > 0) {
				return;
			}
		}

		this->cells_to_unrefine.insert(cell);
	}


	/*!
	As unrefine_completely, but uses the smallest existing cell at given coordinates.
	Does nothing in the same cases as unrefine_completely and additionally if the coordinate is outside of the grid.
	*/
	void unrefine_completely_at(const double x, const double y, const double z)
	{
		const uint64_t cell = this->get_smallest_cell_from_coordinate(x, y, z);
		if (cell == 0) {
			return;
		}

		this->unrefine_completely(cell);
	}


	/*!
	Executes refines / unrefines that have been requested so far.

	Must be called simultaneously on all processes.
	Returns cells that were created by refinement on this process.
	Moves user data of unrefined cells to the process of their parent.
	*/
	std::vector<uint64_t> stop_refining(void)
	{
		this->comm.barrier();
		this->induce_refines();
		this->override_unrefines();
		return this->execute_refines();
	}


	/*!
	Returns cells that were removed by unrefinement whose parent is on this process
	Removed cells data is also on this process, but only until balance_load() is called
	*/
	std::vector<uint64_t> get_removed_cells(void) const
	{
		std::vector<uint64_t> unref_removed_cells;
		unref_removed_cells.reserve(this->unrefined_cell_data.size());

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->unrefined_cell_data.begin(); cell != this->unrefined_cell_data.end(); cell++) {
			unref_removed_cells.push_back(cell->first);
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

		const uint64_t parent = get_cell_from_indices(this->get_x_index(cell), this->get_y_index(cell), this->get_z_index(cell), this->get_refinement_level(cell) - 1);
		if (this->cell_process.count(parent) > 0) {
			return parent;
		} else {
			return cell;
		}
	}

	/*!
	Returns the parent of given cell
	Returns the given cell if its refinement level == 0 or > maximum refinement level
	TODO: replace with get_{,all_}siblings?
	*/
	uint64_t get_parent_for_removed(const uint64_t cell) const
	{
		const int refinement_level = this->get_refinement_level(cell);
		assert(refinement_level >= 0);
		if (refinement_level == 0 || refinement_level > this->max_refinement_level) {
			return cell;
		}

		return get_cell_from_indices(this->get_x_index(cell), this->get_y_index(cell), this->get_z_index(cell), refinement_level - 1);
	}


	/*!
	Returns true if given cells exist and share at least one vertex.
	Returns false otherwise
	FIXME: only works for stencil size 1
	*/
	bool shared_vertex(const uint64_t cell1, const uint64_t cell2) const
	{
		if (this->cell_process.count(cell1) == 0 || this->cell_process.count(cell2) == 0) {
			return false;
		}
		if (this->is_neighbour(cell1, cell2) && this->is_neighbour(cell2, cell1)) {
			return true;
		} else {
			return false;
		}
	}

	/*!
	Returns true if given cells exist and share at least one edge.
	Returns false otherwise
	FIXME: only works for stencil size 1
	*/
	bool shared_edge(const uint64_t cell1, const uint64_t cell2) const
	{
		if (this->cell_process.count(cell1) == 0 || this->cell_process.count(cell2) == 0) {
			return false;
		}
		if (this->overlapping_indices(cell1, cell2) > 0 && this->is_neighbour(cell1, cell2) && this->is_neighbour(cell2, cell1)) {
			return true;
		} else {
			return false;
		}
	}

	/*!
	Returns true if given cells exist and share at least one face.
	Returns false otherwise
	FIXME: only works for stencil size 1
	*/
	bool shared_face(const uint64_t cell1, const uint64_t cell2) const
	{
		if (this->cell_process.count(cell1) == 0 || this->cell_process.count(cell2) == 0) {
			return false;
		}
		if (this->overlapping_indices(cell1, cell2) > 1 && this->is_neighbour(cell1, cell2) && this->is_neighbour(cell2, cell1)) {
			return true;
		} else {
			return false;
		}
	}


	/*!
	Returns the existing neighbours (that don't have children) of given cell even if it is on another process.

	Returns nothing if the following cases:
		-given cell has children
		-given doesn't exist on any process
	Doesn't use existing neighbour lists and hence is slow but works if for example given cell was moved to another process by load balancing.
	*/
	std::vector<uint64_t> find_neighbours_of(const uint64_t cell) const
	{
		std::vector<uint64_t> return_neighbours;

		if (cell == 0
		|| cell > this->max_cell_number
		|| this->cell_process.count(cell) == 0) {
			return return_neighbours;
		}

		if (cell != this->get_child(cell)) {
			return return_neighbours;
		}

		const uint64_t x_index = this->get_x_index(cell);
		const uint64_t y_index = this->get_y_index(cell);
		const uint64_t z_index = this->get_z_index(cell);

		// search neighbours in cells of the same size as the given cell (times neighbourhood size)
		const uint64_t size_in_indices = this->get_cell_size_in_indices(cell);

		const int refinement_level = this->get_refinement_level(cell);
		#ifdef DEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << cell << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << cell << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// can limit search due to maximum refinement level difference of 1 between neighbours
		const int search_min_ref_level = (refinement_level == 0) ? 0 : refinement_level - 1;
		const int search_max_ref_level = (refinement_level == this->max_refinement_level) ? refinement_level : refinement_level + 1;

		// must have some neighbours even if neighbourhood_size == 0
		const int temp_neighbourhood_size = (this->neighbourhood_size > 0) ? this->neighbourhood_size : 1;

		// grid length in indices
		const uint64_t x_length_in_indices = this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level);
		const uint64_t y_length_in_indices = this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level);
		const uint64_t z_length_in_indices = this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level);

		// search neighbourhood_size number of cells (of given cell's size) away from the given cell and not outside of the grid
		const uint64_t outer_min_x = (x_index < size_in_indices * temp_neighbourhood_size) ? 0 : x_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_x = (x_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 < x_length_in_indices) ? x_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 : x_length_in_indices - 1;

		const uint64_t outer_min_y = (y_index < size_in_indices * temp_neighbourhood_size) ? 0 : y_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_y = (y_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 < y_length_in_indices) ? y_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 : y_length_in_indices - 1;

		const uint64_t outer_min_z = (z_index < size_in_indices * temp_neighbourhood_size) ? 0 : z_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_z = (z_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 < z_length_in_indices) ? z_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 : z_length_in_indices - 1;

		// don't search within the given cell
		const uint64_t inner_max_x = x_index + size_in_indices - 1;
		const uint64_t inner_max_y = y_index + size_in_indices - 1;
		const uint64_t inner_max_z = z_index + size_in_indices - 1;

		boost::unordered_set<uint64_t> unique_neighbours = find_cells(
			outer_min_x, outer_min_y, outer_min_z,
			outer_max_x, outer_max_y, outer_max_z,
			x_index, y_index, z_index,
			inner_max_x, inner_max_y, inner_max_z,
			search_min_ref_level, search_max_ref_level);

		return_neighbours.reserve(unique_neighbours.size());
		return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
		return return_neighbours;
	}



	/*!
	Returns cells (which don't have children) that consider given cell as a neighbour but aren't considered neighbours by given cell.

	Returns nothing if the following cases:
		-given cell has children
		-given doesn't exist on any process
	Doesn't use existing neighbour lists and hence is slow but works if for example given cell was moved to another process by load balancing.
	*/
	std::vector<uint64_t> find_neighbours_to(const uint64_t cell) const
	{
		std::vector<uint64_t> return_neighbours;

		if (cell == 0
		|| cell > this->max_cell_number
		|| this->cell_process.count(cell) == 0) {
			return return_neighbours;
		}

		if (cell != this->get_child(cell)) {
			return return_neighbours;
		}

		const int refinement_level = this->get_refinement_level(cell);
		if (refinement_level == 0) {
			// largest cell cannot have neighbours_to
			return return_neighbours;
		}
		#ifdef DEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << cell << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << cell << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// can limit search due to maximum refinement level difference of 1 between neighbours
		const int search_min_ref_level = (refinement_level == 0) ? refinement_level : refinement_level - 1;
		const int search_max_ref_level = (refinement_level == this->max_refinement_level) ? refinement_level : refinement_level + 1;

		// search first for neighbours of the parent, then discard neighbours_of given cell
		const uint64_t parent = this->get_parent(cell);
		const uint64_t x_index = this->get_x_index(parent);
		const uint64_t y_index = this->get_y_index(parent);
		const uint64_t z_index = this->get_z_index(parent);

		// search neighbours in cells of the same size as the given cell's parent (times neighbourhood size)
		const uint64_t size_in_indices = this->get_cell_size_in_indices(parent);

		#ifdef DEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << cell << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << cell << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// must have some neighbours even if neighbourhood_size == 0
		const int temp_neighbourhood_size = (this->neighbourhood_size > 0) ? this->neighbourhood_size : 1;

		// grid length in indices
		const uint64_t x_length_in_indices = this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level);
		const uint64_t y_length_in_indices = this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level);
		const uint64_t z_length_in_indices = this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level);

		const uint64_t outer_min_x = (x_index < size_in_indices * temp_neighbourhood_size) ? 0 : x_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_x = (x_index + size_in_indices * temp_neighbourhood_size < x_length_in_indices) ? x_index + size_in_indices * temp_neighbourhood_size : x_length_in_indices - 1;

		const uint64_t outer_min_y = (y_index < size_in_indices * temp_neighbourhood_size) ? 0 : y_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_y = (y_index + size_in_indices * temp_neighbourhood_size < y_length_in_indices) ? y_index + size_in_indices * temp_neighbourhood_size : y_length_in_indices - 1;

		const uint64_t outer_min_z = (z_index < size_in_indices * temp_neighbourhood_size) ? 0 : z_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_z = (z_index + size_in_indices * temp_neighbourhood_size < z_length_in_indices) ? z_index + size_in_indices * temp_neighbourhood_size : z_length_in_indices - 1;

		// don't search within the given cell's parent
		const uint64_t inner_max_x = x_index + size_in_indices - 1;
		const uint64_t inner_max_y = y_index + size_in_indices - 1;
		const uint64_t inner_max_z = z_index + size_in_indices - 1;

		boost::unordered_set<uint64_t> unique_neighbours = find_cells(
			outer_min_x, outer_min_y, outer_min_z,
			outer_max_x, outer_max_y, outer_max_z,
			x_index, y_index, z_index,
			inner_max_x, inner_max_y, inner_max_z,
			search_min_ref_level, search_max_ref_level);

		return_neighbours.reserve(unique_neighbours.size());
		for (boost::unordered_set<uint64_t>::const_iterator neighbour = unique_neighbours.begin(); neighbour != unique_neighbours.end(); neighbour++) {
			if (!this->is_neighbour(cell, *neighbour)
			&& this->is_neighbour(*neighbour, cell)) {
				return_neighbours.push_back(*neighbour);
			}
		}

		return return_neighbours;
	}


	/*!
	Returns cells between given outer and inner shells (given as indices, inclusive) and refinement levels.

	Ignores cells diagonal to given inner shell if neighbourhood_size == 0.
	*/
	boost::unordered_set<uint64_t> find_cells
	(
		const uint64_t outer_min_x,
		const uint64_t outer_min_y,
		const uint64_t outer_min_z,
		const uint64_t outer_max_x,
		const uint64_t outer_max_y,
		const uint64_t outer_max_z,
		const uint64_t inner_min_x,
		const uint64_t inner_min_y,
		const uint64_t inner_min_z,
		const uint64_t inner_max_x,
		const uint64_t inner_max_y,
		const uint64_t inner_max_z,
		const int minimum_refinement_level,
		const int maximum_refinement_level
	) const
	{
		#ifdef DEBUG
		if (minimum_refinement_level > maximum_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid refinement levels given" << std::endl;
			exit(EXIT_FAILURE);
		}

		// check that outer shell makes sense
		if (outer_min_x > outer_max_x) {
			std::cerr << __FILE__ << ":" << __LINE__ << " outer_min_x > outer_max_x" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (outer_min_y > outer_max_y) {
			std::cerr << __FILE__ << ":" << __LINE__ << " outer_min_y > outer_max_y" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (outer_min_z > outer_max_z) {
			std::cerr << __FILE__ << ":" << __LINE__ << " outer_min_z > outer_max_z" << std::endl;
			exit(EXIT_FAILURE);
		}

		// check that inner shell makes sense
		if (inner_min_x > inner_max_x) {
			std::cerr << __FILE__ << ":" << __LINE__ << " inner_min_x > inner_max_x" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (inner_min_y > inner_max_y) {
			std::cerr << __FILE__ << ":" << __LINE__ << " inner_min_y > inner_max_y" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (inner_min_z > inner_max_z) {
			std::cerr << __FILE__ << ":" << __LINE__ << " inner_min_z > inner_max_z" << std::endl;
			exit(EXIT_FAILURE);
		}

		// check that inner shell is at least partially inside of the outer shell
		if (inner_min_x == outer_min_x
		&& inner_min_y == outer_min_y
		&& inner_min_z == outer_min_z
		&& inner_max_x == outer_max_x
		&& inner_max_y == outer_max_y
		&& inner_max_z == outer_max_z) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Nothing to search" << std::endl;
			exit(EXIT_FAILURE);
		}

		// check that inner shell isn't outside of the outer shell
		if (inner_min_x < outer_min_x
		|| inner_min_y < outer_min_y
		|| inner_min_z < outer_min_z
		|| inner_max_x > outer_max_x
		|| inner_max_y > outer_max_y
		|| inner_max_z > outer_max_z) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Inner shell is outside of the outer shell" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		boost::unordered_set<uint64_t> result;

		const uint64_t index_increase = uint64_t(1) << (this->max_refinement_level - maximum_refinement_level);
		for (uint64_t x = outer_min_x; x <= outer_max_x; x += index_increase) {
			for (uint64_t y = outer_min_y; y <= outer_max_y; y += index_increase) {
				for (uint64_t z = outer_min_z; z <= outer_max_z; z += index_increase) {

					// skip inner shell
					if (x >= inner_min_x
					&& x <= inner_max_x
					&& y >= inner_min_y
					&& y <= inner_max_y
					&& z >= inner_min_z
					&& z <= inner_max_z) {
						continue;
					}

					if (this->neighbourhood_size == 0) {

						// don't search diagonally
						int overlaps = 0;
						if (this->indices_overlap(inner_min_x, 1 + inner_max_x - inner_min_x, x, 1)) {
							overlaps++;
						}
						if (this->indices_overlap(inner_min_y, 1 + inner_max_y - inner_min_y, y, 1)) {
							overlaps++;
						}
						if (this->indices_overlap(inner_min_z, 1 + inner_max_z - inner_min_z, z, 1)) {
							overlaps++;
						}

						if (overlaps < 2) {
							continue;
						}
					}

					const uint64_t neighbour = this->get_cell_from_indices(x, y, z, minimum_refinement_level, maximum_refinement_level);

					#ifdef DEBUG
					if (neighbour == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " No neighbour found between refinement levels [" << minimum_refinement_level << ", " << maximum_refinement_level << "] at indices " << x << " " << y << " " << z << std::endl;
						const uint64_t smallest = this->get_cell_from_indices(x, y, z, minimum_refinement_level, maximum_refinement_level);
						std::cerr << __FILE__ << ":" << __LINE__ << " smallest neighbour there is " << smallest << " with refinement level " << this->get_refinement_level(smallest) << std::endl;
					}

					if (neighbour > this->max_cell_number) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour can't exist" << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif

					/*
					When searching for neighbours_to cells may exist with larger refinement level than given in find_neighbours_to and they don't consider this cell as a neighbour.
					*/
					if (neighbour != this->get_child(neighbour)) {
						continue;
					}

					result.insert(neighbour);
				}
			}
		}

		return result;
	}


	/*!
	Returns the smallest cell at given coordinates or 0 if outside of the grid
	*/
	uint64_t get_cell(const double x, const double y, const double z) const
	{
		#ifdef DCCRG_ARBITRARY_STRETCH
		if (x < this->x_coordinates[0]
		|| x > this->x_coordinates[this->x_coordinates.size() - 1]
		|| y < this->y_coordinates[0]
		|| y > this->y_coordinates[this->y_coordinates.size() - 1]
		|| z < this->z_coordinates[0]
		|| z > this->z_coordinates[this->z_coordinates.size() - 1]) {
			return 0;
		}
		#else
		if (x < this->geometry.get_x_start()
		|| x > this->geometry.get_x_start() + this->geometry.get_x_length() * this->geometry.get_unrefined_cell_x_size()
		|| y < this->geometry.get_y_start()
		|| y > this->geometry.get_y_start() + this->geometry.get_y_length() * this->geometry.get_unrefined_cell_y_size()
		|| z < this->geometry.get_z_start()
		|| z > this->geometry.get_z_start() + this->geometry.get_z_length() * this->geometry.get_unrefined_cell_z_size()) {
			return 0;
		}
		#endif

		return this->get_cell_from_indices(this->geometry.get_x_index(x), this->geometry.get_y_index(y), this->geometry.get_z_index(z), 0, this->geometry.get_maximum_refinement_level());
	}

	#ifdef DCCRG_ARBITRARY_STRETCH
	#else
	/*!
	The following return the start and end corners of the grid
	*/
	double get_x_start(void) const { return this->geometry.get_x_start(); }
	double get_y_start(void) const { return this->geometry.get_y_start(); }
	double get_z_start(void) const { return this->geometry.get_z_start(); }
	double get_x_end(void) const { return this->geometry.get_x_end(); }
	double get_y_end(void) const { return this->geometry.get_y_end(); }
	double get_z_end(void) const { return this->geometry.get_z_end(); }
	#endif

	/*!
	The following return the grid length in unrefined cells
	*/
	uint64_t get_x_length(void) const { return this->geometry.get_x_length(); }
	uint64_t get_y_length(void) const { return this->geometry.get_y_length(); }
	uint64_t get_z_length(void) const { return this->geometry.get_z_length(); }


	/*!
	These return the index of the cell with given id in x, y or z direction of the grid, starting from 0
	For cells that are larger than the smallest possible according to max_refinement_level, the index closest to the grids starting corner is returned
	 */
	uint64_t get_x_index(uint64_t id) const
	{
		assert(id);
		assert(id <= this->max_cell_number);

		// substract ids of larger cells
		int refinement_level = get_refinement_level(id);
		for (int i = 0; i < refinement_level; i++) {
			id -= this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << i * 3);
		}

		// get the index at this cells refinement level
		id -= 1;	// cell numbering starts from 1
		uint64_t this_level_index = id % (this->geometry.get_x_length() * (uint64_t(1) << refinement_level));

		assert(this_level_index * (uint64_t(1) << (this->max_refinement_level - refinement_level)) < this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level));
		return this_level_index * (uint64_t(1) << (max_refinement_level - refinement_level));
	}

	uint64_t get_y_index(uint64_t id) const
	{
		assert(id);
		assert(id <= this->max_cell_number);

		// substract ids of larger cells
		int refinement_level = get_refinement_level(id);
		for (int i = 0; i < refinement_level; i++) {
			id -= this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << i * 3);
		}

		// get the index at this cells refinement level
		id -= 1;	// cell numbering starts from 1
		uint64_t this_level_index =  (id / (this->geometry.get_x_length() * (uint64_t(1) << refinement_level))) % (this->geometry.get_y_length()  * (uint64_t(1) << refinement_level));

		return this_level_index * (uint64_t(1) << (max_refinement_level - refinement_level));
	}

	uint64_t get_z_index(uint64_t id) const
	{
		assert(id);
		assert(id <= this->max_cell_number);

		// substract ids of larger cells
		int refinement_level = get_refinement_level(id);
		for (int i = 0; i < refinement_level; i++) {
			id -= this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << i * 3);
		}

		// get the index at this cells refinement level
		id -= 1;	// cell numbering starts from 1
		uint64_t this_level_index =  id / (this->geometry.get_x_length() * this->geometry.get_y_length() *  (uint64_t(1) << 2 * refinement_level));

		return this_level_index * (uint64_t(1) << (max_refinement_level - refinement_level));
	}


	/*!
	Returns the cell of given refinement level at given indices even if it doesn't exist.
	*/
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int refinement_level) const
	{
		assert(x_index < this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level));
		assert(y_index < this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level));
		assert(z_index < this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level));
		assert(refinement_level <= this->max_refinement_level);

		uint64_t id = 1;

		// add ids of larger cells
		for (int i = 0; i < refinement_level; i++) {
			id += this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << i * 3);
		}

		// convert to indices of this cells refinement level
		uint64_t this_level_x_index = x_index / (uint64_t(1) << (max_refinement_level - refinement_level));
		uint64_t this_level_y_index = y_index / (uint64_t(1) << (max_refinement_level - refinement_level));
		uint64_t this_level_z_index = z_index / (uint64_t(1) << (max_refinement_level - refinement_level));

		// get the size of the grid in terms of cells of this level
		uint64_t this_level_x_length = this->geometry.get_x_length() *  (uint64_t(1) << refinement_level);
		uint64_t this_level_y_length = this->geometry.get_y_length() *  (uint64_t(1) << refinement_level);

		id += this_level_x_index + this_level_y_index * this_level_x_length + this_level_z_index * this_level_x_length * this_level_y_length;

		assert(id > 0);
		assert(id <= this->max_cell_number);
		return id;
	}


	// Returns the lengths of given cell in indices in every direction
	uint64_t get_cell_size_in_indices(const uint64_t cell) const
	{
		assert(cell);
		return uint64_t(1) << (this->max_refinement_level - this->get_refinement_level(cell));
	}


	/*!
	Removes user data of refined and unrefined cells from this process.
	*/
	void clear_refined_unrefined_data(void)
	{
		this->refined_cell_data.clear();
		this->unrefined_cell_data.clear();
	}


	/*!
	Sets the given option for non-hierarchial partitioning.

	Does nothing if option name is one of: RETURN_LISTS, EDGE_WEIGHT_DIM, NUM_GID_ENTRIES, OBJ_WEIGHT_DIM
	Call this with name = LB_METHOD and value = HIER to use hierarchial partitioning and set those options using the other function with this name.
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
	Adds a new level for hierarchial partitioning, with each part of that level having given number of processes.

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
		-option name is one of: RETURN_LISTS, ...
		-given level doesn't exist
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

		for (boost::unordered_map<std::string, std::string>::const_iterator option = this->partitioning_options[hierarchial_partitioning_level].begin(); option != this->partitioning_options[hierarchial_partitioning_level].end(); option++) {
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
	int get_process(const uint64_t cell)
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
		this->pin(cell, this->comm.rank());
	}

	/*!
	Given cell is sent to the given process and kept there during subsequent load balancing.

	Does nothing in the following cases:
		-given cell doesn't exist
		-given cell exists on another process
		-given cell has children
		-given process doesn't exist
	*/
	void pin(const uint64_t cell, const int process)
	{
		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cell_process.at(cell) != this->comm.rank()) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		if (process < 0 || process >= this->comm.size()) {
			return;
		}

		// do nothing if the request already exists
		if (this->pin_requests.count(cell) > 0
		&& this->pin_requests.at(cell) == process) {
			return;
		}

		this->new_pin_requests[cell] = process;
	}

	/*!
	Allows the given cell to be moved to another process during subsequent load balancing.

	Does nothing in the following cases:
		-given cell has children
		-given cell doesn't exist
		-given cell exists on another process
	*/
	void unpin(const uint64_t cell)
	{
		if (this->cell_process.count(cell) == 0) {
			return;
		}

		if (this->cell_process.at(cell) != this->comm.rank()) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		if (this->pin_requests.count(cell) > 0) {
			this->new_pin_requests[cell] = -1;
		} else {
			this->new_pin_requests.erase(cell);
		}
	}

	/*!
	Executes unpin(cell) for all cells on this process.
	*/
	void unpin_local_cells(void)
	{
		#ifdef DEBUG
		// check that all child cells on this process are also in this->cells.
		for (auto i = this->cell_process.cbegin(); i != this->cell_process.cend(); i++) {
			const uint64_t cell = i->first;

			if (this->cell_process.at(cell) != this->comm.rank()) {
				return;
			}

			if (cell == this->get_child(cell)) {
				if (this->cells.count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " should be in this->cells of process " << this->comm.rank() << std::endl;
					abort();
				}
			} else {
				if (this->cells.count(cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " shouldn't be in this->cells of process " << this->comm.rank() << std::endl;
					abort();
				}
			}
		}
		#endif

		for (auto i = this->cells.cbegin(); i != this->cells.cend(); i++) {
			this->unpin(i->first);
		}
	}

	/*!
	Allows all cells of all processes to be moved to another process during subsequent load balancing.

	Must be called simultaneously on all processes.
	*/
	void unpin_all_cells(void)
	{
		this->comm.barrier();
		this->new_pin_requests.clear();
		this->pin_requests.clear();
	}



private:

	bool initialized;

	CellGeometry geometry;
	#ifdef DCCRG_ARBITRARY_STRETCH
	/*!
	The coordinates of unrefined cells in respective directions
	First value is the starting point of the grid, following ith value is the end point of the ith unrefined cell
	*/
	std::vector<double> x_coordinates, y_coordinates, z_coordinates;
	#else
	// length of unrefined cells in all directions
	double cell_x_size, cell_y_size, cell_z_size;
	#endif
	// maximum refinemet level of any cell in the grid, 0 means unrefined
	int max_refinement_level;
	// the id of the last cell in the grid at maximum refinement level
	uint64_t max_cell_number;
	// size of the neighbour stencil of a cells in cells (of the same size as the cell itself)
	unsigned int neighbourhood_size;
	// the grid is distributed between these processes
	boost::mpi::communicator comm;

	// cells and their data on this process
	boost::unordered_map<uint64_t, UserData> cells;

	// cell on this process and its neighbours
	boost::unordered_map<uint64_t, std::vector<uint64_t> > neighbours;

	/*!
	Cell on this process and those cells that aren't neighbours of this cell but whose neighbour this cell is.
	For example with a stencil size of 1 in the following grid:
\verbatim
|-----------|
|     |5 |6 |
|  1  |--|--|
|     |9 |10|
|-----------|
\endverbatim
	neighbours_to[6] = 1 because neighbours[6] = 5, 9, 10 while
	neighbours_to[5] is empty because neighbours[5] = 1, 6, 9, 10
	*/
	boost::unordered_map<uint64_t, std::vector<uint64_t> > neighbours_to;

	// on which process every cell in the grid is
	boost::unordered_map<uint64_t, int> cell_process;

	// cells on this process that have a neighbour on another process or are considered as a neighbour of a cell on another process
	boost::unordered_set<uint64_t> cells_with_remote_neighbours;

	// cells on other processes that have a neighbour on this process or are considered as a neighbour of a cell on this process
	boost::unordered_set<uint64_t> remote_cells_with_local_neighbours;

	// remote neighbours and their data, of cells on this process
	boost::unordered_map<uint64_t, UserData> remote_neighbours;

	#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::unordered_map<int, std::vector<MPI_Request> > send_requests, receive_requests;
	#else
	boost::unordered_map<int, std::vector<boost::mpi::request> > send_requests, receive_requests;
	#endif

	// cells whose data has to be received / sent by this process from the process as the key
	boost::unordered_map<int, std::vector<uint64_t> > cells_to_send, cells_to_receive;

	#ifndef DCCRG_SEND_SINGLE_CELLS
	#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
	// storage for cells' user data that awaits transfer to or from this process
	boost::unordered_map<int, std::vector<UserData> > incoming_data, outgoing_data;
	#endif
	#endif

	// cells to be refined / unrefined after a call to stop_refining()
	boost::unordered_set<uint64_t> cells_to_refine, cells_to_unrefine;

	// stores user data of cells whose children were created while refining
	boost::unordered_map<uint64_t, UserData> refined_cell_data;
	// stores user data of cells that were removed while unrefining
	boost::unordered_map<uint64_t, UserData> unrefined_cell_data;

	// cell that should be kept on a particular process
	boost::unordered_map<uint64_t, int> pin_requests;
	// pin requests given since that last time load was balanced
	boost::unordered_map<uint64_t, int> new_pin_requests;

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


	/*!
	Moves cells between processes due to load balancing or user request.

	Recalculates neighbour lists, etc.
	*/
	void move_cells(
		const boost::unordered_set<uint64_t>& added_cells,
		const boost::unordered_set<uint64_t>& removed_cells)
	{
		// TODO: get rid of added_cells and removed_cells and use cells_to_send and receive instead?
		#ifdef DEBUG
		if (!this->verify_remote_neighbour_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " virhe" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// clear user data which is about to get old
		this->cells_with_remote_neighbours.clear();
		this->remote_cells_with_local_neighbours.clear();
		this->remote_neighbours.clear();
		this->cells_to_refine.clear();
		this->refined_cell_data.clear();
		this->cells_to_unrefine.clear();
		this->unrefined_cell_data.clear();

		/*
		Calculate where cells have migrated to update internal data structures
		Any cell can end up on any process and any neighbour of any cell can end up on yet another process
		*/

		// removed cells on all processes
		std::vector<uint64_t> temp_removed_cells(removed_cells.begin(), removed_cells.end());
		std::vector<std::vector<uint64_t> > all_removed_cells;
		all_gather(this->comm, temp_removed_cells, all_removed_cells);

		// created cells on all processes
		std::vector<uint64_t> temp_added_cells(added_cells.begin(), added_cells.end());
		std::vector<std::vector<uint64_t> > all_added_cells;
		all_gather(this->comm, temp_added_cells, all_added_cells);

		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->cells
		#else
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->cells
		#endif
		#endif
		);

		#ifdef DEBUG
		// check that there are no duplicate adds / removes
		boost::unordered_set<uint64_t> all_adds, all_removes;

		for (std::vector<std::vector<uint64_t> >::const_iterator item = all_removed_cells.begin(); item != all_removed_cells.end(); item++) {
			for (std::vector<uint64_t>::const_iterator removed_cell = item->begin(); removed_cell != item->end(); removed_cell++) {
				if (all_removes.count(*removed_cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << *removed_cell << " was already removed" << std::endl;
					abort();
				}
				all_removes.insert(*removed_cell);
			}
		}

		for (std::vector<std::vector<uint64_t> >::const_iterator item = all_added_cells.begin(); item != all_added_cells.end(); item++) {
			for (std::vector<uint64_t>::const_iterator added_cell = item->begin(); added_cell != item->end(); added_cell++) {
				if (all_adds.count(*added_cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << *added_cell << " was already removed" << std::endl;
					abort();
				}
				all_adds.insert(*added_cell);
			}
		}

		// check that cells were removed by their process
		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
			for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {
				if (this->cell_process.at(*removed_cell) != cell_remover) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << *removed_cell << " doesn't belong to process " << cell_remover << std::endl;
					abort();
				}
			}
		}
		#endif

		// update cell to process mappings
		for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {
				this->cell_process.at(*created_cell) = cell_creator;
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

		// create neighbour lists for cells without children that came to this process
		for (boost::unordered_set<uint64_t>::const_iterator added_cell = added_cells.begin(); added_cell != added_cells.end(); added_cell++) {

			if (*added_cell != this->get_child(*added_cell)) {
				continue;
			}

			this->neighbours[*added_cell] = this->find_neighbours_of(*added_cell);
			this->neighbours_to[*added_cell] = this->find_neighbours_to(*added_cell);
		}

		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->cells
		#endif
		#endif
		);
		this->wait_user_data_transfer_sends();
		this->cells_to_send.clear();
		this->cells_to_receive.clear();

		// free user data and neighbour lists of cells removed from this process
		for (boost::unordered_set<uint64_t>::const_iterator removed_cell = removed_cells.begin(); removed_cell != removed_cells.end(); removed_cell++) {
			this->cells.erase(*removed_cell);
			this->neighbours.erase(*removed_cell);
			this->neighbours_to.erase(*removed_cell);
		}

		this->update_remote_neighbour_info();

		this->recalculate_neighbour_update_send_receive_lists();

		#ifdef DEBUG
		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " The grid is inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_neighbours()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour lists are incorrect" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_remote_neighbour_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " virhe" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}


	/*!
	Updates user pin requests globally based on new_pin_requests.

	Must be called simultaneously on all processes.
	*/
	void update_pin_requests(void)
	{
		std::vector<uint64_t> new_pinned_cells, new_pinned_processes;

		new_pinned_cells.reserve(this->new_pin_requests.size());
		new_pinned_processes.reserve(this->new_pin_requests.size());
		for (auto item = this->new_pin_requests.cbegin(); item != this->new_pin_requests.cend(); item++) {
			new_pinned_cells.push_back(item->first);
			new_pinned_processes.push_back(item->second);
		}

		std::vector<std::vector<uint64_t> > all_new_pinned_cells, all_new_pinned_processes;
		all_gather(this->comm, new_pinned_cells, all_new_pinned_cells);
		all_gather(this->comm, new_pinned_processes, all_new_pinned_processes);

		for (int process = 0; process < int(all_new_pinned_cells.size()); process++) {
			for (unsigned int i = 0; i < all_new_pinned_cells[process].size(); i++) {

				const int requested_process = all_new_pinned_processes[process][i];

				if (requested_process == -1) {
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
	Calculates what to send and where during a remote neighbour data update.

	Assumes up-to-date neighbour lists, clears previous send / receive lists.
	*/
	void recalculate_neighbour_update_send_receive_lists(void)
	{
		// clear previous lists
		this->cells_to_send.clear();
		this->cells_to_receive.clear();

		// only send a cell to a process once
		boost::unordered_map<int, boost::unordered_set<uint64_t> > unique_cells_to_send, unique_cells_to_receive;

		// calculate new lists for neighbour data updates
		for (boost::unordered_set<uint64_t>::const_iterator cell = this->cells_with_remote_neighbours.begin(); cell != this->cells_with_remote_neighbours.end(); cell++) {

			assert(*cell == this->get_child(*cell));

			int current_process = this->comm.rank();

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*cell].begin(); neighbour != this->neighbours[*cell].end(); neighbour++) {
				if (this->cell_process[*neighbour] != current_process) {
					// *neighbours process has to send *neighbours cell data to current_process
					unique_cells_to_receive[this->cell_process[*neighbour]].insert(*neighbour);
					// current process has to send currents cell data to neighbour
					unique_cells_to_send[this->cell_process[*neighbour]].insert(*cell);
				}
			}

			// also cells that have this one as neighbour
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[*cell].begin(); neighbour != this->neighbours_to[*cell].end(); neighbour++) {
				if (this->cell_process[*neighbour] != current_process) {
					// *neighbours process has to send *neighbours cell data to current_process
					unique_cells_to_receive[this->cell_process[*neighbour]].insert(*neighbour);
					// current process has to send currents cell data to neighbour
					unique_cells_to_send[this->cell_process[*neighbour]].insert(*cell);
				}
			}
		}

		// populate final send / receive list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator receiver = unique_cells_to_send.begin(); receiver != unique_cells_to_send.end(); receiver++) {
			this->cells_to_send[receiver->first].insert(this->cells_to_send[receiver->first].begin(), receiver->second.begin(), receiver->second.end());
			sort(this->cells_to_send[receiver->first].begin(), this->cells_to_send[receiver->first].end());
		}

		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator sender = unique_cells_to_receive.begin(); sender != unique_cells_to_receive.end(); sender++) {
			this->cells_to_receive[sender->first].insert(this->cells_to_receive[sender->first].begin(), sender->second.begin(), sender->second.end());
			sort(this->cells_to_receive[sender->first].begin(), this->cells_to_receive[sender->first].end());
		}
	}


	/*!
	Updates given cell's neighbour and neighbour_to lists.

	Does nothing in the following cases:
		-given cell has children
		-given isn't on this process
	Uses existing neighbour list of given cell and hence is fast but gives the wrong result if the cell's neighbour list isn't up-to-date (excluding refined / unrefined neighbours) or doesn't exist.
	Also takes into account that neighbours in given cell's neighbour list could've been refined / urefined.
	*/
	void update_neighbours(const uint64_t cell)
	{
		const uint64_t parent = this->get_parent(cell);

		#ifdef DEBUG
		if (cell == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given: " << cell << std::endl;
			exit(EXIT_FAILURE);
		}

		if (this->get_refinement_level(cell) > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of given cell (" << cell << ") is too large: " << this->get_refinement_level(cell) << std::endl;
			exit(EXIT_FAILURE);
		}

		if (this->get_refinement_level(cell) < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid refinement level for cell " << cell << ": " << this->get_refinement_level(cell) << std::endl;
			exit(EXIT_FAILURE);
		}

		if (parent == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent for cell " << cell << ": " << parent << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		if (this->cell_process.at(cell) != this->comm.rank()) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		bool cell_result_of_refining;
		if (this->cells_to_refine.count(parent) > 0) {
			cell_result_of_refining = true;
		} else {
			cell_result_of_refining = false;
		}

		// choose which neighbour lists to use
		std::vector<uint64_t> old_neighbours;
		// use given cell's parent's neighbour lists
		if (cell_result_of_refining) {

			#ifdef DEBUG
			if (this->neighbours.count(parent) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour list of cell's " << cell << " parent " << this->get_parent(cell) << " doesn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (this->neighbours_to.count(parent) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour_to list of cell's " << cell << " parent " << this->get_parent(cell) << " doesn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			const std::vector<uint64_t> siblings = this->get_all_children(parent);

			// don't insert the given cell itself
			for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {
				// TODO: use boost::phoenix (e.g. http://www.boost.org/doc/libs/1_45_0/libs/spirit/phoenix/example/users_manual/if.cpp) instead of a for loop?
				if (*sibling != cell) {
					old_neighbours.push_back(*sibling);

					#ifdef DEBUG
					if (*sibling == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Invalid sibling" << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif
				}
			}

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(parent).begin(); neighbour != this->neighbours.at(parent).end(); neighbour++) {
				old_neighbours.push_back(*neighbour);

				#ifdef DEBUG
				if (*neighbour == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour for parent of cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

			for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(parent).begin(); neighbour_to != this->neighbours_to.at(parent).end(); neighbour_to++) {
				old_neighbours.push_back(*neighbour_to);

				#ifdef DEBUG
				if (*neighbour_to == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour_to for parent of cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

		// use given cell's neighbour lists
		} else {

			#ifdef DEBUG
			if (this->neighbours.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour list of cell " << cell << " doesn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}

			if (this->neighbours_to.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour_to list of cell " << cell << " doesn't exist" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {
				old_neighbours.push_back(*neighbour);

				#ifdef DEBUG
				if (*neighbour == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour for cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

			for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(cell).begin(); neighbour_to != this->neighbours_to.at(cell).end(); neighbour_to++) {
				old_neighbours.push_back(*neighbour_to);

				#ifdef DEBUG
				if (*neighbour_to == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour_to for cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

		}
		this->neighbours[cell].clear();
		this->neighbours_to[cell].clear();

		// neighbour candidates of given cell based on old neighbour lists
		boost::unordered_set<uint64_t> neighbour_candidates;

		for (std::vector<uint64_t>::const_iterator old_neighbour = old_neighbours.begin(); old_neighbour != old_neighbours.end(); old_neighbour++) {

			#ifdef DEBUG
			if (*old_neighbour == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid old neighbour for cell " << cell << std::endl;
				exit(EXIT_FAILURE);
			}

			if (*old_neighbour == cell) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " has itself as an old neighbour" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			// add the parents of unrefined cells as neighbour candidates instead of the unrefined cells
			bool add_parent = false;

			// check all siblings because only one is recorded even though all have been removed
			std::vector<uint64_t> siblings = this->get_all_children(this->get_parent_for_removed(*old_neighbour));
			for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {
				if (this->cells_to_unrefine.count(*sibling) > 0) {

					#ifdef DEBUG
					if (this->cell_process.count(*sibling) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << *sibling << " shouldn't exist" << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif

					add_parent = true;
					break;
				}
			}

			if (add_parent) {
				#ifdef DEBUG
				if (this->cells_to_refine.count(*old_neighbour) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Old neighbour " << *old_neighbour << " was refined and it or its sibling was also unrefined" << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif

				neighbour_candidates.insert(this->get_parent_for_removed(*old_neighbour));
				continue;
			}

			// add an unchanged old neighbour as a neighbour candidate
			if (this->cells_to_refine.count(*old_neighbour) == 0) {
				neighbour_candidates.insert(*old_neighbour);

				#ifdef DEBUG
				if (this->cell_process.count(*old_neighbour) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Old neighbour " << *old_neighbour << " doesn't exist but was supposed to have been refined" << std::endl;
					exit(EXIT_FAILURE);
				}

				if (*old_neighbour != this->get_child(*old_neighbour)) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Old neighbour " << *old_neighbour << " for cell " << cell << " (child of " << this->get_parent(cell) << ") has children (at least " << this->get_child(*old_neighbour) << ") but wasn't refined" << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif

			// add children of refined cells as neighbour candidates
			} else {

				std::vector<uint64_t> children = this->get_all_children(*old_neighbour);
				neighbour_candidates.insert(children.begin(), children.end());

				#ifdef DEBUG
				for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
					if (this->cell_process.count(*child) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << *child << " doesn't exist" << std::endl;
						exit(EXIT_FAILURE);
					}
				}
				#endif
			}
		}

		// parents of unrefined cells might also become neighbours_to of given cell
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = this->cells_to_unrefine.begin(); unrefined != this->cells_to_unrefine.end(); unrefined++) {
			neighbour_candidates.insert(this->get_parent_for_removed(*unrefined));
		}

		// don't include self as a neighbour candidate
		neighbour_candidates.erase(cell);

		#ifdef DEBUG
		if (neighbour_candidates.count(0) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell in neighbour candidates" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// update neighbour lists from neighbour candidates
		for (boost::unordered_set<uint64_t>::const_iterator candidate = neighbour_candidates.begin(); candidate != neighbour_candidates.end(); candidate++) {
			if (this->is_neighbour(cell, *candidate)) {
				this->neighbours[cell].push_back(*candidate);
			} else if (this->is_neighbour(*candidate, cell)) {
				this->neighbours_to[cell].push_back(*candidate);
			}

			#ifdef DEBUG
			if (*candidate != this->get_child(*candidate)) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour candidate " << *candidate << " of cell " << cell << " has children" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif
		}

		#ifdef DEBUG
		if (!this->verify_neighbours(cell)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour update failed for cell " << cell << " (child of " << this->get_parent(cell) << ")" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}


	/*!
	Updates the remote neighbour info of given cell on this process without children.

	Uses current neighbour lists.
	Does nothing if given cell doesn't exist on this process or has children
	*/
	void update_remote_neighbour_info(const uint64_t cell)
	{
		if (this->cells.count(cell) == 0) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		// TODO: also update remote_cells_with_local_neighbours
		this->cells_with_remote_neighbours.erase(cell);

		#ifdef DEBUG
		if (this->neighbours.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour list for cell " << cell << " doesn't exist" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (this->neighbours_to.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbours_to list for cell " << cell << " doesn't exist" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// neighbours of given cell
		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {
			if (this->cell_process.at(*neighbour) != this->comm.rank()) {
				this->cells_with_remote_neighbours.insert(cell);
				this->remote_cells_with_local_neighbours.insert(*neighbour);
			}
		}
		// cells with given cell as neighbour
		for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(cell).begin(); neighbour_to != this->neighbours_to.at(cell).end(); neighbour_to++) {
			if (this->cell_process.at(*neighbour_to) != this->comm.rank()) {
				this->cells_with_remote_neighbours.insert(cell);
				this->remote_cells_with_local_neighbours.insert(*neighbour_to);
			}
		}

		#ifdef DEBUG
		if (!this->verify_remote_neighbour_info(cell)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info for cell " << cell << " is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}


	/*!
	Updates the remote neighbour info of all cells on this process without children.

	Uses current neighbour lists.
	*/
	void update_remote_neighbour_info(void)
	{
		// TODO this probably can't be optimized without storing neighbour lists also for remote neighbours
		this->cells_with_remote_neighbours.clear();
		this->remote_cells_with_local_neighbours.clear();

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

			if (cell->first != this->get_child(cell->first)) {
				continue;
			}

			this->update_remote_neighbour_info(cell->first);

			#ifdef DEBUG
			if (!this->verify_remote_neighbour_info(cell->first)) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info for cell " << cell->first << " is not consistent" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif
		}

		#ifdef DEBUG
		if (!this->verify_remote_neighbour_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}


	/*!
	Same as find_neighbours_of but for the parent of given cell and different assumptions.
	*/
	std::vector<uint64_t> find_neighbours_of_parent(const uint64_t cell) const
	{
		#ifdef DEBUG
		if (cell == 0
		|| cell > this->max_cell_number) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		std::vector<uint64_t> return_neighbours;

		if (cell != this->get_child(cell)) {
			return return_neighbours;
		}

		const uint64_t parent = this->get_parent(cell);	// TODO just use find_neighbours_of(parent)?
		#ifdef DEBUG
		if (parent == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent for cell " << cell << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		const uint64_t x_index = this->get_x_index(parent);
		const uint64_t y_index = this->get_y_index(parent);
		const uint64_t z_index = this->get_z_index(parent);

		// search neighbours in cells of the same size as the given cell (times neighbourhood size)
		const uint64_t size_in_indices = this->get_cell_size_in_indices(parent);

		const int refinement_level = this->get_refinement_level(parent);
		#ifdef DEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << parent << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << parent << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// can limit search due to maximum refinement level difference of 1 between neighbours
		const int search_min_ref_level = (refinement_level == 0) ? 0 : refinement_level - 1;
		const int search_max_ref_level = (refinement_level == this->max_refinement_level) ? refinement_level : refinement_level + 1;

		// must have some neighbours even if neighbourhood_size == 0
		const int temp_neighbourhood_size = (this->neighbourhood_size > 0) ? this->neighbourhood_size : 1;

		// grid length in indices
		const uint64_t x_length_in_indices = this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level);
		const uint64_t y_length_in_indices = this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level);
		const uint64_t z_length_in_indices = this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level);

		// search neighbourhood_size number of cells (of given cell's parent's size) away from the given cell and not outside of the grid
		const uint64_t outer_min_x = (x_index < size_in_indices * temp_neighbourhood_size) ? 0 : x_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_x = (x_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 < x_length_in_indices) ? x_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 : x_length_in_indices - 1;

		const uint64_t outer_min_y = (y_index < size_in_indices * temp_neighbourhood_size) ? 0 : y_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_y = (y_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 < y_length_in_indices) ? y_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 : y_length_in_indices - 1;

		const uint64_t outer_min_z = (z_index < size_in_indices * temp_neighbourhood_size) ? 0 : z_index - size_in_indices * temp_neighbourhood_size;
		const uint64_t outer_max_z = (z_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 < z_length_in_indices) ? z_index + size_in_indices * (1 + temp_neighbourhood_size) - 1 : z_length_in_indices - 1;

		// don't search within the given cell
		const uint64_t inner_max_x = x_index + size_in_indices - 1;
		const uint64_t inner_max_y = y_index + size_in_indices - 1;
		const uint64_t inner_max_z = z_index + size_in_indices - 1;

		boost::unordered_set<uint64_t> unique_neighbours = find_cells(
			outer_min_x, outer_min_y, outer_min_z,
			outer_max_x, outer_max_y, outer_max_z,
			x_index, y_index, z_index,
			inner_max_x, inner_max_y, inner_max_z,
			search_min_ref_level, search_max_ref_level);

		return_neighbours.reserve(unique_neighbours.size());
		return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
		return return_neighbours;
	}


	/*!
	Returns true if cell1 considers cell2 as a neighbour, even if neither of them exists
	*/
	bool is_neighbour(const uint64_t cell1, const uint64_t cell2) const
	{
		assert(cell1 > 0);
		assert(cell1 <= this->max_cell_number);
		assert(cell2 > 0);
		assert(cell2 <= this->max_cell_number);

		const uint64_t cell1_x_index = this->get_x_index(cell1), cell1_y_index = this->get_y_index(cell1), cell1_z_index = this->get_z_index(cell1);
		const uint64_t cell2_x_index = this->get_x_index(cell2), cell2_y_index = this->get_y_index(cell2), cell2_z_index = this->get_z_index(cell2);
		const uint64_t cell1_size = this->get_cell_size_in_indices(cell1);
		const uint64_t cell2_size = this->get_cell_size_in_indices(cell2);

		// must have some neighbours even if neighbourhood_size == 0
		const unsigned int temp_neighbourhood_size = (this->neighbourhood_size > 0) ? this->neighbourhood_size : 1;

		const uint64_t dindex1 = cell2_size + cell1_size * temp_neighbourhood_size;
		const uint64_t dindex2 = (cell2_size < cell1_size) ? cell1_size * temp_neighbourhood_size + cell2_size : cell1_size * temp_neighbourhood_size;

		if (cell1_x_index < cell2_x_index + dindex1
		&& cell1_y_index < cell2_y_index + dindex1
		&& cell1_z_index < cell2_z_index + dindex1
		&& cell1_x_index + dindex2 >= cell2_x_index
		&& cell1_y_index + dindex2 >= cell2_y_index
		&& cell1_z_index + dindex2 >= cell2_z_index) {
			if (this->neighbourhood_size == 0 && this->overlapping_indices(cell1, cell2) < 2) {
				// diagonal cell isn't a neighbour
				return false;
			}
			return true;
		}
		return false;
	}


	/*!
	Given a cell that exists and has children returns one of the children
	Returns the given cell if it doesn't have children or 0 if the cell doesn't exist
	*/
	uint64_t get_child(const uint64_t cell) const
	{
		if (this->cell_process.count(cell) == 0) {
			return 0;
		}

		const int refinement_level = this->get_refinement_level(cell);

		// given cell cannot have children
		if (refinement_level == this->max_refinement_level) {
			return cell;
		}

		const uint64_t child = get_cell_from_indices(this->get_x_index(cell), this->get_y_index(cell), this->get_z_index(cell), refinement_level + 1);
		if (this->cell_process.count(child) > 0) {
			return child;
		} else {
			return cell;
		}
	}


	/*!
	Adds new cells to cells_to_refine in order to enforce maximum refinement level difference of one between neighbours (also across processes).

	After this function cells_to_refine will contain the refines of all processes.
	*/
	void induce_refines(void)
	{
		std::vector<uint64_t> new_refines(this->cells_to_refine.begin(), this->cells_to_refine.end());
		while (all_reduce(this->comm, new_refines.size(), std::plus<uint64_t>()) > 0) {

			std::vector<std::vector<uint64_t> > all_new_refines;
			all_gather(this->comm, new_refines, all_new_refines);
			new_refines.clear();

			boost::unordered_set<uint64_t> unique_induced_refines;

			// induced refines on this process
			for (std::vector<uint64_t>::const_iterator refined = all_new_refines[this->comm.rank()].begin(); refined != all_new_refines[this->comm.rank()].end(); refined++) {

				// refine local neighbours that are too large
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*refined].begin(); neighbour != this->neighbours[*refined].end(); neighbour++) {

					#ifdef DEBUG
					if (this->cell_process.count(*neighbour) == 0) {
						std::cerr << "Process " << this->comm.rank() << ": Cell " << *refined << " had a non-existing neighbour in neighbour list: " << *neighbour << std::endl;
					}
					#endif

					if (this->cell_process[*neighbour] != this->comm.rank()) {
						continue;
					}

					if (this->get_refinement_level(*neighbour) < this->get_refinement_level(*refined)) {
						if (this->cells_to_refine.count(*neighbour) == 0) {
							unique_induced_refines.insert(*neighbour);
						}
					}
				}
				for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to[*refined].begin(); neighbour_to != this->neighbours_to[*refined].end(); neighbour_to++) {

					#ifdef DEBUG
					if (this->cell_process.count(*neighbour_to) == 0) {
						std::cerr << "Process " << this->comm.rank() << ": Cell " << *refined << " had a non-existing neighbour in neighbour list: " << *neighbour_to << std::endl;
					}
					#endif

					if (this->cell_process[*neighbour_to] != this->comm.rank()) {
						continue;
					}

					if (this->get_refinement_level(*neighbour_to) < this->get_refinement_level(*refined)) {
						if (this->cells_to_refine.count(*neighbour_to) == 0) {
							unique_induced_refines.insert(*neighbour_to);
						}
					}
				}
			}

			// refines induced here by other processes
			for (int process = 0; process < this->comm.size(); process++) {

				if (process == this->comm.rank()) {
					continue;
				}

				for (std::vector<uint64_t>::const_iterator refined = all_new_refines[process].begin(); refined != all_new_refines[process].end(); refined++) {

					if (this->remote_cells_with_local_neighbours.count(*refined) == 0) {
						continue;
					}

					// refine all local cells that are too large and neighbouring the refined cell
					// TODO: probably faster to search for local neighbours of refined cell, even faster would be to also store neighbours lists of remote cells with local neighbours
					for (boost::unordered_set<uint64_t>::const_iterator local = this->cells_with_remote_neighbours.begin(); local != this->cells_with_remote_neighbours.end(); local++) {
						if (this->is_neighbour(*local, *refined)
						&& this->get_refinement_level(*local) < this->get_refinement_level(*refined)
						&& this->cells_to_refine.count(*local) == 0) {
							unique_induced_refines.insert(*local);
						}
					}
				}
			}
			all_new_refines.clear();

			new_refines.insert(new_refines.end(), unique_induced_refines.begin(), unique_induced_refines.end());
			this->cells_to_refine.insert(unique_induced_refines.begin(), unique_induced_refines.end());
			unique_induced_refines.clear();
		}

		// reduce future global communication by adding refines from all processes to cells_to_refine
		std::vector<uint64_t> refines(this->cells_to_refine.begin(), this->cells_to_refine.end());
		std::vector<std::vector<uint64_t> > all_refines;
		all_gather(this->comm, refines, all_refines);

		for (int process = 0; process < this->comm.size(); process++) {
			this->cells_to_refine.insert(all_refines[process].begin(), all_refines[process].end());
		}

		#ifdef DEBUG
		// check that all required refines have been induced
		for (boost::unordered_set<uint64_t>::const_iterator refined = this->cells_to_refine.begin(); refined != this->cells_to_refine.end(); refined++) {

			// neighbours_of
			std::vector<uint64_t> neighbours_of = this->find_neighbours_of(*refined);

			for (std::vector<uint64_t>::const_iterator neighbour_of = neighbours_of.begin(); neighbour_of != neighbours_of.end(); neighbour_of++) {

				if (this->get_refinement_level(*neighbour_of) < this->get_refinement_level(*refined)
				&& this->cells_to_refine.count(*neighbour_of) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour (" << *neighbour_of << ") of cell that will be refined (" << *refined << ", ref lvl " << this->get_refinement_level(*refined) << ") has too small refinement level: " << this->get_refinement_level(*neighbour_of) << std::endl;
					exit(EXIT_FAILURE);
				}
			}

			// neighbours_to
			std::vector<uint64_t> neighbours_to = this->find_neighbours_to(*refined);
			for (std::vector<uint64_t>::const_iterator neighbour_to = neighbours_to.begin(); neighbour_to != neighbours_to.end(); neighbour_to++) {

				if (this->get_refinement_level(*neighbour_to) < this->get_refinement_level(*refined)
				&& this->cells_to_refine.count(*neighbour_to) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour (" << *neighbour_to << ") of cell that will be refined (" << *refined << ", ref lvl " << this->get_refinement_level(*refined) << ") has too small refinement level: " << this->get_refinement_level(*neighbour_to) << std::endl;
					exit(EXIT_FAILURE);
				}
			}
		}

		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Grid isn't consistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}


	/*!
	Removes cells from cells_to_unrefine in order to enforce maximum refinement level difference of one between neighbours (also across processes).

	cells_to_refine must contain the refines to be done by all processes.
	After this function cells_to_unrefine will contain the unrefines of all processes.
	*/
	void override_unrefines(void)
	{
		// unrefines that were not overridden by refines or too small neighbours
		boost::unordered_set<uint64_t> final_unrefines;

		for (boost::unordered_set<uint64_t>::const_iterator unrefined = this->cells_to_unrefine.begin(); unrefined != this->cells_to_unrefine.end(); unrefined++) {

			#ifdef DEBUG
			if (*unrefined == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell to unrefine" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			bool can_unrefine = true;

			// any sibling being refined will override this unrefine
			std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(*unrefined));
			for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {
				if (this->cells_to_refine.count(*sibling) > 0) {
					can_unrefine = false;
					break;
				}
			}

			// TODO improve performance by first using local neighbour lists to override unrefines
			std::vector<uint64_t> neighbours_of_parent = this->find_neighbours_of_parent(*unrefined);
			for (std::vector<uint64_t>::const_iterator neighbour_of_parent = neighbours_of_parent.begin(); neighbour_of_parent != neighbours_of_parent.end(); neighbour_of_parent++) {

				if (this->get_refinement_level(*neighbour_of_parent) < this->get_refinement_level(*unrefined)) {
					can_unrefine = false;
					break;
				}

				if (this->cells_to_refine.count(*neighbour_of_parent) > 0
				&& this->get_refinement_level(*neighbour_of_parent) <= this->get_refinement_level(*unrefined)) {
					can_unrefine = false;
					break;
				}
			}

			if (can_unrefine) {
				final_unrefines.insert(*unrefined);
			}
		}
		this->cells_to_unrefine.clear();

		// reduce future global communication by adding unrefines from all processes to cells_to_unrefine
		std::vector<uint64_t> unrefines(final_unrefines.begin(), final_unrefines.end());
		std::vector<std::vector<uint64_t> > all_unrefines;
		all_gather(this->comm, unrefines, all_unrefines);

		for (int process = 0; process < this->comm.size(); process++) {
			this->cells_to_unrefine.insert(all_unrefines[process].begin(), all_unrefines[process].end());
		}

		#ifdef DEBUG
		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Grid isn't consistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif
	}


	/*!
	Adds refined cells to the grid, removes unrefined cells from the grid.

	cells_to_refine and cells_to_unrefine must contain the cells to refine/unrefine of all processes.
	Returns new cells created on this process by refinement.
	Moves unrefined cell data to the process of their parent.
	*/
	std::vector<uint64_t> execute_refines(void)
	{
		#ifdef DEBUG
		if (!this->verify_remote_neighbour_info()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Remote neighbour info is not consistent" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " User data is inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		std::vector<uint64_t> new_cells;

		this->remote_neighbours.clear();
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
		sort(ordered_cells_to_refine.begin(), ordered_cells_to_refine.end());

		std::vector<std::vector<uint64_t> > all_ordered_cells_to_refine;
		all_gather(this->comm, ordered_cells_to_refine, all_ordered_cells_to_refine);

		for (int process = 0; process < this->comm.size(); process++) {
			if (!std::equal(all_ordered_cells_to_refine[process].begin(), all_ordered_cells_to_refine[process].end(), all_ordered_cells_to_refine[0].begin())) {
				std::cerr << __FILE__ << ":" << __LINE__ << " cells_to_refine differ between processes 0 and " << process << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// check that cells_to_unrefine is identical between processes
		std::vector<uint64_t> ordered_cells_to_unrefine(this->cells_to_unrefine.begin(), this->cells_to_unrefine.end());
		sort(ordered_cells_to_unrefine.begin(), ordered_cells_to_unrefine.end());

		std::vector<std::vector<uint64_t> > all_ordered_cells_to_unrefine;
		all_gather(this->comm, ordered_cells_to_unrefine, all_ordered_cells_to_unrefine);

		for (int process = 0; process < this->comm.size(); process++) {
			if (!std::equal(all_ordered_cells_to_unrefine[process].begin(), all_ordered_cells_to_unrefine[process].end(), all_ordered_cells_to_unrefine[0].begin())) {
				std::cerr << __FILE__ << ":" << __LINE__ << " cells_to_unrefine differ between processes 0 and " << process << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		#endif

		// cells whose neighbour lists have to be updated afterwards
		boost::unordered_set<uint64_t> update_neighbours;

		// a separate neighbourhood update function has to be used for cells whose children were removed by unrefining
		boost::unordered_set<uint64_t> update_neighbours_unrefined;

		// refines
		for (boost::unordered_set<uint64_t>::const_iterator refined = this->cells_to_refine.begin(); refined != this->cells_to_refine.end(); refined++) {

			const int process_of_refined = this->cell_process.at(*refined);

			// move user data of refined cells into refined_cell_data
			if (this->comm.rank() == process_of_refined) {
				// TODO: move data instead of copying, using boost::move or c++0x move?
				this->refined_cell_data[*refined] = this->cells.at(*refined);
				this->cells.erase(*refined);
			}

			// add children of refined cells into the grid
			std::vector<uint64_t> children = this->get_all_children(*refined);
			for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
				this->cell_process[*child] = process_of_refined;

				if (this->comm.rank() == process_of_refined) {
					this->cells[*child];
					new_cells.push_back(*child);
				}
			}

			// children of refined cells inherit their pin request status
			if (this->pin_requests.count(*refined) > 0) {
				for (auto child = children.cbegin(); child != children.cend(); child++) {
					this->pin_requests[*child] = this->pin_requests.at(*refined);
				}
				this->pin_requests.erase(*refined);
			}
			if (this->new_pin_requests.count(*refined) > 0) {
				for (auto child = children.cbegin(); child != children.cend(); child++) {
					this->new_pin_requests[*child] = this->new_pin_requests.at(*refined);
				}
				this->new_pin_requests.erase(*refined);
			}

			// use local neighbour lists to find cells whose neighbour lists have to updated
			if (process_of_refined == this->comm.rank()) {
				// update the neighbour lists of created local cells
				for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
						update_neighbours.insert(*child);

						#ifdef DEBUG
						if (this->neighbours.count(*child) > 0) {
							std::cerr << __FILE__ << ":" << __LINE__ << " Neighbours for cell " << *child << " shouldn't exist yet" << std::endl;
							exit(EXIT_FAILURE);
						}
						#endif
				}

				// update neighbour lists of all the parent's neighbours
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(*refined).begin(); neighbour != this->neighbours.at(*refined).end(); neighbour++) {
					if (this->cell_process.at(*neighbour) == this->comm.rank()) {
						update_neighbours.insert(*neighbour);
					}
				}
				for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(*refined).begin(); neighbour_to != this->neighbours_to.at(*refined).end(); neighbour_to++) {
					if (this->cell_process.at(*neighbour_to) == this->comm.rank()) {
						update_neighbours.insert(*neighbour_to);
					}
				}
			}

			// without using local neighbour lists figure out rest of the neighbour lists that need updating
			for (boost::unordered_set<uint64_t>::const_iterator with_remote_neighbour = this->cells_with_remote_neighbours.begin(); with_remote_neighbour != this->cells_with_remote_neighbours.end(); with_remote_neighbour++) {
				if (this->is_neighbour(*with_remote_neighbour, *refined)
				|| this->is_neighbour(*refined , *with_remote_neighbour)) {
					update_neighbours.insert(*with_remote_neighbour);
				}
			}
		}

		// needed for checking which neighbourhoods to update due to unrefining
		boost::unordered_set<uint64_t> parents_of_unrefined;

		// initially only one sibling is recorded per process when unrefining, insert the rest of them now
		boost::unordered_set<uint64_t> all_to_unrefine;
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = this->cells_to_unrefine.begin(); unrefined != this->cells_to_unrefine.end(); unrefined++) {

			const uint64_t parent_of_unrefined = this->get_parent(*unrefined);
			#ifdef DEBUG
			if (parent_of_unrefined == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent cell" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			parents_of_unrefined.insert(parent_of_unrefined);

			std::vector<uint64_t> siblings = this->get_all_children(parent_of_unrefined);
			all_to_unrefine.insert(siblings.begin(), siblings.end());
		}

		// unrefines
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = all_to_unrefine.begin(); unrefined != all_to_unrefine.end(); unrefined++) {

			const uint64_t parent_of_unrefined = this->get_parent(*unrefined);
			const int process_of_parent = this->cell_process.at(parent_of_unrefined);
			const int process_of_unrefined = this->cell_process.at(*unrefined);

			// remove unrefined cells and their siblings from the grid, but don't remove user data yet
			this->cell_process.erase(*unrefined);
			update_neighbours.erase(*unrefined);
			this->pin_requests.erase(*unrefined);
			this->new_pin_requests.erase(*unrefined);

			// don't send unrefined cells' user data to self
			if (this->comm.rank() == process_of_unrefined
			&& this->comm.rank() == process_of_parent) {
				// TODO move data instead of copying
				this->unrefined_cell_data[*unrefined] = this->cells.at(*unrefined);
				this->cells.erase(*unrefined);

			// send user data of removed cell to the parent's process
			} else if (this->comm.rank() == process_of_unrefined) {
				this->cells_to_send[process_of_parent].push_back(*unrefined);

			// receive user data of removed cell from its process
			} else if (this->comm.rank() == process_of_parent) {
				this->cells_to_receive[process_of_unrefined].push_back(*unrefined);
			}
		}

		this->start_user_data_transfers(
		#ifdef DCCRG_SEND_SINGLE_CELLS
		this->unrefined_cell_data
		#else
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->unrefined_cell_data
		#endif
		#endif
		);

		// update data for parents of unrefined cells
		for (boost::unordered_set<uint64_t>::const_iterator parent = parents_of_unrefined.begin(); parent != parents_of_unrefined.end(); parent++) {

			// find neighbours of unrefined cells' parents, update also those neighbourhoods
			// neighbours_of is sufficient because unrefined cells shouldn't have had a neighbour with ref. lvl. + 2
			std::vector<uint64_t> found_neighbours = find_neighbours_of(*parent);

			if (this->cell_process.at(*parent) == this->comm.rank()) {
				this->neighbours.erase(*parent);
				this->neighbours_to.erase(*parent);

				this->neighbours[*parent] = found_neighbours;
				// neighbours_to should be empty due to max refinement level difference <= 1
				this->neighbours_to[*parent];
			}

			for (std::vector<uint64_t>::const_iterator neighbour = found_neighbours.begin(); neighbour != found_neighbours.end(); neighbour++) {

				if (this->cell_process.at(*neighbour) == this->comm.rank()
				&& parents_of_unrefined.count(*neighbour) == 0) {
					update_neighbours.insert(*neighbour);
				}
			}

			// default construct user data for local parents of unrefined cells
			if (this->cell_process.at(*parent) == this->comm.rank()) {
				this->cells[*parent];
			}
		}

		// update neighbour lists of cells affected by refining / unrefining
		for (boost::unordered_set<uint64_t>::const_iterator cell = update_neighbours.begin(); cell != update_neighbours.end(); cell++) {
			this->update_neighbours(*cell);
		}

		// remove neighbour lists of added cells' parents
		for (boost::unordered_set<uint64_t>::const_iterator refined = this->cells_to_refine.begin(); refined != this->cells_to_refine.end(); refined++) {
			if (this->cell_process.at(*refined) == this->comm.rank()) {

				#ifdef DEBUG
				if (this->neighbours.count(*refined) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour list for cell " << *refined << " doesn't exist" << std::endl;
					exit(EXIT_FAILURE);
				}

				if (this->neighbours_to.count(*refined) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour_to list for cell " << *refined << " doesn't exist" << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif

				this->neighbours.erase(*refined);
				this->neighbours_to.erase(*refined);
				update_neighbours.erase(*refined);
			}
		}

		// remove neighbour lists of removed cells
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = all_to_unrefine.begin(); unrefined != all_to_unrefine.end(); unrefined++) {
			this->neighbours.erase(*unrefined);
			this->neighbours_to.erase(*unrefined);
		}

		this->update_remote_neighbour_info();

		#ifdef DEBUG
		if (!this->verify_neighbours()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour lists are inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->wait_user_data_transfer_receives(
		#ifndef DCCRG_SEND_SINGLE_CELLS
		#ifndef DCCRG_CELL_DATA_SIZE_FROM_USER
		this->unrefined_cell_data
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

		this->recalculate_neighbour_update_send_receive_lists();

		return new_cells;
	}


	/*!
	Starts user data transfers between processes based on cells_to_send and cells_to_receive.

	User data arriving to this process is saved in given destination.
	*/
	void start_user_data_transfers(
	#ifdef DCCRG_SEND_SINGLE_CELLS
	boost::unordered_map<uint64_t, UserData>& destination
	#else
	#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::unordered_map<uint64_t, UserData>& destination
	#endif
	#endif
	)
	{
		#ifdef DCCRG_USER_MPI_DATA_TYPE
		MPI_Datatype user_datatype = UserData::mpi_datatype();
		MPI_Type_commit(&user_datatype);
		#endif

		#ifdef DCCRG_SEND_SINGLE_CELLS

		// post all receives, messages are unique between different senders so just iterate over processes in random order
		for (auto sender = this->cells_to_receive.begin(); sender != this->cells_to_receive.end(); sender++) {

			#ifdef DEBUG
			if (sender->first == this->comm.rank()
			&& sender->second.size() > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Trying to transfer to self" << std::endl;
				abort();
			}
			#endif

			std::sort(sender->second.begin(), sender->second.end());
			for (auto cell = sender->second.cbegin(); cell != sender->second.cend(); cell++) {

				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->receive_requests[sender->first].push_back(MPI_Request());

				// FIXME: make sure message tags between two processes are unique

				MPI_Irecv(
					destination[*cell].at(),
					#ifdef DCCRG_USER_MPI_DATA_TYPE
					1,
					user_datatype,
					#else
					UserData::size(),
					MPI_BYTE,
					#endif
					sender->first,
					*cell % boost::mpi::environment::max_tag(),
					this->comm,
					&(this->receive_requests[sender->first].back())
				);

				#else // ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

				this->receive_requests[sender->first].push_back(
					this->comm.irecv(
						sender->first,
						*cell % boost::mpi::environment::max_tag(),
						destination[*cell]
					)
				);
				#endif
			}
		}

		// post all sends
		for (auto receiver = this->cells_to_send.begin(); receiver != this->cells_to_send.end(); receiver++) {

			#ifdef DEBUG
			if (receiver->first == this->comm.rank()
			&& receiver->second.size() > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Trying to transfer to self" << std::endl;
				abort();
			}
			#endif

			std::sort(receiver->second.begin(), receiver->second.end());
			for (auto cell = receiver->second.cbegin(); cell != receiver->second.cend(); cell++) {
				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->send_requests[receiver->first].push_back(MPI_Request());

				MPI_Isend(
					this->cells.at(*cell).at(),
					#ifdef DCCRG_USER_MPI_DATA_TYPE
					1,
					user_datatype,
					#else
					UserData::size(),
					MPI_BYTE,
					#endif
					receiver->first,
					*cell % boost::mpi::environment::max_tag(),
					this->comm,
					&(this->send_requests[receiver->first].back())
				);

				#else

				this->send_requests[receiver->first].push_back(
					this->comm.isend(
						receiver->first,
						*cell % boost::mpi::environment::max_tag(),
						this->cells.at(*cell)
					)
				);
				#endif
			}
		}

		// all user data is sent using one MPI message / process
		#else	// ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		// receive one MPI datatype per process
		for (auto sender = this->cells_to_receive.begin(); sender != this->cells_to_receive.end(); sender++) {

			std::sort(sender->second.begin(), sender->second.end());

			// reserve space for incoming user data at our end
			for (uint64_t i = 0; i < sender->second.size(); i++) {
				destination[sender->second[i]];
			}

			// get displacements in bytes for incoming user data
			std::vector<MPI_Aint> displacements(sender->second.size(), 0);
			for (uint64_t i = 0; i < sender->second.size(); i++) {
				displacements[i] = (uint8_t*) destination.at(sender->second[i]).at() - (uint8_t*) destination.at(sender->second[0]).at();
			}

			MPI_Datatype receive_datatype;

			#ifdef DCCRG_USER_MPI_DATA_TYPE
			std::vector<int> block_lengths(displacements.size(), 1);
			#else
			std::vector<int> block_lengths(displacements.size(), UserData::size());
			#endif

			MPI_Type_create_hindexed(
				displacements.size(),
				&block_lengths[0],
				&displacements[0],
				#ifdef DCCRG_USER_MPI_DATA_TYPE
				user_datatype,
				#else
				MPI_BYTE,
				#endif
				&receive_datatype
			);

			MPI_Type_commit(&receive_datatype);

			int receive_tag = sender->first * this->comm.size() + this->comm.rank();

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
		}

		// send one MPI datatype per process
		for (auto receiver = this->cells_to_send.begin(); receiver != this->cells_to_send.end(); receiver++) {

			std::sort(receiver->second.begin(), receiver->second.end());

			// get displacements in bytes for outgoing user data
			std::vector<MPI_Aint> displacements(receiver->second.size(), 0);
			for (uint64_t i = 0; i < receiver->second.size(); i++) {
				displacements[i] = (uint8_t*) this->cells.at(receiver->second[i]).at() - (uint8_t*) this->cells.at(receiver->second[0]).at();
			}

			MPI_Datatype send_datatype;

			#ifdef DCCRG_USER_MPI_DATA_TYPE
			std::vector<int> block_lengths(displacements.size(), 1);
			#else
			std::vector<int> block_lengths(displacements.size(), UserData::size());
			#endif

			MPI_Type_create_hindexed(
				displacements.size(),
				&block_lengths[0],
				&displacements[0],
				#ifdef DCCRG_USER_MPI_DATA_TYPE
				user_datatype,
				#else
				MPI_BYTE,
				#endif
				&send_datatype
			);

			MPI_Type_commit(&send_datatype);

			int send_tag = this->comm.rank() * this->comm.size() + receiver->first;

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
		}

		#ifdef DCCRG_USER_MPI_DATA_TYPE
		MPI_Type_free(&user_datatype);
		#endif

		#else	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		// post all receives
		for (int sender = 0; sender < this->comm.size(); sender++) {

			if (sender == this->comm.rank()) {
				continue;
			}

			if (this->cells_to_receive.count(sender) == 0) {
				// no data to send / receive
				continue;
			}

			int receive_tag = sender * this->comm.size() + this->comm.rank();

			this->receive_requests[sender].push_back(
				this->comm.irecv(
					sender,
					receive_tag,
					this->incoming_data[sender]
				)
			);
		}

		// gather all data to send
		for (int receiver = 0; receiver < this->comm.size(); receiver++) {

			if (receiver == this->comm.rank()) {
				// don't send to self
				continue;
			}

			if (this->cells_to_send.count(receiver) == 0) {
				// no data to send / receive
				continue;
			}

			std::sort(this->cells_to_send.at(receiver).begin(), this->cells_to_send.at(receiver).end());
			// construct the outgoing data vector
			for (auto cell = this->cells_to_send[receiver].cbegin(); cell != this->cells_to_send[receiver].cend(); cell++) {
				UserData* user_data = (*this)[*cell];
				assert(user_data != NULL);
				this->outgoing_data[receiver].push_back(*user_data);
			}
			assert(this->outgoing_data[receiver].size() == this->cells_to_send[receiver].size());
		}

		// post all sends
		for (int receiver = 0; receiver < this->comm.size(); receiver++) {

			if (receiver == this->comm.rank()) {
				continue;
			}

			if (this->cells_to_send.count(receiver) == 0) {
				// no data to send / receive
				continue;
			}

			int send_tag = this->comm.rank() * this->comm.size() + receiver;

			this->send_requests[receiver].push_back(
				this->comm.isend(
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
	boost::unordered_map<uint64_t, UserData>& destination
	#endif
	#endif
	)
	{
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (auto process = this->receive_requests.begin(); process != this->receive_requests.end(); process++) {

			std::vector<MPI_Status> statuses;
			statuses.resize(process->second.size());

			if (MPI_Waitall(process->second.size(), &(process->second[0]), &(statuses[0])) != MPI_SUCCESS) {
				for (std::vector<MPI_Status>::const_iterator status = statuses.begin(); status != statuses.end(); status++) {
					if (status->MPI_ERROR != MPI_SUCCESS) {
						std::cerr << "MPI receive failed from process " << status->MPI_SOURCE << " with tag " << status->MPI_TAG << std::endl;
					}
				}
			}
		}

		#else	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (auto process = this->receive_requests.begin(); process != this->receive_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}

		#ifndef DCCRG_SEND_SINGLE_CELLS

		// incorporate received data
		for (auto sender = this->incoming_data.cbegin(); sender != this->incoming_data.cend(); sender++) {

			assert(this->incoming_data.at(sender->first).size() == this->cells_to_receive.at(sender->first).size());

			std::sort(this->cells_to_receive.at(sender->first).begin(), this->cells_to_receive.at(sender->first).end());

			int i = 0;
			for (
				auto cell = this->cells_to_receive.at(sender->first).cbegin();
				cell != this->cells_to_receive.at(sender->first).cend();
				cell++, i++
			) {
				// TODO move data instead of copying
				destination[*cell] = this->incoming_data.at(sender->first)[i];
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
	void wait_user_data_transfer_sends(void)
	{
		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (auto process = this->send_requests.begin(); process != this->send_requests.end(); process++) {

			std::vector<MPI_Status> statuses;
			statuses.resize(process->second.size());

			if (MPI_Waitall(process->second.size(), &(process->second[0]), &(statuses[0])) != MPI_SUCCESS) {
				for (std::vector<MPI_Status>::const_iterator status = statuses.begin(); status != statuses.end(); status++) {
					if (status->MPI_ERROR != MPI_SUCCESS) {
						std::cerr << "MPI receive failed from process " << status->MPI_SOURCE << " with tag " << status->MPI_TAG << std::endl;
					}
				}
			}
		}

		#else	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (auto process = this->send_requests.begin(); process != this->send_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}

		#ifndef DCCRG_SEND_SINGLE_CELLS

		this->outgoing_data.clear();

		#endif	// ifndef DCCRG_SEND_SINGLE_CELLS
		#endif	// ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		this->send_requests.clear();
	}


	/*!
	Returns the smallest existing cell at the given coordinate
	Returns 0 if the coordinate is outside of the grid or the cell is on another process
	*/
	uint64_t get_smallest_cell_from_coordinate(const double x, const double y, const double z) const
	{
		#ifdef DCCRG_ARBITRARY_STRETCH
		if (x < this->x_coordinates[0]
			|| x > this->x_coordinates[this->geometry.get_x_length()]
			|| y < this->y_coordinates[0]
			|| y > this->y_coordinates[this->geometry.get_y_length()]
			|| z < this->z_coordinates[0]
			|| z > this->z_coordinates[this->geometry.get_z_length()]) {
		#else
		if (x < this->geometry.get_x_start()
			|| x > this->geometry.get_x_start() + this->cell_size * this->geometry.get_x_length()
			|| y < this->geometry.get_y_start()
			|| y > this->geometry.get_y_start() + this->cell_size * this->geometry.get_y_length()
			|| z < this->geometry.get_z_start()
			|| z > this->geometry.get_z_start() + this->cell_size * this->geometry.get_z_length()) {
		#endif
			return 0;
		}

		return this->get_cell_from_indices(this->get_x_index(x), this->get_y_index(y), this->get_z_index(z), 0, this->max_refinement_level);
	}


	/*!
	These return the x, y or z index of the given coordinate
	*/
	uint64_t get_x_index(const double x) const
	{
		#ifdef DCCRG_ARBITRARY_STRETCH
		assert((x >= this->x_coordinates[0]) && (x <= this->x_coordinates[this->geometry.get_x_length()]));

		uint64_t x_coord_start_index = 0;
		while (x_coordinates[x_coord_start_index] < x) {
			x_coord_start_index++;
		}
		x_coord_start_index--;

		double x_size_of_index = (this->x_coordinates[x_coord_start_index + 1] - this->x_coordinates[x_coord_start_index]) / this->get_cell_size_in_indices(1);

		uint64_t index_offset = 0;
		while (x_coordinates[x_coord_start_index] + index_offset * x_size_of_index < x) {
			index_offset++;
		}
		index_offset--;

		return x_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;
		#else
		assert((x >= this->geometry.get_x_start()) and (x <= this->geometry.get_x_start() + this->geometry.get_x_length() * this->cell_size));
		return uint64_t((x - this->geometry.get_x_start()) / (this->cell_size / (uint64_t(1) << this->max_refinement_level)));
		#endif
	}

	uint64_t get_y_index(const double y) const
	{
		#ifdef DCCRG_ARBITRARY_STRETCH
		assert((y >= this->y_coordinates[0]) && (y <= this->y_coordinates[this->geometry.get_y_length()]));

		uint64_t y_coord_start_index = 0;
		while (y_coordinates[y_coord_start_index] < y) {
			y_coord_start_index++;
		}
		y_coord_start_index--;

		double y_size_of_index = (this->y_coordinates[y_coord_start_index + 1] - this->y_coordinates[y_coord_start_index]) / this->get_cell_size_in_indices(1);

		uint64_t index_offset = 0;
		while (y_coordinates[y_coord_start_index] + index_offset * y_size_of_index < y) {
			index_offset++;
		}
		index_offset--;

		return y_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;
		#else
		assert((y >= this->geometry.get_y_start()) and (y <= this->geometry.get_y_start() + this->geometry.get_y_length() * this->cell_size));
		return uint64_t((y - this->geometry.get_y_start()) / (this->cell_size / (uint64_t(1) << this->max_refinement_level)));
		#endif
	}

	uint64_t get_z_index(const double z) const
	{
		#ifdef DCCRG_ARBITRARY_STRETCH
		assert((z >= this->z_coordinates[0]) && (z <= this->z_coordinates[this->geometry.get_z_length()]));

		uint64_t z_coord_start_index = 0;
		while (z_coordinates[z_coord_start_index] < z) {
			z_coord_start_index++;
		}
		z_coord_start_index--;

		double z_size_of_index = (this->z_coordinates[z_coord_start_index + 1] - this->z_coordinates[z_coord_start_index]) / this->get_cell_size_in_indices(1);

		uint64_t index_offset = 0;
		while (z_coordinates[z_coord_start_index] + index_offset * z_size_of_index < z) {
			index_offset++;
		}
		index_offset--;

		return z_coord_start_index * this->get_cell_size_in_indices(1) + index_offset;
		#else
		assert((z >= this->geometry.get_z_start()) and (z <= this->geometry.get_z_start() + this->geometry.get_z_length() * this->cell_size));
		return uint64_t((z - this->geometry.get_z_start()) / (this->cell_size / (uint64_t(1) << this->max_refinement_level)));
		#endif
	}


	/*!
	Returns true if cells with given index properties overlap.

	Sizes are also given in indices.
	*/
	bool indices_overlap(const uint64_t index1, const uint64_t size1, const uint64_t index2, const uint64_t size2) const
	{
		#ifdef DEBUG
		if (index1 >= this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level)
		&& index1 >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)
		&& index1 >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid index given" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (index2 >= this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level)
		&& index2 >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)
		&& index2 >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
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
	Returns true if x indices of given cells overlap, even if they don't exist
	*/
	bool x_indices_overlap(const uint64_t cell1, const uint64_t cell2) const
	{
		assert(cell1 > 0);
		assert(cell1 <= this->max_cell_number);
		assert(cell2 > 0);
		assert(cell2 <= this->max_cell_number);

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
		assert(cell1 <= this->max_cell_number);
		assert(cell2 > 0);
		assert(cell2 <= this->max_cell_number);

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
		assert(cell1 <= this->max_cell_number);
		assert(cell2 > 0);
		assert(cell2 <= this->max_cell_number);

		const uint64_t index1 = this->get_z_index(cell1);
		const uint64_t index2 = this->get_z_index(cell2);
		const uint64_t size1 = this->get_cell_size_in_indices(cell1);
		const uint64_t size2 = this->get_cell_size_in_indices(cell2);

		return this->indices_overlap(index1, size1, index2, size2);
	}


	/*!
	Returns the number of directions in which given cells' indices overlap
	Returns 0 if even one of given cell's doesn't exist
	*/
	int overlapping_indices(const uint64_t cell1, const uint64_t cell2) const
	{
		if (this->cell_process.count(cell1) == 0 || this->cell_process.count(cell2) == 0) {
			return 0;
		}

		int ret = 0;
		if (this->x_indices_overlap(cell1, cell2)) {
			ret++;
		}
		if (this->y_indices_overlap(cell1, cell2)) {
			ret++;
		}
		if (this->z_indices_overlap(cell1, cell2)) {
			ret++;
		}
		return ret;
	}



	/*!
	Returns the smallest existing cell at given indices between given refinement levels inclusive.

	Returns 0 if no cell between given refinement ranges exists or an index is outside of the grid or minimum_refinement_level > maximum_refinement_level.
	*/
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int minimum_refinement_level, const int maximum_refinement_level) const
	{
		if (x_index >= this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level)) {
			#ifdef DEBUG
			std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell for x index " << x_index << std::endl;
			#endif
			return 0;
		}

		if (y_index >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)) {
			#ifdef DEBUG
			std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell for y index " << y_index << std::endl;
			#endif
			return 0;
		}

		if (z_index >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
			#ifdef DEBUG
			std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell for z index " << z_index << std::endl;
			#endif
			return 0;
		}

		if (minimum_refinement_level > maximum_refinement_level) {
			#ifdef DEBUG
			std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell because of refinement levels: " << minimum_refinement_level << " > " << maximum_refinement_level << std::endl;
			#endif
			return 0;
		}

		int average_refinement_level = (maximum_refinement_level + minimum_refinement_level) / 2;
		uint64_t id = this->get_cell_from_indices(x_index, y_index, z_index, average_refinement_level);

		// use binary search recursively (assumes that a cells refine to 8 children)
		if (this->cell_process.count(id) == 0) {
			// doesn't exist, search the bin of smaller refinement_level values
			if (average_refinement_level > minimum_refinement_level) {

				uint64_t smaller_refinement_value_cell = this->get_cell_from_indices(x_index, y_index, z_index, minimum_refinement_level, average_refinement_level - 1);

				#ifdef DEBUG
				if (this->cell_process.count(smaller_refinement_value_cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell: " << smaller_refinement_value_cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif

				return smaller_refinement_value_cell;

			} else {

				// nothing left to search
				return 0;

			}
		} else {
			// does exist, search the bin of larger refinement_level values
			if (average_refinement_level < maximum_refinement_level) {
				uint64_t larger_refinement_value_cell = this->get_cell_from_indices(x_index, y_index, z_index, average_refinement_level + 1, maximum_refinement_level);

				if (larger_refinement_value_cell > 0) {

					#ifdef DEBUG
					if (this->cell_process.count(larger_refinement_value_cell) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell: " << larger_refinement_value_cell << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif

					return larger_refinement_value_cell;

				} else {

					// current cell has the largest refinement value at given indices

					#ifdef DEBUG
					if (this->cell_process.count(id) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell: " << id << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif

					return id;
				}
			} else {
				// nothing left to search

				#ifdef DEBUG
				if (this->cell_process.count(id) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell: " << id << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif

				return id;
			}
		}
	}


	/*!
	Returns all children of given cell regardless of whether they exist
	Returns no cells if childrens' refinement level would exceed max_refinement_level
	 */
	std::vector<uint64_t> get_all_children(const uint64_t id) const
	{
		#ifdef DEBUG
		if (id == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (this->cell_process.count(id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Given cell doesn't exist: " << id << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		std::vector<uint64_t> children;
		children.reserve(8);

		// given cell cannot have children
		int refinement_level = this->get_refinement_level(id);
		if (refinement_level >= this->max_refinement_level) {
			return children;
		}

		uint64_t x_index = this->get_x_index(id);
		uint64_t y_index = this->get_y_index(id);
		uint64_t z_index = this->get_z_index(id);

		// get indices of next refinement level within this cell
		refinement_level++;
		uint64_t index_offset = (uint64_t(1) << (max_refinement_level - refinement_level));
		for (uint64_t x_index_offset = 0; x_index_offset < 2 * index_offset; x_index_offset += index_offset) {
			for (uint64_t y_index_offset = 0; y_index_offset < 2 * index_offset; y_index_offset += index_offset) {
				for (uint64_t z_index_offset = 0; z_index_offset < 2 * index_offset; z_index_offset += index_offset) {
					children.push_back(get_cell_from_indices(x_index + x_index_offset, y_index + y_index_offset, z_index + z_index_offset, refinement_level));
				}
			}
		}

		return children;
	}


	/*!
	Returns the number of values needed to represent the coordinate of a cell
	*/
	static int get_grid_dimensionality(void* /*data*/, int* error)
	{
		*error = ZOLTAN_OK;
		return 3;
	}


	// TODO: Zoltan assumes global ids are integers, which works as long as there are less than 2^32 cells in the grid

	/*!
	Fills geom_vec with the coordinates of cells given in global_id
	*/
	static void fill_with_cell_coordinates(void *data, int /*global_id_size*/, int /*local_id_size*/, int number_of_cells, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int /*number_of_dimensions*/, double *geom_vec, int *error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		for (int i = 0; i < number_of_cells; i++) {
			uint64_t cell = uint64_t(global_ids[i]);
			if (dccrg_instance->cells.count(cell) == 0) {
				*error = ZOLTAN_FATAL;
				std::cerr << "Process " << dccrg_instance->comm.rank() << ": Zoltan wanted the coordinates of a non-existing cell " << cell << std::endl;
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
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;
		return dccrg_instance->cells.size();
	}


	/*!
	Writes all cell ids on this process to the global_ids array
	*/
	static void fill_cell_list(void* data, int /*global_id_size*/, int /*local_id_size*/, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int number_of_weights_per_object, float* object_weights, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		int i = 0;
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = dccrg_instance->cells.begin(); cell != dccrg_instance->cells.end(); cell++, i++) {

			#ifdef DEBUG
			if (cell->first == 0) {
				std::cerr << "User data exist for an illegal cell" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif

			global_ids[i] = cell->first;

			if (number_of_weights_per_object > 0) {
				object_weights[i] = 1e-10;
			}
		}
	}


	/*!
	Writes the number of neighbours into number_of_neighbours for all cells given in global_ids with length number_of_cells
	*/
	static void fill_number_of_neighbours_for_cells(void* data, int /*global_id_size*/, int /*local_id_size*/, int number_of_cells, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int* number_of_neighbours, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		for (int i = 0; i < number_of_cells; i++) {
			uint64_t cell = uint64_t(global_ids[i]);
			if (dccrg_instance->cells.count(cell) == 0) {
				*error = ZOLTAN_FATAL;
				std::cerr << "Process " << dccrg_instance->comm.rank() << ": Zoltan wanted the number of neighbours of a non-existing cell " << cell << std::endl;
				return;
			}

			number_of_neighbours[i] = dccrg_instance->neighbours[cell].size();
		}
	}


	/*!
	Writes neighbour lists of given cells into neighbours, etc.
	*/
	static void fill_neighbour_lists(void* data, int /*global_id_size*/, int /*local_id_size*/, int number_of_cells, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int* number_of_neighbours, ZOLTAN_ID_PTR neighbours, int* processes_of_neighbours, int number_of_weights_per_edge, float* edge_weights, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		int current_neighbour_number = 0;
		for (int i = 0; i < number_of_cells; i++) {
			uint64_t cell = uint64_t(global_ids[i]);
			if (dccrg_instance->cells.count(cell) == 0) {
				*error = ZOLTAN_FATAL;
				std::cerr << "Process " << dccrg_instance->comm.rank() << ": Zoltan wanted neighbour list of a non-existing cell " << cell << std::endl;
				return;
			}

			number_of_neighbours[i] = dccrg_instance->neighbours[cell].size();

			for (std::vector<uint64_t>::const_iterator neighbour = dccrg_instance->neighbours[cell].begin(); neighbour != dccrg_instance->neighbours[cell].end(); neighbour++, current_neighbour_number++) {

				neighbours[current_neighbour_number] = *neighbour;
				processes_of_neighbours[current_neighbour_number] = dccrg_instance->cell_process[*neighbour];

				// weight of edge from cell to *neighbour
				if (number_of_weights_per_edge > 0) {
					edge_weights[current_neighbour_number] = 1.0;
				}
			}
		}
	}


	/*!
	Writes the number of hyperedges (one per existing cell) in the grid on this process.
	*/
	static void fill_number_of_hyperedges(void* data, int* number_of_hyperedges, int* number_of_connections, int* format, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		*number_of_hyperedges = dccrg_instance->cells.size();
		*format = ZOLTAN_COMPRESSED_EDGE;

		*number_of_connections = 0;
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = dccrg_instance->cells.begin(); cell != dccrg_instance->cells.end(); cell++) {
			*number_of_connections += 1 + dccrg_instance->neighbours[cell->first].size();
		}
	}


	/*!
	Writes the hypergraph in compressed edge format
	*/
	static void fill_hyperedge_lists(void* data, int /*global_id_size*/, int number_of_hyperedges, int number_of_connections, int format, ZOLTAN_ID_PTR hyperedges, int* hyperedge_connection_offsets, ZOLTAN_ID_PTR connections, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		if (format != ZOLTAN_COMPRESSED_EDGE) {
			std::cerr << "Only compressed edge format supported for hypergraph partitioning" << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		if ((unsigned int) number_of_hyperedges != dccrg_instance->cells.size()) {
			std::cerr << "Zoltan is expecting wrong number of hyperedges: " << number_of_hyperedges << " instead of " << dccrg_instance->cells.size() << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		int i = 0;
		int connection_number = 0;
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = dccrg_instance->cells.begin(); cell != dccrg_instance->cells.end(); cell++, i++) {

			hyperedges[i] = cell->first;
			hyperedge_connection_offsets[i] = connection_number;

			// add a connection to the cell itself from its hyperedge
			connections[connection_number++] = cell->first;

			for (std::vector<uint64_t>::const_iterator neighbour = dccrg_instance->neighbours[cell->first].begin(); neighbour != dccrg_instance->neighbours[cell->first].end(); neighbour++, connection_number++) {
				connections[connection_number] = *neighbour;
			}
		}

		if (connection_number != number_of_connections) {
			std::cerr << "Zoltan is expecting wrong number of connections from hyperedges: " << number_of_connections << " instead of " << connection_number << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}
	}


	/*!
	Writes the number of hyperedge weights (one per hyperedge) on this process
	*/
	static void fill_number_of_edge_weights(void* data, int* number_of_edge_weights, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		*number_of_edge_weights = dccrg_instance->cells.size();
		return;
	}


	/*!
	Writes hyperedge weights (one per hyperedge) on this process
	*/
	static void fill_edge_weights(void* data, int /*global_id_size*/, int /*local_id_size*/, int number_of_hyperedges, int number_of_weights_per_hyperedge, ZOLTAN_ID_PTR hyperedges, ZOLTAN_ID_PTR /*hyperedges_local_ids*/, float* hyperedge_weights, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;

		if ((unsigned int) number_of_hyperedges != dccrg_instance->cells.size()) {
			std::cerr << "Zoltan is expecting wrong number of hyperedges: " << number_of_hyperedges << " instead of " << dccrg_instance->cells.size() << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		int i = 0;
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = dccrg_instance->cells.begin(); cell != dccrg_instance->cells.end(); cell++, i++) {
			hyperedges[i] = cell->first;

			if (number_of_weights_per_hyperedge > 0) {
				hyperedge_weights[i] = 1.0 * dccrg_instance->neighbours[cell->first].size();
			}
		}
	}


	/*!
	Returns the number of hierarchies to use for load balancing.
	*/
	static int get_number_of_load_balancing_hierarchies(void* data, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;
		return dccrg_instance->processes_per_part.size();
	}


	/*!
	Returns the part number of this process on given hierarchy level of load balancing.
	*/
	static int get_part_number(void* data, int level, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);

		if (level < 0 || level >= int(dccrg_instance->processes_per_part.size())) {
			std::cerr << "Zoltan wanted a part number for an invalid hierarchy level (level should be between 0 and " << dccrg_instance->processes_per_part.size() - 1 << " inclusive): " << level << std::endl;
			*error = ZOLTAN_FATAL;
			return -1;
		} else {
			*error = ZOLTAN_OK;
		}

		int process = dccrg_instance->comm.rank();
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

		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);

		if (level < 0 || level >= int(dccrg_instance->processes_per_part.size())) {
			std::cerr << "Zoltan wanted partitioning options for an invalid hierarchy level (level should be between 0 and " << dccrg_instance->processes_per_part.size() - 1 << " inclusive): " << level << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		} else {
			*error = ZOLTAN_OK;
		}

		for (boost::unordered_map<std::string, std::string>::const_iterator option = dccrg_instance->partitioning_options[level].begin(); option != dccrg_instance->partitioning_options[level].end(); option++) {
			Zoltan_Set_Param(zz, option->first.c_str(), option->second.c_str());
		}
	}



	#ifdef DEBUG
	/*!
	Returns false if the same cells don't exist on the same process for all processes.
	*/
	bool is_consistent(void)
	{
		// sort existing cells from this process
		std::vector<uint64_t> local_cells;
		local_cells.reserve(this->cell_process.size());
		for (typename boost::unordered_map<uint64_t, int>::const_iterator cell = this->cell_process.begin(); cell != this->cell_process.end(); cell++) {
			local_cells.push_back(cell->first);
		}
		sort(local_cells.begin(), local_cells.end());

		// processes of existing cells from this process
		std::vector<int> local_processes;
		local_processes.reserve(this->cell_process.size());
		for (std::vector<uint64_t>::const_iterator cell = local_cells.begin(); cell != local_cells.end(); cell++) {
			local_processes.push_back(this->cell_process[*cell]);
		}

		// compare the above between processes
		std::vector<std::vector<uint64_t> > all_cells;
		all_gather(this->comm, local_cells, all_cells);
		std::vector<std::vector<int> > all_processes;
		all_gather(this->comm, local_processes, all_processes);

		for (int process = 0; process < this->comm.size(); process++) {
			if (!std::equal(all_cells[process].begin(), all_cells[process].end(), all_cells[0].begin())) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Grid has different cells between processes 0 and " << process << std::endl;
				return false;
			}

			if (!std::equal(all_processes[process].begin(), all_processes[process].end(), all_processes[0].begin())) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Grid's cells have different processes between processes 0 and " << process << std::endl;
				return false;
			}
		}

		return true;
	}


	/*!
	Return false if neighbours lists of the given cell aren't consistent
	*/
	bool verify_neighbours(const uint64_t cell)
	{
		if (cell == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid cell given" << std::endl;
			return false;
		}

		if (cell > this->max_cell_number) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " shouldn't exist" << std::endl;
			return false;
		}

		if (this->cell_process.count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << cell << " doesn't exist" << std::endl;
			return false;
		}

		if (cell == this->get_child(cell)) {

			if (this->neighbours.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " No neighbour list for cell " << cell << std::endl;
				return false;
			}

			if (this->neighbours_to.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " No neighbour_to list for cell " << cell << std::endl;
				return false;
			}

		} else {

			if (this->neighbours.count(cell) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour list for cell " << cell << " shouldn't exist" << std::endl;
				return false;
			}

			if (this->neighbours_to.count(cell) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour_to list for cell " << cell << " shouldn't exist" << std::endl;
				return false;
			}

		}

		if (cell != this->get_child(cell)) {
			return true;
		}

		// neighbours
		sort(this->neighbours.at(cell).begin(), this->neighbours.at(cell).end());
		std::vector<uint64_t> compare_neighbours = this->find_neighbours_of(cell);
		sort(compare_neighbours.begin(), compare_neighbours.end());

		if ((
			this->neighbours.at(cell).size() != compare_neighbours.size())
		|| (
			this->neighbours.at(cell).size() > 0
			&& compare_neighbours.size() > 0
			&& !std::equal(this->neighbours.at(cell).begin(), this->neighbours.at(cell).end(), compare_neighbours.begin()))
		) {
			std::cerr << "Process " << this->comm.rank() << " neighbour counts for cell " << cell << " (child of " << this->get_parent(cell) << ") don't match " << this->neighbours.at(cell).size() << ": ";
			for (std::vector<uint64_t>::const_iterator c = this->neighbours.at(cell).begin(); c != this->neighbours.at(cell).end(); c++) {
				std::cerr << *c << " ";
			}
			std::cerr << ", should be (+ child of) " << compare_neighbours.size() << ": ";
			for (std::vector<uint64_t>::const_iterator c = compare_neighbours.begin(); c != compare_neighbours.end(); c++) {
				std::cerr << *c << "(" << this->get_parent(*c) << ") ";
			}
			std::cerr << std::endl;
			return false;
		}

		// neighbours_to
		if (this->neighbours_to.at(cell).size() > 0) {
			sort(this->neighbours_to.at(cell).begin(), this->neighbours_to.at(cell).end());
		}
		std::vector<uint64_t> compare_neighbours_to = this->find_neighbours_to(cell);
		if (compare_neighbours_to.size() > 0) {
			sort(compare_neighbours_to.begin(), compare_neighbours_to.end());
		}

		if (!std::equal(this->neighbours_to.at(cell).begin(), this->neighbours_to.at(cell).end(), compare_neighbours_to.begin())) {
			std::cerr << "Process " << this->comm.rank() << " neighbour_to counts for cell " << cell << " (child of " << this->get_parent(cell) << ") don't match: " << this->neighbours_to.at(cell).size() << " (";
			for (std::vector<uint64_t>::const_iterator c = this->neighbours_to.at(cell).begin(); c != this->neighbours_to.at(cell).end(); c++) {
				std::cerr << *c;
				if (*c != this->get_child(*c)) {
					std::cerr << " [has a child " << this->get_child(*c) << "], ";
				} else {
					std::cerr << ", ";
				}
			}
			std::cerr << ") should be " << compare_neighbours_to.size() << " (";
			for (std::vector<uint64_t>::const_iterator c = compare_neighbours_to.begin(); c != compare_neighbours_to.end(); c++) {
				std::cerr << *c << ", ";
			}
			std::cerr << ")" << std::endl;
			return false;
		}

		return true;
	}


	/*!
	Returns false if neighbour lists on this process aren't consistent
	*/
	bool verify_neighbours(void)
	{
		for (typename boost::unordered_map<uint64_t, int>::const_iterator cell = this->cell_process.begin(); cell != this->cell_process.end(); cell++) {

			if (cell->second != this->comm.rank()) {
				continue;
			}

			if (!this->verify_neighbours(cell->first)) {
				return false;
			}
		}

		return true;
	}


	/*!
	Returns false if remote neighbour info for given cell is inconsistent.

	Remote neighbour info consists of cells_with_remote_neighbours and remote_cells_with_local_neighbours.
	*/
	bool verify_remote_neighbour_info(const uint64_t cell)
	{
		if (!this->verify_neighbours(cell)) {
			return false;
		}

		if (cell != this->get_child(cell)) {
			return true;
		}

		std::vector<uint64_t> all_neighbours(this->neighbours.at(cell).cbegin(), this->neighbours.at(cell).cend());
		all_neighbours.insert(all_neighbours.end(), this->neighbours_to.at(cell).cbegin(), this->neighbours_to.at(cell).cend());

		for (auto neighbour = all_neighbours.cbegin(); neighbour != all_neighbours.cend(); neighbour++) {
			if (this->cell_process.at(*neighbour) != this->comm.rank()) {

				if (this->cells_with_remote_neighbours.count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Local cell " << cell << " should be in cells_with_remote_neighbours" << std::endl;
					return false;
				}

				if (this->remote_cells_with_local_neighbours.count(*neighbour) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Remote cell " << *neighbour << " should be in remote_cells_with_local_neighbours" << std::endl;
					return false;
				}
			}
		}

		return true;
	}


	/*!
	Returns false if remote neighbour info on this process is inconsistent.

	Remote neighbour info consists of cells_with_remote_neighbours and remote_cells_with_local_neighbours.
	*/
	bool verify_remote_neighbour_info(void)
	{
		for (boost::unordered_map<uint64_t, int>::const_iterator item = this->cell_process.begin(); item != this->cell_process.end(); item++) {

			if (item->first != this->get_child(item->first)) {
				continue;
			}

			// check whether this cell should be in remote_cells_with_local_neighbours
			if (item->second != this->comm.rank()) {

				bool should_be_in_remote_cells = false;

				for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

					if (cell->first != this->get_child(cell->first)) {
						continue;
					}

					if (this->is_neighbour(item->first, cell->first)
					|| this->is_neighbour(cell->first, item->first)) {
						should_be_in_remote_cells = true;
					}
				}

				if (should_be_in_remote_cells) {

					if (this->remote_cells_with_local_neighbours.count(item->first) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Remote cell " << item->first << " should be in remote_cells_with_local_neighbours because:" << std::endl;

						for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
							if (this->is_neighbour(item->first, cell->first)
							|| this->is_neighbour(cell->first, item->first)) {
								std::cout << "\tremote cell " << item->first << " has a local neighbour " << cell->first << std::endl;
							}
						}
						return false;
					}

				} else {

					if (this->remote_cells_with_local_neighbours.count(item->first) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Remote cell " << item->first << " should not be in remote_cells_with_local_neighbours" << std::endl;
						return false;
					}
				}

			// check whether this cell should be in cells_with_remote_neighbour
			} else {

				bool no_remote_neighbour = true;

				// search in neighbours_of
				std::vector<uint64_t> neighbours_of = this->find_neighbours_of(item->first);
				for (std::vector<uint64_t>::const_iterator neighbour = neighbours_of.begin(); neighbour != neighbours_of.end(); neighbour++) {
					if (this->cell_process.at(*neighbour) != this->comm.rank()) {
						no_remote_neighbour = false;
					}

					if (!this->is_neighbour(item->first, *neighbour)) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << *neighbour << " should be a neighbour of cell " << item->first << std::endl;
						exit(EXIT_FAILURE);
					}
				}

				// search in neighbours_to
				std::vector<uint64_t> neighbours_to = this->find_neighbours_to(item->first);
				for (std::vector<uint64_t>::const_iterator neighbour = neighbours_to.begin(); neighbour != neighbours_to.end(); neighbour++) {
					if (this->cell_process.at(*neighbour) != this->comm.rank()) {
						no_remote_neighbour = false;
					}

					if (!this->is_neighbour(*neighbour, item->first)) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << item->first << " should be a neighbour of cell " << *neighbour << std::endl;
						exit(EXIT_FAILURE);
					}
				}

				if (no_remote_neighbour) {
					if (this->cells_with_remote_neighbours.count(item->first) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Local cell " << item->first << " should not be in cells_with_remote_neighbours" << std::endl;
						return false;
					}
				} else {
					if (this->cells_with_remote_neighbours.count(item->first) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Local cell " << item->first << " should be in cells_with_remote_neighbours" << std::endl;
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
	bool verify_user_data(void)
	{
		for (auto item = this->cell_process.cbegin(); item != this->cell_process.cend(); item++) {
			if (item->second == this->comm.rank()
			&& item->first == this->get_child(item->first)
			&& this->cells.count(item->first) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " User data for local cell " << item->first << " should exist" << std::endl;
				return false;
			}
			if (item->second != this->comm.rank()
			&& this->cells.count(item->first) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " User data for local cell " << item->first << " shouldn't exist" << std::endl;
				return false;
			}
		}

		return true;
	}


	/*!
	Returns true if all cells are where pin reqests should have placed them.
	*/
	bool pin_requests_succeeded(void)
	{
		for (auto pin_request = this->pin_requests.cbegin(); pin_request != this->pin_requests.cend(); pin_request++) {
			if (this->cell_process.at(pin_request->first) != pin_request->second) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Cell " << pin_request->first << " not at requested process " << pin_request->second << " but at " << this->cell_process.at(pin_request->first) << std::endl;
				return false;
			}
		}

		return true;
	}
	#endif

};

#endif
