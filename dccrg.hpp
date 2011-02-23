/*
A distributed cartesian cell-refinable grid

Copyright 2009, 2010, 2011 Ilja Honkonen

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
Cells of arbitrary size in x, y and z directions can be created by defining DCCRG_ARBITRARY_STRETCH
DCCRG_CONSTANT_STRETCH is not supported at the moment
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
	-At the moment neighbour data updates are supported only one cell at a time, hence DCCRG_SEND_SINGLE_CELLS must also be defined.
*/
#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
	#ifndef DCCRG_SEND_SINGLE_CELLS
		#error Using cells with user defined size of data requires that DCCRG_SEND_SINGLE_CELLS also be defined
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

		// set reserved options
		Zoltan_Set_Param(this->zoltan, "EDGE_WEIGHT_DIM", "0");	// 0 because Zoltan crashes with larger values
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
		this->geometry.set_coordinates(x_coordinates, y_coordinates, z_coordinates);
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
		#ifndef NDEBUG
		if (!this->verify_neighbours()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour lists are inconsistent" << std::endl;
			// TODO: throw an exception instead when debugging?
			exit(EXIT_FAILURE);
		}
		#endif

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
			this->update_remote_neighbour_info(cell->first);
		}

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

			uint64_t child = this->get_child(cell->first);
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
	Load balances the grid's cells among processes.

	Must be called simultaneously on all processes.	Does not update remote neighbour data between processes. Discards refines / unrefines.
	*/
	void balance_load(void)
	{
		this->comm.barrier();

		int partition_changed, global_id_size, local_id_size, number_to_receive, number_to_send;
		ZOLTAN_ID_PTR global_ids_to_receive, local_ids_to_receive, global_ids_to_send, local_ids_to_send;
		int *sender_processes, *receiver_processes;

		// here these record where cells have migrated
		boost::unordered_set<uint64_t> added_cells, removed_cells;

		if (Zoltan_LB_Balance(this->zoltan, &partition_changed, &global_id_size, &local_id_size, &number_to_receive, &global_ids_to_receive, &local_ids_to_receive, &sender_processes, &number_to_send, &global_ids_to_send, &local_ids_to_send, &receiver_processes) != ZOLTAN_OK) {
			if (!this->no_load_balancing) {
				std::cerr << "Zoltan_LB_Partition failed" << std::endl;
				Zoltan_Destroy(&this->zoltan);
				// TODO: throw an exception instead
				exit(EXIT_FAILURE);
			}
		}

		if (partition_changed == 0) {
			return;
		}

		// clear user data which is about to get old
		this->cells_with_remote_neighbours.clear();
		this->remote_neighbours.clear();
		this->cells_to_refine.clear();
		this->refined_cell_data.clear();
		this->cells_to_unrefine.clear();
		this->unrefined_cell_data.clear();
		this->remote_cells_with_local_neighbours.clear();

		// clear send / receive lists, here they mean cells that will be moved between processes
		this->cells_to_receive.clear();
		this->cells_to_send.clear();

		// TODO: move identical code from here and start_neighbour_data_update into a separate function

		// processes and the cells for which data has to be received by this process
		for (int i = 0; i < number_to_receive; i++) {
			assert(this->cell_process[global_ids_to_receive[i]] == sender_processes[i]);

			this->cells_to_receive[sender_processes[i]].push_back(global_ids_to_receive[i]);
			added_cells.insert(global_ids_to_receive[i]);
		}

		// post all receives
		for (int sender = 0; sender < this->comm.size(); sender++) {

			// Zoltan returns cell migration info also for cells that don't chage their process
			if (sender == this->comm.rank()) {
				continue;
			}

			if (this->cells_to_receive.count(sender) == 0) {
				continue;
			}

			sort(this->cells_to_receive[sender].begin(), this->cells_to_receive[sender].end());

			int send_receive_tag = sender * this->comm.size() + this->comm.rank();

			// don't actually send one cell at a time or anything, just use the correct request structures
			#ifdef DCCRG_SEND_SINGLE_CELLS

			#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
			this->load_balance_requests[sender].push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#else
			this->receive_requests[sender].push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#endif


			#else
			this->receive_requests.push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#endif
		}

		// processes and the cells for which data has to be sent by this process
		for (int i = 0; i < number_to_send; i++) {
			if (this->cell_process[global_ids_to_send[i]] != this->comm.rank()) {
				std::cout << "Process " << this->comm.rank() << ", requested to send cell " << global_ids_to_send[i] << " but that cell is currently on process " << this->cell_process[global_ids_to_send[i]] << std::endl;
			}
			assert(this->cell_process[global_ids_to_send[i]] == this->comm.rank());
			assert(this->cells.count(global_ids_to_send[i]) > 0);

			this->cells_to_send[receiver_processes[i]].push_back(global_ids_to_send[i]);
			removed_cells.insert(global_ids_to_send[i]);
		}

		Zoltan_LB_Free_Data(&global_ids_to_receive, &local_ids_to_receive, &sender_processes, &global_ids_to_send, &local_ids_to_send, &receiver_processes);

		// gather data to send
		for (int receiver = 0; receiver < this->comm.size(); receiver++) {

			// Zoltan returns cell migration info also for cells that don't change their process
			if (receiver == this->comm.rank()) {
				continue;
			}

			if (this->cells_to_send.count(receiver) == 0) {
				continue;
			}

			sort(this->cells_to_send[receiver].begin(), this->cells_to_send[receiver].end());

			// construct the outgoing data vector
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_send[receiver].begin(); cell != this->cells_to_send[receiver].end(); cell++) {
				UserData* user_data = (*this)[*cell];
				assert(user_data != NULL);
				this->outgoing_data[receiver].push_back(*user_data);
				this->cells.erase(*cell);
				this->cells_with_remote_neighbours.erase(*cell);
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

			int send_receive_tag = this->comm.rank() * this->comm.size() + receiver;

			// don't actually receive one cell at a time or anything, just use the correct request structures
			#ifdef DCCRG_SEND_SINGLE_CELLS

			#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
			this->load_balance_requests[receiver].push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#else
			this->send_requests[receiver].push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#endif

			#else
			this->send_requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#endif
		}

		// wait for transfers to complete
		#ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->load_balance_requests.begin(); process != this->load_balance_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
			process->second.clear();
		}
		this->load_balance_requests.clear();

		#else

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->receive_requests.begin(); process != this->receive_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
			process->second.clear();
		}
		this->receive_requests.clear();
		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->send_requests.begin(); process != this->send_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
			process->second.clear();
		}
		this->send_requests.clear();

		#endif

		#else

		boost::mpi::wait_all(this->receive_requests.begin(), this->receive_requests.end());
		this->receive_requests.clear();
		boost::mpi::wait_all(this->send_requests.begin(), this->send_requests.end());
		this->send_requests.clear();

		#endif

		this->outgoing_data.clear();

		// incorporate received data
		for (typename boost::unordered_map<int, std::vector<UserData> >::const_iterator sender = this->incoming_data.begin(); sender != this->incoming_data.end(); sender++) {

			int i = 0;
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_receive[sender->first].begin(); cell != this->cells_to_receive[sender->first].end(); cell++, i++) {
				this->cells[*cell] = this->incoming_data[sender->first][i];
			}
		}
		this->cells_to_receive.clear();
		this->incoming_data.clear();

		/*
		Calculate where cells have migrated to update internal data structures
		Any cell can end up on any process and any neighbour of any cell can end up on yet another process
		*/
		// removed cells on all processes
		std::vector<uint64_t> temp_removed_cells(removed_cells.begin(), removed_cells.end());
		std::vector<std::vector<uint64_t> > all_removed_cells;
		all_gather(this->comm, temp_removed_cells, all_removed_cells);
		removed_cells.clear();
		// created cells on all processes
		std::vector<uint64_t> temp_added_cells(added_cells.begin(), added_cells.end());
		std::vector<std::vector<uint64_t> > all_added_cells;
		all_gather(this->comm, temp_added_cells, all_added_cells);
		added_cells.clear();

		// check that cells were removed by their process
		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
			for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {
				assert(this->cell_process[*removed_cell] == cell_remover);
			}
		}

		// update cell to process mappings
		for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {
				this->cell_process[*created_cell] = cell_creator;
			}
		}

		// if a child cell left this process, update remote neighbour info of nearby cells
		for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[this->comm.rank()].begin(); removed_cell != all_removed_cells[this->comm.rank()].end(); removed_cell++) {

			if (*removed_cell != this->get_child(*removed_cell)) {
				continue;
			}

			// neighbours of removed child
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*removed_cell].begin(); neighbour != this->neighbours[*removed_cell].end(); neighbour++) {
				if (this->cells.count(*neighbour) > 0) {
					this->update_remote_neighbour_info(*neighbour);
				}
			}
			this->neighbours.erase(*removed_cell);
			// cells with removed child as neighbour
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[*removed_cell].begin(); neighbour != this->neighbours_to[*removed_cell].end(); neighbour++) {
				if (this->cells.count(*neighbour) > 0) {
					this->update_remote_neighbour_info(*neighbour);
				}
			}
			this->neighbours_to.erase(*removed_cell);
		}

		// create neighbour lists for children that came to this process
		for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {

			if (*created_cell != this->get_child(*created_cell)) {
				continue;
			}

			this->neighbours[*created_cell] = this->find_neighbours_of(*created_cell);
			this->neighbours_to[*created_cell] = this->find_neighbours_to(*created_cell);
		}

		// if a child cell came to this process, update remote neighbour info of nearby cells
		for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {

			if (*created_cell != this->get_child(*created_cell)) {
				continue;
			}

			this->update_remote_neighbour_info(*created_cell);

			// neighbours of added child
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*created_cell].begin(); neighbour != this->neighbours[*created_cell].end(); neighbour++) {
				if (this->cells.count(*neighbour) > 0) {
					this->update_remote_neighbour_info(*neighbour);
				}
			}

			// cells with added child as neighbour
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[*created_cell].begin(); neighbour != this->neighbours_to[*created_cell].end(); neighbour++) {
				if (this->cells.count(*neighbour) > 0) {
					this->update_remote_neighbour_info(*neighbour);
				}
			}
		}

		// if a child cell changed its process and is neighbouring a local cell update the local cells neighbourhood
		for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {

				if (this->comm.rank() == cell_creator) {
					continue;
				}

				if (*created_cell != this->get_child(*created_cell)) {
					continue;
				}

				if (this->remote_cells_with_local_neighbours.count(*created_cell) == 0) {
					continue;
				}

				std::vector<uint64_t> temp_neighbours = this->find_neighbours_to(*created_cell);
				for (std::vector<uint64_t>::const_iterator neighbour = temp_neighbours.begin(); neighbour != temp_neighbours.end(); neighbour++) {
					if (this->cells.count(*neighbour)) {
						this->update_remote_neighbour_info(*neighbour);
					}
				}
			}
		}

		// free user data from migrated cells
		for (boost::unordered_set<uint64_t>::const_iterator removed_cell = removed_cells.begin(); removed_cell != removed_cells.end(); removed_cell++) {
			this->cells.erase(*removed_cell);
		}

		this->recalculate_neighbour_update_send_receive_lists();
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

		/*
		TODO: Find out why setting the message tags to zero here leads to this in wait_neighbour_data_update(), at least when using OpenMPI:
		terminate called after throwing an instance of 'boost::exception_detail::clone_impl<boost::exception_detail::error_info_injector<boost::mpi::exception> >'
		  what():  MPI_Wait: MPI_ERR_TRUNCATE: message truncated
		*/

		// user data is sent to another process one cell at a time
		#ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_USER_MPI_DATA_TYPE
		MPI_Datatype data_type = UserData::mpi_data_type();
		MPI_Type_commit(&data_type);
		#endif

		// post all receives, messages are unique between different senders so just iterate over processes in random order
		for (boost::unordered_map<int, std::vector<uint64_t> >::const_iterator sender = this->cells_to_receive.begin(); sender != this->cells_to_receive.end(); sender++) {

			for (std::vector<uint64_t>::const_iterator cell = sender->second.begin(); cell != sender->second.end(); cell++) {

				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->receive_requests[sender->first].push_back(MPI_Request());

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Irecv(this->remote_neighbours[*cell].at(), 1, data_type, sender->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->receive_requests[sender->first].back()));
				#else
				MPI_Irecv(this->remote_neighbours[*cell].at(), UserData::size(), MPI_BYTE, sender->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->receive_requests[sender->first].back()));
				#endif

				#else

				this->receive_requests[sender->first].push_back(this->comm.irecv(sender->first, *cell % boost::mpi::environment::max_tag(), this->remote_neighbours[*cell]));	// FIXME: make sure message tags between two processes are unique

				#endif
			}
		}

		// post all sends
		for (boost::unordered_map<int, std::vector<uint64_t> >::const_iterator receiver = this->cells_to_send.begin(); receiver != this->cells_to_send.end(); receiver++) {

			for (std::vector<uint64_t>::const_iterator cell = receiver->second.begin(); cell != receiver->second.end(); cell++) {
				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->send_requests[receiver->first].push_back(MPI_Request());

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Isend(this->cells[*cell].at(), 1, data_type, receiver->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->send_requests[receiver->first].back()));
				#else
				MPI_Isend(this->cells[*cell].at(), UserData::size(), MPI_BYTE, receiver->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->send_requests[receiver->first].back()));
				#endif

				#else

				this->send_requests[receiver->first].push_back(this->comm.isend(receiver->first, *cell % boost::mpi::environment::max_tag(), this->cells.at(*cell)));

				#endif
			}
		}

		#ifdef DCCRG_USER_MPI_DATA_TYPE
		MPI_Type_free(&data_type);
		#endif

		// user data is packed into a vector which is sent to another process
		#else

		// post all receives
		for (int sender = 0; sender < this->comm.size(); sender++) {

			if (sender == this->comm.rank()) {
				continue;
			}

			if (this->cells_to_receive.count(sender) == 0) {
				// no data to send / receive
				continue;
			}

			int send_receive_tag = sender * this->comm.size() + this->comm.rank();

			this->receive_requests.push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
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

			// construct the outgoing data vector
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_send[receiver].begin(); cell != this->cells_to_send[receiver].end(); cell++) {
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

			int send_receive_tag = this->comm.rank() * this->comm.size() + receiver;

			this->send_requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
		}
		#endif // ifdef DCCRG_SEND_SINGLE_CELLS
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
		#ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator process = this->send_requests.begin(); process != this->send_requests.end(); process++) {

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
		this->send_requests.clear();

		#else

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->send_requests.begin(); process != this->send_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}
		this->send_requests.clear();

		#endif

		#else

		boost::mpi::wait_all(this->send_requests.begin(), this->send_requests.end());
		this->send_requests.clear();
		this->outgoing_data.clear();

		#endif
	}


	/*!
	Waits until all receives associated with neighbour data update transfers between processes have completed and incorporates that data.
	Must be called simultaneously on all processes and probably must be called before wait...update_sends(void).
	*/
	void wait_neighbour_data_update_receives(void)
	{
		#ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator process = this->receive_requests.begin(); process != this->receive_requests.end(); process++) {

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
		this->receive_requests.clear();

		#else

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->receive_requests.begin(); process != this->receive_requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}
		this->receive_requests.clear();

		#endif

		#else

		boost::mpi::wait_all(this->receive_requests.begin(), this->receive_requests.end());
		this->receive_requests.clear();

		// incorporate received data
		for (typename boost::unordered_map<int, std::vector<UserData> >::const_iterator sender = this->incoming_data.begin(); sender != this->incoming_data.end(); sender++) {

			assert(this->incoming_data[sender->first].size() == this->cells_to_receive[sender->first].size());

			int i = 0;
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_receive[sender->first].begin(); cell != this->cells_to_receive[sender->first].end(); cell++, i++) {
				this->remote_neighbours[*cell] = this->incoming_data[sender->first][i];
			}
		}
		this->incoming_data.clear();

		#endif
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
	int get_refinement_level(uint64_t cell) const
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

		const uint64_t x_index = this->get_x_index(cell), y_index = this->get_y_index(cell), z_index = this->get_z_index(cell);

		// search neighbours in cells of the same size as the given cell (times neighbourhood size)
		const uint64_t size_in_indices = this->get_cell_size_in_indices(cell);

		const int refinement_level = this->get_refinement_level(cell);
		const int search_min_ref_level = (refinement_level == 0) ? 0 : refinement_level - 1;
		const int search_max_ref_level = (refinement_level == this->max_refinement_level) ? refinement_level : refinement_level + 1;

		#ifndef NDEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << cell << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << cell << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// the refinement level difference between neighbours is <= 1 so search index must increase by at most half the size of given cell
		const uint64_t index_increase = (size_in_indices < 2) ? 1 : size_in_indices / 2;

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours;

		// if neighbour_size == 0 just check the volume inside cells of the same size and that share a face with the given cell
		if (this->neighbourhood_size == 0) {

			// -x direction
			if (x_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index - size_in_indices; current_x_index < x_index; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +x direction
			if (x_index < this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index + size_in_indices; current_x_index < x_index + 2 * size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// -y direction
			if (y_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index - size_in_indices; current_y_index < y_index; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +y direction
			if (y_index < this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index + size_in_indices; current_y_index < y_index + 2 * size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// -z direction
			if (z_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index - size_in_indices; current_z_index < z_index; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +z direction
			if (z_index < this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index + size_in_indices; current_z_index < z_index + 2 * size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
			return return_neighbours;
		}


		// don't start searching outside of the grid
		uint64_t current_x_index;
		if (x_index < size_in_indices * this->neighbourhood_size) {
			current_x_index = 0;
		} else {
			current_x_index = x_index - size_in_indices * this->neighbourhood_size;
		}

		// search for neighbours in cells that share a vertex with the given cell
		for (; current_x_index < x_index + size_in_indices * (1 + this->neighbourhood_size); current_x_index += index_increase) {

			// don't search outside of the grid
			if (current_x_index >= this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level)) {
				continue;
			}

			// don't start searching outside of the grid
			uint64_t current_y_index;
			if (y_index < size_in_indices * this->neighbourhood_size) {
				current_y_index = 0;
			} else {
				current_y_index = y_index - size_in_indices * this->neighbourhood_size;
			}

			for (; current_y_index < y_index + size_in_indices * (1 + this->neighbourhood_size); current_y_index += index_increase) {

				if (current_y_index >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)) {
					continue;
				}

				// don't start searching outside of the grid
				uint64_t current_z_index;
				if (z_index < size_in_indices * this->neighbourhood_size) {
					current_z_index = 0;
				} else {
					current_z_index = z_index - size_in_indices * this->neighbourhood_size;
				}

				for (; current_z_index < z_index + size_in_indices * (1 + this->neighbourhood_size); current_z_index += index_increase) {

					if (current_z_index >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
						continue;
					}

					// don't search inside of the given cell
					if (current_x_index >= x_index
					    && current_y_index >= y_index
					    && current_z_index >= z_index
					    && current_x_index < x_index + size_in_indices
					    && current_y_index < y_index + size_in_indices
					    && current_z_index < z_index + size_in_indices) {
						continue;
					}

					const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
					assert(neighbour > 0);
					assert(neighbour <= this->max_cell_number);
					assert(neighbour == this->get_child(neighbour));
					unique_neighbours.insert(neighbour);
				}
			}
		}

		// a cell isn't a neighbour of itself
		unique_neighbours.erase(cell);

		return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
		return return_neighbours;
	}


	/* TODO: use this in both neighbour finding functions?
	std::vector<uint64_t> find_cells(x_index, y_index, z_index, min_distance, max_distance, min_refinement_level, max_refinement_level)
	{
	}*/


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

		const int search_min_ref_level = (refinement_level == 0) ? 0 : refinement_level - 1;
		const int search_max_ref_level = (refinement_level == this->max_refinement_level) ? refinement_level : refinement_level + 1;

		const uint64_t x_index = this->get_x_index(cell), y_index = this->get_y_index(cell), z_index = this->get_z_index(cell);

		// can be a neighbour to larger cells further away than own size * neighbourhood
		const uint64_t size_in_indices = 2 * this->get_cell_size_in_indices(cell);

		#ifndef NDEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << cell << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << cell << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// the refinement level difference between neighbours is <= 1 so search index must increase by at most half the size of given cell
		const uint64_t index_increase = (size_in_indices < 4) ? 1 : size_in_indices / 4;

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours;

		// if neighbour_size == 0 just check the volume inside cells of the same size and that share a face with the given cell
		if (this->neighbourhood_size == 0) {

			// -x direction
			if (x_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index - size_in_indices; current_x_index < x_index; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +x direction
			if (x_index < this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index + size_in_indices; current_x_index < x_index + 2 * size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// -y direction
			if (y_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index - size_in_indices; current_y_index < y_index; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +y direction
			if (y_index < this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index + size_in_indices; current_y_index < y_index + 2 * size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// -z direction
			if (z_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index - size_in_indices; current_z_index < z_index; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +z direction
			if (z_index < this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index + size_in_indices; current_z_index < z_index + 2 * size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
			return return_neighbours;
		}


		// don't start searching outside of the grid
		uint64_t current_x_index;
		if (x_index < size_in_indices * this->neighbourhood_size) {
			current_x_index = 0;
		} else {
			current_x_index = x_index - size_in_indices * this->neighbourhood_size;
		}

		// search for neighbours in cells that share a vertex with the given cell
		for (; current_x_index < x_index + size_in_indices * (1 + this->neighbourhood_size); current_x_index += index_increase) {

			// don't search outside of the grid
			if (current_x_index >= this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level)) {
				continue;
			}

			// don't start searching outside of the grid
			uint64_t current_y_index;
			if (y_index < size_in_indices * this->neighbourhood_size) {
				current_y_index = 0;
			} else {
				current_y_index = y_index - size_in_indices * this->neighbourhood_size;
			}

			for (; current_y_index < y_index + size_in_indices * (1 + this->neighbourhood_size); current_y_index += index_increase) {

				if (current_y_index >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)) {
					continue;
				}

				// don't start searching outside of the grid
				uint64_t current_z_index;
				if (z_index < size_in_indices * this->neighbourhood_size) {
					current_z_index = 0;
				} else {
					current_z_index = z_index - size_in_indices * this->neighbourhood_size;
				}

				for (; current_z_index < z_index + size_in_indices * (1 + this->neighbourhood_size); current_z_index += index_increase) {

					if (current_z_index >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
						continue;
					}

					// don't search inside of the given cell
					if (current_x_index >= x_index
					    && current_y_index >= y_index
					    && current_z_index >= z_index
					    && current_x_index < x_index + size_in_indices
					    && current_y_index < y_index + size_in_indices
					    && current_z_index < z_index + size_in_indices) {
						continue;
					}

					const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
					assert(neighbour > 0);
					assert(neighbour <= this->max_cell_number);
					#ifndef NDEBUG
					if (neighbour != this->get_child(neighbour)) {
						std::cerr << __FILE__ << ":" << __LINE__ << " While searching for neighbours_to of cell " << cell << " (ref lvl " << this->get_refinement_level(cell) << ") between refinement levels " << search_min_ref_level << " and " << search_max_ref_level << ": returned neighbour (" << neighbour << ", ref lvl " << this->get_refinement_level(neighbour) << ") has a child: " << this->get_child(neighbour) << " (ref lvl " << this->get_refinement_level(this->get_child(neighbour)) << ")" << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif
					unique_neighbours.insert(neighbour);
				}
			}
		}

		// a cell isn't a neighbour of itself
		unique_neighbours.erase(cell);

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
			#ifndef NDEBUG
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
			#ifndef NDEBUG
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
			#ifndef NDEBUG
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

	#ifdef DCCRG_SEND_SINGLE_CELLS	// user data is sent to another process one cell at a time
	// list of pending transfers between this process and the process as the key
	#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::unordered_map<int, std::vector<MPI_Request> > send_requests;
	boost::unordered_map<int, std::vector<MPI_Request> > receive_requests;
	boost::unordered_map<int, std::vector<boost::mpi::request> > load_balance_requests;	// TODO: get rid of this
	#else
	boost::unordered_map<int, std::vector<boost::mpi::request> > send_requests;
	boost::unordered_map<int, std::vector<boost::mpi::request> > receive_requests;
	#endif
	#else	// user data is packed into a vector which is sent to another process
	// pending neighbour data requests for this process
	std::vector<boost::mpi::request> send_requests;
	std::vector<boost::mpi::request> receive_requests;
	#endif

	// cells whose data has to be received / sent by this process from the process as the key
	boost::unordered_map<int, std::vector<uint64_t> > cells_to_receive;
	boost::unordered_map<int, std::vector<uint64_t> > cells_to_send;

	// storage for cells' user data that awaits transfer to or from this process
	boost::unordered_map<int, std::vector<UserData> > incoming_data, outgoing_data;

	// cells to be refined / unrefined after a call to stop_refining()
	boost::unordered_set<uint64_t> cells_to_refine, cells_to_unrefine;

	// stores user data of cells whose children were created while refining
	boost::unordered_map<uint64_t, UserData> refined_cell_data;
	// stores user data of cells that were removed while unrefining
	boost::unordered_map<uint64_t, UserData> unrefined_cell_data;

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

		#ifndef NDEBUG
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

		if (this->cell_process[cell] != this->comm.rank()) {
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

			#ifndef NDEBUG
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

					#ifndef NDEBUG
					if (*sibling == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Invalid sibling" << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif
				}
			}

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(parent).begin(); neighbour != this->neighbours.at(parent).end(); neighbour++) {
				old_neighbours.push_back(*neighbour);

				#ifndef NDEBUG
				if (*neighbour == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour for parent of cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

			for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(parent).begin(); neighbour_to != this->neighbours_to.at(parent).end(); neighbour_to++) {
				old_neighbours.push_back(*neighbour_to);

				#ifndef NDEBUG
				if (*neighbour_to == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour_to for parent of cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

		// use given cell's neighbour lists
		} else {

			#ifndef NDEBUG
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

				#ifndef NDEBUG
				if (*neighbour == 0) {
					std::cerr << __FILE__ << ":" << __LINE__ << " Invalid neighbour for cell " << cell << std::endl;
					exit(EXIT_FAILURE);
				}
				#endif
			}

			for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(cell).begin(); neighbour_to != this->neighbours_to.at(cell).end(); neighbour_to++) {
				old_neighbours.push_back(*neighbour_to);

				#ifndef NDEBUG
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

			#ifndef NDEBUG
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

					#ifndef NDEBUG
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
				#ifndef NDEBUG
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

				#ifndef NDEBUG
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

				#ifndef NDEBUG
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

		#ifndef NDEBUG
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

			#ifndef NDEBUG
			if (*candidate != this->get_child(*candidate)) {
				std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour candidate " << *candidate << " of cell " << cell << " has children" << std::endl;
				exit(EXIT_FAILURE);
			}
			#endif
		}

		#ifndef NDEBUG
		// verify the result for neighbours_of
		sort(this->neighbours[cell].begin(), this->neighbours[cell].end());
		std::vector<uint64_t> compare_neighbours = this->find_neighbours_of(cell);
		sort(compare_neighbours.begin(), compare_neighbours.end());

		if (!std::equal(this->neighbours[cell].begin(), this->neighbours[cell].end(), compare_neighbours.begin())) {
			std::cerr << "Process " << this->comm.rank() << " neighbour counts for cell " << cell << " don't match: " << this->neighbours[cell].size() << " (";
			for (std::vector<uint64_t>::const_iterator c = this->neighbours[cell].begin(); c != this->neighbours[cell].end(); c++) {
				std::cerr << *c << " ";
			}
			std::cerr << ") should be " << compare_neighbours.size() << " (";
			for (std::vector<uint64_t>::const_iterator c = compare_neighbours.begin(); c != compare_neighbours.end(); c++) {
				std::cerr << *c << " ";
			}
			std::cerr << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		// verify the result for neighbours_to
		sort(this->neighbours_to[cell].begin(), this->neighbours_to[cell].end());
		std::vector<uint64_t> compare_neighbours_to = this->find_neighbours_to(cell);
		sort(compare_neighbours_to.begin(), compare_neighbours_to.end());

		if (!std::equal(this->neighbours_to[cell].begin(), this->neighbours_to[cell].end(), compare_neighbours_to.begin())) {
			std::cerr << "Process " << this->comm.rank() << " neighbour_to counts for cell " << cell << " at indices " << this->get_x_index(cell) << ", " << this->get_y_index(cell) << ", " << this->get_z_index(cell) << " (child of " << this->get_parent(cell) << " at indices " << this->get_x_index(this->get_parent(cell)) << ", " << this->get_y_index(this->get_parent(cell)) << ", " << this->get_z_index(this->get_parent(cell)) << ") don't match: " << this->neighbours_to[cell].size() << " (";
			for (std::vector<uint64_t>::const_iterator c = this->neighbours_to[cell].begin(); c != this->neighbours_to[cell].end(); c++) {
				std::cerr << *c << " ";
			}
			std::cerr << ") should be " << compare_neighbours_to.size() << " (";
			for (std::vector<uint64_t>::const_iterator c = compare_neighbours_to.begin(); c != compare_neighbours_to.end(); c++) {
				std::cerr << *c << " ";
			}
			std::cerr << ")" << std::endl;

			std::cerr << "Neighbour candidates of cell " << cell << ": ";
			for (boost::unordered_set<uint64_t>::const_iterator candidate = neighbour_candidates.begin(); candidate != neighbour_candidates.end(); candidate++) {
				std::cerr << *candidate << " ";
			}
			std::cerr << std::endl;

			if (cell_result_of_refining) {
				std::cerr << "Neighbours of parent " << parent << ": ";
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(parent).begin(); neighbour != this->neighbours.at(parent).end(); neighbour++) {
					std::cerr << *neighbour << " ";
				}
				std::cerr << std::endl;
			}
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

		// neighbours of given cell
		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[cell].begin(); neighbour != this->neighbours[cell].end(); neighbour++) {
			if (this->cell_process[*neighbour] != this->comm.rank()) {
				this->cells_with_remote_neighbours.insert(cell);
				this->remote_cells_with_local_neighbours.insert(*neighbour);
			}
		}
		// cells with given cell as neighbour
		for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to[cell].begin(); neighbour_to != this->neighbours_to[cell].end(); neighbour_to++) {
			if (this->cell_process[*neighbour_to] != this->comm.rank()) {
				this->cells_with_remote_neighbours.insert(cell);
				this->remote_cells_with_local_neighbours.insert(*neighbour_to);
			}
		}
	}


	/*!
	Returns the existing neighbours of given cell's parent or existing neighbours of given cell if it is an unrefined cell (= without a parent)

	Assumes that given cell will be unrefined so the refinement level of possible neighbours of the parent is <= n and >= n - 2, where n is the refinement level of given cell
	*/
	std::vector<uint64_t> find_neighbours_of_parent(const uint64_t cell) const
	{

		std::vector<uint64_t> return_neighbours;

		const uint64_t parent = this->get_parent(cell);
		#ifndef NDEBUG
		if (parent == 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid parent for cell " << cell << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		if (cell == 0
		|| cell > this->max_cell_number) {
			return return_neighbours;
		}

		if (cell != this->get_child(cell)) {
			return return_neighbours;
		}

		const uint64_t x_index = this->get_x_index(parent), y_index = this->get_y_index(parent), z_index = this->get_z_index(parent);

		// search neighbours in cells of the same size as the given cell (times neighbourhood size)
		const uint64_t size_in_indices = this->get_cell_size_in_indices(parent);

		const int refinement_level = this->get_refinement_level(parent);
		const int search_min_ref_level = (refinement_level == 0) ? 0 : refinement_level - 1;
		const int search_max_ref_level = (refinement_level == this->max_refinement_level) ? refinement_level : refinement_level + 1;

		#ifndef NDEBUG
		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level (" << refinement_level << ") of cell " << parent << " exceeds maximum refinement level of the grid (" << this->max_refinement_level << ")" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Refinement level of cell " << parent << " is less than 0: " << refinement_level << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// the refinement level difference between neighbours is <= 1 so search index must increase by at most half the size of given cell
		const uint64_t index_increase = (size_in_indices < 2) ? 1 : size_in_indices / 2;

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours;

		// if neighbour_size == 0 just check the volume inside cells of the same size and that share a face with the given cell
		if (this->neighbourhood_size == 0) {

			// -x direction
			if (x_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index - size_in_indices; current_x_index < x_index; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +x direction
			if (x_index < this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index + size_in_indices; current_x_index < x_index + 2 * size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// -y direction
			if (y_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index - size_in_indices; current_y_index < y_index; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +y direction
			if (y_index < this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index + size_in_indices; current_y_index < y_index + 2 * size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index; current_z_index < z_index + size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// -z direction
			if (z_index >= size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index - size_in_indices; current_z_index < z_index; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			// +z direction
			if (z_index < this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level) - size_in_indices) {
			for (uint64_t current_x_index = x_index; current_x_index < x_index + size_in_indices; current_x_index += index_increase) {
			for (uint64_t current_y_index = y_index; current_y_index < y_index + size_in_indices; current_y_index += index_increase) {
			for (uint64_t current_z_index = z_index + size_in_indices; current_z_index < z_index + 2 * size_in_indices; current_z_index += index_increase) {
				const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
				assert(neighbour > 0);
				assert(neighbour <= this->max_cell_number);
				assert(neighbour == this->get_child(neighbour));
				unique_neighbours.insert(neighbour);
			}}}}

			return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
			return return_neighbours;
		}


		// don't start searching outside of the grid
		uint64_t current_x_index;
		if (x_index < size_in_indices * this->neighbourhood_size) {
			current_x_index = 0;
		} else {
			current_x_index = x_index - size_in_indices * this->neighbourhood_size;
		}

		// search for neighbours in cells that share a vertex with the given cell
		for (; current_x_index < x_index + size_in_indices * (1 + this->neighbourhood_size); current_x_index += index_increase) {

			// don't search outside of the grid
			if (current_x_index >= this->geometry.get_x_length() * (uint64_t(1) << this->max_refinement_level)) {
				continue;
			}

			// don't start searching outside of the grid
			uint64_t current_y_index;
			if (y_index < size_in_indices * this->neighbourhood_size) {
				current_y_index = 0;
			} else {
				current_y_index = y_index - size_in_indices * this->neighbourhood_size;
			}

			for (; current_y_index < y_index + size_in_indices * (1 + this->neighbourhood_size); current_y_index += index_increase) {

				if (current_y_index >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)) {
					continue;
				}

				// don't start searching outside of the grid
				uint64_t current_z_index;
				if (z_index < size_in_indices * this->neighbourhood_size) {
					current_z_index = 0;
				} else {
					current_z_index = z_index - size_in_indices * this->neighbourhood_size;
				}

				for (; current_z_index < z_index + size_in_indices * (1 + this->neighbourhood_size); current_z_index += index_increase) {

					if (current_z_index >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
						continue;
					}

					// don't search inside of the given cell's parent
					if (current_x_index >= x_index
					    && current_y_index >= y_index
					    && current_z_index >= z_index
					    && current_x_index < x_index + size_in_indices
					    && current_y_index < y_index + size_in_indices
					    && current_z_index < z_index + size_in_indices) {
						continue;
					}

					const uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
					assert(neighbour > 0);
					assert(neighbour <= this->max_cell_number);
					assert(neighbour == this->get_child(neighbour));
					unique_neighbours.insert(neighbour);
				}
			}
		}

		// a cell isn't a neighbour to itself
		unique_neighbours.erase(parent);

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

		if (this->neighbourhood_size == 0) {
			// for cells to share a face 2 indices must overlap
			if (this->overlapping_indices(cell1, cell2) < 2) {
				return false;
			}

			// cells are close enough in x direction
			if (cell1_x_index + cell1_size >= cell2_x_index
			&& cell1_x_index <= cell2_x_index + cell2_size
			&& !this->x_indices_overlap(cell1, cell2)) {
				return true;
			}

			// cells are close enough in y direction
			if (cell1_y_index + cell1_size >= cell2_y_index
			&& cell1_y_index <= cell2_y_index + cell2_size
			&& !this->y_indices_overlap(cell1, cell2)) {
				return true;
			}

			// cells are close enough in z direction
			if (cell1_z_index + cell1_size >= cell2_z_index
			&& cell1_z_index <= cell2_z_index + cell2_size
			&& !this->z_indices_overlap(cell1, cell2)) {
				return true;
			}

			return false;
		}

		const uint64_t dindex1 = cell2_size + cell1_size * this->neighbourhood_size;
		uint64_t dindex2 = cell1_size * this->neighbourhood_size;	// TODO: make this also const
		if (cell2_size < cell1_size) {
			dindex2 += cell2_size;
		}

		if (cell1_x_index < cell2_x_index + dindex1
		    && cell1_y_index < cell2_y_index + dindex1
		    && cell1_z_index < cell2_z_index + dindex1
		    && cell1_x_index + dindex2 >= cell2_x_index
		    && cell1_y_index + dindex2 >= cell2_y_index
		    && cell1_z_index + dindex2 >= cell2_z_index) {
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

					#ifndef NDEBUG
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

					#ifndef NDEBUG
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

		#ifndef NDEBUG
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

			#ifndef NDEBUG
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

		#ifndef NDEBUG
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
		std::vector<uint64_t> new_cells;

		this->remote_neighbours.clear();
		this->cells_to_send.clear();
		this->cells_to_receive.clear();
		this->incoming_data.clear();
		this->outgoing_data.clear();
		this->refined_cell_data.clear();
		this->unrefined_cell_data.clear();

		#ifndef NDEBUG
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

			// use local neighbour lists to find cells whose neighbour lists have to updated
			if (process_of_refined == this->comm.rank()) {
				// update the neighbour lists of created local cells
				for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
						update_neighbours.insert(*child);

						#ifndef NDEBUG
						if (this->neighbours.count(*child) > 0) {
							std::cerr << __FILE__ << ":" << __LINE__ << " Neighbours for cell " << *child << " shouldn't exist yet" << std::endl;
							exit(EXIT_FAILURE);
						}
						#endif
				}

				// update neighbour lists of all the parent's neighbours
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(*refined).begin(); neighbour != this->neighbours.at(*refined).end(); neighbour++) {
					if (this->cell_process[*neighbour] == this->comm.rank()) {
						update_neighbours.insert(*neighbour);
					}
				}
				for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to.at(*refined).begin(); neighbour_to != this->neighbours_to.at(*refined).end(); neighbour_to++) {
					if (this->cell_process[*neighbour_to] == this->comm.rank()) {
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

			this->cells_with_remote_neighbours.erase(*refined);
			this->remote_cells_with_local_neighbours.erase(*refined);
		}

		// needed for checking which neighbourhoods to update due to unrefining
		boost::unordered_set<uint64_t> parents_of_unrefined;

		// initially only one sibling is recorded per process when unrefining, insert the rest of them now
		boost::unordered_set<uint64_t> all_to_unrefine;
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = this->cells_to_unrefine.begin(); unrefined != this->cells_to_unrefine.end(); unrefined++) {

			const uint64_t parent_of_unrefined = this->get_parent(*unrefined);
			#ifndef NDEBUG
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
			this->cells_with_remote_neighbours.erase(*unrefined);
			this->remote_cells_with_local_neighbours.erase(*unrefined);
			update_neighbours.erase(*unrefined);

			// don't send unrefined cells' user data to self
			if (this->comm.rank() == process_of_unrefined
			&& this->comm.rank() == process_of_parent) {
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

		// start sending unrefined cells' user data to the process of their parent
		// TODO: use the same function for all data transfer stuff
		// post all receives
		for (int sender = 0; sender < this->comm.size(); sender++) {

			if (sender == this->comm.rank()) {
				continue;
			}

			if (this->cells_to_receive.count(sender) == 0) {
				// no data to send / receive
				continue;
			}

			int send_receive_tag = sender * this->comm.size() + this->comm.rank();

			this->receive_requests.push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
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

			// construct the outgoing data vector
			sort(this->cells_to_send[receiver].begin(), this->cells_to_send[receiver].end());
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_send[receiver].begin(); cell != this->cells_to_send[receiver].end(); cell++) {
				UserData* user_data = (*this)[*cell];
				assert(user_data != NULL);
				this->outgoing_data[receiver].push_back(*user_data);
				// remove local data of unrefined cell
				this->cells.erase(*cell);
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

			int send_receive_tag = this->comm.rank() * this->comm.size() + receiver;

			this->send_requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
		}

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

				this->update_remote_neighbour_info(*parent);
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
			this->update_remote_neighbour_info(*cell);
		}

		// remove neighbour lists of refined cells' parents
		for (boost::unordered_set<uint64_t>::const_iterator refined = this->cells_to_refine.begin(); refined != this->cells_to_refine.end(); refined++) {
			if (this->cell_process[*refined] == this->comm.rank()) {

				#ifndef NDEBUG
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
			}
		}

		// remove neighbour lists of removed cells
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = this->cells_to_unrefine.begin(); unrefined != this->cells_to_unrefine.end(); unrefined++) {

			// only one sibling of cells to unrefine is stored, but were removed
			std::vector<uint64_t> siblings = this->get_all_children(this->get_parent_for_removed(*unrefined));
			for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {
				this->neighbours.erase(*sibling);
				this->neighbours_to.erase(*sibling);
			}
		}

		#ifndef NDEBUG
		if (!this->verify_neighbours()) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Neighbour lists are inconsistent" << std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		// wait for receives to complete
		boost::mpi::wait_all(this->receive_requests.begin(), this->receive_requests.end());
		this->receive_requests.clear();

		// incorporate received data
		for (typename boost::unordered_map<int, std::vector<UserData> >::const_iterator sender = this->incoming_data.begin(); sender != this->incoming_data.end(); sender++) {

			assert(this->incoming_data[sender->first].size() == this->cells_to_receive[sender->first].size());

			int i = 0;
			sort(this->cells_to_receive[sender->first].begin(), this->cells_to_receive[sender->first].end());
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_receive[sender->first].begin(); cell != this->cells_to_receive[sender->first].end(); cell++, i++) {
				this->unrefined_cell_data[*cell] = this->incoming_data[sender->first][i];
			}
		}
		this->incoming_data.clear();

		// remove user data of unrefined cells from this->cells
		for (boost::unordered_set<uint64_t>::const_iterator unrefined = this->cells_to_unrefine.begin(); unrefined != this->cells_to_unrefine.end(); unrefined++) {

			std::vector<uint64_t> siblings = this->get_all_children(this->get_parent_for_removed(*unrefined));
			for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {
				this->cells.erase(*sibling);
			}
		}

		// wait for sends to complete
		boost::mpi::wait_all(this->send_requests.begin(), this->send_requests.end());
		this->send_requests.clear();
		this->outgoing_data.clear();

		this->cells_to_refine.clear();
		this->cells_to_unrefine.clear();

		this->recalculate_neighbour_update_send_receive_lists();

		return new_cells;
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

		if (index2 + size2 > index1 && index2 < index1 + size1) {
			return true;
		} else {
			return false;
		}
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

		if (index2 + size2 > index1 && index2 < index1 + size1) {
			return true;
		} else {
			return false;
		}
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

		if (index2 + size2 > index1 && index2 < index1 + size1) {
			return true;
		} else {
			return false;
		}
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
			return 0;
		}

		if (y_index >= this->geometry.get_y_length() * (uint64_t(1) << this->max_refinement_level)) {
			return 0;
		}

		if (z_index >= this->geometry.get_z_length() * (uint64_t(1) << this->max_refinement_level)) {
			return 0;
		}

		if (minimum_refinement_level > maximum_refinement_level) {
			return 0;
		}

		int average_refinement_level = (maximum_refinement_level + minimum_refinement_level) / 2;
		uint64_t id = this->get_cell_from_indices(x_index, y_index, z_index, average_refinement_level);

		// use binary search recursively (assumes that a cells refine to 8 children)
		if (this->cell_process.count(id) == 0) {
			// doesn't exist, search the bin of smaller refinement_level values
			if (average_refinement_level > minimum_refinement_level) {

				uint64_t smaller_refinement_value_cell = this->get_cell_from_indices(x_index, y_index, z_index, minimum_refinement_level, average_refinement_level - 1);

				#ifndef NDEBUG
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

					#ifndef NDEBUG
					if (this->cell_process.count(larger_refinement_value_cell) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell: " << larger_refinement_value_cell << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif

					return larger_refinement_value_cell;

				} else {

					// current cell has the largest refinement value at given indices

					#ifndef NDEBUG
					if (this->cell_process.count(id) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__ << " Returning non-existing cell: " << id << std::endl;
						exit(EXIT_FAILURE);
					}
					#endif

					return id;
				}
			} else {
				// nothing left to search

				#ifndef NDEBUG
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
		#ifndef NDEBUG
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



	#ifndef NDEBUG
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
	Returns false if neighbour lists on this process aren't consistent
	*/
	bool verify_neighbours(void)
	{
		for (typename boost::unordered_map<uint64_t, int>::const_iterator cell = this->cell_process.begin(); cell != this->cell_process.end(); cell++) {

			if (cell->second != this->comm.rank()) {
				continue;
			}

			if (cell->first != this->get_child(cell->first)) {
				continue;
			}

			if (this->neighbours.count(cell->first) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << " No neighbour list for cell " << cell->first << std::endl;
				return false;
			}

			// neighbours
			sort(this->neighbours[cell->first].begin(), this->neighbours[cell->first].end());
			std::vector<uint64_t> compare_neighbours = this->find_neighbours_of(cell->first);
			sort(compare_neighbours.begin(), compare_neighbours.end());

			if (!std::equal(this->neighbours[cell->first].begin(), this->neighbours[cell->first].end(), compare_neighbours.begin())) {
				std::cerr << "Process " << this->comm.rank() << " neighbour counts for cell " << cell->first << " don't match: " << this->neighbours[cell->first].size() << " (";
				for (std::vector<uint64_t>::const_iterator c = this->neighbours[cell->first].begin(); c != this->neighbours[cell->first].end(); c++) {
					std::cerr << *c << " ";
				}
				std::cerr << ") should be " << compare_neighbours.size() << " (";
				for (std::vector<uint64_t>::const_iterator c = compare_neighbours.begin(); c != compare_neighbours.end(); c++) {
					std::cerr << *c << " ";
				}
				std::cerr << ")" << std::endl;
				return false;
			}

			// neighbours_to
			sort(this->neighbours_to[cell->first].begin(), this->neighbours_to[cell->first].end());
			std::vector<uint64_t> compare_neighbours_to = this->find_neighbours_to(cell->first);
			sort(compare_neighbours_to.begin(), compare_neighbours_to.end());

			if (!std::equal(this->neighbours_to[cell->first].begin(), this->neighbours_to[cell->first].end(), compare_neighbours_to.begin())) {
				std::cerr << "Process " << this->comm.rank() << " neighbour_to counts for cell " << cell->first << " (child of " << this->get_parent(cell->first) << ") don't match: " << this->neighbours_to[cell->first].size() << " (";
				for (std::vector<uint64_t>::const_iterator c = this->neighbours_to[cell->first].begin(); c != this->neighbours_to[cell->first].end(); c++) {
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
		}

		return true;
	}
	#endif

};

#endif
