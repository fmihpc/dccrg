/*
A distributed cartesian cell-refinable grid

Copyright 2009, 2010 Ilja Honkonen

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
	Creates a new grid where each cell stores one instance of given UserData
	Zoltan_Initialize must be called before calling this constructor

	comm: the grid will span all the processes in the communicator comm

	load_balancing_method:
		The method that Zoltan will use for load balancing given as a string
		Currently supported methods: NONE, BLOCK, RANDOM, RCB, RIB, HSFC, GRAPH and HYPERGRAPH

	#ifdef DCCRG_ARBITRARY_STRETCH
	x, y and z_coordinates:
		The coordinates of unrefined cells in the respective direction
		First coordinate is the starting point of the grid, the following ith value is the endpoint of the ith unrefined cell
	#else
	x_start, y_start, z_start:
		the starting corner of the grid
	cell_size:
		the size of each unrefined cell in every direction
	x_length, y_length, z_length:
		the number of cells in the grid in x, y and z direction
	#endif

	neighbourhood_size:
		Determines which cells are considered neighbours.
		When calculating the neighbours of a given cell a cube of length neighbourhood_size + 1 in every direction is considered, centered at the cell for which neighbours are being calculated.
		The unit lenght of the cube is the cell for which neighbours are being calculated.
		If neighbourhood_size == 0, only cells within the volume of cells that share a face are considered.

	maximum_refinement_level:
		The maximum number of times an unrefined cell can be refined (replacing it with 8 smaller cells)
		Optional: if not given it is maximized based on the grids initial size
	 */
	#ifdef DCCRG_ARBITRARY_STRETCH
	dccrg(boost::mpi::communicator comm, const char* load_balancing_method, const std::vector<double> x_coordinates, const std::vector<double> y_coordinates, const std::vector<double> z_coordinates, const unsigned int neighbourhood_size, const int maximum_refinement_level = -1)
	#else
	dccrg(boost::mpi::communicator comm, const char* load_balancing_method, const double x_start, const double y_start, const double z_start, const double cell_x_size, const double cell_y_size, const double cell_z_size, const uint64_t x_length, const uint64_t y_length, const uint64_t z_length, const unsigned int neighbourhood_size, const int maximum_refinement_level = -1)
	#endif
	{
		this->comm = comm;

		// Setup Zoltan
		this->zoltan = Zoltan_Create(this->comm);
		if (this->zoltan == NULL) {
			std::cerr << "Zoltan_Create failed"  << std::endl;
			exit(EXIT_FAILURE);
		}

		// Set Zoltan parameters
		// http://www.cs.sandia.gov/Zoltan/ug_html/ug_alg.html#LB_METHOD
		Zoltan_Set_Param(this->zoltan, "LB_METHOD", load_balancing_method);

		// check whether Zoltan_LB_Partition is expected to fail
		if (strncmp(load_balancing_method, "NONE", sizeof("NONE")) == 0) {
			this->no_load_balancing = true;
		} else {
			this->no_load_balancing = false;
		}

		Zoltan_Set_Param(this->zoltan, "DEBUG_LEVEL", "0");
		// size of cells id in unsigned ints, but has to be 1 even when global id is uint64_t, for some reason
		/*char global_id_length_string[10];
		snprintf(global_id_length_string, 10, "%0i", int(sizeof(uint64_t) / sizeof(unsigned int)));*/
		Zoltan_Set_Param(this->zoltan, "NUM_GID_ENTRIES", "1");
		// no object weights
		Zoltan_Set_Param(this->zoltan, "OBJ_WEIGHT_DIM", "1");
		Zoltan_Set_Param(this->zoltan, "EDGE_WEIGHT_DIM", "1");
		//Zoltan_Set_Param(this->zoltan, "PHG_REFINEMENT_QUALITY", "1.5");
		// try to minimize moving of data between processes
		Zoltan_Set_Param(this->zoltan, "REMAP", "1");
		// when load balancing return only cells whose process changed
		Zoltan_Set_Param(this->zoltan, "RETURN_LISTS", "ALL");
		Zoltan_Set_Param(this->zoltan, "IMBALANCE_TOL", "1.1");

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
			exit(EXIT_FAILURE);
		}
		this->geometry.set_x_length(x_length);

		if (y_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the y direction" << std::endl;
			exit(EXIT_FAILURE);
		}
		this->geometry.set_y_length(y_length);

		if (z_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the z direction" << std::endl;
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
			exit(EXIT_FAILURE);
		}


		if (maximum_refinement_level > refinement_level) {

			std::cerr << "Given max_refinement_level (" << maximum_refinement_level << ") is too large: " << "x_length * this->geometry.get_y_length() * this->geometry.get_z_length() * 8^max_refinement_level / (2^64 - 1) >= " << this->geometry.get_x_length() *  this->geometry.get_y_length() * this->geometry.get_z_length() * (uint64_t(1) << maximum_refinement_level * 3) / max_id << " but must be < 1" << std::endl;
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
			this->update_neighbours(cell->first);
		}

		// update neighbour_to lists of created cells
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
			this->update_remote_neighbour_info(cell->first);
		}

		this->recalculate_neighbour_update_send_receive_lists();
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
		} else if (this->removed_cell_data.count(cell) > 0) {
			return (UserData*) &(this->removed_cell_data.at(cell));
		} else {
			return NULL;
		}
	}


	/*!
	Load balances the grids cells among processes
	Must be called simultaneously on all processes
	Does not update remote neighbour data of processes afterwards
	*/
	void balance_load(void)
	{
		this->comm.barrier();

		int partition_changed, global_id_size, local_id_size, number_to_receive, number_to_send;
		ZOLTAN_ID_PTR global_ids_to_receive, local_ids_to_receive, global_ids_to_send, local_ids_to_send;
		int *sender_processes, *receiver_processes;

		// here these record where cells have migrated
		assert(this->added_cells.size() == 0);
		assert(this->removed_cells.size() == 0);

		if (Zoltan_LB_Balance(this->zoltan, &partition_changed, &global_id_size, &local_id_size, &number_to_receive, &global_ids_to_receive, &local_ids_to_receive, &sender_processes, &number_to_send, &global_ids_to_send, &local_ids_to_send, &receiver_processes) != ZOLTAN_OK) {
			if (!this->no_load_balancing) {
				std::cerr << "Zoltan_LB_Partition failed" << std::endl;
				Zoltan_Destroy(&this->zoltan);
				exit(EXIT_FAILURE);
			}
		}

		if (partition_changed == 0) {
			return;
		}
		this->removed_cell_data.clear();

		// clear send / receive lists, here they mean cells that will be moved between processes
		this->cells_to_receive.clear();
		this->cells_to_send.clear();

		// TODO: move identical code from here and start_neighbour_data_update into a separate function

		// processes and the cells for which data has to be received by this process
		for (int i = 0; i < number_to_receive; i++) {
			assert(this->cell_process[global_ids_to_receive[i]] == sender_processes[i]);

			this->cells_to_receive[sender_processes[i]].push_back(global_ids_to_receive[i]);
			this->added_cells.insert(global_ids_to_receive[i]);
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
			this->requests[sender].push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#endif


			#else
			this->requests.push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
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
			this->removed_cells.insert(global_ids_to_send[i]);
		}

		Zoltan_LB_Free_Data(&global_ids_to_receive, &local_ids_to_receive, &sender_processes, &global_ids_to_send, &local_ids_to_send, &receiver_processes);

		// gather data to send
		for (int receiver = 0; receiver < this->comm.size(); receiver++) {

			// Zoltan returns cell migration info also for cells that don't chage their process
			if (receiver == this->comm.rank()) {
				continue;
			}

			if (cells_to_send.count(receiver) == 0) {
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

			if (cells_to_send.count(receiver) == 0) {
				// no data to send / receive
				continue;
			}

			int send_receive_tag = this->comm.rank() * this->comm.size() + receiver;

			// don't actually receive one cell at a time or anything, just use the correct request structures
			#ifdef DCCRG_SEND_SINGLE_CELLS

			#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
			this->load_balance_requests[receiver].push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#else
			this->requests[receiver].push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#endif

			#else
			this->requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
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

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->requests.begin(); process != this->requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
			process->second.clear();
		}
		this->requests.clear();
		#endif

		#else

		boost::mpi::wait_all(this->requests.begin(), this->requests.end());
		this->requests.clear();

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
		std::vector<uint64_t> temp_removed_cells(this->removed_cells.begin(), this->removed_cells.end());
		std::vector<std::vector<uint64_t> > all_removed_cells;
		all_gather(this->comm, temp_removed_cells, all_removed_cells);
		this->removed_cells.clear();
		// created cells on all processes
		std::vector<uint64_t> temp_added_cells(this->added_cells.begin(), this->added_cells.end());
		std::vector<std::vector<uint64_t> > all_added_cells;
		all_gather(this->comm, temp_added_cells, all_added_cells);
		this->added_cells.clear();

		// TODO: only delete those that aren't remote neighbours anymore
		this->remote_neighbours.clear();
		this->remote_cells_with_local_neighbours.clear();

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

		// create neighbour list for children that came to this process
		for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {

			if (*created_cell != this->get_child(*created_cell)) {
				continue;
			}

			this->update_neighbours(*created_cell);
		}

		// create neighbour_to list for children that came to this process
		for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {

			if (*created_cell != this->get_child(*created_cell)) {
				continue;
			}

			this->update_neighbours_to(*created_cell);
		}

		// if a child cell came to this process, update remote neighbour info of nearby cells
		for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {

			if (*created_cell != this->get_child(*created_cell)) {
				continue;
			}

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

				boost::unordered_set<uint64_t> temp_neighbours = this->get_neighbours_to(*created_cell);
				for (boost::unordered_set<uint64_t>::const_iterator neighbour = temp_neighbours.begin(); neighbour != temp_neighbours.end(); neighbour++) {
					if (this->cells.count(*neighbour)) {
						this->update_remote_neighbour_info(*neighbour);
					}
				}
			}
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
				this->requests[sender->first].push_back(MPI_Request());

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Irecv(this->remote_neighbours[*cell].at(), 1, data_type, sender->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->requests[sender->first].back()));
				#else
				MPI_Irecv(this->remote_neighbours[*cell].at(), UserData::size(), MPI_BYTE, sender->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->requests[sender->first].back()));
				#endif

				#else

				this->requests[sender->first].push_back(this->comm.irecv(sender->first, *cell % boost::mpi::environment::max_tag(), this->remote_neighbours[*cell]));	// FIXME: make sure message tags between two processes are unique

				#endif
			}
		}

		// post all sends
		for (boost::unordered_map<int, std::vector<uint64_t> >::const_iterator receiver = this->cells_to_send.begin(); receiver != this->cells_to_send.end(); receiver++) {

			for (std::vector<uint64_t>::const_iterator cell = receiver->second.begin(); cell != receiver->second.end(); cell++) {
				#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
				this->requests[receiver->first].push_back(MPI_Request());

				#ifdef DCCRG_USER_MPI_DATA_TYPE
				MPI_Isend(this->cells[*cell].at(), 1, data_type, receiver->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->requests[receiver->first].back()));
				#else
				MPI_Isend(this->cells[*cell].at(), UserData::size(), MPI_BYTE, receiver->first, *cell % boost::mpi::environment::max_tag(), this->comm, &(this->requests[receiver->first].back()));
				#endif

				#else

				this->requests[receiver->first].push_back(this->comm.isend(receiver->first, *cell % boost::mpi::environment::max_tag(), this->cells[*cell]));

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

			this->requests.push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
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

			this->requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
		}
		#endif // ifdef DCCRG_SEND_SINGLE_CELLS
	}


	/*!
	Waits until all neighbour data update transfers between processes have completed and incorporates that data
	Must be called simultaneously on all processes
	*/
	void wait_neighbour_data_update(void)
	{
		#ifdef DCCRG_SEND_SINGLE_CELLS

		#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator process = this->requests.begin(); process != this->requests.end(); process++) {

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
		this->requests.clear();

		#else

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->requests.begin(); process != this->requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}
		this->requests.clear();

		#endif

		#else

		boost::mpi::wait_all(this->requests.begin(), this->requests.end());
		this->requests.clear();
		this->outgoing_data.clear();

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
			return (std::vector<uint64_t>*) &(this->neighbours.at(cell));
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
			return &(this->neighbours_to[cell]);
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


	bool is_local(const uint64_t cell)
	{
		if (this->cell_process[cell] == this->comm.rank()) {
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
			exit(EXIT_FAILURE);
		}

		outfile.close();
	}


	/*!
	Creates all children of given cell (and possibly more due to induced refinement), overriding any unrefines of this cell (and possibly of its neighbours) before or after this call
	After refining / unrefining even one cell on any process stop_refining() must be called before doing anything else with the grid, except refining / unrefining
	Does nothing in any of the following cases:
		-given cell has already been refined (including induced refinement) and stop_refining() has not been called afterwards
		-given cell doesn't exist
		-given cell exists on another process
		-given cells children already exist
		-the created childrens' refinement level would exceed max_refinement_level
	The children are created on their parents process
	 */
	void refine_completely(const uint64_t cell)
	{
		if (this->cells.count(cell) == 0) {
			return;
		}

		if (this->get_refinement_level(cell) == this->max_refinement_level) {
			return;
		}

		if (cell != this->get_child(cell)) {
			// cell already has children
			return;
		}

		std::vector<uint64_t> children = this->get_all_children(cell);

		// don't refine the same cell (and possibly its neighbours) again
		if (this->added_cells.count(children[0]) > 0) {
			return;
		}

		this->added_cells.insert(children.begin(), children.end());

		// don't unrefine given cells or its siblings
		std::vector<uint64_t> peers = this->get_all_children(this->get_parent(cell));
		for (std::vector<uint64_t>::const_iterator peer = peers.begin(); peer != peers.end(); peer++) {
			this->removed_cells.erase(*peer);
		}

		// refine given cell's neighbours on this process if they are too large
		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[cell].begin(); neighbour != this->neighbours[cell].end(); neighbour++) {

			if (this->cells.count(*neighbour) == 0) {
				continue;
			}

			// stop unrefining the neighbour (and its parent's other children) if it is too large
			if (this->get_refinement_level(*neighbour) <= this->get_refinement_level(cell)) {
				std::vector<uint64_t> peers = this->get_all_children(this->get_parent(*neighbour));
				for (std::vector<uint64_t>::const_iterator peer = peers.begin(); peer != peers.end(); peer++) {
					this->removed_cells.erase(*peer);
				}
			}

			// induced refinement
			if (this->get_refinement_level(*neighbour) < this->get_refinement_level(cell)) {
				this->refine_completely(*neighbour);
			}
		}

		// refine cells on this process that consider given cell as neighbour if they are too large
		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[cell].begin(); neighbour != this->neighbours_to[cell].end(); neighbour++) {

			if (this->cells.count(*neighbour) == 0) {
				continue;
			}

			// stop unrefining the neighbour (and its parent's other children) if it is too large
			if (this->get_refinement_level(*neighbour) <= this->get_refinement_level(cell)) {
				std::vector<uint64_t> peers = this->get_all_children(this->get_parent(*neighbour));
				for (std::vector<uint64_t>::const_iterator peer = peers.begin(); peer != peers.end(); peer++) {
					this->removed_cells.erase(*peer);
				}
			}

			// induced refinement
			if (this->get_refinement_level(*neighbour) < this->get_refinement_level(cell)) {
				this->refine_completely(*neighbour);
			}
		}
	}

	/*!
	As refine_completely, but uses the smallest existing cell at given coordinates
	Does nothing in the same cases as refine_completely and additionally if the coordinate is outside of the grid
	*/
	void refine_completely_at(const double x, const double y, const double z)
	{
		uint64_t cell = this->get_smallest_cell_from_coordinate(x, y, z);
		if (cell == 0) {
			return;
		}

		this->refine_completely(cell);
	}


	/*!
	Removes the given cell and its siblings from the grid if none of their neighbours are smaller than the given cell
	After refining / unrefining even one cell on any process stop_refining() must be called before doing anything else with the grid, except refining / unrefining
	In case a cell would be both refined and unrefined it will be refined
	Does nothing in any of the following cases:
		-given cell has already been unrefined and stop_refining() has not been called yet
		-given cell doesn't exist
		-given cell exists on another process
		-given cell has children
		-given cells refinement level is 0
	*/
	void unrefine_completely(const uint64_t cell)
	{
		if (this->cells.count(cell) == 0) {
			return;
		}

		// don't unrefine the same cell again
		if (this->removed_cells.count(cell) > 0) {
			return;
		}

		if (this->get_refinement_level(cell) == 0) {
			return;
		}

		if (cell != this->get_child(cell)) {
			// cell already has children
			return;
		}

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[cell].begin(); neighbour != this->neighbours[cell].end(); neighbour++) {
			if (this->get_refinement_level(cell) < this->get_refinement_level(*neighbour)) {
				return;
			}
		}

		this->removed_cells.insert(cell);
	}


	/*!
	As unrefine_completely, but uses the smallest existing cell at given coordinates
	Does nothing in the same cases as unrefine_completely and additionally if the coordinate is outside of the grid
	*/
	void unrefine_completely_at(const double x, const double y, const double z)
	{
		uint64_t cell = this->get_smallest_cell_from_coordinate(x, y, z);
		if (cell == 0) {
			return;
		}

		this->unrefine_completely(cell);
	}


	/*!
	Informs other processes about cells created / deleted on this process by refining / unrefining
	Must be called simultaneously on all processes if even one of the has refined or unrefined even one cell, before doing anything else with the grid except refining / unrefining
	Returns all cells that were created on this process since the last call to this function
	*/
	std::vector<uint64_t> stop_refining(void)
	{
		this->comm.barrier();

		// tells the user which cells were created eventually on this process
		std::vector<uint64_t> new_cells;

		// first refine cells globally so that unrefines can be overridden locally
		while (true) {

			std::vector<std::vector<uint64_t> > all_added_cells;
			std::vector<uint64_t> temp_added_cells(this->added_cells.begin(), this->added_cells.end());
			all_gather(this->comm, temp_added_cells, all_added_cells);
			this->added_cells.clear();

			// continue until induced refinement across processes has stopped
			bool cells_created = false;
			for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
				if (all_added_cells[cell_creator].size() > 0) {
					cells_created = true;
					break;
				}
			}
			if (!cells_created) {
				break;
			}


			// add new cells to the cell to process mappings
			for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
				for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {

					this->cell_process[*created_cell] = cell_creator;

					// cells should be created on the same process as their parent
					assert(this->cell_process[this->get_parent(*created_cell)] == cell_creator);

					if (this->comm.rank() == cell_creator) {
						this->cells[*created_cell];
						// tell the user which cells were created on this process
						new_cells.push_back(*created_cell);
					}
				}
			}

			// only children should be created by refinement
			for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
				for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {
					assert(*created_cell != this->get_parent(*created_cell));
				}
			}

			// all cells whose neighbour and _to list will be updated
			boost::unordered_set<uint64_t> neighbour_lists_to_update;

			// locally created cells' neighbour lists
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {
				neighbour_lists_to_update.insert(*created_cell);
			}

			// locally created cells' parent's neighbours' neighbour lists and those that considered these their neighbours
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {

				uint64_t parent_of_created = this->get_parent(*created_cell);

				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[parent_of_created].begin(); neighbour != this->neighbours[parent_of_created].end(); neighbour++) {
					neighbour_lists_to_update.insert(*neighbour);
				}

				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[parent_of_created].begin(); neighbour != this->neighbours_to[parent_of_created].end(); neighbour++) {
					neighbour_lists_to_update.insert(*neighbour);
				}
			}

			// local cells that should be refined because their neighbours on other processes are too small
			boost::unordered_set<uint64_t> induced_refines;

			// local neighbours of remotely created cells
			for (boost::unordered_set<uint64_t>::const_iterator cell = this->cells_with_remote_neighbours.begin(); cell != this->cells_with_remote_neighbours.end(); cell++) {

				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*cell].begin(); neighbour != this->neighbours[*cell].end(); neighbour++) {

					if (*neighbour == this->get_child(*neighbour)) {
						// neighbour wasn't refined
						continue;
					}

					if (*cell == this->get_child(*cell)) {
						neighbour_lists_to_update.insert(*cell);
					// if cell was refined update its childrens' neighbour lists instead
					} else {
						std::vector<uint64_t> children = this->get_all_children(*cell);
						neighbour_lists_to_update.insert(children.begin(), children.end());
					}

					// remote refined cell doesn't count as a neighbour anymore
					this->remote_cells_with_local_neighbours.erase(*neighbour);
					this->remote_neighbours.erase(*neighbour);

					// induced refinement
					if (this->get_refinement_level(*neighbour) > this->get_refinement_level(*cell)) {
						induced_refines.insert(*cell);
					}
				}

				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[*cell].begin(); neighbour != this->neighbours_to[*cell].end(); neighbour++) {

					if (*neighbour == this->get_child(*neighbour)) {
						// neighbour wasn't refined
						continue;
					}

					if (*cell == this->get_child(*cell)) {
						neighbour_lists_to_update.insert(*cell);
					// if cell was refined update its childrens' neighbour lists instead
					} else {
						std::vector<uint64_t> children = this->get_all_children(*cell);
						neighbour_lists_to_update.insert(children.begin(), children.end());
					}

					// remote refined cell doesn't count as a neighbour anymore
					this->remote_cells_with_local_neighbours.erase(*neighbour);
					this->remote_neighbours.erase(*neighbour);

					// induced refinement
					if (this->get_refinement_level(*neighbour) > this->get_refinement_level(*cell)) {
						induced_refines.insert(*cell);
					}
				}
			}

			// update neighbour lists
			for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {

				if (this->cells.count(*cell) == 0) {
					continue;
				}

				if (*cell != this->get_child(*cell)) {
					continue;
				}

				// cell was just created if parent's neighbour list still exists
				if (this->neighbours.count(this->get_parent(*cell)) > 0) {
					this->neighbours[*cell] = this->get_neighbours_of(*cell, &(this->neighbours[this->get_parent(*cell)]));
				// cell's neighbourhood changed, start from its old neighbour list
				} else {
					this->neighbours[*cell] = this->get_neighbours_of(*cell, &(this->neighbours[*cell]));
				}
			}

			// update neighbour_to lists
			for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {
				if (this->cells.count(*cell) == 0) {
					continue;
				}

				if (*cell != this->get_child(*cell)) {
					continue;
				}

				// cell was just created if parent's neighbour_to list still exists
				if (this->neighbours_to.count(this->get_parent(*cell)) > 0) {
					this->neighbours_to[*cell] = this->get_neighbours_to(*cell, &(this->neighbours[this->get_parent(*cell)]), &(this->neighbours_to[this->get_parent(*cell)]));
				// cell's neighbourhood changed, start from its old neighbour_to lists
				} else {
					this->neighbours_to[*cell] = this->get_neighbours_to(*cell, &(this->neighbours[*cell]), &(this->neighbours_to[*cell]));
				}
				this->update_remote_neighbour_info(*cell);
			}

			// override obvious unrefines locally
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {
				this->removed_cells.erase(*created_cell);
				this->removed_cells.erase(this->get_parent(*created_cell));
			}

			// remove locally created cells' parent's neighbour info
			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[this->comm.rank()].begin(); created_cell != all_added_cells[this->comm.rank()].end(); created_cell++) {
				uint64_t parent_of_created = this->get_parent(*created_cell);
				this->neighbours.erase(parent_of_created);
				this->neighbours_to.erase(parent_of_created);
				this->cells_with_remote_neighbours.erase(parent_of_created);
			}

			// induce refines
			for (boost::unordered_set<uint64_t>::const_iterator cell = induced_refines.begin(); cell != induced_refines.end(); cell++) {
				this->refine_completely(*cell);
			}
		}

		// since the grid has been refined globally, finally really override local unrefines
		boost::unordered_set<uint64_t> final_removed_cells;
		for (boost::unordered_set<uint64_t>::const_iterator removed_cell = this->removed_cells.begin(); removed_cell != this->removed_cells.end(); removed_cell++) {

			bool unrefine = true;

			// override if any of the neighbours of the removed cell would be too small
			if (this->cells.count(*removed_cell) > 0) {

				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*removed_cell].begin(); neighbour != this->neighbours[*removed_cell].end(); neighbour++) {
					if (this->get_refinement_level(*removed_cell) < this->get_refinement_level(*neighbour)) {
						unrefine = false;
						break;
					}
				}

			} else {

				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(*removed_cell);
				for (std::vector<uint64_t>::const_iterator neighbour = temp_neighbours.begin(); neighbour != temp_neighbours.end(); neighbour++) {
					if (this->get_refinement_level(*removed_cell) < this->get_refinement_level(*neighbour)) {
						unrefine = false;
						break;
					}
				}
			}
			if (!unrefine) {
				continue;
			}

			// override if any of the neighbours of the removed cell's parent would be too small
			boost::unordered_set<uint64_t> neighbours_of_parent = this->get_neighbours_of_parent(*removed_cell);
			for (boost::unordered_set<uint64_t>::const_iterator neighbour_of_parent = neighbours_of_parent.begin(); neighbour_of_parent != neighbours_of_parent.end(); neighbour_of_parent++) {

				if (this->get_refinement_level(*removed_cell) < this->get_refinement_level(*neighbour_of_parent)) {
					unrefine = false;
					break;
				}
			}
			if (!unrefine) {
				continue;
			}

			final_removed_cells.insert(*removed_cell);
		}
		this->removed_cells.clear();


		/*
		Unrefine
		*/
		std::vector<std::vector<uint64_t> > all_removed_cells;
		std::vector<uint64_t> temp_removed_cells(final_removed_cells.begin(), final_removed_cells.end());
		all_gather(this->comm, temp_removed_cells, all_removed_cells);

		this->cells_to_receive.clear();
		this->cells_to_send.clear();

		// send user data of removed cells to their parent cell's process
		this->removed_cell_data.clear();

		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
			for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {

				int process_of_parent = this->cell_process[this->get_parent(*removed_cell)];

				std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(*removed_cell));
				for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {

					if (process_of_parent == this->cell_process[*sibling] && process_of_parent == this->comm.rank()) {
						// don't send to self, but still save user data
						this->removed_cell_data[*sibling] = this->cells[*sibling];
						continue;
					}

					if (this->comm.rank() == process_of_parent) {
						this->cells_to_receive[this->cell_process[*sibling]].push_back(*sibling);
					}

					if (this->comm.rank() == this->cell_process[*sibling]) {
						this->cells_to_send[process_of_parent].push_back(*sibling);
					}
				}
			}
		}

		// post all receives
		for (int sender = 0; sender < this->comm.size(); sender++) {

			if (this->cells_to_receive.count(sender) == 0) {
				continue;
			}

			int send_receive_tag = sender * this->comm.size() + this->comm.rank();

			// receive data in known order
			sort(this->cells_to_receive[sender].begin(), this->cells_to_receive[sender].end());

			#ifdef DCCRG_SEND_SINGLE_CELLS

			#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
			this->load_balance_requests[sender].push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#else
			this->requests[sender].push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#endif

			#else
			this->requests.push_back(this->comm.irecv(sender, send_receive_tag, this->incoming_data[sender]));
			#endif
		}

		// gather data to send
		for (int receiver = 0; receiver < this->comm.size(); receiver++) {

			if (this->cells_to_send.count(receiver) == 0) {
				continue;
			}

			// send data in known order
			sort(this->cells_to_send[receiver].begin(), this->cells_to_send[receiver].end());

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

			if (this->cells_to_send.count(receiver) == 0) {
				continue;
			}

			int send_receive_tag = this->comm.rank() * this->comm.size() + receiver;

			#ifdef DCCRG_SEND_SINGLE_CELLS

			#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
			this->load_balance_requests[receiver].push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#else
			this->requests[receiver].push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
			#endif

			#else
			this->requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
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
		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator process = this->requests.begin(); process != this->requests.end(); process++) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
			process->second.clear();
		}
		this->requests.clear();
		#endif

		#else

		boost::mpi::wait_all(this->requests.begin(), this->requests.end());
		this->requests.clear();

		#endif

		this->cells_to_send.clear();
		this->outgoing_data.clear();

		// incorporate received data
		for (typename boost::unordered_map<int, std::vector<UserData> >::const_iterator sender = this->incoming_data.begin(); sender != this->incoming_data.end(); sender++) {

			assert(this->incoming_data[sender->first].size() == this->cells_to_receive[sender->first].size());

			int i = 0;
			for (std::vector<uint64_t>::const_iterator cell = this->cells_to_receive[sender->first].begin(); cell != this->cells_to_receive[sender->first].end(); cell++, i++) {
				this->removed_cell_data[*cell] = this->incoming_data[sender->first][i];
			}
		}
		this->cells_to_receive.clear();
		this->incoming_data.clear();


		// all cells whose neighbour and _to list will be updated
		boost::unordered_set<uint64_t> neighbour_lists_to_update;

		// local cells whose data structures will be updated
		boost::unordered_set<uint64_t> all_locally_removed_cells;

		// remove cells from cell to process mappings
		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
			for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {

				// initial cells can't be removed
				assert(this->get_refinement_level(*removed_cell) > 0);

				if (this->cell_process.count(*removed_cell) == 0) {
					continue;
				}

				neighbour_lists_to_update.insert(this->get_parent(*removed_cell));

				// update neighbour and _to lists of all cells that are close enough to a removed cells parent...
				boost::unordered_set<uint64_t> neighbours_of_parent = this->get_neighbours_of_parent(*removed_cell);


				std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(*removed_cell));
				for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {

					if (this->cell_process[*sibling] == this->comm.rank()) {
						this->cells.erase(*sibling);
						all_locally_removed_cells.insert(*sibling);
					}

					this->cell_process.erase(*sibling);
				}

				// ...already here because next removed cell might have had these as a neighbour
				for (boost::unordered_set<uint64_t>::const_iterator neighbour_of_parent = neighbours_of_parent.begin(); neighbour_of_parent != neighbours_of_parent.end(); neighbour_of_parent++) {
					if (*neighbour_of_parent != this->get_child(*neighbour_of_parent)) {
						continue;
					}
					if (this->cells.count(*neighbour_of_parent) > 0) {
						this->update_neighbours(*neighbour_of_parent);
					}
				}
				for (boost::unordered_set<uint64_t>::const_iterator neighbour_of_parent = neighbours_of_parent.begin(); neighbour_of_parent != neighbours_of_parent.end(); neighbour_of_parent++) {
					if (*neighbour_of_parent != this->get_child(*neighbour_of_parent)) {
						continue;
					}
					if (this->cells.count(*neighbour_of_parent) > 0) {
						this->update_neighbours_to(*neighbour_of_parent);
					}
				}
			}
		}

		// update local neighbour and _to lists
		for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {

			if (this->cells.count(*cell) == 0) {
				continue;
			}

			if (*cell != this->get_child(*cell)) {
				continue;
			}

			this->update_neighbours(*cell);
		}
		for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {

			if (this->cells.count(*cell) == 0) {
				continue;
			}

			if (*cell != this->get_child(*cell)) {
				continue;
			}

			this->update_neighbours_to(*cell);
		}

		// update locally removed cells' data structures
		for (boost::unordered_set<uint64_t>::const_iterator removed_cell = all_locally_removed_cells.begin(); removed_cell != all_locally_removed_cells.end(); removed_cell++) {
			this->neighbours.erase(*removed_cell);
			this->neighbours_to.erase(*removed_cell);
			this->cells_with_remote_neighbours.erase(*removed_cell);
		}

		this->recalculate_neighbour_update_send_receive_lists();
		return new_cells;
	}


	/*!
	Returns cells that were removed by unrefinement whose parent is on this process
	Removed cells data is also on this process, but only until balance_load() is called
	*/
	std::vector<uint64_t> get_removed_cells(void) const
	{
		std::vector<uint64_t> unref_removed_cells;
		unref_removed_cells.reserve(this->removed_cell_data.size());

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->removed_cell_data.begin(); cell != this->removed_cell_data.end(); cell++) {
			unref_removed_cells.push_back(cell->first);
		}

		return unref_removed_cells;
	}


	/*!
	Given a cell that exists and has a parent returns the parent cell
	Returns the given cell if it doesn't have a parent or 0 if the cell doesn't exist
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

		uint64_t parent = get_cell_from_indices(this->get_x_index(cell), this->get_y_index(cell), this->get_z_index(cell), this->get_refinement_level(cell) - 1);
		if (this->cell_process.count(parent) > 0) {
			return parent;
		} else {
			return cell;
		}
	}

	/*!
	Returns the parent of given cell
	Returns the given cell if its refinement level == 0 or > maximum refinement level
	*/
	uint64_t get_parent_for_removed(const uint64_t cell) const
	{
		int refinement_level = this->get_refinement_level(cell);
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
	Returns the existing neighbours without children of given cell
	Returns nothing if the given cell has children
	*/
	std::vector<uint64_t> get_neighbours_of(const uint64_t id) const
	{
		std::vector<uint64_t> return_neighbours;

		if (id == 0
		|| id > this->max_cell_number
		|| this->cell_process.count(id) == 0) {
			return return_neighbours;
		}

		if (id != this->get_child(id)) {
			return return_neighbours;
		}

		const uint64_t x_index = this->get_x_index(id), y_index = this->get_y_index(id), z_index = this->get_z_index(id);

		// search neighbours in cells of the same size as the given cell
		const uint64_t size_in_indices = this->get_cell_size_in_indices(id);

		// the refinement level difference between neighbours is <= 2, only search for those cells
		const int refinement_level = this->get_refinement_level(id);
		int search_min_ref_level = refinement_level - 2, search_max_ref_level = refinement_level + 2;	// TODO: make these also const
		if (search_min_ref_level < 0) {
			search_min_ref_level = 0;
		}
		if (search_max_ref_level > this->max_refinement_level) {
			search_max_ref_level = this->max_refinement_level;
		}

		// the refinement level difference between neighbours is <= 2 because induced refines haven't propagated across processes when neighbours are updated, search index can increase by this amount at most
		uint64_t index_increase = size_in_indices / 4;
		if (index_increase == 0) {
			index_increase = 1;
		}

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours;

		// if neighbour_size == 0 just check the volume inside cells of the same size and that share a face with the given one
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
		unique_neighbours.erase(id);

		return_neighbours.insert(return_neighbours.end(), unique_neighbours.begin(), unique_neighbours.end());
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


	// Returns the cell of given refinement level at given indices even if it doesn't exist
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



private:

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

	// bookkeeping for Zoltan
	Zoltan_Struct* zoltan;
	// record whether Zoltan_LB_Partition is expected to fail (when the user selects NONE as the load balancing algorithm)
	bool no_load_balancing;

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

	// cells added to or removed from the grid on this process that haven't been communicated to other processes yet
	boost::unordered_set<uint64_t> added_cells, removed_cells;

	#ifdef DCCRG_SEND_SINGLE_CELLS	// user data is sent to another process one cell at a time
	// list of pending transfers between this process and the process as the key
	#ifdef DCCRG_CELL_DATA_SIZE_FROM_USER
	boost::unordered_map<int, std::vector<MPI_Request> > requests;
	boost::unordered_map<int, std::vector<boost::mpi::request> > load_balance_requests;	// TODO: get rid of this
	#else
	boost::unordered_map<int, std::vector<boost::mpi::request> > requests;
	#endif
	#else	// user data is packed into a vector which is sent to another process
	// pending neighbour data requests for this process
	std::vector<boost::mpi::request> requests;
	#endif

	// cells whose data has to be received / sent by this process from the process as the key
	boost::unordered_map<int, std::vector<uint64_t> > cells_to_receive;
	boost::unordered_map<int, std::vector<uint64_t> > cells_to_send;

	// storage for cells' user data that awaits transfer to or from this process
	boost::unordered_map<int, std::vector<UserData> > incoming_data, outgoing_data;

	// stores user data of cells that were removed when unrefining (possibly on another process)
	boost::unordered_map<uint64_t, UserData> removed_cell_data;


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
	Updates the neighbour list of given cell without children on this process
	*/
	void update_neighbours(const uint64_t cell)
	{
		assert(cell > 0 && cell <= this->max_cell_number);
		assert(this->cells.count(cell) > 0);
		assert(cell == this->get_child(cell));

		this->neighbours[cell] = this->get_neighbours_of(cell);
	}


	/*!
	Returns the existing neighbours without children of given cell from the given list of cells
	Returns nothing if the given cell has children
	Checks also for given cells siblings, in case it has been refined so they aren't in given cells
	*/
	std::vector<uint64_t> get_neighbours_of(const uint64_t cell, const std::vector<uint64_t>* cells) const
	{
		assert(cell > 0 && cell <= this->max_cell_number);
		assert(this->cell_process.count(cell) > 0);

		std::vector<uint64_t> return_neighbours;
		boost::unordered_set<uint64_t> unique_return_neighbours;

		if (cell != this->get_child(cell)) {
			return return_neighbours;
		}

		// if given cell has been refined, its siblings are not in given cells
		std::vector<uint64_t> siblings = this->get_all_children(this->get_parent(cell));
		for (std::vector<uint64_t>::const_iterator sibling = siblings.begin(); sibling != siblings.end(); sibling++) {

			if (cell == *sibling) {
				continue;
			}

			if (*sibling == this->get_child(*sibling)) {
				unique_return_neighbours.insert(*sibling);
			} else {
				// if sibling has children this cell wasn't refined and sibling is included in given cells
			}
		}

		// add neighbours from given list of cells
		for (std::vector<uint64_t>::const_iterator neighbour_candidate = cells->begin(); neighbour_candidate != cells->end(); neighbour_candidate++) {

			if (!this->is_neighbour(cell, *neighbour_candidate)) {
				continue;
			}

			if (*neighbour_candidate == this->get_child(*neighbour_candidate)) {
				unique_return_neighbours.insert(*neighbour_candidate);
			// check the children in case the cell's been refined
			} else {
				std::vector<uint64_t> children = this->get_all_children(*neighbour_candidate);
				assert(children.size() > 0);
				for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
					if (this->is_neighbour(cell, *child)) {
						unique_return_neighbours.insert(*child);
					}
				}
			}
		}

		return_neighbours.insert(return_neighbours.end(), unique_return_neighbours.begin(), unique_return_neighbours.end());

		#ifndef NDEBUG
		sort(return_neighbours.begin(), return_neighbours.end());
		std::vector<uint64_t> compare_neighbours = this->get_neighbours_of(cell);
		sort(compare_neighbours.begin(), compare_neighbours.end());

		if (return_neighbours.size() != compare_neighbours.size()) {
			std::cout << "Process " << this->comm.rank() << " neighbour counts for cell " << cell << " don't match: " << return_neighbours.size() << " (";
			for (std::vector<uint64_t>::const_iterator c = return_neighbours.begin(); c != return_neighbours.end(); c++) {
				std::cout << *c << " ";
			}
			std::cout << ") should be " << compare_neighbours.size() << " (";
			for (std::vector<uint64_t>::const_iterator c = compare_neighbours.begin(); c != compare_neighbours.end(); c++) {
				std::cout << *c << " ";
			}
			std::cout << ")" << std::endl;
			assert(0);
		}
		#endif

		return return_neighbours;
	}


	/*!
	Updates the neighbours_to structure of given cell on this process without children and remote neighbour info
	*/
	void update_neighbours_to(const uint64_t cell)
	{
		assert(this->cells.count(cell) > 0);
		assert(cell == this->get_child(cell));

		// get cells which consider given cell as neighbour, but aren't neighbours of given cell
		this->neighbours_to[cell].clear();
		boost::unordered_set<uint64_t> temp_neighbours_to = this->get_neighbours_to(cell);

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[cell].begin(); neighbour != this->neighbours[cell].end(); neighbour++) {
			temp_neighbours_to.erase(*neighbour);
		}

		for (boost::unordered_set<uint64_t>::const_iterator neighbour_to = temp_neighbours_to.begin(); neighbour_to != temp_neighbours_to.end(); neighbour_to++) {
			this->neighbours_to[cell].push_back(*neighbour_to);
		}

		this->update_remote_neighbour_info(cell);
	}


	/*!
	Returns the existing cells without children that have the given cell as a neighbour
	Returns nothing if the given cell has children
	Doesn't update neighbour lists
	*/
	boost::unordered_set<uint64_t> get_neighbours_to(const uint64_t id) const
	{
		assert(this->cell_process.count(id) > 0);

		boost::unordered_set<uint64_t> return_neighbours;

		if (id != this->get_child(id)) {
			return return_neighbours;
		}

		/*
		The largest cell that considers this as a neighbour has a refinement level of no less than this cells refinement level -2 if this cell has been refined
		All of those can be found by searching the neighbours of this cell's neighbours' neighbours
		*/

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours_to;

		// first round of neighbour search, use existing neighbour lists if available
		boost::unordered_set<uint64_t> neighbours1, neighbours2;
		if (this->cells.count(id) > 0) {

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(id).begin(); neighbour != this->neighbours.at(id).end(); neighbour++) {
				assert(*neighbour == this->get_child(*neighbour));
				neighbours1.insert(*neighbour);
				unique_neighbours_to.insert(*neighbour);
			}

		} else {
			std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(id);
			neighbours1.insert(temp_neighbours.begin(), temp_neighbours.end());
			unique_neighbours_to.insert(temp_neighbours.begin(), temp_neighbours.end());
		}

		// second round of neighbour search...
		for (boost::unordered_set<uint64_t>::const_iterator neighbour1 = neighbours1.begin(); neighbour1 != neighbours1.end(); neighbour1++) {

			if (this->cells.count(*neighbour1) > 0) {
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(*neighbour1).begin(); neighbour != this->neighbours.at(*neighbour1).end(); neighbour++) {
					assert(*neighbour == this->get_child(*neighbour));
					neighbours2.insert(*neighbour);
					unique_neighbours_to.insert(*neighbour);
				}
			} else {
				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(*neighbour1);
				neighbours2.insert(temp_neighbours.begin(), temp_neighbours.end());
				unique_neighbours_to.insert(temp_neighbours.begin(), temp_neighbours.end());
			}
		}

		// third round of neighbour search...
		for (boost::unordered_set<uint64_t>::const_iterator neighbour2 = neighbours2.begin(); neighbour2 != neighbours2.end(); neighbour2++) {

			if (this->cells.count(*neighbour2) > 0) {
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(*neighbour2).begin(); neighbour != this->neighbours.at(*neighbour2).end(); neighbour++) {
					assert(*neighbour == this->get_child(*neighbour));
					unique_neighbours_to.insert(*neighbour);
				}
			} else {
				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(*neighbour2);
				unique_neighbours_to.insert(temp_neighbours.begin(), temp_neighbours.end());
			}
		}
		unique_neighbours_to.erase(id);

		// return only cells which consider given cell as a neighbour
		for (boost::unordered_set<uint64_t>::const_iterator neighbour = unique_neighbours_to.begin(); neighbour != unique_neighbours_to.end(); neighbour++) {
			if (this->is_neighbour(*neighbour, id)) {
				return_neighbours.insert(*neighbour);
			}
		}

		#ifndef NDEBUG
		// check that wrong cells weren't excluded / included
		uint64_t x_index = this->get_x_index(id), y_index = this->get_y_index(id), z_index = this->get_z_index(id);
		for (boost::unordered_set<uint64_t>::const_iterator neighbour = unique_neighbours_to.begin(); neighbour != unique_neighbours_to.end(); neighbour++) {

			uint64_t neighbour_x_index = this->get_x_index(*neighbour), neighbour_y_index = this->get_y_index(*neighbour), neighbour_z_index = this->get_z_index(*neighbour);
			uint64_t neighbour_size = this->get_cell_size_in_indices(*neighbour);

			std::vector<uint64_t> temp_neighbours1 = this->get_neighbours_of(*neighbour);
			boost::unordered_set<uint64_t> temp_neighbours2;
			temp_neighbours2.insert(temp_neighbours1.begin(), temp_neighbours1.end());

			if (temp_neighbours2.count(id) == 0 && return_neighbours.count(*neighbour) > 0) {
				std::cout << "Process " << this->comm.rank() << ": cell " << *neighbour << " considers cell " << id << " as neighbour but it shouldn't" << std::endl;
				std::cout << "Cell " << id << " size " << this->get_cell_size_in_indices(id) << ", indices: " << x_index << " " << y_index << " " << z_index << "; Cell " << *neighbour << " indices: " << neighbour_x_index << " " << neighbour_y_index << " " << neighbour_z_index << ", size in indices " << neighbour_size << std::endl;
				assert(false);
			}
			assert(!(temp_neighbours2.count(id) == 0 && return_neighbours.count(*neighbour) > 0));
			if (temp_neighbours2.count(id) > 0 && return_neighbours.count(*neighbour) == 0) {
				std::cout << "Process " << this->comm.rank() << ": cell " << *neighbour << " doesn't consider cell " << id << " as neighbour but it should" << std::endl;
				std::cout << "Cell " << id << " size " << this->get_cell_size_in_indices(id) << ",  indices: " << x_index << " " << y_index << " " << z_index << "; Cell " << *neighbour << " indices: " << neighbour_x_index << " " << neighbour_y_index << " " << neighbour_z_index << ", size in indices " << neighbour_size << std::endl;
				assert(false);
			}
			assert(!(temp_neighbours2.count(id) > 0 && return_neighbours.count(*neighbour) == 0));
		}
		#endif

		return return_neighbours;
	}


	/*!
	Returns the existing cells without children from the given lists that have the given cell as a neighbour
	Returns nothing if the given cell has children
	Doesn't update neighbour lists
	*/
	std::vector<uint64_t> get_neighbours_to(const uint64_t cell, const std::vector<uint64_t>* cells1, const std::vector<uint64_t>* cells2) const
	{
		assert(this->cell_process.count(cell) > 0);

		std::vector<uint64_t> return_neighbours;

		if (cell != this->get_child(cell)) {
			return return_neighbours;
		}

		for (std::vector<uint64_t>::const_iterator candidate = cells1->begin(); candidate != cells1->end(); candidate++) {
			if (this->is_neighbour(cell, *candidate)) {
				continue;
			}

			if (!this->is_neighbour(*candidate, cell)) {
				continue;
			}

			if (*candidate == this->get_child(*candidate)) {
				return_neighbours.push_back(*candidate);
			// check the children in case the cell's been refined
			} else {
				std::vector<uint64_t> children = this->get_all_children(*candidate);
				assert(children.size() > 0);
				for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
					if (this->is_neighbour(*child, cell) && !this->is_neighbour(cell, *child)) {
						return_neighbours.push_back(*child);
					}
				}
			}
		}

		for (std::vector<uint64_t>::const_iterator candidate = cells2->begin(); candidate != cells2->end(); candidate++) {
			if (this->is_neighbour(cell, *candidate)) {
				continue;
			}

			if (!this->is_neighbour(*candidate, cell)) {
				continue;
			}

			if (*candidate == this->get_child(*candidate)) {
				return_neighbours.push_back(*candidate);
			// check the children in case the cell's been refined
			} else {
				std::vector<uint64_t> children = this->get_all_children(*candidate);
				assert(children.size() > 0);
				for (std::vector<uint64_t>::const_iterator child = children.begin(); child != children.end(); child++) {
					if (this->is_neighbour(*child, cell) && !this->is_neighbour(cell, *child)) {
						return_neighbours.push_back(*child);
					}
				}
			}
		}

		#ifndef NDEBUG
		boost::unordered_set<uint64_t> compare_neighbours = this->get_neighbours_to(cell);
		// remove neighbours of given cell from compare neighbours
		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(cell).begin(); neighbour != this->neighbours.at(cell).end(); neighbour++) {
			compare_neighbours.erase(*neighbour);
		}

		if (return_neighbours.size() != compare_neighbours.size()) {
			std::cout << "Process " << this->comm.rank() << " neighbours_to counts for cell " << cell << " don't match: " << return_neighbours.size() << " (";
			for (std::vector<uint64_t>::const_iterator c = return_neighbours.begin(); c != return_neighbours.end(); c++) {
				std::cout << *c << " ";
			}
			std::cout << ") should be " << compare_neighbours.size() << " (";
			for (boost::unordered_set<uint64_t>::const_iterator c = compare_neighbours.begin(); c != compare_neighbours.end(); c++) {
				std::cout << *c << " ";
			}
			std::cout << ")" << std::endl;
			assert(0);
		}
		#endif

		return return_neighbours;
	}


	/*!
	Updates the remote neighbour info of given cell on this process without children
	Doesn't update neighbour or neighbour_to lists
	*/
	void update_remote_neighbour_info(const uint64_t cell)
	{
		assert(this->cells.count(cell) > 0);
		assert(cell == this->get_child(cell));

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
	*/
	boost::unordered_set<uint64_t> get_neighbours_of_parent(const uint64_t cell) const
	{
		assert(this->cell_process.count(cell) > 0);

		boost::unordered_set<uint64_t> return_neighbours;

		if (this->get_refinement_level(cell) == 0) {
			// given cells doesn't have a parent
			if (this->cells.count(cell) > 0) {
				return_neighbours.insert(this->neighbours.at(cell).begin(), this->neighbours.at(cell).end());
			} else {
				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(cell);
				return_neighbours.insert(temp_neighbours.begin(), temp_neighbours.end());
			}
			return return_neighbours;
		}

		uint64_t parent = this->get_parent(cell);

		// search for given cell's siblings' neighbours of neighbours of... until a neighbour isn't really a neighbour of given cell's parent
		boost::unordered_set<uint64_t> new_neighbours;
		std::vector<uint64_t> peers = this->get_all_children(parent);
		for (std::vector<uint64_t>::const_iterator peer = peers.begin(); peer != peers.end(); peer++) {

			// try to use existing neighbour list
			if (this->cells.count(*peer) > 0) {
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours.at(*peer).begin(); neighbour != this->neighbours.at(*peer).end(); neighbour++) {
					if (this->is_neighbour(parent, *neighbour)) {

						if (this->cell_process.count(*neighbour) == 0) {
							std::cout << "outoo! " << cell << " " << parent << " " << *peer << " ei ole: " << *neighbour << std::endl;
							exit(EXIT_FAILURE);
						}

						new_neighbours.insert(*neighbour);
					}
				}
			} else {
				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(*peer);
				for (std::vector<uint64_t>::const_iterator neighbour = temp_neighbours.begin(); neighbour != temp_neighbours.end(); neighbour++) {
					if (this->is_neighbour(parent, *neighbour)) {
						new_neighbours.insert(*neighbour);
					}
				}
			}
		}

		// continue searching until all neighbours have been found
		while (new_neighbours.size() > 0) {

			// stop when a neighbour isn't really a neighbour to given cell's parent
			boost::unordered_set<uint64_t> search_neighbours;
			for (boost::unordered_set<uint64_t>::const_iterator new_neighbour = new_neighbours.begin(); new_neighbour != new_neighbours.end(); new_neighbour++) {

				if (!this->is_neighbour(parent, *new_neighbour)) {
					continue;
				}

				// stop also if neighbour has already been considered
				if (return_neighbours.count(*new_neighbour) > 0) {
					continue;
				}

				search_neighbours.insert(*new_neighbour);
			}
			new_neighbours.clear();

			// get new neighbours
			for (boost::unordered_set<uint64_t>::const_iterator search_neighbour = search_neighbours.begin(); search_neighbour != search_neighbours.end(); search_neighbour++) {

				return_neighbours.insert(*search_neighbour);

				// try to use existing neighbour list
				if (this->cells.count(*search_neighbour) > 0) {
					for (std::vector<uint64_t>::const_iterator neigh_of_new_neigh = this->neighbours.at(*search_neighbour).begin(); neigh_of_new_neigh != this->neighbours.at(*search_neighbour).end(); neigh_of_new_neigh++)  {
						if (this->is_neighbour(parent, *neigh_of_new_neigh) && return_neighbours.count(*neigh_of_new_neigh) == 0) {
							new_neighbours.insert(*neigh_of_new_neigh);
						}
					}
				} else {
					std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(*search_neighbour);
					for (std::vector<uint64_t>::const_iterator neigh_of_new_neigh = temp_neighbours.begin(); neigh_of_new_neigh != temp_neighbours.end(); neigh_of_new_neigh++)  {
						if (this->is_neighbour(parent, *neigh_of_new_neigh) && return_neighbours.count(*neigh_of_new_neigh) == 0) {
							new_neighbours.insert(*neigh_of_new_neigh);
						}
					}
				}
			}
		}

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

		// given cell cannot have children
		if (get_refinement_level(cell) == this->max_refinement_level) {
			return cell;
		}

		uint64_t child = get_cell_from_indices(this->get_x_index(cell), this->get_y_index(cell), this->get_z_index(cell), this->get_refinement_level(cell) + 1);
		if (this->cell_process.count(child) > 0) {
			return child;
		} else {
			return cell;
		}
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
	Returns the smallest cell at given indices between given refinement levels inclusive.
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

		// use binary search recursively (assumes that a cell refines to max. 8 children)
		if (this->cell_process.count(id) == 0) {
			// doesn't exist, search the bin of smaller refinement_level values
			if (average_refinement_level > minimum_refinement_level) {
				return this->get_cell_from_indices(x_index, y_index, z_index, minimum_refinement_level, average_refinement_level - 1);
			} else {
				// nothing left to search
				return 0;
			}
		} else {
			// does exist, search the bin of larger refinement_level values
			if (average_refinement_level < maximum_refinement_level) {
				uint64_t larger_refinement_value_cell = this->get_cell_from_indices(x_index, y_index, z_index, average_refinement_level + 1, maximum_refinement_level);
				if (larger_refinement_value_cell > 0) {
					return larger_refinement_value_cell;
				} else {
					// current cell has the largest refinement value at given indices
					return id;
				}
			} else {
				// nothing left to search
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
		assert(id > 0);
		assert(this->cell_process.count(id) > 0);

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
			for (std::vector<uint64_t>::const_iterator neighbour = dccrg_instance->neighbours[cell->first].begin(); neighbour != dccrg_instance->neighbours[cell->first].end(); neighbour++) {
				(*number_of_connections)++;
			}
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

};

#endif
