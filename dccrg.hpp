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


template <class UserData> class dccrg
{

public:

	/*
	Creates a new grid where each cell stores one instance of given UserData
	Zoltan_Initialize must be called before calling this constructor

	comm: the grid will span all the processes in the communicator comm

	load_balancing_method:
		The method that Zoltan will use for load balancing given as a string
		Currently supported methods: NONE, BLOCK, RANDOM, RCB, RIB and HSFC

	x_start, y_start, z_start: the starting corner of the grid
	cell_size: the size of each unrefined cell in every direction

	x_length, y_length, z_length: the number of cells in the grid in x, y and z direction

	neighbourhood_size:
		Determines which cells are considered neighbours.
		When calculating the neighbours of a given cell a cube of length neighbourhood_size + 1 in every direction is considered, centered at the cell for which neighbours are being calculated.
		The unit lenght of the cube is the cell for which neighbours are being calculated.
		TODO: If neighbourhood_size == 0, only cells that share a face are considered.

	maximum_refinement_level:
		The maximum number of times an unrefined cell can be refined (replacing it with 8 smaller cells)
		Optional: if not given it is maximized based on the grids initial size
	 */
	dccrg(boost::mpi::communicator comm, const char* load_balancing_method, const double x_start, const double y_start, const double z_start, const double cell_size, const unsigned int x_length, const unsigned int y_length, const unsigned int z_length, const unsigned int neighbourhood_size, const int maximum_refinement_level = -1)
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
		Zoltan_Set_Param(this->zoltan, "OBJ_WEIGHT_DIM", "0");
		Zoltan_Set_Param(this->zoltan, "AUTO_MIGRATE", "1");
		// try to minimize moving of data between processes
		Zoltan_Set_Param(this->zoltan, "REMAP", "1");
		// RTFM...
		Zoltan_Set_Param(this->zoltan, "MIGRATE_ONLY_PROC_CHANGES", "1");

		// set the grids callback functions in Zoltan
		Zoltan_Set_Num_Obj_Fn(this->zoltan, &dccrg<UserData>::get_number_of_cells, this);
		Zoltan_Set_Obj_List_Fn(this->zoltan, &dccrg<UserData>::get_cell_list, this);
		Zoltan_Set_Obj_Size_Fn(this->zoltan, &dccrg<UserData>::user_data_size, NULL);
		Zoltan_Set_Pack_Obj_Fn(this->zoltan, &dccrg<UserData>::pack_user_data, this);
		Zoltan_Set_Unpack_Obj_Fn(this->zoltan, &dccrg<UserData>::unpack_user_data, this);
		Zoltan_Set_Num_Geom_Fn(this->zoltan, &dccrg<UserData>::get_grid_dimensionality, NULL);
		Zoltan_Set_Geom_Fn(this->zoltan, &dccrg<UserData>::fill_with_cell_coordinates, this);


		// Set grid parameters
		this->x_start = x_start;
		this->y_start = y_start;
		this->z_start = z_start;

		if (cell_size <= 0) {
			std::cerr << "Cell size must be > 0" << std::endl;
			exit(EXIT_FAILURE);
		}
		this->cell_size = cell_size;

		if (x_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the x direction" << std::endl;
			exit(EXIT_FAILURE);
		}
		this->x_length = x_length;

		if (y_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the y direction" << std::endl;
			exit(EXIT_FAILURE);
		}
		this->y_length = y_length;

		if (z_length == 0) {
			std::cerr << "Length of the grid in cells must be > 0 in the z direction" << std::endl;
			exit(EXIT_FAILURE);
		}
		this->z_length = z_length;

		if (neighbourhood_size == 0) {
			std::cerr << "Neighbour stencil size has to be > 0" << std::endl;
			exit(EXIT_FAILURE);
		}
		this->neighbourhood_size = neighbourhood_size;

		// get the maximum refinement level based on the size of the grid when using uint64_t for cell ids
		double max_id = pow(2, 64) - 1, last_id = x_length * y_length * z_length;
		int refinement_level = 0;
		while (last_id / max_id < 1) {
			refinement_level++;
			last_id += x_length * y_length * z_length * pow(8, refinement_level);
		}
		refinement_level--;

		// grid is too large even without refinement
		if (refinement_level < 0) {
			std::cerr << "Given grid would contain more than 2^64 - 1 unrefined cells" << std::endl;
			exit(EXIT_FAILURE);
		}


		if (maximum_refinement_level > refinement_level) {

			std::cerr << "Given max_refinement_level (" << maximum_refinement_level << ") is too large: " << "x_length * y_length * z_length * 8^max_refinement_level / (2^64 - 1) == " << x_length * y_length * z_length * pow(8, maximum_refinement_level) / max_id << " but must be < 1" << std::endl;
			exit(EXIT_FAILURE);

		} else if (maximum_refinement_level < 0) {
			this->max_refinement_level = refinement_level;
		} else {
			this->max_refinement_level = maximum_refinement_level;
		}

		// the number of the last cell at maximum refinement level
		uint64_t id = 0;
		for (refinement_level = 0; refinement_level <= this->max_refinement_level; refinement_level++) {
			id += x_length * y_length * z_length * pow(8, refinement_level);
		}
		this->max_cell_number = id;

		// create unrefined cells on process 0
		for (uint64_t id = 1; id <= x_length * y_length * z_length; id++) {
			if (this->comm.rank() == 0) {
				this->cells[id];
			}
			this->cell_process[id] = 0;
		}

		// update neighbour lists of created cells
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
			this->update_neighbours(cell->first);
		}

		// update neighbour_to lists of created cells
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
			this->update_remote_neighbour_info(cell->first);
		}
	}


	/*
	Returns all cells on this process that don't have children (e.g. leaf cells)
	*/
	std::vector<uint64_t> get_cells(void)
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


	/*
	Returns all cells on this process that don't have children (e.g. leaf cells) and don't have neighbours on other processes
	*/
	std::vector<uint64_t> get_cells_with_local_neighbours(void)
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
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[cell->first].begin(); neighbour != this->neighbours[cell->first].end(); neighbour++) {
				if (this->cell_process[*neighbour] != this->comm.rank()) {
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


	/*
	Returns all cells on this process that don't have children (e.g. leaf cells) and have at least one neighbour on another processes
	*/
	std::vector<uint64_t> get_cells_with_remote_neighbour(void)
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
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[cell->first].begin(); neighbour != this->neighbours[cell->first].end(); neighbour++) {
				if (this->cell_process[*neighbour] != this->comm.rank()) {
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


	/*
	Returns a pointer to the user supplied data of given cell
	Return NULL if the given cell isn't on this process and if the given cell isn't a neighbour of any cell on this process
	*/
	UserData* operator [] (const uint64_t cell)
	{
		if (this->cells.count(cell) > 0) {
			return &(this->cells[cell]);
		} else if (this->remote_neighbours.count(cell) > 0) {
			return &(this->remote_neighbours[cell]);
		} else {
			return NULL;
		}
	}


	/*
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

		// update cell to process mappings
		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
			for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {
				assert(this->cell_process[*removed_cell] == cell_remover);
			}
		}
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

				if (this->remote_neighbours.count(*created_cell) == 0) {
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
	}


	/*
	Updates the user data of those cells that have at least one neighbour or are considered as a neighbour of a cell on another process
	Must be called simultaneously on all processes
	*/
	void update_remote_neighbour_data(void)
	{
		this->start_remote_neighbour_data_update();
		this->wait_neighbour_data_update();
	}


	/*
	Starts the update of neighbour data between processes and returns before (probably) it has completed
	Must be called simultaneously on all processes
	*/
	void start_remote_neighbour_data_update(void)
	{
		this->comm.barrier();

		// processes and the cells for which data has to be sent by this process
		boost::unordered_map<int, boost::unordered_set<uint64_t> > cells_to_send;

		// calculate where this process has to send and from where to receive cell data
		for (boost::unordered_set<uint64_t>::const_iterator cell = this->cells_with_remote_neighbours.begin(); cell != this->cells_with_remote_neighbours.end(); cell++) {

			assert(*cell == this->get_child(*cell));

			int current_process = this->comm.rank();

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*cell].begin(); neighbour != this->neighbours[*cell].end(); neighbour++) {
				if (this->cell_process[*neighbour] != current_process) {
					// *neighbours process has to send *neighbours cell data to current_process
					this->cells_to_receive[this->cell_process[*neighbour]].insert(*neighbour);
					// current process has to send currents cell data to neighbour
					cells_to_send[this->cell_process[*neighbour]].insert(*cell);
				}
			}

			// also cells that have this one as neighbour
			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours_to[*cell].begin(); neighbour != this->neighbours_to[*cell].end(); neighbour++) {
				if (this->cell_process[*neighbour] != current_process) {
					// *neighbours process has to send *neighbours cell data to current_process
					this->cells_to_receive[this->cell_process[*neighbour]].insert(*neighbour);
					// current process has to send currents cell data to neighbour
					cells_to_send[this->cell_process[*neighbour]].insert(*cell);
				}
			}
		}

		// post all receives
		for (int sender = 0; sender < this->comm.size(); sender++) {

			if (sender == this->comm.rank()) {
				continue;
			}

			if (cells_to_receive.count(sender) == 0) {
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

			if (cells_to_send.count(receiver) == 0) {
				// no data to send / receive
				continue;
			}

			// send and receive cell data in an order known in advance
			std::vector<uint64_t> current_cells_to_send(cells_to_send[receiver].begin(), cells_to_send[receiver].end());
			sort(current_cells_to_send.begin(), current_cells_to_send.end());

			// construct the outgoing data vector
			for (std::vector<uint64_t>::const_iterator cell = current_cells_to_send.begin(); cell != current_cells_to_send.end(); cell++) {
				UserData* user_data = (*this)[*cell];
				assert(user_data != NULL);
				this->outgoing_data[receiver].push_back(*user_data);
			}
			assert(this->outgoing_data[receiver].size() == current_cells_to_send.size());
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

			this->requests.push_back(this->comm.isend(receiver, send_receive_tag, this->outgoing_data[receiver]));
		}
	}


	/*
	Waits until all neighbour data update transfers between processes have completed and incorporates that data
	Must be called simultaneously on all processes
	*/
	void wait_neighbour_data_update(void)
	{
		boost::mpi::wait_all(this->requests.begin(), this->requests.end());
		this->outgoing_data.clear();
		this->requests.clear();

		// incorporate received data
		for (typename boost::unordered_map<int, std::vector<UserData> >::const_iterator sender = incoming_data.begin(); sender != incoming_data.end(); sender++) {

			std::vector<uint64_t> current_cells_to_receive(this->cells_to_receive[sender->first].begin(), this->cells_to_receive[sender->first].end());
			sort(current_cells_to_receive.begin(), current_cells_to_receive.end());
			assert(this->incoming_data[sender->first].size() == current_cells_to_receive.size());

			int i = 0;
			for (std::vector<uint64_t>::const_iterator cell = current_cells_to_receive.begin(); cell != current_cells_to_receive.end(); cell++, i++) {
				this->remote_neighbours[*cell] = this->incoming_data[sender->first][i];
			}
		}
		this->cells_to_receive.clear();
		this->incoming_data.clear();
	}


	/*
	Returns a pointer to the neighbours (some of which might be on another process) of given cell
	Returns NULL if given cell doesn't exist or is on another process
	*/
	const std::vector<uint64_t>* get_neighbours(const uint64_t cell)
	{
		if (this->cells.count(cell) > 0) {
			return &(this->neighbours[cell]);
		} else {
			return NULL;
		}
	}

	/*
	Returns a pointer to the cells that consider given cell as a neighbour but that are not neighbours of given cell
	Some cell might be on another process or the structure might be empty
	Returns NULL if given cell doesn't exist or is on another process
	*/
	const std::vector<uint64_t>* get_neighbours2(const uint64_t cell)
	{
		if (this->cells.count(cell) > 0) {
			return &(this->neighbours_to[cell]);
		} else {
			return NULL;
		}
	}


	/*
	Returns the given cells neighbours that are on another process
	Returns nothing if given cell doesn't exist or is on another process or doesn't have remote neighbours
	*/
	std::vector<uint64_t> get_remote_neighbours(const uint64_t id)
	{
		std::vector<uint64_t> remote_neighbours;

		if (this->cells.count(id) == 0) {
			return remote_neighbours;
		}

		for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[id].begin(); neighbour != this->neighbours[id].end(); neighbour++) {
			if (this->cell_process[*neighbour] != this->comm.rank()) {
				remote_neighbours.push_back(*neighbour);
			}
		}

		return remote_neighbours;
	}


	/*
	Returns the maximum possible refinement level of any cell in the grid (0 means unrefined)
	*/
	int get_max_refinement_level(void)
	{
		return this->max_refinement_level;
	}


	// The following return the x, y or z coordinate of given cell regardless of whether it exists
	double get_cell_x(const uint64_t cell)
	{
		return this->x_start + this->get_x_index(cell) * this->cell_size / int(pow(2, this->max_refinement_level)) + this->get_cell_size(cell) / 2;
	}
	double get_cell_y(const uint64_t cell)
	{
		return this->y_start + this->get_y_index(cell) * this->cell_size / int(pow(2, this->max_refinement_level)) + this->get_cell_size(cell) / 2;
	}
	double get_cell_z(const uint64_t cell)
	{
		return this->z_start + this->get_z_index(cell) * this->cell_size / int(pow(2, this->max_refinement_level)) + this->get_cell_size(cell) / 2;
	}


	// Returns the length of given cell regardless of whether it exists
	double get_cell_size(const uint64_t cell)
	{
		return this->cell_size / int(pow(2, this->get_refinement_level(cell)));
	}


	// Returns the refinement level of given cell, even if it doesn't exist
	int get_refinement_level(uint64_t id)
	{
		assert(id > 0);
		assert(id <= this->max_cell_number);

		int refinement_level;
		for (refinement_level = 0; refinement_level < this->max_refinement_level; refinement_level++) {
			if (id <= (this->x_length * this->y_length * this->z_length * int(pow(8, refinement_level)))) {
				break;
			}
			// substract ids of larger cells
			id -= this->x_length * this->y_length * this->z_length * int(pow(8, refinement_level));
		}

		return refinement_level;
	}


	/*
	Writes the cells on this process into a vtk file with given name in ASCII format
	The cells are written in ascending order
	Must be called simultaneously on all processes
	*/
	void write_vtk_file(const char* file_name)
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

			double x = get_cell_x(leaf_cells[i]), y = get_cell_y(leaf_cells[i]), z = get_cell_z(leaf_cells[i]);
			// make the cells a little smaller than in theory so that the data at every face is well defined and VisIt plots it correctly
			for (double z_offset = -get_cell_size(leaf_cells[i]) / 2; z_offset < get_cell_size(leaf_cells[i]); z_offset += get_cell_size(leaf_cells[i])) {
				for (double y_offset = -get_cell_size(leaf_cells[i]) / 2; y_offset < get_cell_size(leaf_cells[i]); y_offset += get_cell_size(leaf_cells[i])) {
					for (double x_offset = -get_cell_size(leaf_cells[i]) / 2; x_offset < get_cell_size(leaf_cells[i]); x_offset += get_cell_size(leaf_cells[i])) {
						outfile << x + x_offset << " " << y + y_offset << " " << z + z_offset << std::endl;
					}
				}
			}
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
			exit(1);
		}

		outfile.close();
	}


	/*
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

		// don't unrefine this or other children of this cell's parent
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

	/*
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


	/*
	Removes the given cell and all other children of the given cells parent (and possibly more cells due to induced refinement)
	After refining / unrefining even one cell on any process stop_refining() must be called before doing anything else with the grid, except refining / unrefining
	In case any cell would be both refined and unrefined it will be refined
	Does nothing in any of the following cases:
		-given cell or any other child of given cells parent has already been unrefined and stop_refining() has not been called afterwards
		-given cell doesn't exist
		-given cell exists on another process
		-given cell has children
		-given cells refinement level is 0
		-given cell (or cells close enough due to induced refinement) has or will be refined before calling stop_refining()
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

		this->removed_cells.insert(cell);
	}


	/*
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


	/*
	Informs other processes about cells created / deleted on this process by refining / unrefining
	Must be called simultaneously on all processes if even one of the has refined or unrefined even one cell, before doing anything else with the grid except refining / unrefining
	Returns all cells that were created on this process since the last call to this function
	*/
	std::vector<uint64_t> stop_refining(void)
	{
		this->comm.barrier();

		// tells the user which cells were created eventually
		std::vector<uint64_t> new_cells;

		// first refine cells globally so that unrefines can be overridden locally
		while (true) {

			std::vector<std::vector<uint64_t> > all_added_cells;
			std::vector<uint64_t> temp_added_cells(this->added_cells.begin(), this->added_cells.end());
			all_gather(this->comm, temp_added_cells, all_added_cells);
			this->added_cells.clear();

			// continue until induced refinement has stopped
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

					neighbour_lists_to_update.insert(*cell);

					// remote refined cell doesn't count as a neighbour anymore
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

					neighbour_lists_to_update.insert(*cell);

					// remote refined cell doesn't count as a neighbour anymore
					this->remote_neighbours.erase(*neighbour);

					// induced refinement
					if (this->get_refinement_level(*neighbour) > this->get_refinement_level(*cell)) {
						induced_refines.insert(*cell);
					}
				}
			}

			// update neighbour and _to lists
			for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {
				if (this->cells.count(*cell) > 0 && *cell == this->get_child(*cell)) {
					this->update_neighbours(*cell);
				}
			}
			for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {
				if (this->cells.count(*cell) > 0 && *cell == this->get_child(*cell)) {
					this->update_neighbours_to(*cell);
				}
			}

			// don't unrefine if another process has refined a small enough remote neighbour
			for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {
				for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {

					if (cell_creator == this->comm.rank()) {
						continue;
					}

					uint64_t parent_of_created = this->get_parent(*created_cell);

					if (this->remote_neighbours.count(parent_of_created) == 0) {
						continue;
					}

					std::vector<uint64_t> local_neighbours = this->get_neighbours_of(parent_of_created);
					for (std::vector<uint64_t>::const_iterator local_neighbour = local_neighbours.begin(); local_neighbour != local_neighbours.end(); local_neighbour++) {

						if (this->cells.count(*local_neighbour) == 0) {
							continue;
						}

						if (this->get_refinement_level(*local_neighbour) > this->get_refinement_level(parent_of_created)) {
							continue;
						}

						// remote refined neighbour is too small
						std::vector<uint64_t> peers = this->get_all_children(this->get_parent(*local_neighbour));
						for (std::vector<uint64_t>::const_iterator peer = peers.begin(); peer != peers.end(); peer++) {
							this->removed_cells.erase(*peer);
						}
					}
				}
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

		// unrefine
		while (true) {

			std::vector<std::vector<uint64_t> > all_removed_cells;
			std::vector<uint64_t> temp_removed_cells(this->removed_cells.begin(), this->removed_cells.end());
			all_gather(this->comm, temp_removed_cells, all_removed_cells);
			this->removed_cells.clear();

			// continue until induced unrefinement has stopped
			bool cells_removed = false;
			for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
				if (all_removed_cells[cell_remover].size() > 0) {
					cells_removed = true;
					break;
				}
			}
			if (!cells_removed) {
				break;
			}

			// all cells whose neighbour and _to list will be updated
			boost::unordered_set<uint64_t> neighbour_lists_to_update;

			// induced unrefinement from other processes
			for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {

				// start unrefining from smallest cells and unrefine one cell at a time so neighbour searching doesn't get confused
				sort(all_removed_cells[cell_remover].begin(), all_removed_cells[cell_remover].end());
				reverse(all_removed_cells[cell_remover].begin(), all_removed_cells[cell_remover].end());

				for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {

					// initial cells shouldn't be removed
					assert(this->get_refinement_level(*removed_cell) > 0);

					// update neighbour lists of removed cells' local parent
					if (cell_remover == this->comm.rank()) {
						neighbour_lists_to_update.erase(*removed_cell);
						uint64_t parent_of_removed = this->get_parent(*removed_cell);
						if (this->cells.count(parent_of_removed) > 0) {
							neighbour_lists_to_update.insert(parent_of_removed);
						}
					}

					// induce unrefinement if remotely unrefined cell is a neighbour to any local cells
					bool induce_refinement = false;
					for (boost::unordered_set<uint64_t>::const_iterator cell = this->cells_with_remote_neighbours.begin(); cell != this->cells_with_remote_neighbours.end(); cell++) {
						if (this->is_neighbour(*cell, *removed_cell)) {
							induce_refinement = true;
						}
					}

					if (!induce_refinement) {
						continue;
					}

					// unrefine local neighbours that are too small
					boost::unordered_set<uint64_t> neighbours_of_parent = this->get_neighbours_of_parent(*removed_cell);
					for (boost::unordered_set<uint64_t>::const_iterator neighbour_of_parent = neighbours_of_parent.begin(); neighbour_of_parent != neighbours_of_parent.end(); neighbour_of_parent++) {

						if (this->cells.count(*neighbour_of_parent) == 0) {
							continue;
						}

						if (this->get_refinement_level(*neighbour_of_parent) > this->get_refinement_level(*removed_cell)) {
							this->unrefine_completely(*neighbour_of_parent);
						}
					}
				}
			}

			// remove cells from cell to process mappings
			for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {
				for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {

					// initial cells shouldn't be removed
					assert(this->get_refinement_level(*removed_cell) > 0);

					this->cell_process.erase(*removed_cell);

					if (cell_remover == this->comm.rank()) {
						// update local cells' and their parents' data structures
						/*TODO: tell the user what cell were removed, retain their data for interpolating their parent's new value and move that data to the process of the parent
							this->removed_cells_data[*removed_cell] = this->cells[*removed_cell];*/
						this->cells.erase(*removed_cell);
						this->neighbours.erase(*removed_cell);
						this->neighbours_to.erase(*removed_cell);
						this->cells_with_remote_neighbours.erase(*removed_cell);
					}
				}
			}

			// update local neighbour and _to lists
			for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {
				this->update_neighbours(*cell);
			}
			for (boost::unordered_set<uint64_t>::const_iterator cell = neighbour_lists_to_update.begin(); cell != neighbour_lists_to_update.end(); cell++) {
				this->update_neighbours_to(*cell);
			}
		}

		return new_cells;
	}



private:

	// starting corner coordinate of the grid
	double x_start, y_start, z_start;
	// length of unrefined cells in all directions
	double cell_size;
	// size of the grid in unrefined cells
	uint64_t x_length, y_length, z_length;
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

	// cell on this process and its neithis->cells.count(current_parent) > 0ghbours
	boost::unordered_map<uint64_t, std::vector<uint64_t> > neighbours;

	/*
	Cell on this process and those cells that aren't neighbours of this cell but whose neighbour this cell is.
	For example with a stencil size of 1 in the following grid:

	|-----------|
	|     |5 |6 |
	|  1  |--|--|
	|     |9 |10|
	|-----------|

	neighbours_to[6] = 1 because neighbours[6] = 5, 9, 10 while
	neighbours_to[5] is empty because neighbours[5] = 1, 6, 9, 10
	*/
	boost::unordered_map<uint64_t, std::vector<uint64_t> > neighbours_to;

	// on which process every cell in the grid is
	boost::unordered_map<uint64_t, int> cell_process;

	// cells on this process that have a neighbour on another process or are considered as a neighbour of a cell on another process
	boost::unordered_set<uint64_t> cells_with_remote_neighbours;

	// remote neighbours and their data, of cells on this process
	boost::unordered_map<uint64_t, UserData> remote_neighbours;

	// cells added to or removed from the grid on this process that haven't been communicated to other processes yet
	boost::unordered_set<uint64_t> added_cells, removed_cells;

	// pending neighbour data requests for this process
	std::vector<boost::mpi::request> requests;

	// processes and their cells whose data has to be sent to this process during neighbour data updates
	boost::unordered_map<int, boost::unordered_set<uint64_t> > cells_to_receive;

	// storage for user cell data that awaits transfer to or from this process
	boost::unordered_map<int, std::vector<UserData> > incoming_data, outgoing_data;


	/*
	Unrefines local cells whose unrefining is induced by unrefining given cell
	Doesn't override unrefines with refines
	In contrast to refine_* functions, removes cells directly from cell to process mappings
	*/
	void unrefine_completely_internal(const uint64_t cell)
	{

		assert(this->cell_process.count(cell) > 0);
		assert(this->get_refinement_level(cell) > 0);
		assert(cell == this->get_child(cell));

		std::vector<uint64_t> children = this->get_all_children(cell);

		// don't unrefine the same cell again
		if (this->removed_cells.count(cell) > 0) {
			return;
		}

		// don't unrefine if any neighbour of given cell's parent is small enough and has been refined
		boost::unordered_set<uint64_t> parents_neighbours = this->get_neighbours_of_parent(cell);
		for (boost::unordered_set<uint64_t>::const_iterator parents_neighbour = parents_neighbours.begin(); parents_neighbour != parents_neighbours.end(); parents_neighbour++) {

			// neighbour is large enough even if its being refined
			if (this->get_refinement_level(*parents_neighbour) < this->get_refinement_level(cell)) {
				continue;
			}

			std::vector<uint64_t> neighbours_children = this->get_all_children(*parents_neighbour);
			if (neighbours_children.size() > 0 && this->added_cells.count(neighbours_children[0]) > 0) {
				return;
			}
		}

		this->removed_cells.insert(cell);

		if (this->cells.count(cell) > 0) {
			std::vector<uint64_t> peers = this->get_all_children(this->get_parent(cell));
			this->removed_cells.insert(peers.begin(), peers.end());
			std::cout << "Process " << this->comm.rank() << " unrefining cell " << cell << ", induced unrefines of cells: ";
		}

		// unrefine
		for (boost::unordered_set<uint64_t>::const_iterator parents_neighbour = parents_neighbours.begin(); parents_neighbour != parents_neighbours.end(); parents_neighbour++) {

			// not remote cells
			if (this->cells.count(*parents_neighbour) == 0) {
				continue;
			}

			if (this->get_refinement_level(*parents_neighbour) > this->get_refinement_level(cell)) {
				this->removed_cells.insert(*parents_neighbour);
				std::cout << *parents_neighbour << " ";
			}

			// induced unrefinement, unrefine also neighbours' parents until they are large enough
			uint64_t current_parent = this->get_parent(*parents_neighbour);
			while (this->get_refinement_level(current_parent) > this->get_refinement_level(cell)) {
				if (this->cells.count(current_parent) > 0) {
					this->removed_cells.insert(current_parent);
					std::cout << current_parent << " ";
				}
				current_parent = this->get_parent(current_parent);
			}
		}
		std::cout << std::endl;
	}


	/*
	Updates the neighbour list of given cell without children on this process
	*/
	void update_neighbours(const uint64_t cell)
	{
		assert(this->cells.count(cell) > 0);
		assert(cell == this->get_child(cell));

		this->neighbours[cell] = this->get_neighbours_of(cell);
	}


	/*
	Returns the existing neighbours without children of given cell
	Returns nothing if the given cell has children
	*/
	std::vector<uint64_t> get_neighbours_of(const uint64_t id)
	{
		assert(this->cell_process.count(id) > 0);

		std::vector<uint64_t> return_neighbours;

		if (id != this->get_child(id)) {
			return return_neighbours;
		}

		uint64_t x_index = this->get_x_index(id), y_index = this->get_y_index(id), z_index = this->get_z_index(id);

		// search neighbours in cells of the same size as the given cell
		uint64_t size_in_indices = this->get_cell_size_in_indices(id);

		// don't start searching outside of the grid
		uint64_t current_x_index;
		if (x_index < size_in_indices * this->neighbourhood_size) {
			current_x_index = 0;
		} else {
			current_x_index = x_index - size_in_indices * this->neighbourhood_size;
		}

		// the refinement level difference between neighbours is <= 2 because induced refines haven't propagated across processes when neighbours are updated, search index can increase by this amount at most
		uint64_t index_increase = size_in_indices / 4;
		if (index_increase == 0) {
			index_increase = 1;
		}

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours;

		// search for neighbours in cells that share a vertex with the given cell
		for (; current_x_index < x_index + size_in_indices * (1 + this->neighbourhood_size); current_x_index += index_increase) {

			// don't search outside of the grid
			if (current_x_index >= this->x_length * int(pow(2, this->max_refinement_level))) {
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

				if (current_y_index >= this->y_length * int(pow(2, this->max_refinement_level))) {
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

					if (current_z_index >= this->z_length * int(pow(2, this->max_refinement_level))) {
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

					// the refinement level difference between neighbours is <= 2, only search for those cells
					int search_min_ref_level = this->get_refinement_level(id) - 2, search_max_ref_level = this->get_refinement_level(id) + 2;
					if (search_min_ref_level < 0) {
						search_min_ref_level = 0;
					}
					if (search_max_ref_level > this->max_refinement_level) {
						search_max_ref_level = this->max_refinement_level;
					}

					uint64_t neighbour = this->get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level);
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


	/*
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


	/*
	Returns the existing cells without children that have the given cell as a neighbour
	Returns nothing if the given cell has children
	Doesn't update neighbour lists
	*/
	boost::unordered_set<uint64_t> get_neighbours_to(const uint64_t id)
	{
		assert(this->cell_process.count(id) > 0);

		boost::unordered_set<uint64_t> return_neighbours;

		if (id != this->get_child(id)) {
			return return_neighbours;
		}

		/*
		The largest cell that considers this as a neighbour has a refinement level of no less than this cells refinement level - 1
		All of those can be found by searching the neighbours of this cell's neighbours
		*/

		// don't add the same neighbour more than once
		boost::unordered_set<uint64_t> unique_neighbours_to;

		// first round of neighbour search, use existing neighbour lists if available
		boost::unordered_set<uint64_t> neighbours1;
		if (this->cells.count(id) > 0) {

			for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[id].begin(); neighbour != this->neighbours[id].end(); neighbour++) {
				assert(*neighbour == this->get_child(*neighbour));
				neighbours1.insert(*neighbour);
				unique_neighbours_to.insert(*neighbour);
			}

		} else {
			std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(id);
			neighbours1.insert(temp_neighbours.begin(), temp_neighbours.end());
			unique_neighbours_to.insert(temp_neighbours.begin(), temp_neighbours.end());
		}

		// second round of neighbour search, try to use existing neighbour lists
		for (boost::unordered_set<uint64_t>::const_iterator neighbour1 = neighbours1.begin(); neighbour1 != neighbours1.end(); neighbour1++) {

			if (this->cells.count(*neighbour1) > 0) {
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*neighbour1].begin(); neighbour != this->neighbours[*neighbour1].end(); neighbour++) {
					assert(*neighbour == this->get_child(*neighbour));
					unique_neighbours_to.insert(*neighbour);
				}
			} else {
				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(*neighbour1);
				unique_neighbours_to.insert(temp_neighbours.begin(), temp_neighbours.end());
			}
		}
		unique_neighbours_to.erase(id);

		// return only cells which consider given cell as a neighbour
		//uint64_t cell_size = this->get_cell_size_in_indices(id);
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


	/*
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
				this->remote_neighbours[*neighbour];
			}
		}
		// cells with given cell as neighbour
		for (std::vector<uint64_t>::const_iterator neighbour_to = this->neighbours_to[cell].begin(); neighbour_to != this->neighbours_to[cell].end(); neighbour_to++) {
			if (this->cell_process[*neighbour_to] != this->comm.rank()) {
				this->cells_with_remote_neighbours.insert(cell);
				this->remote_neighbours[*neighbour_to];
			}
		}
	}


	/*
	Returns the existing neighbours of given cell's parent or existing neighbours of given cell if it is an unrefined cell (= without a parent)
	*/
	boost::unordered_set<uint64_t> get_neighbours_of_parent(const uint64_t cell)
	{
		assert(this->cell_process.count(cell) > 0);

		boost::unordered_set<uint64_t> return_neighbours;

		// given cells doesn't have a parent
		if (this->get_refinement_level(cell) == 0) {
			if (this->cells.count(cell) > 0) {
				return_neighbours.insert(this->neighbours[cell].begin(), this->neighbours[cell].end());
			} else {
				std::vector<uint64_t> temp_neighbours = this->get_neighbours_of(cell);
				return_neighbours.insert(temp_neighbours.begin(), temp_neighbours.end());
			}
			return return_neighbours;
		}

		uint64_t parent = this->get_parent(cell);

		// search for given cell's peers' neighbours of neighbours of... until a neighbour isn't really a neighbour of given cell's parent
		boost::unordered_set<uint64_t> new_neighbours;
		std::vector<uint64_t> peers = this->get_all_children(parent);
		for (std::vector<uint64_t>::const_iterator peer = peers.begin(); peer != peers.end(); peer++) {

			// try to use existing neighbour list
			if (this->cells.count(*peer) > 0) {
				for (std::vector<uint64_t>::const_iterator neighbour = this->neighbours[*peer].begin(); neighbour != this->neighbours[*peer].end(); neighbour++) {
					if (this->is_neighbour(parent, *neighbour)) {
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
					for (std::vector<uint64_t>::const_iterator neigh_of_new_neigh = this->neighbours[*search_neighbour].begin(); neigh_of_new_neigh != this->neighbours[*search_neighbour].end(); neigh_of_new_neigh++)  {
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


	// Returns true if cell1 considers cell2 as a neighbour, even if neither of them exists
	bool is_neighbour(const uint64_t cell1, const uint64_t cell2)
	{
		assert(cell1 > 0);
		assert(cell1 <= this->max_cell_number);
		assert(cell2 > 0);
		assert(cell2 <= this->max_cell_number);

		uint64_t cell1_x_index = this->get_x_index(cell1), cell1_y_index = this->get_y_index(cell1), cell1_z_index = this->get_z_index(cell1);
		uint64_t cell2_x_index = this->get_x_index(cell2), cell2_y_index = this->get_y_index(cell2), cell2_z_index = this->get_z_index(cell2);
		uint64_t cell1_size = this->get_cell_size_in_indices(cell1);
		uint64_t cell2_size = this->get_cell_size_in_indices(cell2);

		uint64_t dindex1 = cell2_size + cell1_size * this->neighbourhood_size;
		uint64_t dindex2 = cell1_size * this->neighbourhood_size;
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


	/*
	Given a cell that exists and has children returns one of the children
	Returns the given cell if it doesn't have children or 0 if the cell doesn't exist
	*/
	uint64_t get_child(const uint64_t cell)
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


	/*
	Given a cell that exists and has a parent returns the parent cell
	Returns the given cell if it doesn't have a parent or 0 if the cell doesn't exist
	*/
	uint64_t get_parent(const uint64_t cell)
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


	/*
	Returns the smallest existing cell at the given coordinate
	Returns 0 if the coordinate is outside of the grid or the cell is on another process
	*/
	uint64_t get_smallest_cell_from_coordinate(const double x, const double y, const double z)
	{
		if (x < this->x_start
			|| x > this->x_start + this->cell_size * this->x_length
			|| y < this->y_start
			|| y > this->y_start + this->cell_size * this->y_length
			|| z < this->z_start
			|| z > this->z_start + this->cell_size * this->z_length)
		{
			return 0;
		}

		return this->get_cell_from_indices(this->get_x_index(x), this->get_y_index(y), this->get_z_index(z), 0, this->max_refinement_level);
	}


	/*
	These return the index of the cell with given id in x, y or z direction of the grid, starting from 0
	For cells that are larger than the smallest possible according to max_refinement_level, the index closest to the grids starting corner is returned
	 */
	uint64_t get_x_index(uint64_t id)
	{
		assert(id);
		assert(id <= this->max_cell_number);

		// substract ids of larger cells
		int refinement_level = get_refinement_level(id);
		for (int i = 0; i < refinement_level; i++) {
			id -= x_length * y_length * z_length * int(pow(8, i));
		}

		// get the index at this cells refinement level
		id -= 1;	// cell numbering starts from 1
		int this_level_index = id % (x_length * int(pow(2, refinement_level)));

		return this_level_index * int(pow(2, max_refinement_level - refinement_level));
	}

	uint64_t get_y_index(uint64_t id)
	{
		assert(id);
		assert(id <= this->max_cell_number);

		// substract ids of larger cells
		int refinement_level = get_refinement_level(id);
		for (int i = 0; i < refinement_level; i++) {
			id -= x_length * y_length * z_length * int(pow(8, i));
		}

		// get the index at this cells refinement level
		id -= 1;	// cell numbering starts from 1
		int this_level_index =  int(id / (x_length * int(pow(2, refinement_level)))) % (y_length  * int(pow(2, refinement_level)));

		return this_level_index * int(pow(2, max_refinement_level - refinement_level));
	}

	uint64_t get_z_index(uint64_t id)
	{
		assert(id);
		assert(id <= this->max_cell_number);

		// substract ids of larger cells
		int refinement_level = get_refinement_level(id);
		for (int i = 0; i < refinement_level; i++) {
			id -= x_length * y_length * z_length * int(pow(8, i));
		}

		// get the index at this cells refinement level
		id -= 1;	// cell numbering starts from 1
		int this_level_index =  int(id / (x_length * y_length * int(pow(2, 2 * refinement_level))));

		return this_level_index * int(pow(2, max_refinement_level - refinement_level));
	}

	// These return the x, y or z index of the given coordinate
	uint64_t get_x_index(const double x)
	{
		assert((x >= this->x_start) and (x <= this->x_start + this->x_length * this->cell_size));

		return uint64_t((x - this->x_start) / (this->cell_size / int(pow(2, this->max_refinement_level))));
	}

	uint64_t get_y_index(const double y)
	{
		assert((y >= this->y_start) and (y <= this->y_start + this->y_length * this->cell_size));

		return uint64_t((y - this->y_start) / (this->cell_size / int(pow(2, this->max_refinement_level))));
	}

	uint64_t get_z_index(const double z)
	{
		assert((z >= this->z_start) and (z <= this->z_start + this->z_length * this->cell_size));

		return uint64_t((z - this->z_start) / (this->cell_size / int(pow(2, this->max_refinement_level))));
	}


	/*
	Returns the smallest cell at given indices between given refinement ranges inclusive
	Returns 0 if no cell between given refinement ranges exists
	*/
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int minimum_refinement_level, const int maximum_refinement_level)
	{
		assert(x_index < this->x_length * int(pow(2, this->max_refinement_level)));
		assert(y_index < this->y_length * int(pow(2, this->max_refinement_level)));
		assert(z_index < this->z_length * int(pow(2, this->max_refinement_level)));
		assert(minimum_refinement_level <= maximum_refinement_level);

		int average_refinement_level = (maximum_refinement_level + minimum_refinement_level) / 2;
		uint64_t id = this->get_cell_from_indices(x_index, y_index, z_index, average_refinement_level);

		// use binary search recursively (assumes that max_refinement_level_difference == 1)
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


	// Returns the cell of given refinement level at given indices even if it doesn't exist
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int refinement_level)
	{
		assert(x_index < this->x_length * int(pow(2, this->max_refinement_level)));
		assert(y_index < this->y_length * int(pow(2, this->max_refinement_level)));
		assert(z_index < this->z_length * int(pow(2, this->max_refinement_level)));
		assert(refinement_level <= this->max_refinement_level);

		uint64_t id = 1;

		// add ids of larger cells
		for (int i = 0; i < refinement_level; i++) {
			id += x_length * y_length * z_length * int(pow(8, i));
		}

		// convert to indices of this cells refinement level
		int this_level_x_index = x_index / int(pow(2, max_refinement_level - refinement_level));
		int this_level_y_index = y_index / int(pow(2, max_refinement_level - refinement_level));
		int this_level_z_index = z_index / int(pow(2, max_refinement_level - refinement_level));

		// get the size of the grid in terms of cells of this level
		uint64_t this_level_x_length = x_length * int(pow(2, refinement_level));
		uint64_t this_level_y_length = y_length * int(pow(2, refinement_level));

		id += this_level_x_index + this_level_y_index * this_level_x_length + this_level_z_index * this_level_x_length * this_level_y_length;

		return id;
	}


	// Returns the lengths of given cell in indices in every direction
	uint64_t get_cell_size_in_indices(const uint64_t id)
	{
		assert(id);
		return uint64_t(pow(2, max_refinement_level - get_refinement_level(id)));
	}


	/*
	Returns all children of given cell regardless of whether they exist
	Returns no cells if childrens' refinement level would exceed max_refinement_level
	 */
	std::vector<uint64_t> get_all_children(const uint64_t id)
	{
		assert(id);
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
		int index_offset = pow(2, max_refinement_level - refinement_level);
		for (int x_index_offset = 0; x_index_offset < 2 * index_offset; x_index_offset += index_offset) {
			for (int y_index_offset = 0; y_index_offset < 2 * index_offset; y_index_offset += index_offset) {
				for (int z_index_offset = 0; z_index_offset < 2 * index_offset; z_index_offset += index_offset) {
					children.push_back(get_cell_from_indices(x_index + x_index_offset, y_index + y_index_offset, z_index + z_index_offset, refinement_level));
				}
			}
		}

		return children;
	}


	/*
	Returns the size of user data in given cell (global_id) in bytes
	*/
	static int user_data_size(void* /*data*/, int /*global_id_size*/, int /*local_id_size*/, ZOLTAN_ID_PTR /*global_id*/, ZOLTAN_ID_PTR /*local_id*/, int* error)
	{
		*error = ZOLTAN_OK;
		return sizeof(UserData);
	}


	/*
	Returns the number of values needed to represent the coordinate of a cell
	*/
	static int get_grid_dimensionality(void* /*data*/, int* error)
	{
		*error = ZOLTAN_OK;
		return 3;
	}


	/*
	Fills geom_vec with the coordinate of the given cell (global_id)
	*/
	static void fill_with_cell_coordinates(void *data, int /*global_id_size*/, int /*local_id_size*/, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR /*local_id*/, double *geom_vec, int *error)
	{
		uint64_t cell = uint64_t(*global_id);
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);

		if (dccrg_instance->cells.count(cell) == 0) {
			*error = ZOLTAN_FATAL;
			std::cerr << "Process " << dccrg_instance->comm.rank() << ": Zoltan wanted the coordinate of a non-existing cell " << cell << std::endl;
			return;
		} else {
			*error = ZOLTAN_OK;
		}

		geom_vec[0] = dccrg_instance->get_cell_x(cell);
		geom_vec[1] = dccrg_instance->get_cell_y(cell);
		geom_vec[2] = dccrg_instance->get_cell_z(cell);
	}


	/*
	Copies the data of given cell (global_id) to a buffer created by Zoltan
	*/
	static void pack_user_data(void* data, int /*global_id_size*/, int /*local_id_size*/, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR /*local_id*/, int /*destination_part*/, int /*buffer_size*/, char* buffer, int* error)
	{
		uint64_t cell = uint64_t(*global_id);
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);

		if (dccrg_instance->cells.count(cell) == 0) {
			*error = ZOLTAN_FATAL;
			std::cerr << "Process " << dccrg_instance->comm.rank() << ": Zoltan wanted to pack a non-existing cell " << cell << std::endl;
			return;
		} else {
			*error = ZOLTAN_OK;
		}

		UserData* user_data = &(dccrg_instance->cells[cell]);
		memcpy(buffer, user_data, sizeof(UserData));

		// remove the cell from this process' data structures
		dccrg_instance->cells.erase(cell);
		dccrg_instance->cells_with_remote_neighbours.erase(cell);

		// tell also others that the cell doesn't exist on this process anymore
		dccrg_instance->removed_cells.insert(cell);
	}


	/*
	Copies the data of given cell (global_id) from a buffer created by Zoltan
	*/
	static void unpack_user_data(void* data, int /*global_id_size*/, ZOLTAN_ID_PTR global_id, int /*buffer_size*/, char* buffer, int* error)
	{
		uint64_t cell = uint64_t(*global_id);
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);

		if (dccrg_instance->get_refinement_level(cell) < 0 || dccrg_instance->get_refinement_level(cell) > dccrg_instance->max_refinement_level) {
			*error = ZOLTAN_FATAL;
			std::cerr << "Process " << dccrg_instance->comm.rank() << ": Zoltan gave a non-existing cell " << cell << " for unpacking" << std::endl;
			return;
		} else {
			*error = ZOLTAN_OK;
		}

		// add the cell to this process
		dccrg_instance->cells[cell];
		memcpy(&(dccrg_instance->cells[cell]), buffer, sizeof(UserData));

		// record that the cell exist now on this process
		dccrg_instance->added_cells.insert(cell);
	}


	// Returns the number of cells on this process
	static int get_number_of_cells(void* data, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;
		return dccrg_instance->cells.size();
	}


	/*
	Writes all cell ids on this process to the global_ids array
	*/
	static void get_cell_list(void* data, int /*global_id_size*/, int /*local_id_size*/, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR /*local_ids*/, int /*weight_dimension*/, float* /*object_weights*/, int* error)
	{
		dccrg<UserData>* dccrg_instance = reinterpret_cast<dccrg<UserData> *>(data);
		*error = ZOLTAN_OK;
		int i = 0;
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = dccrg_instance->cells.begin(); cell != dccrg_instance->cells.end(); cell++, i++) {
			global_ids[i] = cell->first;
		}
	}

};

#endif
