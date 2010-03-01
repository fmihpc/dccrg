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


#include "boost/mpi.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cassert"
#include "cstdio"
#include "cstdlib"
#include "cstring"
#include "fstream"
#include "limits"
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

	x_length, y_length, z_length: the number of cell in the grid in x, y and z direction

	neighbourhood_size:
		Determines which cells are considered neighbours.
		When calculating the neighbours of a given cell a cube of length neighbourhood_size + 1 in every direction is considered, centered at the cell for which neighbours are being calculated.
		The unit lenght of the cube is the cell for neighbours are being calculated.
		TODO: If neighbourhood_size == 0, only cells that share a face are considered.

	maximum_refinement_level:
		The maximum number of times an unrefined cell can be refined (replacing it with 8 smaller cells)
		Optional: if not given it is maximized based on the grids size
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
		if (comm.rank() == 0) {
			for (uint64_t id = 1; id <= x_length * y_length * z_length; id++) {
				this->added_cells.insert(id);
			}
		}
		// and tell that to the other processes
		stop_refining();
	}


	/*
	Returns all cells on the same process that don't have children
	*/
	std::vector<uint64_t> get_cells(void)
	{
		std::vector<uint64_t> all_cells;
		all_cells.reserve(this->cells.size());

		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

			uint64_t child = get_child(cell->first);
			assert(child > 0);

			if (child == cell->first) {
				all_cells.push_back(cell->first);
			}
		}

		return all_cells;
	}


	/*
	Returns a pointer to the user supplied data of given cell
	Return NULL if the given cell isn't on this process and if the given cell isn't a neighbour of any cell on this process
	*/
	UserData* operator [] (uint64_t cell)
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
	Informs other processes about cells created / removed on this process
	Must be called simultaneously on all processes, after even one of them has refined / recoarsened its cells, before doing anything else with the grid
	*/
	void stop_refining(void)
	{

		std::vector<uint64_t> temp_added_cells(this->added_cells.begin(), this->added_cells.end());
		std::vector<std::vector<uint64_t> > all_added_cells;

		all_gather(this->comm, temp_added_cells, all_added_cells);

		// update created cells from other processes
		for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {

			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {

				this->cell_process[*created_cell] = cell_creator;

				// cells created on this process also store user data
				if (this->comm.rank() == cell_creator) {
					this->cells[*created_cell];
				}
			}
		}

		// TODO: also redistribute removed cells and handle the case when the same cell would be refined and recoarsened (the cells children added but the cell itself removed)

		// handle removed cells that were actually only migrated to another process

		// update neighbour lists of all cells on this process
		// TODO: only update for cells whose neighbourhood changed
		this->cells_with_remote_neighbours.clear();
		for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {

			this->update_neighbours(cell->first);

			// check whether any of the neighbours are on another process
			for (boost::unordered_set<uint64_t>::const_iterator neighbour = this->neighbours[cell->first].begin(); neighbour != neighbours[cell->first].end(); neighbour++) {
				if (this->cell_process[*neighbour] != this->comm.rank()) {
					this->cells_with_remote_neighbours.insert(cell->first);
				}
			}
		}

		this->added_cells.clear();
		this->comm.barrier();
	}


	/*
	Load balances the grids cells among processes
	Must be called simultaneously on all processes
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

		// calculate where cells have migrated to update cell_process[]
		// TODO: could probably be combined with stop_refining()
		std::vector<uint64_t> temp_removed_cells(this->removed_cells.begin(), this->removed_cells.end());
		std::vector<std::vector<uint64_t> > all_removed_cells;
		all_gather(this->comm, temp_removed_cells, all_removed_cells);
		// removed cells
		for (int cell_remover = 0; cell_remover < int(all_removed_cells.size()); cell_remover++) {

			for (std::vector<uint64_t>::const_iterator removed_cell = all_removed_cells[cell_remover].begin(); removed_cell != all_removed_cells[cell_remover].end(); removed_cell++) {

				this->cell_process.erase(*removed_cell);
			}
		}
		this->removed_cells.clear();

		// created cells
		std::vector<uint64_t> temp_added_cells(this->added_cells.begin(), this->added_cells.end());
		std::vector<std::vector<uint64_t> > all_added_cells;
		all_gather(this->comm, temp_added_cells, all_added_cells);
		for (int cell_creator = 0; cell_creator < int(all_added_cells.size()); cell_creator++) {

			for (std::vector<uint64_t>::const_iterator created_cell = all_added_cells[cell_creator].begin(); created_cell != all_added_cells[cell_creator].end(); created_cell++) {

				this->cell_process[*created_cell] = cell_creator;
			}
		}
		this->added_cells.clear();
	}


	/*
	Updates the user data of those cells that have at least one neighbour on another process
	Must be called simultaneously on all processes
	*/
	void update_remote_neighbour_data(void)
	{
		this->comm.barrier();

		// a matrix of sets where each sets holds the cells whose data has to be sent from one process to another
		// outermost index is the process from which to send and innermost index is the process that will receive the cells
		std::vector<std::vector<boost::unordered_set<uint64_t> > > send_receive_matrix;

		// fill the matrix with empty sets
		for (int i = 0; i < this->comm.size(); i++) {
			std::vector<boost::unordered_set<uint64_t> > send_receive_vector;
			boost::unordered_set<uint64_t> temp_set;
			for (int j = 0; j < this->comm.size(); j++) {
				send_receive_vector.push_back(temp_set);
			}
			send_receive_matrix.push_back(send_receive_vector);
		}

		// go through all cells and all their neighbours
		for (boost::unordered_map<uint64_t, int>::const_iterator cell = this->cell_process.begin(); cell != this->cell_process.end(); cell++) {

			uint64_t current_cell = cell->first;
			int current_process = cell->second;

			boost::unordered_set<uint64_t> tmp_neighbours = this->get_neighbours_internal(current_cell);
			for (boost::unordered_set<uint64_t>::const_iterator neighbour = tmp_neighbours.begin(); neighbour != tmp_neighbours.end(); neighbour++) {

				if (this->cell_process[*neighbour] != current_process) {
					// *neighbours process has to send *neighbours cell data to current_process
					send_receive_matrix[this->cell_process[*neighbour]][current_process].insert(*neighbour);
					// current process has to send currents cell data to neighbour
					send_receive_matrix[current_process][this->cell_process[*neighbour]].insert(current_cell);
				}
			}
		}

		#ifndef NDEBUG
		// print the send_receive_matrix of every process
		for (int process = 0; process < this->comm.size(); process++) {
			this->comm.barrier();

			if (process != this->comm.rank()) {
				continue;
			}

			for (int sender = 0; sender < this->comm.size(); sender++) {
				for (int receiver = 0; receiver < this->comm.size(); receiver++) {

					if (sender == receiver) {
						continue;
					}

					std::cout << "Process " << sender << " sending these cells' data to process " << receiver << ": ";
					for (boost::unordered_set<uint64_t>::const_iterator cell = send_receive_matrix[sender][receiver].begin(); cell != send_receive_matrix[sender][receiver].end(); cell++) {
						std::cout << *cell << " ";
					}
					std::cout << std::endl;
				}
			}
			std::cout.flush();
		}
		#endif

		// spread the neighbour data
		for (int sender = 0; sender < this->comm.size(); sender++) {
			for (int receiver = 0; receiver < this->comm.size(); receiver++) {

				if (sender == receiver) {
					// don't send to self
					continue;
				}

				if (sender != this->comm.rank() && receiver != this->comm.rank()) {
					// this process is neither sending nor receiving
					continue;
				}

				if (send_receive_matrix[sender][receiver].size() == 0) {
					// no data to send / receive
					continue;
				}

				int send_receive_tag = sender * this->comm.size() + receiver;

				// send cell data in the same order
				std::vector<uint64_t> send_receive_cells(send_receive_matrix[sender][receiver].begin(), send_receive_matrix[sender][receiver].end());
				sort(send_receive_cells.begin(), send_receive_cells.end());

				std::vector<UserData> incoming_data, outgoing_data;
				if (this->comm.rank() == sender) {

					// construct the outgoing data vector
					for (std::vector<uint64_t>::const_iterator cell = send_receive_cells.begin(); cell != send_receive_cells.end(); cell++) {
						UserData* user_data = (*this)[*cell];
						assert(user_data != NULL);
						outgoing_data.push_back(*user_data);
					}

					this->comm.send(receiver, send_receive_tag, outgoing_data);

				} else {

					this->comm.recv(sender, send_receive_tag, incoming_data);
					// incorporate received data
					int i = 0;
					for (std::vector<uint64_t>::const_iterator cell = send_receive_cells.begin(); cell != send_receive_cells.end(); cell++, i++) {
						this->remote_neighbours[*cell] = incoming_data[i];
					}

				}
			}
		}

		this->comm.barrier();
	}


	/*
	Returns the neighbours (some of which might be on another process) of given cell
	Returns nothing if given cell doesn't exist or is on another process
	*/
	boost::unordered_set<uint64_t> get_neighbours(const uint64_t id)
	{
		if (this->cells.count(id) > 0) {
			return this->neighbours[id];
		} else {
			boost::unordered_set<uint64_t> empty_set;
			return empty_set;
		}
	}


	/*
	Given a cell that exists and has children returns one of the children
	Returns the given cell if it doesnt have children or 0 if the cell doesn't exist
	*/
	uint64_t get_child(const uint64_t id)
	{
		if (this->cell_process.count(id) == 0) {
			return 0;
		}

		// given cell cannot have children
		if (get_refinement_level(id) == this->max_refinement_level) {
			return id;
		}

		uint64_t child = get_cell_from_indices(get_x_index(id), get_y_index(id), get_z_index(id), get_refinement_level(id) + 1);
		if (this->cell_process.count(child) > 0) {
			return child;
		} else {
			return id;
		}
	}


	/*
	Returns the maximum possible refinement level of any cell in the grid (0 means unrefined)
	*/
	int get_max_refinement_level(void)
	{
		return this->max_refinement_level;
	}


	// The following return the x, y or z coordinate of given cells center or NaN if the cell doesn't exist or exists on another process
	double get_cell_x(const uint64_t id)
	{
		if (this->cells.count(id) == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		return x_start + get_x_index(id) * cell_size / int(pow(2, max_refinement_level)) + get_cell_size(id) / 2;
	}

	double get_cell_y(const uint64_t id)
	{
		if (this->cells.count(id) == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		return y_start + get_y_index(id) * cell_size / int(pow(2, max_refinement_level)) + get_cell_size(id) / 2;
	}

	double get_cell_z(const uint64_t id)
	{
		if (this->cells.count(id) == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		return z_start + get_z_index(id) * cell_size / int(pow(2, max_refinement_level)) + get_cell_size(id) / 2;
	}


	// Returns the length of given cell or NaN if the cell doesn't exist or exitst on another process
	double get_cell_size(const uint64_t id)
	{
		if (this->cells.count(id) == 0) {
			return std::numeric_limits<double>::quiet_NaN();
		}
		return cell_size / int(pow(2, get_refinement_level(id)));
	}


	// Returns the refinement level of given cell
	int get_refinement_level(uint64_t id)
	{
		if (this->cell_process.count(id) == 0) {
			return -1;
		}

		int refinement_level;
		for (refinement_level = 0; refinement_level < max_refinement_level; refinement_level++) {
			if (id <= (x_length * y_length * z_length * int(pow(8, refinement_level)))) {
				break;
			}
			// substract ids of larger cells
			id -= x_length * y_length * z_length * int(pow(8, refinement_level));
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

		std::vector<uint64_t> leaf_cells = get_cells();
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
	Returns the smallest existing cell at the given coordinate
	Returns 0 if the coordinate is outside of the grid
	*/
	/*uint64_t get_cell_from_coordinate(const double x, const double y, const double z);
	{
		if (x < x_start
			|| x > x_start + cell_size * x_length
			|| y < y_start
			|| y > y_start + cell_size * y_length
			|| z < z_start
			|| z > z_start + cell_size * z_length)
		{
			return 0;
		}

		return get_cell_from_indices(get_x_index(x), get_y_index(y), get_z_index(z), 0, max_refinement_level);
	}*/

	/*
	Creates all children of given cell
	Does nothing in any of the following cases:
		-given cell doesn't exist
		-given cell exists on another process
		-given cells children already exist
		-the created childrens' refinement level would exceed max_refinement_level
	 */
	/*void refine_completely(const uint64_t id)
	{
		assert(id);

		std::vector<uint64_t> children;

		if (cells.count(id) == 0) {
			return children;
		}

		children = get_all_children(id);

		for (unsigned int i = 0; i < children.size(); i++) {
			cells.insert(children[i]);
		}

		return children;
	}*/


	/*
	Creates all children of the smallest existing cell at given location
	Does nothing in any of the following cases:
		-the coordinate is outside of the grid
		-the smallest cell at given location exists on another process
		-the created childrens' refinement level would exceed max_refinement_level
	*/
	/*void refine_completely(const double x, const double y, const double z);
	{
		uint64_t id = get_cell_from_coordinate(x, y, z);
		if (id == 0) {
			return;
		}

		int refinement_level = get_refinement_level(id);
		if (refinement_level >= max_refinement_level) {
			return;
		}
		refinement_level++;

		uint64_t new_id = get_cell_from_indices(get_x_index(x), get_y_index(y), get_z_index(z), refinement_level);
		assert(new_id);

		cells.insert(new_id);
		return;
	}*/


	/*
	Recursively removes all the children of given cell
	Does nothing if the given cell doesn't exist or in on another process
	*/
	/*void unrefine(const uint64_t id);
	{
		assert(id);

		if (cells.count(id) == 0) {
			return;
		}

		std::vector<uint64_t> children = get_all_children(id);

		for (unsigned int i = 0; i < children.size(); i++) {
			if (cells.count(children[i]) == 0) {
				break;
			}
			unrefine(children[i]);
			cells.erase(children[i]);
		}
	}*/


	/*
	Debugging functions
	Must be called simultaneously on all processes
	*/
	void print_existing_cells(void)
	{
		for (int process = 0; process < this->comm.size(); process++) {
			this->comm.barrier();
			if (process != this->comm.rank()) {
				continue;
			}

			std::cout << "Process " << process << " cells: ";
			for (typename boost::unordered_map<uint64_t, UserData>::const_iterator cell = this->cells.begin(); cell != this->cells.end(); cell++) {
				std::cout << cell->first << " ";
			}
			std::cout << std::endl;
			std::cout.flush();
		}
	}

	void print_cell_to_process_mappings(void)
	{
		for (int process = 0; process < this->comm.size(); process++) {
			this->comm.barrier();
			if (process != this->comm.rank()) {
				continue;
			}

			std::cout << "Process " << process << " cells and processes: ";
			for (typename boost::unordered_map<uint64_t, int>::const_iterator cell = this->cell_process.begin(); cell != this->cell_process.end(); cell++) {
				std::cout << cell->first << ": " << cell->second << ", ";
			}
			std::cout << std::endl;
			std::cout.flush();
		}
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
	// size of the neighbour stencil
	unsigned int neighbourhood_size;
	// the grid is distributed between these processes
	boost::mpi::communicator comm;

	// bookkeeping for Zoltan
	Zoltan_Struct* zoltan;
	// record whether Zoltan_LB_Partition is expected to fail (when the user selects NONE as the load balancing algorithm)
	bool no_load_balancing;

	// cells and their data on this process
	boost::unordered_map<uint64_t, UserData> cells;

	// cells on this process and their neighbours
	boost::unordered_map<uint64_t, boost::unordered_set<uint64_t> > neighbours;

	// on which process every cell in the grid is
	boost::unordered_map<uint64_t, int> cell_process;

	// cells on this process that have a neighbour on another process
	boost::unordered_set<uint64_t> cells_with_remote_neighbours;

	// remote neighbours and their data, of cells on this process
	boost::unordered_map<uint64_t, UserData> remote_neighbours;

	// cells added to or removed from the grid on this process that haven't been communicated to other processes yet
	boost::unordered_set<uint64_t> added_cells, removed_cells;


	/*
	Updates the neighbour list of given cell that is on this process
	*/
	void update_neighbours(const uint64_t id)
	{
		assert(id);
		assert(this->cells.count(id) > 0);

		this->neighbours[id] = get_neighbours_internal(id);
	}


	/*
	Returns the neighbours of given cell, if it doesn't have children
	*/
	boost::unordered_set<uint64_t> get_neighbours_internal(const uint64_t id)
	{
		assert(id);
		assert(this->cell_process.count(id) > 0);

		boost::unordered_set<uint64_t> return_neighbours;

		// return nothing if given cell has children
		uint64_t child = get_child(id);
		if (child != id) {
			return return_neighbours;
		}

		uint64_t x_index = get_x_index(id), y_index = get_y_index(id), z_index = get_z_index(id);

		// search neighbours in cells of the same size as the given cell
		uint64_t size_in_indices = get_cell_size_in_indices(id);

		// don't start searching outside of the grid
		uint64_t current_x_index;
		if (x_index < size_in_indices * this->neighbourhood_size) {
			current_x_index = 0;
		} else {
			current_x_index = x_index - size_in_indices * this->neighbourhood_size;
		}

		// the refinement level difference between neighbours is <= 1, search index can increase by this amount
		uint64_t index_increase = size_in_indices / 2;
		if (index_increase == 0) {
			index_increase = 1;
		}

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

			// the refinement level difference between neighbours is <= 1
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

				// the refinement level difference between neighbours is <= 1
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

					// the refinement level difference between neighbours is <= 1, only search for those cells
					int search_min_ref_level = get_refinement_level(id) - 1, search_max_ref_level = get_refinement_level(id) + 1;
					if (search_min_ref_level < 0) {
						search_min_ref_level = 0;
					}
					if (search_max_ref_level > this->max_refinement_level) {
						search_max_ref_level = this->max_refinement_level;
					}

					return_neighbours.insert(get_cell_from_indices(current_x_index, current_y_index, current_z_index, search_min_ref_level, search_max_ref_level));
				}
			}
		}

		return_neighbours.erase(id);	// a cell isn't a neighbour of itself

		return return_neighbours;
	}


	/*
	These return the index of the cell with given id in x, y or z direction of the grid, starting from 0
	For cells that are larger than the smallest possible according to max_refinement_level, the index closest to the grids starting corner is returned
	 */
	uint64_t get_x_index(uint64_t id)
	{
		assert(id);
		assert(id <= max_cell_number);

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

		if (id > max_cell_number) {
			assert(0);
			exit(EXIT_FAILURE);
		}

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

		if (id > max_cell_number) {
			assert(0);
			exit(EXIT_FAILURE);
		}

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
		if ((x < x_start) or (x > x_start + x_length * cell_size)) {
			assert(0);
			exit(EXIT_FAILURE);
		}

		return uint64_t((x - x_start) / (cell_size / int(pow(2, max_refinement_level))));
	}

	uint64_t get_y_index(const double y)
	{
		if ((y < y_start) or (y > y_start + y_length * cell_size)) {
			assert(0);
			exit(EXIT_FAILURE);
		}

		return uint64_t((y - y_start) / (cell_size / int(pow(2, max_refinement_level))));
	}

	uint64_t get_z_index(const double z)
	{
		if ((z < z_start) or (z > z_start + z_length * cell_size)) {
			assert(0);
			exit(EXIT_FAILURE);
		}

		return uint64_t((z - z_start) / (cell_size / int(pow(2, max_refinement_level))));
	}


	/*
	Returns the smallest cell at given indices between given refinement ranges inclusive
	Returns 0 if no cell between given refinement ranges exists
	*/
	uint64_t get_cell_from_indices(const uint64_t x_index, const uint64_t y_index, const uint64_t z_index, const int minimum_refinement_level, const int maximum_refinement_level)
	{
		assert(x_index < x_length * int(pow(2, max_refinement_level)));
		assert(y_index < y_length * int(pow(2, max_refinement_level)));
		assert(z_index < z_length * int(pow(2, max_refinement_level)));
		assert(minimum_refinement_level <= maximum_refinement_level);

		int average_refinement_level = (maximum_refinement_level + minimum_refinement_level) / 2;
		uint64_t id = get_cell_from_indices(x_index, y_index, z_index, average_refinement_level);

		// use binary search recursively (assumes that max_refinement_level_difference == 1)
		if (cell_process.count(id) == 0) {
			// doesn't exist, search the bin of smaller refinement_level values
			if (average_refinement_level > minimum_refinement_level) {
				return get_cell_from_indices(x_index, y_index, z_index, minimum_refinement_level, average_refinement_level - 1);
			} else {
				// nothing left to search
				return 0;
			}
		} else {
			// does exist, search the bin of larger refinement_level values
			if (average_refinement_level < maximum_refinement_level) {
				uint64_t larger_refinement_value_cell = get_cell_from_indices(x_index, y_index, z_index, average_refinement_level + 1, maximum_refinement_level);
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
		assert(x_index < x_length * int(pow(2, max_refinement_level)));
		assert(y_index < y_length * int(pow(2, max_refinement_level)));
		assert(z_index < z_length * int(pow(2, max_refinement_level)));
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
	/*std::vector<uint64_t> get_all_children(const uint64_t id)
	{
		assert(id);
		assert(this->cell_process.count(id) > 0);

		std::vector<uint64_t> children;
		children.reserve(8);

		// given cell cannot have children
		int refinement_level = get_refinement_level(id);
		if (refinement_level >= max_refinement_level) {
			return children;
		}

		uint64_t x_index = get_x_index(id);
		uint64_t y_index = get_y_index(id);
		uint64_t z_index = get_z_index(id);

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
	}*/


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

		// remove the cell from this process
		dccrg_instance->cells.erase(cell);
		dccrg_instance->neighbours.erase(cell);

		// record that the cell doesn't exist on this process anymore
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
		dccrg_instance->neighbours[cell] = dccrg_instance->get_neighbours_internal(cell);

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
