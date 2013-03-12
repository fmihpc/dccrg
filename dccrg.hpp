/*
A distributed cartesian cell-refinable grid.

Copyright 2009, 2010, 2011, 2012, 2013 Finnish Meteorological Institute

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


/*!
\mainpage Distributed Cartesian Cell-Refinable Grid.

\section intro_sec Introduction
dccrg is a grid library for simulations using the finite volume method.
See the examples directory for some simple examples and the tests directory
for more advanced usage of dccrg.
*/


#include "algorithm"
#include "boost/array.hpp"
#include "boost/foreach.hpp"
#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
#include "boost/mpi.hpp"
#endif
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "cstdio"
#include "cstdlib"
#include "fstream"
#include "functional"
#include "limits"
#include "mpi.h"
#include "stdint.h"
#include "utility"
#include "vector"
#include "zoltan.h"

#ifdef USE_SFC
#include "sfc++.hpp"
#endif


/*
If compilation fails with a complaint about
MPI_UNSIGNED_LONG_LONG try replacing it with
MPI_UNSIGNED_LONG in the following:
*/
#ifndef MPI_UINT64_T
#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
#endif


#include "dccrg_index.hpp"
#include "dccrg_mpi_support.hpp"
#include "dccrg_types.hpp"
#include "dccrg_cartesian_geometry.hpp"


/*!
Namespace where all dccrg classes, functions, etc are defined.
*/
namespace dccrg
{

static const int
	/*! @var
	*/

	/*!
	Id of the default neighborhood created when dccrg is initialized
	\see Dccrg::initialize() Dccrg::add_neighborhood()
	*/
	default_neighborhood_id = -0xdcc,

	/*!
	This bit is set for a cell that does not consider any cell as a
	neighbor and is not considered as a neighbor by any cell
	\see Dccrg::get_cells()
	*/
	has_no_neighbor = 0,

	/*!
	This bit is set for a cell that considers a cell on this
	process as a neighbor
	\see Dccrg::get_cells()
	*/
	has_local_neighbor_of = (1 << 0),

	/*!
	This bit is set for a cell that is considered as a neighbor
	by a cell on this process
	\see Dccrg::get_cells()
	*/
	has_local_neighbor_to = (1 << 1),

	/*!
	This bit is set for a cell that considers a cell on another
	process as a neighbor
	\see Dccrg::get_cells()
	*/
	has_remote_neighbor_of = (1 << 2),

	/*!
	This bit is set for a cell that c is considered as a neighbor
	by a cell on another process
	\see Dccrg::get_cells()
	*/
	has_remote_neighbor_to = (1 << 3),

	/*!
	This bit is set for a cell which is both dccrg::has_local_neighbor_of
	and dccrg::has_local_neighbor_to
	\see Dccrg::get_cells()
	*/
	has_local_neighbor_both = has_local_neighbor_of | has_local_neighbor_to,

	/*!
	This bit is set for a cell which is both dccrg::has_remote_neighbor_of
	and dccrg::has_remote_neighbor_to
	\see Dccrg::get_cells()
	*/
	has_remote_neighbor_both = has_remote_neighbor_of | has_remote_neighbor_to;



template <
	class Cell_Data,
	class Geometry = Cartesian_Geometry
> class Dccrg : public Geometry
{

public:

	/*!
	Helper type for iterating over local cells and their data using BOOST_FOREACH
	\see
	operator[]()
	begin()
	get_cells()
	*/
	typedef typename std::pair<const uint64_t&, const Cell_Data&> cell_and_data_pair_t;

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
	Creates an instance of the grid based on another instance of the grid.

	Call this with all processes unless you know what you are doing.
	Must not be used while the instance being copied is updating remote
	neighbor data or balancing the load.
	The Cell_Data class can differ between the two grids but the geometry
	must be the same.
	The following data is not included in the new dccrg instance:
		- refined/unrefined cells and their Cell_Data
		- Cell_Data of copies of remote neighbors
		- everything related to load balancing or remote neighbor updates
	*/
	template<class Other_Cell_Data> Dccrg(const Dccrg<Other_Cell_Data, Geometry>& other) :
		initialized(other.get_initialized()),
		neighborhood_length(other.get_neighborhood_length()),
		max_tag(other.get_max_tag()),
		max_ref_lvl_diff(other.get_max_ref_lvl_diff()),
		send_single_cells(other.get_send_single_cells()),
		comm(other.get_communicator()),
		rank(other.get_rank()),
		comm_size(other.get_comm_size()),
        #ifdef DCCRG_TRANSFER_USING_BOOST_MPI
        boost_comm(other.get_boost_comm()),
        #endif
		neighbors(other.get_cell_neighbor()),
		neighborhood_of(other.get_neighborhood_of()),
		neighborhood_to(other.get_neighborhood_to()),
		user_hood_of(other.get_user_hood_of()),
		user_hood_to(other.get_user_hood_to()),
		neighbors_to(other.get_all_neighbors_to()),
		user_neigh_of(other.get_all_user_neigh_of()),
		user_neigh_to(other.get_all_user_neigh_to()),
		cell_process(other.get_cell_process()),
		local_cells_on_process_boundary(other.get_local_cells_on_process_boundary_internal()),
		remote_cells_on_process_boundary(other.get_remote_cells_on_process_boundary_internal()),
		user_local_cells_on_process_boundary(other.get_user_local_cells_on_process_boundary()),
		user_remote_cells_on_process_boundary(other.get_user_remote_cells_on_process_boundary()),
		cells_to_send(other.get_cells_to_send()),
		cells_to_receive(other.get_cells_to_receive()),
		user_neigh_cells_to_send(other.get_user_neigh_cells_to_send()),
		user_neigh_cells_to_receive(other.get_user_neigh_cells_to_receive()),
		pin_requests(other.get_pin_requests()),
		new_pin_requests(other.get_new_pin_requests()),
		processes_per_part(other.get_processes_per_part()),
		partitioning_options(other.get_partitioning_options()),
		no_load_balancing(other.get_no_load_balancing()),
		reserved_options(other.get_reserved_options()),
		cell_weights(other.get_cell_weights()),
		neighbor_processes(other.get_neighbor_processes()),
		balancing_load(other.get_balancing_load())
	{
		if (other.get_balancing_load()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Copy constructor called while the instance being copied is balancing load"
				<< std::endl;
			abort();
		}

		// copy grid geometry
		// FIXME: support geometries other than constant
		if (!this->set_geometry(
			other.get_length_x(),
			other.get_length_y(),
			other.get_length_z(),
			other.get_start_x(),
			other.get_start_y(),
			other.get_start_z(),
			other.get_unrefined_cell_length_x(),
			other.get_unrefined_cell_length_y(),
			other.get_unrefined_cell_length_z()
		)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't copy geometry when copy constructing"
				<< std::endl;
			abort();
		}

		// maximum refinement level
		if (!this->set_maximum_refinement_level(other.get_maximum_refinement_level())) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't set maximum refinement level when copy constructing"
				<< std::endl;
			abort();
		}

		// periodicity
		for (size_t i = 0; i < 3; i++) {
			this->set_periodicity(i, other.is_periodic(i));
		}

		// zoltan
		this->zoltan = Zoltan_Copy(other.get_zoltan());

		// default construct Other_Cell_Data of local cells
		for (typename boost::unordered_map<uint64_t, Other_Cell_Data>::const_iterator
			cell_item = other.get_cell_data().begin();
			cell_item != other.get_cell_data().end();
			cell_item++
		) {
			this->cells[cell_item->first];
		}
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

	neighborhood_length:
		- Determines which cells are considered neighbors.
		- When calculating the neighbors of a given cell a cube of length
		  2 * neighborhood_length + 1 in every direction is considered, centered
		  at the cell for which neighbors are being calculated.
		- The unit lenght of the cube is the cell for which neighbors are being calculated.
		- If neighborhood_length == 0, only cells (or children within the volume of
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
		const unsigned int given_neighborhood_length,
		const int maximum_refinement_level = -1,
		const bool periodic_in_x = false,
		const bool periodic_in_y = false,
		const bool periodic_in_z = false,
		const uint64_t sfc_caching_batches = 1
	) {
		if (this->initialized) {
			std::cerr << "Initialize function called for an already initialized dccrg" << std::endl;
			// TODO: throw an exception instead?
			abort();
		}

		this->balancing_load = false;
		this->send_single_cells = false;

		if (sfc_caching_batches == 0) {
			std::cerr << "sfc_caching_batches must be > 0" << std::endl;
			abort();
		}

		int ret_val = -1;

		if (MPI_Comm_dup(given_comm, &this->comm) != MPI_SUCCESS) {
			std::cerr << "Couldn't duplicate given communicator" << std::endl;
			abort();
		}

		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
		this->boost_comm = boost::mpi::communicator(this->comm, boost::mpi::comm_attach);
		#endif

		int temp_size = 0;
		if (MPI_Comm_size(this->comm, &temp_size) != MPI_SUCCESS) {
			std::cerr << "Couldn't get size of communicator" << std::endl;
			abort();
		}
		if (temp_size < 0) {
			std::cerr << "Negative MPI comm size not supported: " << temp_size << std::endl;
			abort();
		}
		this->comm_size = (uint64_t) temp_size;

		int temp_rank = 0;
		if (MPI_Comm_rank(this->comm, &temp_rank) != MPI_SUCCESS) {
			std::cerr << "Couldn't get rank for communicator" << std::endl;
			abort();
		}
		if (temp_rank < 0) {
			std::cerr << "Negative MPI rank not supported: " << temp_rank << std::endl;
			abort();
		}
		this->rank = (uint64_t) temp_rank;

		// get maximum tag value
		int attr_flag = -1, *attr = NULL;
		ret_val = MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &attr, &attr_flag);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't get MPI_TAG_UB: " << Error_String()(ret_val)
				<< std::endl;
			abort();
		}
		if (attr == NULL) {
			// guaranteed by MPI
			this->max_tag = 32767;
		} else {
			this->max_tag = (unsigned int) *attr;
		}

		this->max_ref_lvl_diff = 1;

		/*
		Setup Zoltan
		*/
		MPI_Comm temp; // give a separate comminucator to zoltan
		if (MPI_Comm_dup(this->comm, &temp) != MPI_SUCCESS) {
			std::cerr << "Couldn't duplicate communicator for Zoltan" << std::endl;
			abort();
		}
		this->zoltan = Zoltan_Create(temp);
		if (this->zoltan == NULL) {
			std::cerr << "Zoltan_Create failed"  << std::endl;
			abort();
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

		/*
		Set reserved options
		*/
		// 0 because Zoltan crashes in hierarchial with larger values
		Zoltan_Set_Param(this->zoltan, "EDGE_WEIGHT_DIM", "0");
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

		// set the grids callback functions in Zoltan
		Zoltan_Set_Num_Obj_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::get_number_of_cells,
			this
		);

		Zoltan_Set_Obj_List_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_cell_list,
			this
		);

		Zoltan_Set_Num_Geom_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::get_grid_dimensionality,
			NULL);

		Zoltan_Set_Geom_Multi_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_with_cell_coordinates,
			this
		);

		Zoltan_Set_Num_Edges_Multi_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_number_of_neighbors_for_cells,
			this
		);

		Zoltan_Set_Edge_List_Multi_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_neighbor_lists,
			this
		);

		Zoltan_Set_HG_Size_CS_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_number_of_hyperedges,
			this
		);

		Zoltan_Set_HG_CS_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_hyperedge_lists,
			this
		);

		Zoltan_Set_HG_Size_Edge_Wts_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_number_of_edge_weights,
			this
		);

		Zoltan_Set_HG_Edge_Wts_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::fill_edge_weights,
			this
		);

		Zoltan_Set_Hier_Num_Levels_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::get_number_of_load_balancing_hierarchies,
			this
		);

		Zoltan_Set_Hier_Part_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::get_part_number,
			this
		);

		Zoltan_Set_Hier_Method_Fn(
			this->zoltan,
			&Dccrg<Cell_Data, Geometry>::set_partitioning_options,
			this
		);


		/*
		Set grid parameters
		*/

		this->set_periodicity(0, periodic_in_x);
		this->set_periodicity(1, periodic_in_y);
		this->set_periodicity(2, periodic_in_z);

		// set / check neighborhood_of
		this->neighborhood_length = given_neighborhood_length;
		if (this->neighborhood_length == 0) {

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

			for (int
				z = -this->neighborhood_length;
				(unsigned int) abs(z) < this->neighborhood_length + 1;
				z++
			)
			for (int
				y = -this->neighborhood_length;
				(unsigned int) abs(y) < this->neighborhood_length + 1;
				y++
			)
			for (int
				x = -this->neighborhood_length;
				(unsigned int) abs(x) < this->neighborhood_length + 1;
				x++
			) {
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Couldn't set maximum refinement level to " << maximum_refinement_level
				<< std::endl;
			abort();
		}

		// TODO: check that the last index in the grid in every direction is less than error_index


		// create unrefined cells
		const uint64_t grid_length = this->length_x * this->length_y * this->length_z;
		uint64_t cells_per_process = 0;
		if (grid_length < this->comm_size) {
			cells_per_process = 1;
		} else if (grid_length % this->comm_size > 0) {
			cells_per_process = grid_length / this->comm_size + 1;
		} else {
			cells_per_process = grid_length / this->comm_size;
		}

		// some processes get fewer cells if grid size not divisible by this->comm_size
		uint64_t procs_with_fewer = cells_per_process * this->comm_size - grid_length;

		#ifndef USE_SFC

		uint64_t cell_to_create = 1;
		for (uint64_t process = 0; process < this->comm_size; process++) {

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

		#ifdef DEBUG
		if (cell_to_create != grid_length + 1) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of cells created: " << cell_to_create - 1
				<< ", should be " << grid_length
				<< std::endl;
			abort();
		}
		#endif

		#else

		const dccrg::Types<3>::indices_t length = {{
			this->length_x,
			this->length_y,
			this->length_z
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

		if (sfc_index != grid_length) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Process " << this->rank
				<< ": Incorrect number of cells created: " << sfc_index
				<< ", should be " << grid_length
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Neighbor lists are inconsistent"
				<< std::endl;
			abort();
		}
		#endif

		this->update_remote_neighbor_info();
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
	Returns true if given cell's neighbor types match given criterion, false otherwise.

	Returns false if:
		- given neighborhood doesn't exist
		- given cell doesn't exist
		- given cell is on another process

	\see get_cells()
	*/
	bool is_neighbor_type_match(
		const uint64_t cell,
		const std::vector<int>& criteria,
		const bool exact_match,
		const int neighborhood_id
	) const {

		if (cell == error_cell) {
			return false;
		}

		if (this->cell_process.count(cell) == 0) {
			return false;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return false;
		}

		if (neighborhood_id != default_neighborhood_id
		&& this->user_hood_of.count(neighborhood_id) == 0) {
			return false;
		}


		int neighbor_types = 0;

		const std::vector<uint64_t>& neighs_of
			= (neighborhood_id == default_neighborhood_id)
			? this->neighbors.at(cell)
			: this->user_neigh_of.at(neighborhood_id).at(cell);

		BOOST_FOREACH(const uint64_t neighbor, neighs_of) {
			if (neighbor == error_cell) {
				continue;
			}

			if (this->is_local(neighbor)) {
				neighbor_types |= has_local_neighbor_of;
			} else {
				neighbor_types |= has_remote_neighbor_of;
			}
		}

		const std::vector<uint64_t>& neighs_to
			= (neighborhood_id == default_neighborhood_id)
			? this->neighbors_to.at(cell)
			: this->user_neigh_to.at(neighborhood_id).at(cell);

		BOOST_FOREACH(const uint64_t neighbor, neighs_to) {
			if (neighbor == error_cell) {
				continue;
			}

			if (this->is_local(neighbor)) {
				neighbor_types |= has_local_neighbor_to;
			} else {
				neighbor_types |= has_remote_neighbor_to;
			}
		}


		if (exact_match) {
			BOOST_FOREACH(const int criterion, criteria) {
				if (neighbor_types == criterion) {
					return true;
				}
			}
		} else {
			// with inexact matching all criteria can be merged into one
			int merged_criteria = 0;
			BOOST_FOREACH(const int criterion, criteria) {
				merged_criteria |= criterion;
			}

			if ((neighbor_types & merged_criteria) > 0) {
				return true;
			}
		}

		return false;
	}


	/*!
	Returns cells without children on this process fulfilling given criteria.

	By default returns all local cells.	Otherwise only those local cells are
	returned which match one or more of the given criteria.

	A list of criteria can be constructed in-place with boost::assign::list_of
	when calling this function, see further down for examples.

	Criteria represent which type of neighbors a cell must (and possibly must not)
	have in order to be returned. Each criteria is a bitmask of possibly several
	neighbor types.

	If exact_match = false then a cell is returned as long as it has at least
	one neighbor of the type given by any of the given criteria. If
	exact_match == true a cell is returned only if it has neighbors of the
	type(s) in any given criteria and no other types of neighbors in that
	criteria.

	For example to only get cells which:
		- are not on the process boundary, i.e. don't have neighbors or only
		  consider cells on this process as a neighbor or are only considered
		  as a neighbor by cells on this process, give exact_match = true and
		  \verbatim
		  list_of
		  	(has_no_neighbor)
		  	(has_local_neighbor_of)
		  	(has_local_neighbor_to)
		  	(has_local_neighbor_both)
		  \endverbatim
		- are on the process boundary, i.e. consider a cell on another process
		  as a neighbor or are considered as a neighbor by a cell on another
		  process, give exact_match = false and
		  \verbatim
		  list_of(has_remote_neighbor_both)
		  \endverbatim
		- are considered as a neighbor by at least one local and at least one
		  remote cell but do not consider any local or remote cells as their
		  neighbors, give exact_match = true and
		  \verbatim
		  list_of(has_local_neighor_to | has_remote_neighbor_to)
		  \endverbatim
		- are considered as a neighbor by a local cell or by a remote cell
		  but do not consider any local of remote cells as their neighbors, give
		  exact_match = true and
		  \verbatim
		  list_of(has_local_neighor_to)(has_remote_neighbor_to)
		  \endverbatim

	Cells' neighbors from the given neighborhood are used when checking for
	neighbor types.

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort.

	Returns nothing if:
		- this process doesn't have any cells
		- given neighborhood doesn't exist

	\see
	get_local_cells_on_process_boundary()
	get_local_cells_not_on_process_boundary()
	get_remote_cells_on_process_boundary()
	*/
	std::vector<uint64_t> get_cells(
		const std::vector<int>& criteria = std::vector<int>(),
		const bool exact_match = false,
		const int neighborhood_id = default_neighborhood_id,
		const bool sorted = false
	) const {
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " get_cells() must not be called while balancing load"
				<< std::endl;
			abort();
		}

		std::vector<uint64_t> ret_val;

		if (neighborhood_id != default_neighborhood_id
		&& this->user_hood_of.count(neighborhood_id) == 0) {
			return ret_val;
		}

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			const uint64_t cell = item.first;

			#ifdef DEBUG
			if (this->cell_process.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell
					<< " shouldn't exist"
					<< std::endl;
				abort();
			}

			if (this->cell_process.at(cell) != this->rank) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Cell " << cell
					<< " should be on process " << this->cell_process.at(cell)
					<< std::endl;
				abort();
			}

			const uint64_t child = this->get_child(cell);
			if (child == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Child == 0 for cell " << cell
					<< std::endl;
				abort();
			}

			if (child != cell) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< ": Cell " << cell
					<< " has a child"
					<< std::endl;
				abort();
			}
			#endif

			if (criteria.size() == 0) {
				ret_val.push_back(cell);
			}

			if (this->is_neighbor_type_match(
				cell,
				criteria,
				exact_match,
				neighborhood_id
			)) {
				ret_val.push_back(cell);
			}
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
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
	Data is written in native endian format.
	Requires at least as much additional memory as local cell data.
	Does nothing if DCCRG_TRANSFER_USING_BOOST_MPI was defined when
	compiling and returns false.

	During this function the receiving process given to the cells'
	mpi_datatype function is -1 and receiving == false.
	*/
	bool write_grid(const std::string& name, MPI_Offset& start_offset)
	{
		// TODO: use nonblocking versions of ...write_at_all
		/*
		File format:

		uint8_t * start_offset, data skipped by this function
		uint64_t  number of cells
		uint64_t  id of 1st cell
		uint64_t  start of data of 1st cell in bytes
		uint64_t  id of 2nd cell
		uint64_t  start of data...
		...
		uint64_t  start of data of last cell in bytes
		uint8_t*N data of 1st cell
		uint8_t*M data of 2nd cell
		...
		*/

		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		return false;

		#else

		int ret_val = -1;

		// MPI_File_open wants a non-constant string
		char* name_c_string = new char [name.size() + 1];
		strncpy(name_c_string, name.c_str(), name.size() + 1);

		MPI_File outfile;

		ret_val = MPI_File_open(
			this->comm,
			name_c_string,
			MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL,
			&outfile
		);

		delete [] name_c_string;

		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't open file " << name
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		// TODO: write an endianness check here
		// ...0x0123456789101213...

		uint64_t number_of_cells = this->cells.size();

		// write the total number of cells that will be written
		const uint64_t total_number_of_cells
			= All_Reduce()(number_of_cells, this->comm);

		if (this->rank == 0) {
			ret_val = MPI_File_write_at(
				outfile,
				start_offset,
				(void*) &total_number_of_cells,
				1,
				MPI_UINT64_T,
				MPI_STATUS_IGNORE
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " Couldn't write cell list to file " << name
					<< ": " << Error_String()(ret_val)
					<< std::endl;
				return false;
			}
		}

		// contiguous memory version of cell list needed by MPI_Write_...
		std::vector<uint64_t> cells_to_write = this->get_cells();
		if (cells_to_write.size() != number_of_cells) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Number of cells is inconsistent: " << number_of_cells
				<< ", should be " << cells_to_write.size()
				<< std::endl;
			abort();
		}

		/*
		Get datatypes etc. of local cell data in memory so byte offsets
		in output file can be calculated in order to write cell list(s)
		*/

		std::vector<void*> addresses(number_of_cells, NULL);
		std::vector<int> counts(number_of_cells, -1);
		std::vector<MPI_Datatype> datatypes(number_of_cells, MPI_DATATYPE_NULL);

		for (size_t i = 0; i < number_of_cells; i++) {
			const uint64_t cell = cells_to_write[i];
			this->cells.at(cell).mpi_datatype(
				addresses[i],
				counts[i],
				datatypes[i],
				cell,
				(int) this->rank,
				-1,
				false
			);
		}

		// displacements of cell data in memory
		std::vector<MPI_Aint> memory_displacements(number_of_cells, 0);
		for (size_t i = 0; i < number_of_cells; i++) {
			memory_displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
		}

		// create datatype representing all local cell data in memory
		MPI_Datatype memory_datatype;
		if (number_of_cells == 0) {

			memory_datatype = MPI_BYTE;

		} else {

			ret_val = MPI_Type_create_struct(
				number_of_cells,
				&counts[0],
				&memory_displacements[0],
				&datatypes[0],
				&memory_datatype
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " Couldn't create final datatype: "
					<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}

			ret_val = MPI_Type_commit(&memory_datatype);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Process " << this->rank
					<< " Couldn't commit final datatype: "
					<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}

		}


		/*
		Calculate where each local cells' data starts in the file
		*/
		uint64_t current_byte_offset = 0;

		std::vector<MPI_Aint> file_displacements(number_of_cells, 0);
		for (size_t i = 0; i < number_of_cells; i++) {

			file_displacements[i] += current_byte_offset;

			int current_bytes;
			MPI_Type_size(datatypes[i], &current_bytes);
			if (current_bytes < 0) {
				std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
				abort();
			}

			current_byte_offset += (uint64_t) current_bytes;
		}

		// tell everyone how many bytes everyone will write
		std::vector<uint64_t> all_number_of_bytes(this->comm_size, 0);
		ret_val = MPI_Allgather(
			&current_byte_offset,
			1,
			MPI_UINT64_T,
			&(all_number_of_bytes[0]),
			1,
			MPI_UINT64_T,
			comm
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Process " << this->rank
				<< " MPI_Allgather failed: "
				<< Error_String()(ret_val)
				<< std::endl;
			abort();
		}

		// calculate where local cell data starts in file
		uint64_t cell_data_start
			= (uint64_t) start_offset
			+ sizeof(uint64_t)
			+ 2 * total_number_of_cells * sizeof(uint64_t);

		for (size_t i = 0; i < (size_t) this->rank; i++) {
			cell_data_start += all_number_of_bytes[i];
		}

		// make file displacements relative to start of file
		for (size_t i = 0; i < number_of_cells; i++) {
			file_displacements[i] += cell_data_start;
		}


		/*
		Write cell list(s)
		*/

		// tell all processes how many cells each one will write
		std::vector<uint64_t> all_number_of_cells(this->comm_size, 0);
		ret_val = MPI_Allgather(
			&number_of_cells,
			1,
			MPI_UINT64_T,
			&(all_number_of_cells[0]),
			1,
			MPI_UINT64_T,
			comm
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " MPI_Allgather failed: " << Error_String()(ret_val)
				<< std::endl;
			abort();
		}

		// calculate where local cell list will begin in output file
		uint64_t cell_list_start = (uint64_t) start_offset + sizeof(uint64_t);
		for (size_t i = 0; i < (size_t) this->rank; i++) {
			cell_list_start += all_number_of_cells[i] * 2 * sizeof(uint64_t);
		}

		// write cell + data displacement list
		std::vector<uint64_t> cells_and_data_displacements(2 * number_of_cells, 0);
		for (size_t i = 0; i < number_of_cells; i++) {
			cells_and_data_displacements[2 * i] = cells_to_write[i];
			cells_and_data_displacements[2 * i + 1] = file_displacements[i];
		}

		// give a valid buffer to ...write_at_all even if no cells to write
		if (number_of_cells == 0) {
			cells_and_data_displacements.push_back(error_cell);
		}

		ret_val = MPI_File_write_at_all(
			outfile,
			(MPI_Offset) cell_list_start,
			(void*) &cells_and_data_displacements[0],
			2 * number_of_cells,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't write cell and displacement list to file " << name
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		/*
		Write cell data
		*/

		// set as file view the datatype representing local cell data in file
		MPI_Datatype file_datatype;

		// wihtout local cells create an empty view and datatype
		if (number_of_cells == 0) {

			ret_val = MPI_Type_contiguous(0, MPI_BYTE, &file_datatype);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Process " << this->rank
					<< " Couldn't create an empty datatype for file view: "
					<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}

			ret_val = MPI_Type_commit(&file_datatype);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Process " << this->rank
					<< " Couldn't commit an empty datatype for file view: "
					<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}

		} else {

			ret_val = MPI_Type_create_struct(
				number_of_cells,
				&counts[0],
				&file_displacements[0],
				&datatypes[0],
				&file_datatype
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " Couldn't create final datatype: "
					<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}

			ret_val = MPI_Type_commit(&file_datatype);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< "Process " << this->rank
					<< " Couldn't commit final datatype: "
					<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}

		}

		ret_val = MPI_File_set_view(
			outfile,
			0,
			MPI_BYTE,
			file_datatype,
			const_cast<char*>("native"),
			MPI_INFO_NULL
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Process " << this->rank
				<< " MPI_File_set_view failed for cell data: "
				<< Error_String()(ret_val)
				<< std::endl;
			abort();
		}

		// write cell data
		ret_val = MPI_File_write_at_all(
			outfile,
			0,
			addresses[0],
			1,
			memory_datatype,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "Process " << this->rank
				<< "MPI_File_write_at_all failed when writing cell data: "
				<< Error_String()(ret_val)
				<< std::endl;
			abort();
		}


		/*
		Deallocate datatypes
		*/

		ret_val = MPI_Type_free(&memory_datatype);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't free datatype for cell data in memory: "
				<< Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		ret_val = MPI_Type_free(&file_datatype);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't free datatype for cell data in file: "
				<< Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		BOOST_FOREACH(MPI_Datatype& datatype, datatypes) {
			if (!Is_Named_Datatype()(datatype)) {
				ret_val = MPI_Type_free(&datatype);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< " Couldn't free user defined datatype: "
						<< Error_String()(ret_val)
						<< std::endl;
					return false;
				}
			}
		}

		ret_val = MPI_File_close(&outfile);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't close file " << name
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		return true;
		#endif
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

	Does nothing and returns false if DCCRG_TRANSFER_USING_BOOST_MPI
	was defined when compiling.

	During this function the sending process given to the cells'
	mpi_datatype function is -1 and receiving == true.

	TODO: Reads at most number_of_cells number of cell data offsets at a time,
	give a smaller number if all cell ids won't fit	into memory at once.
	*/
	bool read_grid(
		const std::string& name,
		const MPI_Offset start_offset,
		const uint64_t /*number_of_cells*/ = ~uint64_t(0)
	) {
		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		return false;

		#else

		int ret_val = -1;

		// MPI_File_open wants a non-constant string
		char* name_c_string = new char [name.size() + 1];
		strncpy(name_c_string, name.c_str(), name.size() + 1);

		MPI_File infile;

		ret_val = MPI_File_open(
			this->comm,
			name_c_string,
			MPI_MODE_RDONLY,
			MPI_INFO_NULL,
			&infile
		);

		delete [] name_c_string;

		if (ret_val != MPI_SUCCESS) {
			std::cerr << "Couldn't open file " << name_c_string
				<< ": " << Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		// TODO: put an endianness check here
		// ...0x0123456789101213...

		// read the total number of cells in the file
		uint64_t total_number_of_cells = 0;
		ret_val = MPI_File_read_at_all(
			infile,
			start_offset,
			&total_number_of_cells,
			1,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << "Couldn't read total number of cells" << std::endl;
			return false;
		}

		// read cells and data displacements
		std::vector<uint64_t> all_cells_and_data_displacements(2 * total_number_of_cells, error_cell);
		ret_val = MPI_File_read_at_all(
			infile,
			start_offset + sizeof(uint64_t),
			&all_cells_and_data_displacements[0],
			2 * total_number_of_cells,
			MPI_UINT64_T,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << "Couldn't read number of cells" << std::endl;
			return false;
		}

		#ifdef DEBUG
		// check that proper cells were read and properly
		for (size_t i = 0; i < all_cells_and_data_displacements.size(); i += 2) {
			const uint64_t cell = all_cells_and_data_displacements[i];
			if (cell == error_cell) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid cell in cell list at index " << i << ": " << cell
					<< std::endl;
				abort();
			}
		}
		#endif

		// remove all but local cell data displacements
		std::vector<std::pair<uint64_t, uint64_t> > cells_and_data_displacements;
		cells_and_data_displacements.reserve(this->cells.size());

		for (size_t i = 0; i < all_cells_and_data_displacements.size(); i += 2) {
			const uint64_t
				cell = all_cells_and_data_displacements[i],
				offset = all_cells_and_data_displacements[i + 1];

			if (this->cell_overlaps_local(cell)) {
				cells_and_data_displacements.push_back(std::make_pair(cell, offset));
			}
		}
		all_cells_and_data_displacements.clear();

		// refine the grid to create cells that exist in the file
		std::vector<uint64_t> final_cells;
		final_cells.reserve(cells_and_data_displacements.size());
		for (std::vector<std::pair<uint64_t, uint64_t> >::const_iterator
			item = cells_and_data_displacements.begin();
			item != cells_and_data_displacements.end();
			item++
		) {
			final_cells.push_back(item->first);
		}

		if (!this->load(final_cells)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't load grid"
				<< std::endl;
			abort();
		}
		final_cells.clear();

		const uint64_t number_of_cells = cells_and_data_displacements.size();

		#ifdef DEBUG
		if (number_of_cells != cells_and_data_displacements.size()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Incorrect number of cell data displacements: "
				<< cells_and_data_displacements.size()
				<< ", should be " << number_of_cells
				<< std::endl;
			abort();
		}
		#endif


		// datatypes etc. of local cell data in memory
		std::vector<void*> addresses(number_of_cells, NULL);
		std::vector<int> counts(number_of_cells, -1);
		std::vector<MPI_Datatype> datatypes(number_of_cells, MPI_DATATYPE_NULL);
		std::vector<MPI_Aint>
			memory_displacements(number_of_cells, 0),
			file_displacements(number_of_cells, 0);

		// set file view representing local cell data
		MPI_Datatype file_datatype;

		if (number_of_cells == 0) {

			MPI_Type_contiguous(0, MPI_BYTE, &file_datatype);

		} else {

			// get datatype info from local cells in memory
			for (uint64_t i = 0; i < number_of_cells; i++) {
				const uint64_t cell = cells_and_data_displacements[i].first;
				this->cells.at(cell).mpi_datatype(
					addresses[i],
					counts[i],
					datatypes[i],
					cell,
					-1,
					(int) this->rank,
					true
				);
			}

			// displacements for cell data in memory are relative to first cell's data
			for (size_t i = 0; i < number_of_cells; i++) {
				memory_displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
			}

			// displacements for cell data in file are relative to start of file
			for (uint64_t i = 0; i < number_of_cells; i++) {
				file_displacements[i] = cells_and_data_displacements[i].second;
			}

			ret_val = MPI_Type_create_struct(
				number_of_cells,
				&counts[0],
				&file_displacements[0],
				&datatypes[0],
				&file_datatype
			);
			if (ret_val != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " Couldn't create datatype for file view: "<< Error_String()(ret_val)
					<< std::endl;
				abort();
			}
		}

		ret_val = MPI_Type_commit(&file_datatype);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't commit datatype for file view: " << Error_String()(ret_val)
				<< std::endl;
			abort();
		}

		ret_val = MPI_File_set_view(
			infile,
			0,
			MPI_BYTE,
			file_datatype,
			const_cast<char*>("native"),
			MPI_INFO_NULL
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't set file view for cell data: "
				<< Error_String()(ret_val)
				<< std::endl;
			abort();
		}

		// create a datatype representing local cell data in memory
		MPI_Datatype memory_datatype;

		if (number_of_cells == 0) {

			MPI_Type_contiguous(0, MPI_BYTE, &memory_datatype);
			if (MPI_Type_commit(&memory_datatype) != MPI_SUCCESS) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " Couldn't commit datatype for file view"
					<< std::endl;
				abort();
			}

		} else {

			if (MPI_Type_create_struct(
				number_of_cells,
				&counts[0],
				&memory_displacements[0],
				&datatypes[0],
				&memory_datatype) != MPI_SUCCESS
			) {
				std::cerr << "Process " << this->rank
					<< " Couldn't create datatype for local cells"
					<< std::endl;
				abort();
			}

			if (MPI_Type_commit(&memory_datatype) != MPI_SUCCESS) {
				std::cerr << "Process " << this->rank
					<< " Couldn't commit datatype for file view"
					<< std::endl;
				abort();
			}
		}

		// give a valid buffer to ...read_at_all even if no cells to read
		if (number_of_cells == 0) {
			addresses.push_back((void*) &number_of_cells);
		}

		ret_val = MPI_File_read_at_all(
			infile,
			0,
			addresses[0],
			1,
			memory_datatype,
			MPI_STATUS_IGNORE
		);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " couldn't read local cell data: "
				<< Error_String()(ret_val)
				<< std::endl;
			abort();
		}


		/*
		Deallocate datatypes
		*/

		ret_val = MPI_Type_free(&memory_datatype);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't free datatype for cell data in memory: "
				<< Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		ret_val = MPI_Type_free(&file_datatype);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't free datatype for cell data in file: "
				<< Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		BOOST_FOREACH(MPI_Datatype& datatype, datatypes) {
			if (!Is_Named_Datatype()(datatype)) {
				ret_val = MPI_Type_free(&datatype);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< " Couldn't free user defined datatype: "
						<< Error_String()(ret_val)
						<< std::endl;
					return false;
				}
			}
		}

		ret_val = MPI_File_close(&infile);
		if (ret_val != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Process " << this->rank
				<< " Couldn't close input file: "
				<< Error_String()(ret_val)
				<< std::endl;
			return false;
		}

		return true;
		#endif
	}


	/*!
	Returns a begin const_iterator to the internal storage of local cells and their data.
	*/
	typename boost::unordered_map<uint64_t, Cell_Data>::const_iterator begin() const
	{
		return this->cells.begin();
	}

	/*!
	Returns an end const_iterator to the internal storage of local cells and their data.
	*/
	typename boost::unordered_map<uint64_t, Cell_Data>::const_iterator end() const
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
	Returns all cells in the grid that don't have children (e.g. leaf cells).

	Only those cells are returned which this process knows about, this
	might not include all cells that exist on all processes.

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.
	*/
	std::vector<uint64_t> get_all_cells(const bool sorted = false) const
	{
		std::vector<uint64_t> ret_val;
		ret_val.reserve(this->cell_process.size());

		for (boost::unordered_map<uint64_t, uint64_t>::const_iterator
			item = this->cell_process.begin();
			item != this->cell_process.end();
			item++
		) {

			const uint64_t child = this->get_child(item->first);

			if (child == item->first) {
				ret_val.push_back(item->first);
			}
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
	}


	/*!
	Returns a pointer to the user supplied data of given cell.

	The data of local cells is always available, including refined cells
	before the next call to stop_refining.
	The data of cells which are on other processes can also be available if:
		- the cells are neighbors to a local cell and remote neighbor data has been updated
		- the cells were unrefined and their parent is now a local cell
	*/
	Cell_Data* operator [] (const uint64_t cell) const
	{
		if (this->cells.count(cell) > 0) {
			return (Cell_Data*) &(this->cells.at(cell));
		} else if (this->remote_neighbors.count(cell) > 0) {
			return (Cell_Data*) &(this->remote_neighbors.at(cell));
		} else if (this->refined_cell_data.count(cell) > 0) {
			return (Cell_Data*) &(this->refined_cell_data.at(cell));
		} else if (this->unrefined_cell_data.count(cell) > 0) {
			return (Cell_Data*) &(this->unrefined_cell_data.at(cell));
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
	bool cell_overlaps_local(const uint64_t cell) const
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
	boost::unordered_set<uint64_t> get_cells_overlapping_local(const std::vector<uint64_t>& given_cells) const
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
		// get the global maximum refinement level of cells to be loaded
		boost::unordered_set<uint64_t> overlapping
			= this->get_cells_overlapping_local(given_cells);

		int local_max_ref_lvl_of_overlapping = 0;
		BOOST_FOREACH(const uint64_t cell, overlapping) {
			const int refinement_level = this->get_refinement_level(cell);
			local_max_ref_lvl_of_overlapping
				= std::max(refinement_level, local_max_ref_lvl_of_overlapping);
		}

		int max_ref_lvl_of_overlapping = 0, result;
		result = MPI_Allreduce(
			&local_max_ref_lvl_of_overlapping,
			&max_ref_lvl_of_overlapping,
			1,
			MPI_INT,
			MPI_MAX,
			this->comm
		);

		if (result != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< "MPI_Allreduce failed: " << Error_String()(result)
				<< std::endl;
			abort();
		}

		/*
		Starting from refinement level 0 refine each local cell
		that has a child in given_cells
		*/
		std::vector<boost::unordered_set<uint64_t> > cells_and_parents(max_ref_lvl_of_overlapping);

		// refine local cells recursively until all given_cells are created
		BOOST_FOREACH(const uint64_t cell, overlapping) {
			uint64_t current_child = cell;
			const int refinement_level = this->get_refinement_level(current_child);
			for (int i = refinement_level - 1; i >= 0; i--) {
				const uint64_t parent = this->get_parent_for_removed(current_child);
				cells_and_parents[i].insert(parent);
				current_child = parent;
			}
		}

		for (int
			current_ref_lvl = 0;
			current_ref_lvl < max_ref_lvl_of_overlapping;
			current_ref_lvl++
		) {
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

	Must be called by all processes.
	Creates a new cell partition, in other words decides which
	cells should be moved to which processes, and then moves
	the cells' data accordingly.

	If use_zoltan == true Zoltan will be used to create the new partition,
	otherwise only pin requests will move local cells to other processes.
	If Zoltan is used pin requests override decisions made by Zoltan.

	The following items are discarded after a call to this function:
		- cell weights
		- refines/unrefines after the last call to stop_refining()
		- the data of local copies of remote neighbors of local cells
	*/
	void balance_load(const bool use_zoltan = true)
	{
		this->initialize_balance_load(use_zoltan);
		this->continue_balance_load();
		this->finish_balance_load();
	}

	/*!
	Starts the procedure of moving cells between processes.

	Default constructs arriving cells.
	Does not transfer any cell data, use continue_balance_load() for that,
	or if cell data doesn't have to be transferred in several steps use
	balance_load() instead of this function.

	After calling this the functions get_cells_to_send() and
	get_cells_to_receive() can be used to query changes to the
	current partition.

	The next function to be called after this one (from those that
	must be called by all processes) must be either continue_balance_load()
	or finish_balance_load(). Information related to the other functions
	must not be queried after calling this, for example the remote neighbors
	of local cells, local cells, etc.
	*/
	void initialize_balance_load(const bool use_zoltan)
	{
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " initialize_balance_load(...) called the second time "
				<< "before calling finish_balance_load() first"
				<< std::endl;
			abort();
		}

		this->balancing_load = true;

		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info is not consistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " User data not consistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->make_new_partition(use_zoltan);

		// default construct user data of arriving cells
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			sender_item = this->cells_to_receive.begin();
			sender_item != this->cells_to_receive.end();
			sender_item++
		) {
			for (std::vector<std::pair<uint64_t, int> >::const_iterator
				cell_item = sender_item->second.begin();
				cell_item != sender_item->second.end();
				cell_item++
			) {
				this->cells[cell_item->first];
			}
		}

		// clear data related to remote neighbor updates, adaptation, etc.
		this->local_cells_on_process_boundary.clear();
		this->remote_cells_on_process_boundary.clear();
		this->user_local_cells_on_process_boundary.clear();
		this->user_remote_cells_on_process_boundary.clear();
		this->user_neigh_cells_to_send.clear();
		this->user_neigh_cells_to_receive.clear();
		this->remote_neighbors.clear();
		this->cells_to_refine.clear();
		this->refined_cell_data.clear();
		this->cells_to_unrefine.clear();
		this->unrefined_cell_data.clear();
		this->cells_not_to_unrefine.clear();
		this->cell_weights.clear();

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
	}

	/*!
	Transfers cell data between processes based on the new partition.

	Must be called by all processes and not before initialize_balance_load(...)
	has been called.

	The next function to be called after this one (from those that
	must be called by all processes) must be either this again or
	finish_balance_load().
	*/
	void continue_balance_load()
	{
		if (!this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " continue_balance_load() called without "
				<< "calling initialize_balance_load(...) first."
				<< std::endl;
			abort();
		}

		this->start_user_data_transfers(
			this->cells,
			this->cells_to_receive,
			this->cells_to_send
		);

		this->wait_user_data_transfer_receives(
		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
		this->cells, this->cells_to_receive
		#endif
		);

		this->wait_user_data_transfer_sends();
	}

	/*!
	Finishes the procedure of moving cells between processes.

	Must be called by all processes and not before initialize_balance_load(...)
	has been called.
	*/
	void finish_balance_load()
	{
		if (!this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " finish_balance_load() called without "
				<< "calling initialize_balance_load(...) first."
				<< std::endl;
			abort();
		}

		this->cells_to_send.clear();
		this->cells_to_receive.clear();

		/*
		Calculate where cells have migrated to update internal data structures
		*/

		// removed cells on all processes
		std::vector<uint64_t> temp_removed_cells(
			this->removed_cells.begin(),
			this->removed_cells.end()
		);

		std::vector<std::vector<uint64_t> > all_removed_cells;
		All_Gather()(temp_removed_cells, all_removed_cells, this->comm);

		// created cells on all processes
		std::vector<uint64_t> temp_added_cells(
			this->added_cells.begin(),
			this->added_cells.end()
		);

		std::vector<std::vector<uint64_t> > all_added_cells;
		All_Gather()(temp_added_cells, all_added_cells, this->comm);

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
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Grid is not consistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->pin_requests_succeeded()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Pin requests didn't succeed"
				<< std::endl;
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
		// also remote neighbor info of user neighborhoods
		for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
			item = this->user_hood_of.begin();
			item != this->user_hood_of.end();
			item++
		) {
			this->update_user_remote_neighbor_info(item->first);
		}

		this->recalculate_neighbor_update_send_receive_lists();

		#ifdef DEBUG
		if (!this->is_consistent()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " The grid is inconsistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_neighbors()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Neighbor lists are incorrect"
				<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_remote_neighbor_info()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info is not consistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}

		if (!this->verify_user_data()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " User data not consistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->added_cells.clear();
		this->removed_cells.clear();
		this->balancing_load = false;
	}


	/*!
	Updates the cell data of neighboring cells between processes.

	Cell data of any local cell that a cell on another process
	considers as a neighbor is sent to that process.
	Cell data of any cell on another process that a local cell
	considers as a neighbor is received from that process.
	Afterwards a copy of the remote cells' data is available
	through operator[].

	Data of any cell is only exchanged between any two processes
	once, even if a cell the neighbor of more than one cell on
	another process.

	The decision of which cells are neighbors is controlled by
	the neighborhood. By default the neighborhood with which
	this instance of dccrg was initialized is used.

	Returns true if successful and false otherwise (on one or more
	processes), for example if a neighborhood with the given id
	has not been set.

	Must be called simultaneously on all processes and with
	identical neigbhorhood_id.
	Must not be called while load balancing is underway.

	\see
	start_remote_neighbor_copy_updates()
	add_neighborhood()
	get_remote_cells_on_process_boundary()
	set_send_single_cells()
	*/
	bool update_copies_of_remote_neighbors(
		const int neighborhood_id = default_neighborhood_id
	) {
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " update_copies_of_remote_neighbors(...) called while balancing load"
				<< std::endl;
			abort();
		}

		bool ret_val = true;

		if (this->user_hood_of.count(neighborhood_id) == 0) {
			ret_val = false;
		}

		if (!this->start_remote_neighbor_copy_updates(neighborhood_id)) {
			ret_val = false;
		}

		if (!this->wait_remote_neighbor_copy_updates(neighborhood_id)) {
			ret_val = false;
		}

		return ret_val;
	}


	/*!
	An asynchronous version of update_copies_of_remote_neighbors().

	Starts remote neighbor data updates and returns immediately.

	\see
	update_copies_of_remote_neighbors()
	wait_remote_neighbor_copy_updates()
	*/
	bool start_remote_neighbor_copy_updates(
		const int neighborhood_id = default_neighborhood_id
	) {
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " start_remote_neighbor_data_update(...) called while balancing load"
				<< std::endl;
			abort();
		}

		bool ret_val = true;

		if (neighborhood_id == default_neighborhood_id) {
			return this->start_user_data_transfers(
				this->remote_neighbors,
				this->cells_to_receive,
				this->cells_to_send
			);
		}

		if (this->user_hood_of.count(neighborhood_id) == 0) {

			#ifdef DEBUG
			if (this->user_hood_to.count(neighborhood_id) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should not have id " << neighborhood_id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_of.count(neighborhood_id) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should not have id " << neighborhood_id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_to.count(neighborhood_id) > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should not have id " << neighborhood_id
					<< std::endl;
				abort();
			}
			#endif

			ret_val = false;
		}

		#ifdef DEBUG
		if (this->user_hood_to.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_of.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_to.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_send.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_receive.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should have id " << neighborhood_id
				<< std::endl;
			abort();
		}
		#endif

		if (!this->start_user_data_transfers(
			this->remote_neighbors,
			this->user_neigh_cells_to_receive.at(neighborhood_id),
			this->user_neigh_cells_to_send.at(neighborhood_id)
		)) {
			ret_val = false;
		}

		return ret_val;
	}


	/*!
	Finishes what start_remote_neighbor_data_update() started.

	\see
	start_remote_neighbor_copy_updates()
	*/
	bool wait_remote_neighbor_copy_updates(
		const int neighborhood_id = default_neighborhood_id
	) {
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " wait_remote_neighbor_copy_updates(...) called while balancing load"
				<< std::endl;
			abort();
		}

		bool ret_val = true;

		if (!this->wait_remote_neighbor_copy_update_receives(neighborhood_id)) {
			ret_val = false;
		}
		if (!this->wait_remote_neighbor_copy_update_sends()) {
			ret_val = false;
		}

		return ret_val;
	}


	/*!
	Waits for sends started by start_remote_neighbor_copy_updates().

	\see
	start_remote_neighbor_copy_updates()
	*/
	bool wait_remote_neighbor_copy_update_sends()
	{
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " wait_remote_neighbor_copy_update_sends() called while balancing load"
				<< std::endl;
			abort();
		}

		return this->wait_user_data_transfer_sends();
	}


	/*!
	Waits for receives started by start_remote_neighbor_copy_updates().

	\see
	start_remote_neighbor_copy_updates()
	*/
	bool wait_remote_neighbor_copy_update_receives(
		const int neighborhood_id = default_neighborhood_id
	) {
		if (this->balancing_load) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " wait_remote_neighbor_copy_update_receives(...) called while balancing load"
				<< std::endl;
			abort();
		}

		bool ret_val = true;

		if (neighborhood_id == default_neighborhood_id) {
			return this->wait_user_data_transfer_receives(
				#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
				this->remote_neighbors,
				this->cells_to_receive
				#endif
			);
		}

		if (this->user_hood_of.count(neighborhood_id) == 0) {
			ret_val = false;
		}

		if (!this->wait_user_data_transfer_receives(
			#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
			this->remote_neighbors,
			this->user_neigh_cells_to_receive.at(neighborhood_id)
			#endif
		)) {
			ret_val = false;
		}

		return ret_val;
	}


	/*!
	Returns number of cells that will be sent in a remote neighbor data update.

	The total amount of cells to be sent is returned so if one cell's data
	is sent to N processes it is counted N times.

	Returns maximum uint64_t if given neighborhood id doesn't exist.

	\see
	update_copies_of_remote_neighbors()
	add_neighborhood()
	*/
	uint64_t get_number_of_update_send_cells(
		const int neighborhood_id = default_neighborhood_id
	) const {
		uint64_t ret_val = 0;

		if (neighborhood_id == default_neighborhood_id) {
			for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
				receiver = cells_to_send.begin();
				receiver != cells_to_send.end();
				receiver++
			) {
				ret_val += receiver->second.size();
			}
			return ret_val;
		}

		if (this->user_hood_to.count(neighborhood_id) == 0) {
			return std::numeric_limits<uint64_t>::max();
		}

		#ifdef DEBUG
		if (this->user_neigh_cells_to_send.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighborhood with id " << neighborhood_id
				<< std::endl;
			abort();
		}
		#endif

		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			receiver_item = this->user_neigh_cells_to_send.at(neighborhood_id).begin();
			receiver_item != this->user_neigh_cells_to_send.at(neighborhood_id).end();
			receiver_item++
		) {
			ret_val += receiver_item->second.size();
		}

		return ret_val;
	}


	/*!
	Returns the number of cells whose data this process has to receive during a neighbor data update.


	Same as get_number_of_update_receive_cells() but for given neighborhood id.

	Returns maximum uint64_t if given neighborhood id doesn't exist.

	\see
	get_number_of_update_send_cells()
	*/
	uint64_t get_number_of_update_receive_cells(
		const int neighborhood_id = default_neighborhood_id
	) const {
		uint64_t ret_val = 0;

		if (neighborhood_id == default_neighborhood_id) {
			for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
				sender = this->cells_to_receive.begin();
				sender != this->cells_to_receive.end();
				sender++
			) {
				ret_val += sender->second.size();
			}
			return ret_val;
		}

		if (this->user_hood_of.count(neighborhood_id) == 0) {
			return std::numeric_limits<uint64_t>::max();
		}

		#ifdef DEBUG
		if (this->user_neigh_cells_to_receive.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighborhood with id " << neighborhood_id
				<< std::endl;
			abort();
		}
		#endif

		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			sender_item = this->user_neigh_cells_to_receive.at(neighborhood_id).begin();
			sender_item != this->user_neigh_cells_to_receive.at(neighborhood_id).end();
			sender_item++
		) {
			ret_val += sender_item->second.size();
		}

		return ret_val;
	}


	/*!
	Returns a pointer to the neighbors of given cell.

	In case the grid is not periodic in one or more directions,
	neighbors that would be outside of the grid are error_cell.
	Some neighbors might be on another process, but have a copy of their data on this process.
	The local copy of remote neighbors' data is updated, for example, by calling
	update_copies_of_remote_neighbors().

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

	If given a non-default neighborhood neighbors are in the same order as
	the offsets in the given neighborhood.

	Returns NULL if:
		- neighborhood with given id doesn't exist
		- given cell doesn't exist
		- given cell is on another process

	\see
	get_neighbors_to()
	get_neighbors_of_at_offset()
	add_neighborhood()
	*/
	const std::vector<uint64_t>* get_neighbors_of(
		const uint64_t cell,
		const int neighborhood_id = default_neighborhood_id
	) const {
		if (this->cells.count(cell) > 0) {
			if (neighborhood_id == default_neighborhood_id) {
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

			} else if (this->user_hood_of.count(neighborhood_id) > 0) {

				#ifdef DEBUG
				if (this->user_neigh_of.at(neighborhood_id).count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Neighbor list for cell " << cell
						<< " doesn't exist for neighborhood id " << neighborhood_id
						<< std::endl;
					abort();
				}
				#endif

				return &(this->user_neigh_of.at(neighborhood_id).at(cell));
			}
		}

		return NULL;
	}


	/*!
	Returns a pointer to the cells that consider given cell as a neighbor.

	This list doesn't include 0s even if the grid isn't periodic in some direction.
	Returns NULL if given cell doesn't exist or is on another process.
	Neighbors returned by this function are in no particular order.

	Returns NULL if neighborhood with given id doesn't exist.
	*/
	const std::vector<uint64_t>* get_neighbors_to(
		const uint64_t cell,
		const int neighborhood_id = default_neighborhood_id
	) const {
		if (this->cells.count(cell) > 0) {

			if (neighborhood_id == default_neighborhood_id) {
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

			} else if (this->user_hood_of.count(neighborhood_id) > 0) {

				#ifdef DEBUG
				if (this->user_neigh_to.at(neighborhood_id).count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighbors_to list for cell " << cell
						<< " doesn't exist for neighborhood id " << neighborhood_id
						<< std::endl;
					abort();
				}
				#endif
				return &(this->user_neigh_to.at(neighborhood_id).at(cell));
			}
		}

		return NULL;
	}


	/*!
	Returns cells which share a face with the given cell.

	Only those cells are returned which are considered as neighbors
	by given cell.
	Uses the default neighborhood.
	Does not return error_cell as a face neighbor.
	Returns nothing in the same cases as get_neighbors_of().

	uint64_t == neighbor id, int == neighbor direction.
	Directions are +N or -N where N is the Nth dimension and
	+ means positive direction and - negative in that dimension,
	e.g. +1 is positive x direction, -3 negative z.

	TODO:
	By default uses neighborhood with which this dccrg was initialized,
	\see
	default_neighborhood_id()
	get_neighbors_of()
	*/
	std::vector<std::pair<uint64_t, int> > get_face_neighbors_of(
		const uint64_t cell/*,
		const int neighborhood_id = default_neighborhood_id*/
	) const {
		std::vector<std::pair<uint64_t, int> > ret_val;

		if (this->cells.count(cell) == 0) {
			return ret_val;
		}

		// get location of face neighbors' offsets in neighborhood_of
		boost::array<size_t, 2 * 3> neighborhood_of_indices = {{0, 0, 0, 0, 0, 0}};
		for (int direction = -1; direction <= 1; direction += 2)
		for (size_t dimension = 0; dimension < 3; dimension++) {

			// neigh_of_indices[n] == negative direction in dimension n,
			// n + 1 positive direction
			const size_t neigh_of_indices_index
				= 2 * dimension + ((direction > 0) ? 1 : 0);

			for (size_t i = 0; i <= this->neighborhood_of.size(); i++) {
				if (i == this->neighborhood_of.size()) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Neighborhood_of offsets not found for face neighbors in dimension: "
						<< dimension << ", direction: " << direction
						<< std::endl;
					abort();
				}

				bool found = true;
				for (size_t other_dims = 0; other_dims < 3; other_dims++) {
					if (other_dims != dimension
					&& this->neighborhood_of[i][other_dims] != 0) {
						found = false;
						break;
					}
				}

				if (this->neighborhood_of[i][dimension] != direction) {
					found = false;
				}

				if (found) {
					neighborhood_of_indices[neigh_of_indices_index] = i;
					break;
				}
			}
		}

		// gather cells in given cell's neighbor_of list at indices found above
		boost::array<size_t, 2 * 3> current_index = {{0, 0, 0, 0, 0, 0}};
		const int refinement_level = this->get_refinement_level(cell);

		for (size_t
			neighbor_i = 0;
			neighbor_i < this->neighbors.at(cell).size();
			neighbor_i++
		) {

			const uint64_t neighbor = this->neighbors.at(cell)[neighbor_i];

			if (neighbor == error_cell) {
				for (size_t i = 0; i < current_index.size(); i++) {
					current_index[i]++;
				}
				continue;
			}

			const int neigh_ref_lvl = this->get_refinement_level(neighbor);

			for (int direction = -1; direction <= 1; direction += 2)
			for (size_t dimension = 0; dimension < 3; dimension++) {

				// neigh_of_indices[n] == negative direction, n + 1 positive
				const size_t neigh_of_indices_index
					= 2 * dimension + ((direction > 0) ? 1 : 0);

				// at correct index in neighbors_of for current dim & dir
				if (current_index[neigh_of_indices_index]
				== neighborhood_of_indices[neigh_of_indices_index]) {

					int final_dir = int(dimension) + 1;
					if (direction < 0) {
						final_dir *= -1;
					}

					// add one neighbor not smaller than given cell
					if (neigh_ref_lvl <= refinement_level) {

						ret_val.push_back(std::make_pair(neighbor, final_dir));

					// add only face neighbors in current dim & dir
					} else {

						const uint64_t neighs_in_offset = uint64_t(1) << 3;

						#ifdef DEBUG
						if (this->neighbors.at(cell).size() < neighbor_i + neighs_in_offset) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " Invalid number of neighbors for cell " << cell
								<< " while processing dimension " << dimension
								<< " and direction " << direction
								<< " starting at index " << neighbor_i
								<< std::endl;
							abort();
						}
						#endif

						// see find_neighbors_of(...) for the order of these
						const std::vector<uint64_t> dir_neighs(
							this->neighbors.at(cell).begin() + neighbor_i,
							this->neighbors.at(cell).begin() + neighbor_i + neighs_in_offset
						);

						// neighbor at offset    0, 1, 2, 3, 4, 5, 6, 7 is
						// face neighbor of given cell when neighbors are in
						// dim = 0, dir = -1      , y,  , y,  , y,  , y
						// dim = 0, dir = +1     y,  , y,  , y,  , y,
						// dim = 1, dir = -1      ,  , y, y,  ,  , y, y
						// dim = 1, dir = +1     y, y,  ,  , y, y,  ,
						// dim = 2, dir = -1      ,  ,  ,  , y, y, y, y
						// dim = 2, dir = +1     y, y, y, y,  ,  ,  ,
						const size_t
							batch_size = size_t(1) << dimension,
							mod_target = (direction < 0) ? 1 : 0;

						for (size_t i = 0; i < neighs_in_offset; i++) {

							#ifdef DEBUG
							if (dir_neighs[i] == error_cell) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " Invalid neighbor of cell " << cell
								<< " at index " << neighbor_i + i
								<< std::endl;
							abort();
							}
							#endif

							if ((i / batch_size) % 2 == mod_target) {
								ret_val.push_back(std::make_pair(dir_neighs[i], final_dir));
							}
						}
					}
				}

				current_index[neigh_of_indices_index]++;
			}

			// skip all cells in this neighborhood offset
			if (neigh_ref_lvl > refinement_level) {
				neighbor_i += 7;
			}
		}

		return ret_val;
	}


	/*!
	Returns the size of cells' neihgbourhood in every direction.
	*/
	unsigned int get_neighborhood_length() const
	{
		return this->neighborhood_length;
	}


	/*!
	Returns all neighbors of given cell that are at given offset from it.

	Offset is in units of size of the given cell
	Returns nothing in the following cases:
		- given cell doesn't exist
		- given cell is on another process
		- any of given offsets is larger in absolute value than the neighborhood
		  size or larger than 1 if neihgborhood size == 0
		- i == 0 && j == 0 && k == 0

	\see
	get_neighbors_of()
	*/
	std::vector<uint64_t> get_neighbors_of_at_offset(
		const uint64_t cell,
		const int i,
		const int j,
		const int k
	) const {
		std::vector<uint64_t> return_neighbors;
		if (this->cell_process.count(cell) == 0
		|| this->cell_process.at(cell) != this->rank
		|| (i == 0 && j == 0 && k == 0)) {
			return return_neighbors;
		}

		const int refinement_level = this->get_refinement_level(cell);

		// find cell(s) at given indices in the stored neighbor list
		const int last_offset
			= (this->neighborhood_length > 0)
			? int(this->neighborhood_length)
			: 1;

		int index = 0;
		for (int
			current_k = (this->neighborhood_length > 0) ? -int(this->neighborhood_length) : -1;
			current_k <= last_offset;
			current_k++
		)
		for (int
			current_j = (this->neighborhood_length > 0) ? -int(this->neighborhood_length) : -1;
			current_j <= last_offset;
			current_j++
		)
		for (int
			current_i = (this->neighborhood_length > 0) ? -int(this->neighborhood_length) : -1;
			current_i <= last_offset;
			current_i++
		) {
			if (current_i == 0 && current_j == 0 && current_k == 0) {
				continue;
			}

			if (this->neighborhood_length == 0) {
				// skip diagonal offsets
				const int zero_offsets_in_current =
					  ((current_i == 0) ? 1 : 0)
					+ ((current_j == 0) ? 1 : 0)
					+ ((current_k == 0) ? 1 : 0);
				if (zero_offsets_in_current != 2) {
					continue;
				}
			}

			const int current_refinement_level
				= this->get_refinement_level(this->neighbors.at(cell)[index]);

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
	Returns neighbors of given local cell that are on another process.

	Returns nothing if:
		- given cell doesn't exist
		- given cell is on another process
		- given cell doesn't have remote neighbors
		- given neighborhood doesn't exist

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.
	*/
	std::vector<uint64_t> get_remote_neighbors_of(
		const uint64_t cell,
		const int neighborhood_id = default_neighborhood_id,
		const bool sorted = false
	) const {
		std::vector<uint64_t> ret_val;

		if (this->cell_process.count(cell) == 0) {
			return ret_val;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return ret_val;
		}

		if (neighborhood_id != default_neighborhood_id
		&& this->user_hood_of.count(neighborhood_id) == 0) {
			return ret_val;
		}

		const std::vector<uint64_t>& neighbors_ref
			= (neighborhood_id == default_neighborhood_id)
			? this->neighbors.at(cell)
			: this->user_neigh_of.at(neighborhood_id).at(cell);

		BOOST_FOREACH(const uint64_t neighbor, neighbors_ref) {

			if (neighbor == error_cell) {
				continue;
			}

			if (this->cell_process.at(neighbor) != this->rank) {
				ret_val.push_back(neighbor);
			}
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
	}

	/*!
	Returns remote cells that consider given local cell as a neighbor.

	Returns nothing if given cell doesn't exist or is on another process
	or doesn't have remote neighbors.

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.
	*/
	std::vector<uint64_t> get_remote_neighbors_to(
		const uint64_t cell,
		const int neighborhood_id = default_neighborhood_id,
		const bool sorted = false
	) const {
		std::vector<uint64_t> ret_val;

		if (this->cell_process.count(cell) == 0) {
			return ret_val;
		}

		if (this->cell_process.at(cell) != this->rank) {
			return ret_val;
		}

		if (neighborhood_id != default_neighborhood_id
		&& this->user_hood_of.count(neighborhood_id) == 0) {
			return ret_val;
		}

		const std::vector<uint64_t>& neighbors_ref
			= (neighborhood_id == default_neighborhood_id)
			? this->neighbors_to.at(cell)
			: this->user_neigh_to.at(neighborhood_id).at(cell);

		BOOST_FOREACH(const uint64_t neighbor, neighbors_ref) {

			if (neighbor == error_cell) {
				continue;
			}

			if (this->cell_process.at(neighbor) != this->rank) {
				ret_val.push_back(neighbor);
			}
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
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
				const std::vector<uint64_t> neighbor_siblings
					= this->get_all_children(this->get_parent(neighbor));

				BOOST_FOREACH(const uint64_t& sibling, neighbor_siblings) {
					this->cells_to_unrefine.erase(sibling);
				}
			}
		}

		BOOST_FOREACH(const uint64_t& neighbor, this->neighbors_to.at(cell)) {

			if (this->get_refinement_level(neighbor) <= refinement_level) {
				const std::vector<uint64_t> neighbor_siblings
					= this->get_all_children(this->get_parent(neighbor));

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
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid refinement level for parent"
				<< std::endl;
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

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.
	*/
	std::vector<uint64_t> stop_refining(const bool sorted = false)
	{
		this->induce_refines();

		// update dont_refines between processes
		this->all_to_all_set(this->cells_not_to_unrefine);

		this->override_unrefines();
		this->cells_not_to_unrefine.clear();

		std::vector<uint64_t> ret_val = this->execute_refines();

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
	}


	/*!
	Returns cells that were removed by unrefinement whose parent is on this process
	Removed cells data is also on this process, but only until balance_load() is called

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.
	*/
	std::vector<uint64_t> get_removed_cells(const bool sorted = false) const
	{
		std::vector<uint64_t> ret_val;
		ret_val.reserve(this->unrefined_cell_data.size());

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->unrefined_cell_data) {
			ret_val.push_back(item.first);
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
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

	Neighborhood is returned in units of length_in_indices.
	If grid is not periodic in some direction then indices which would fall outside
	of the grid in that direction are returned as error_index.
	*/
	std::vector<Types<3>::indices_t> indices_from_neighborhood(
		const Types<3>::indices_t indices,
		const uint64_t length_in_indices,
		const std::vector<Types<3>::neighborhood_item_t>& neighborhood
	) const {
		// TODO: make neighborhood a const reference
		std::vector<Types<3>::indices_t> return_indices;
		return_indices.reserve(neighborhood.size());

		// grid length in indices
		const uint64_t grid_length[3] = {
			this->get_length_x() * (uint64_t(1) << this->max_refinement_level),
			this->get_length_y() * (uint64_t(1) << this->max_refinement_level),
			this->get_length_z() * (uint64_t(1) << this->max_refinement_level)
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
							if (temp_indices[dimension] < length_in_indices - 1
							&& temp_indices[dimension] > 0) {
								std::cerr << __FILE__ << ":" << __LINE__
									<< " Cells aren't supposed to wrap around the grid."
									<< std::endl;
								abort();
							}
							#endif

							if (temp_indices[dimension] >= length_in_indices) {
								temp_indices[dimension] -= length_in_indices;
							} else {
								temp_indices[dimension] = grid_length[dimension] - length_in_indices;
							}
						}
					// use error_indices to signal that this neighborhood item is outside of the grid
					} else {
						if (indices[dimension] < abs(offsets[dimension]) * length_in_indices) {
							temp_indices[0] = error_index;
							temp_indices[1] = error_index;
							temp_indices[2] = error_index;
							break;
						}

						temp_indices[dimension] += offsets[dimension] * length_in_indices;
					}

				} else {

					if (this->is_periodic(dimension)) {
						for (int i = 0; i < offsets[dimension]; i++) {

							#ifdef DEBUG
							if (temp_indices[dimension] > grid_length[dimension] - length_in_indices) {
								std::cerr << __FILE__ << ":" << __LINE__
									<< " Cells aren't supposed to wrap around the grid."
									<< std::endl;
								abort();
							}
							#endif

							if (temp_indices[dimension] < grid_length[dimension] - length_in_indices) {
								temp_indices[dimension] += length_in_indices;
							} else {
								temp_indices[dimension] = 0;
							}
						}
					} else {
						if (
							indices[dimension] + offsets[dimension] * length_in_indices
							>= grid_length[dimension]
						) {
							temp_indices[0] = error_index;
							temp_indices[1] = error_index;
							temp_indices[2] = error_index;
							break;
						}

						temp_indices[dimension] += offsets[dimension] * length_in_indices;
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

	Cells smaller than given one at any offset are returned in the order in which
	neighbors' index increases first in x, then in y, ...
	In other words the first small cell within an neighborhood offset is
	closest to the origin of the grid, the second one is in direction +1
	from the first, 3rd is in dir +2 from 1st, 4th is in +1 from 3rd, +2 from 2nd,
	8th is +3 from 4th, +2 from 6th, +1 from 7th, assuming max_ref_lvl_diff == 1
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< " max_diff must not be negative"
				<< std::endl;
			abort();
		}

		if (cell == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid cell given: " << cell
				<< std::endl;
			abort();
		}

		if (refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Refinement level of given cell (" << cell
				<< ") is too large: " << refinement_level
				<< std::endl;
			abort();
		}

		if (refinement_level < 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid refinement level for cell " << cell
				<< ": " << refinement_level
				<< std::endl;
			abort();
		}
		#endif

		if (this->cell_process.count(cell) == 0) {
			return return_neighbors;
		}

		if (!has_children && cell != this->get_child(cell)) {
			return return_neighbors;
		}

		const uint64_t cell_length = this->get_cell_length_in_indices(cell);

		const std::vector<Types<3>::indices_t> indices_of = this->indices_from_neighborhood(
			this->get_indices(cell),
			cell_length,
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
				const uint64_t smallest
					= this->get_existing_cell(index_of, 0, this->max_refinement_level);

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
					index_of[0] + cell_length - 1,
					index_of[1] + cell_length - 1,
					index_of[2] + cell_length - 1
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
						<< " of size " << cell_length
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
							<< " doesn't exist between refinement levels "
							<< refinement_level - max_diff
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Refinement level (" << refinement_level
				<< ") of cell " << cell
				<< " exceeds maximum refinement level of the grid ("
				<< this->max_refinement_level << ")"
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
			const uint64_t length_in_indices = this->get_cell_length_in_indices(parent);

			std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
				indices,
				length_in_indices,
				neighborhood
			);

			BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

				if (search_index[0] == error_index) {
					continue;
				}

				const uint64_t found
					= this->get_cell_from_indices(search_index, refinement_level - 1);

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
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Got no children for cell " << cell
					<< std::endl;
				abort();
			}
			#endif

			const uint64_t length_in_indices = this->get_cell_length_in_indices(children[0]);

			BOOST_FOREACH(const uint64_t& child, children) {

				const Types<3>::indices_t indices = this->get_indices(child);

				std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
					indices,
					length_in_indices,
					neighborhood
				);

				BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

					if (search_index[0] == error_index) {
						continue;
					}

					const uint64_t found
						= this->get_cell_from_indices(search_index, refinement_level + 1);

					if (found == this->get_child(found)) {
						unique_neighbors.insert(found);
					}
				}
			}
		}

		// neighbors_to of the same size as given cell
		const Types<3>::indices_t indices = this->get_indices(cell);
		const uint64_t length_in_indices = this->get_cell_length_in_indices(cell);

		std::vector<Types<3>::indices_t> search_indices = this->indices_from_neighborhood(
			indices,
			length_in_indices,
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
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid parent for cell " << cell
					<< std::endl;
				abort();
			}
			#endif

			const Types<3>::indices_t indices = this->get_indices(parent);
			const uint64_t length_in_indices = this->get_cell_length_in_indices(parent);

			const std::vector<Types<3>::indices_t> search_indices
				= this->indices_from_neighborhood(
					indices,
					length_in_indices,
					this->neighborhood_to
				);

			BOOST_FOREACH(const Types<3>::indices_t& search_index, search_indices) {

				if (search_index[0] == error_index) {
					continue;
				}

				const uint64_t found
					= this->get_cell_from_indices(search_index, refinement_level - 1);

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
	std::vector<uint64_t> find_cells(
		const Types<3>::indices_t indices_min,
		const Types<3>::indices_t indices_max,
		const int minimum_refinement_level,
		const int maximum_refinement_level
	) const {
		// size of cells in indices of given maximum_refinement_level
		const uint64_t index_increase
			= uint64_t(1) << (this->max_refinement_level - maximum_refinement_level);

		#ifdef DEBUG
		if (minimum_refinement_level > maximum_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid refinement levels given"
				<< std::endl;
			abort();
		}

		if (maximum_refinement_level > this->max_refinement_level) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Invalid maximum refinement level given"
				<< std::endl;
			abort();
		}

		// check that outer shell makes sense
		if (indices_min[0] > indices_max[0]) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " minimum x index > maximum x index"
				<< std::endl;
			abort();
		}

		if (indices_min[1] > indices_max[1]) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " minimum y index > maximum y index"
				<< std::endl;
			abort();
		}

		if (indices_min[2] > indices_max[2]) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " minimum z index > maximum z index"
				<< std::endl;
			abort();
		}
		#endif

		std::vector<uint64_t> result;
		boost::unordered_set<uint64_t> uniques;

		Types<3>::indices_t indices = {{0, 0, 0}};
		for (indices[2] = indices_min[2]; indices[2] <= indices_max[2]; indices[2] += index_increase)
		for (indices[1] = indices_min[1]; indices[1] <= indices_max[1]; indices[1] += index_increase)
		for (indices[0] = indices_min[0]; indices[0] <= indices_max[0]; indices[0] += index_increase) {

			const uint64_t cell
				= this->get_existing_cell(
					indices,
					minimum_refinement_level,
					maximum_refinement_level
				);

			#ifdef DEBUG
			if (cell == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No cell found between refinement levels [" << minimum_refinement_level
					<< ", " << maximum_refinement_level
					<< "] at indices " << indices[0]
					<< " " << indices[1]
					<< " " << indices[2]
					<< std::endl;

				const uint64_t smallest
					= this->get_existing_cell(indices, 0, this->max_refinement_level);

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
	Removes Cell_Data of refined and unrefined cells from this process.

	\see
	refine_completely()
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< "User tried to set an option reserved for dccrg (" << name
				<< ": " << value << ")"
				<< std::endl;
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< "User tried to assign " << processes
				<< " processes per part for a new hierarchial partitioning level"
				<< std::endl;
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

		this->processes_per_part.erase(
			this->processes_per_part.begin() + hierarchial_partitioning_level
		);
		this->partitioning_options.erase(
			this->partitioning_options.begin() + hierarchial_partitioning_level
		);
	}


	/*!
	Adds (or overwrites) the given option and its value for
	hierarchial partitioning of given level.

	Level numbering starts from 0.
	Does nothing in the following cases:
		- option name is one of: RETURN_LISTS, ...
		- given level doesn't exist
	*/
	void add_partitioning_option(
		const int hierarchial_partitioning_level,
		const std::string name,
		const std::string value
	) {
		if (hierarchial_partitioning_level < 0
		|| hierarchial_partitioning_level >= int(this->processes_per_part.size())) {
			return;
		}

		if (this->reserved_options.count(name) > 0) {
			#ifdef DEBUG
			std::cerr << __FILE__ << ":" << __LINE__
				<< "User tried to set an option reserved for dccrg (" << name
				<< ": " << value
				<< ") for level " << hierarchial_partitioning_level
				<< " of hierarchial partitioning"
				<< std::endl;
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
	void remove_partitioning_option(
		const int hierarchial_partitioning_level,
		const std::string name
	) {
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
	std::string get_partitioning_option_value(
		const int hierarchial_partitioning_level,
		const std::string name
	) const {
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

	TODO: only execute for cells that are pinned
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
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << cell
						<< " should be in this->cells of process " << this->rank
						<< std::endl;
					abort();
				}
			} else {
				if (this->cells.count(cell) > 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cell " << cell
						<< " shouldn't be in this->cells of process " << this->rank
						<< std::endl;
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
	Returns cells that are on the local process boundary.

	Returns cells for which either of the following is true:
		- a cell on another process is considered as a neighbor
		- are considered as a neighbor by a cell on another process

	Given neighborhood is used when considering neighbor relations.

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.

	Returns nothing if neighborhood with given id doesn't exist.

	\see get_cells() is more general than this but slower.
	*/
	std::vector<uint64_t> get_local_cells_on_process_boundary(
		const int neighborhood_id = default_neighborhood_id,
		const bool sorted = false
	) const {
		std::vector<uint64_t> ret_val;

		if (neighborhood_id == default_neighborhood_id) {
			ret_val.insert(
				ret_val.end(),
				this->local_cells_on_process_boundary.begin(),
				this->local_cells_on_process_boundary.end()
			);
		} else if (this->user_hood_of.count(neighborhood_id) > 0) {
			ret_val.insert(
				ret_val.end(),
				this->user_local_cells_on_process_boundary.at(neighborhood_id).begin(),
				this->user_local_cells_on_process_boundary.at(neighborhood_id).end()
			);
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
	}


	/*!
	Returns cells that are not on the local process boundary.

	Returns cells which do not consider a cell on another process
	as a neighbor and which are not considered as a neighbor
	by a cell on another process.

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.

	Returns nothing if neighborhood with given id doesn't exist.

	\see get_cells() is more general than this but a bit slower.
	*/
	std::vector<uint64_t> get_local_cells_not_on_process_boundary(
		const int neighborhood_id = default_neighborhood_id,
		const bool sorted = false
	) const {
		std::vector<uint64_t> ret_val;

		if (neighborhood_id == default_neighborhood_id) {

			BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
				const uint64_t cell = item.first;
				if (this->local_cells_on_process_boundary.count(cell) == 0) {
					ret_val.push_back(cell);
				}
			}

		} else if (this->user_hood_of.count(neighborhood_id) > 0) {

			BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
				const uint64_t cell = item.first;
				if (this->user_local_cells_on_process_boundary.at(neighborhood_id).count(cell) == 0) {
					ret_val.push_back(cell);
				}
			}
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
	}


	/*!
	Returns remote cells that are on the local process boundary.

	Returns cells for which either of the following is true:
		- a local cell is considered as a neighbor
		- are considered as a neighbor by a local cell

	By default returned cells are in random order but if sorted == true
	they are sorted using std::sort before returning.

	Returns nothing if neighborhood with given id doesn't exist.

	\see get_cells() is more general than this but a bit slower.
	*/
	std::vector<uint64_t> get_remote_cells_on_process_boundary(
		const int neighborhood_id = default_neighborhood_id,
		const bool sorted = false
	) const {
		std::vector<uint64_t> ret_val;

		if (neighborhood_id == default_neighborhood_id) {
			ret_val.insert(
				ret_val.end(),
				this->remote_cells_on_process_boundary.begin(),
				this->remote_cells_on_process_boundary.end()
			);
		} else if (this->user_hood_of.count(neighborhood_id) > 0) {
			ret_val.insert(ret_val.end(),
				this->user_remote_cells_on_process_boundary.at(neighborhood_id).begin(),
				this->user_remote_cells_on_process_boundary.at(neighborhood_id).end()
			);
		}

		if (sorted && ret_val.size() > 0) {
			std::sort(ret_val.begin(), ret_val.end());
		}

		return ret_val;
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
	Returns the cells that will be added to this process by load balancing.
	*/
	const boost::unordered_set<uint64_t>& get_cells_added_by_balance_load() const
	{
		return this->added_cells;
	}

	/*!
	Returns the cells that will be removed from this process by load balancing.
	*/
	const boost::unordered_set<uint64_t>& get_cells_removed_by_balance_load() const
	{
		return this->removed_cells;
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
		const uint64_t index_offset
			= (uint64_t(1) << (this->max_refinement_level - refinement_level));

		for (uint64_t
			z_index_offset = 0;
			z_index_offset < 2 * index_offset;
			z_index_offset += index_offset
		)
		for (uint64_t
			y_index_offset = 0;
			y_index_offset < 2 * index_offset;
			y_index_offset += index_offset
		)
		for (uint64_t
			x_index_offset = 0;
			x_index_offset < 2 * index_offset;
			x_index_offset += index_offset
		) {
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


	/*!
	Returns cell weights that have been set.
	*/
	const boost::unordered_map<uint64_t, double>& get_cell_weights() const
	{
		return this->cell_weights;
	}


	/*!
	Adds a new neighborhood for updating Cell_Data between neighbors on different processes.

	Must be called with identical parameters on all processes.
	No neighborhood_item_t should have all offsets equal to 0.
	Returns true on success and false in any of the following cases:
		- given id already exists, use remove_neighborhood() before calling this
		- (part of) given neighborhood is outside of initial neighborhood size

	Can be used to reduce the amount of data transferred between
	processes if Cell_Data from the full neighborhood isn't needed.

	\see
	initialize()
	remove_neighborhood()
	*/
	bool add_neighborhood(
		const int neighborhood_id,
		const std::vector<Types<3>::neighborhood_item_t>& given_neigh
	) {
		if (neighborhood_id == default_neighborhood_id) {
			return false;
		}

		if (this->user_hood_of.count(neighborhood_id) > 0) {

			#ifdef DEBUG
			if (this->user_hood_to.count(neighborhood_id) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should have id " << neighborhood_id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_of.count(neighborhood_id) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should have id " << neighborhood_id
					<< std::endl;
				abort();
			}

			if (this->user_neigh_to.count(neighborhood_id) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Should have id " << neighborhood_id
					<< std::endl;
				abort();
			}
			#endif

			return false;
		}

		#ifdef DEBUG
		if (this->user_hood_to.count(neighborhood_id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_of.count(neighborhood_id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_to.count(neighborhood_id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_send.count(neighborhood_id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_cells_to_receive.count(neighborhood_id) > 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Should not have id " << neighborhood_id
				<< std::endl;
			abort();
		}
		#endif

		if (this->neighborhood_length > 0) {

			BOOST_FOREACH(const Types<3>::neighborhood_item_t& neigh_item, given_neigh) {
				for (size_t i = 0; i < 3; i++) {
					if ((unsigned int) abs(neigh_item[i]) > this->neighborhood_length) {
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
		this->user_hood_of[neighborhood_id] = given_neigh;
		this->user_hood_to[neighborhood_id].clear();
		BOOST_FOREACH(const Types<3>::neighborhood_item_t& neigh_item, given_neigh) {
			const Types<3>::neighborhood_item_t neigh_item_to = {{
				-neigh_item[0],
				-neigh_item[1],
				-neigh_item[2]
			}};
			this->user_hood_to.at(neighborhood_id).push_back(neigh_item_to);
		}

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {
			this->update_user_neighbors(item.first, neighborhood_id);
		}

		this->update_user_remote_neighbor_info(neighborhood_id);

		this->recalculate_neighbor_update_send_receive_lists(neighborhood_id);

		return true;
	}


	/*!
	Removes the given neighborhood.

	Must be called with identical id on all processes.
	Frees local neighbor lists and other resources associated
	with given neighborhood.

	\see
	add_neighborhood()
	*/
	void remove_neighborhood(const int neighborhood_id)
	{
		if (neighborhood_id == default_neighborhood_id) {
			return;
		}

		this->user_hood_of.erase(neighborhood_id);
		this->user_hood_to.erase(neighborhood_id);
		this->user_neigh_of.erase(neighborhood_id);
		this->user_neigh_to.erase(neighborhood_id);
		this->user_neigh_cells_to_send.erase(neighborhood_id);
		this->user_neigh_cells_to_receive.erase(neighborhood_id);
		this->user_local_cells_on_process_boundary.erase(neighborhood_id);
		this->user_remote_cells_on_process_boundary.erase(neighborhood_id);
	}


	/*!
	Returns whether this dccrg instance has been initialized.
	*/
	bool get_initialized() const
	{
		return this->initialized;
	}

	/*!
	Returns the maximum value an MPI tag can have.
	*/
	unsigned int get_max_tag() const
	{
		return this->max_tag;
	}

	/*!
	Returns the maximum allowed difference in refinement level between neighboring cells.

	This difference in enforced both between a cell and its neighbors_of and its neighbors_to.
	*/
	int get_max_ref_lvl_diff() const
	{
		return this->max_ref_lvl_diff;
	}

	/*!
	Returns whether cell data is transferred using one message or more.

	If true then one MPI message per cell is used, otherwise all cells
	are transferred in one message.
	*/
	bool get_send_single_cells() const
	{
		return this->send_single_cells;
	}

	/*!
	Sets whether cell data is transferred using one message or more.

	Do not switch sending type while data transfers are going on,
	for example with start_remote_neighbor_data_update(...).
	See get_send_single_cells() for more info.
	*/
	void set_send_single_cells(const bool given)
	{
		this->send_single_cells = given;
	}

	/*!
	Returns dccrg's communicator.

	Returns a duplicate of the MPI communicator of dccrg
	which itself is a duplicate of the communicator that
	was given to dccrg's initialize function.
	*/
	MPI_Comm get_communicator() const
	{
		MPI_Comm ret_val;
		int result = MPI_Comm_dup(this->comm, &ret_val);
		if (result != MPI_SUCCESS) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Couldn't duplicate communicator: " << Error_String()(result)
				<< std::endl;
		}
		return ret_val;
	}

	/*!
	Returns the MPI process number, i.e. rank, of this process.
	*/
	uint64_t get_rank() const
	{
		return this->rank;
	}

	/*!
	Returns the number of MPI processes participating in this dccrg instance.
	*/
	uint64_t get_comm_size() const
	{
		return this->comm_size;
	}

	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
	/*!
	Returns the boost communicator of this dccrg instance.
	*/
	boost::mpi::communicator get_boost_comm() const
	{
		return this->boost_comm;
	}
	#endif

	/*!
	Returns the storage of mostly local cell ids and their data.
	*/
	const boost::unordered_map<uint64_t, Cell_Data>& get_cell_data() const
	{
		return this->cells;
	}

	/*!
	Returns neighbor lists of local cells.
	At some point also neighbor lists of some other processes' cells might be included.
	*/
	const boost::unordered_map<uint64_t, std::vector<uint64_t> >& get_cell_neighbor() const
	{
		return this->neighbors;
	}

	/*!
	Returns the neighborhood of cells.

	Returns offsets in which cells are considered as neighbors of a cell.
	*/
	const std::vector<Types<3>::neighborhood_item_t>& get_neighborhood_of() const
	{
		return this->neighborhood_of;
	}

	/*!
	Returns the opposite of neighborhood of cells.

	Returns offsets in which cells consider a cell as neighbor.
	*/
	const std::vector<Types<3>::neighborhood_item_t>& get_neighborhood_to() const
	{
		return this->neighborhood_to;
	}

	/*!
	Returns the offsets of user defined neighborhoods.
	*/
	const boost::unordered_map<
		int,
		std::vector<Types<3>::neighborhood_item_t>
	>& get_user_hood_of() const
	{
		return this->user_hood_of;
	}

	/*!
	Returns the opposite of offsets of user defined neighborhoods
	*/
	const boost::unordered_map<
		int,
		std::vector<Types<3>::neighborhood_item_t>
	>& get_user_hood_to() const
	{
		return this->user_hood_to;
	}

	/*!
	Returns cells (2nd value) which consider a cell (1st value) as neighbor
	*/
	const boost::unordered_map<uint64_t, std::vector<uint64_t> >& get_all_neighbors_to() const
	{
		return this->neighbors_to;
	}

	/*!
	User defined neighborhood version of get_cell_neighbor.
	*/
	const boost::unordered_map<
		int,
		boost::unordered_map<uint64_t, std::vector<uint64_t> >
	>& get_all_user_neigh_of() const
	{
		return this->user_neigh_of;
	}

	/*!
	User defined neighborhood version of get_all_neighbors_to().
	*/
	const boost::unordered_map<
		int,
		boost::unordered_map<uint64_t, std::vector<uint64_t> >
	>& get_all_user_neigh_to() const
	{
		return this->user_neigh_to;
	}

	/*!
	Returns the process (2nd value) of each cell (1st value).
	*/
	const boost::unordered_map<uint64_t, uint64_t>& get_cell_process() const
	{
		return this->cell_process;
	}

	/*!
	Returns cells which have a remote neighbor.
	*/
	const boost::unordered_set<uint64_t>& get_local_cells_on_process_boundary_internal() const
	{
		return this->local_cells_on_process_boundary;
	}

	/*!
	Returns remote cells which have a local neighbor.
	*/
	const boost::unordered_set<uint64_t>& get_remote_cells_on_process_boundary_internal() const
	{
		return this->remote_cells_on_process_boundary;
	}

	/*!
	Returns cells which have a remote neighbor for each used defined neighborhood.
	*/
	const boost::unordered_map<
		int,
		boost::unordered_set<uint64_t>
	>& get_user_local_cells_on_process_boundary() const
	{
		return this->user_local_cells_on_process_boundary;
	}

	/*!
	Returns remote cells which have a local neighbor for each used defined neighborhood.
	*/
	const boost::unordered_map<
		int,
		boost::unordered_set<uint64_t>
	>& get_user_remote_cells_on_process_boundary() const
	{
		return this->user_remote_cells_on_process_boundary;
	}

	/*!
	Returns cells that will be sent to other processes.

	First int is the target process, uint64_t is cell id and
	last int is the message tag when dccrg sends one cell / MPI_Isend.
	*/
	const boost::unordered_map<
		int,
		std::vector<std::pair<uint64_t, int> >
	>& get_cells_to_send() const
	{
		return this->cells_to_send;
	}

	/*!
	Returns cells that will be received from other processes.

	First int is the sending process, uint64_t is cell id and
	last int is the message tag when dccrg receives one cell / MPI_Irecv.
	*/
	const boost::unordered_map<
		int,
		std::vector<std::pair<uint64_t, int> >
	>& get_cells_to_receive() const
	{
		return this->cells_to_receive;
	}

	/*!
	Returns cells which will be sent for each user defined neighborhood.
	*/
	const boost::unordered_map<
		int,
		boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >
	>& get_user_neigh_cells_to_send() const
	{
		return this->user_neigh_cells_to_send;
	}

	/*!
	Returns cells which will be received for each user defined neighborhood.
	*/
	const boost::unordered_map<
		int,
		boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >
	>& get_user_neigh_cells_to_receive() const
	{
		return this->user_neigh_cells_to_receive;
	}

	/*!
	Returns pin requests currently in force.
	*/
	const boost::unordered_map<uint64_t, uint64_t>& get_pin_requests() const
	{
		return this->pin_requests;
	}

	/*!
	Returns pin requests of this process.

	These have not been told to other processes yet.
	*/
	const boost::unordered_map<uint64_t, uint64_t>& get_new_pin_requests() const
	{
		return this->new_pin_requests;
	}

	/*!
	Returns the pointer to this instances Zoltan data structure.
	*/
	const Zoltan_Struct* get_zoltan() const
	{
		return this->zoltan;
	}

	/*!
	Returns the number of processes per Zoltan partition.

	First partition has vec[0] number of processes, second vec[1], etc.
	*/
	const std::vector<unsigned int>& get_processes_per_part() const
	{
		return this->processes_per_part;
	}

	/*!
	Returns the options that are given to Zoltan when partitioning.
	*/
	const std::vector<
		boost::unordered_map<std::string, std::string>
	> get_partitioning_options() const
	{
		return this->partitioning_options;
	}

	/*!
	Returns whether load balancing by Zoltan is supposed to fail.
	*/
	bool get_no_load_balancing() const
	{
		return this->no_load_balancing;
	}

	/*!
	Returns options to Zoltan which cannot be set by the user.
	*/
	const boost::unordered_set<std::string>& get_reserved_options() const
	{
		return this->reserved_options;
	}

	/*!
	Returns processes which have cells close enough to any cell of this process.

	Close enough is number of refinement level 0 cells * size of default neighborhood.
	FIXME not implemented yet so doesn't return anything at the moment.
	*/
	const boost::unordered_set<uint64_t>& get_neighbor_processes() const
	{
		return this->neighbor_processes;
	}

	/*!
	Returns whether this instance is currently balancing the load.
	*/
	bool get_balancing_load() const
	{
		return this->balancing_load;
	}



private:

	bool initialized;

	// size of the neighbor stencil of a cells in cells (of the same size as the cell itself)
	unsigned int neighborhood_length;

	// maximum value an MPI tag can have
	unsigned int max_tag;

	// maximum difference in refinement level between neighbors
	int max_ref_lvl_diff;

	// whether to send user's cell data between processes
	// with only one message or one message / cell
	bool send_single_cells;

	// the grid is distributed between these processes
	MPI_Comm comm;
	uint64_t rank, comm_size;
	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
	boost::mpi::communicator boost_comm;
	#endif

	// cells and their data on this process
	boost::unordered_map<uint64_t, Cell_Data> cells;

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
	Cell on this process and those cells that aren't neighbors of
	this cell but whose neighbor this cell is.
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

	// cells on this process that have a neighbor on another
	// process or are considered as a neighbor of a cell on another process
	boost::unordered_set<uint64_t> local_cells_on_process_boundary;

	// cells on other processes that have a neighbor on this process
	// or are considered as a neighbor of a cell on this process
	boost::unordered_set<uint64_t> remote_cells_on_process_boundary;

	/*
	User defined versions of local_cells_on_process_boundary and remote_cells_on_process_boundary
	*/
	boost::unordered_map<
		int, // user defined id of neighbor lists
		boost::unordered_set<uint64_t>
	> user_local_cells_on_process_boundary, user_remote_cells_on_process_boundary;

	// remote neighbors and their data, of cells on this process
	boost::unordered_map<uint64_t, Cell_Data> remote_neighbors;

	boost::unordered_map<
		int,
		std::vector<
			#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
			boost::mpi::request
			#else
			MPI_Request
			#endif
		>
	> send_requests, receive_requests;

	// cells whose data has to be received / sent by this process from the process as the key
	boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >
		cells_to_send, cells_to_receive;

	/*
	User defined neighborhood versions of cells_to_send and _receive.
	*/
	boost::unordered_map<
		int, // user defined id of neighbor lists
		boost::unordered_map<
			int, // process to send to / receive from
			std::vector<std::pair<uint64_t, int> >
		>
	> user_neigh_cells_to_send, user_neigh_cells_to_receive;

	// cells added to / removed from this process by load balancing
	boost::unordered_set<uint64_t> added_cells, removed_cells;

	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
	// storage for cells' user data that awaits transfer to or from this process
	boost::unordered_map<int, std::vector<Cell_Data> > incoming_data, outgoing_data;
	#endif

	// cells to be refined / unrefined after a call to stop_refining()
	boost::unordered_set<uint64_t> cells_to_refine, cells_to_unrefine;

	// cells whose siblings shouldn't be unrefined
	boost::unordered_set<uint64_t> cells_not_to_unrefine;

	// stores user data of cells whose children were created while refining
	boost::unordered_map<uint64_t, Cell_Data> refined_cell_data;
	// stores user data of cells that were removed while unrefining
	boost::unordered_map<uint64_t, Cell_Data> unrefined_cell_data;

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
	// record whether Zoltan_LB_Partition is expected to fail
	// (when the user selects NONE as the load balancing algorithm)
	bool no_load_balancing;
	// reserved options that the user cannot change
	boost::unordered_set<std::string> reserved_options;

	// optional user-given weights of cells on this process
	boost::unordered_map<uint64_t, double> cell_weights;

	// processes which have cells close enough from cells of this process
	boost::unordered_set<uint64_t> neighbor_processes;

	bool balancing_load;



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
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << process
						<< " tried pin cell " << all_new_pinned_cells[process][i]
						<< std::endl;
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

		int
			partition_changed,
			global_id_size,
			local_id_size,
			number_to_receive,
			number_to_send,
			*sender_processes,
			*receiver_processes;

		ZOLTAN_ID_PTR
			global_ids_to_receive,
			local_ids_to_receive,
			global_ids_to_send,
			local_ids_to_send;

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
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Cannot send cell " << global_ids_to_send[i]
						<< " to process " << receiver_processes[i]
						<< std::endl;
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

		this->added_cells.clear();
		this->removed_cells.clear();
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
				this->cells_to_receive[int(current_process_of_cell)].push_back(
					std::make_pair(pin_request->first, -1)
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
					std::make_pair(global_ids_to_receive[i], -1)
				);

				#ifdef DEBUG
				if (this->added_cells.count(global_ids_to_receive[i]) > 0) {
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
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			sender = this->cells_to_receive.begin();
			sender != this->cells_to_receive.end();
			sender++
		) {
			std::sort(sender->second.begin(), sender->second.end());

			for (unsigned int i = 0; i < sender->second.size(); i++) {
				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for receiving cell " << sender->second[i].first
						<< " from process " << sender->first
						<< std::endl;
					abort();
				}
				sender->second[i].second = tag;
			}
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
				this->cells_to_send[int(destination_process)].push_back(
					std::make_pair(pin_request->first, -1)
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
					std::make_pair(global_ids_to_send[i], -1)
				);

				#ifdef DEBUG
				if (this->removed_cells.count(global_ids_to_send[i]) > 0) {
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
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			receiver = this->cells_to_send.begin();
			receiver != this->cells_to_send.end();
			receiver++
		) {
			std::sort(receiver->second.begin(), receiver->second.end());
			for (unsigned int i = 0; i < receiver->second.size(); i++) {
				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for sending cell " << receiver->second[i].first
						<< " to process " << receiver->first
						<< std::endl;
					abort();
				}
				receiver->second[i].second = tag;
			}
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
		BOOST_FOREACH(const uint64_t& cell, this->local_cells_on_process_boundary) {

			#ifdef DEBUG
			if (cell != this->get_child(cell)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell << " has children"
					<< std::endl;
				abort();
			}

			if (this->neighbors.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No neighbor_of list for cell " << cell
					<< std::endl;
				abort();
			}

			if (this->neighbors_to.count(cell) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " No neighbor_to list for cell " << cell
					<< std::endl;
				abort();
			}
			#endif

			const int current_process = int(this->rank);

			// data must be received from neighbors_of
			BOOST_FOREACH(const uint64_t& neighbor, this->neighbors.at(cell)) {

				if (neighbor == error_cell) {
					continue;
				}

				const int other_process = int(this->cell_process.at(neighbor));

				if (other_process != current_process) {
					unique_cells_to_receive[other_process].insert(neighbor);
				}
			}

			// data must be sent to neighbors_to
			BOOST_FOREACH(const uint64_t& neighbor, this->neighbors_to.at(cell)) {

				if (neighbor == error_cell) {
					continue;
				}

				const int other_process = int(this->cell_process.at(neighbor));

				if (other_process != current_process) {
					unique_cells_to_send[other_process].insert(cell);
				}
			}
		}

		// populate final send list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			receiver = unique_cells_to_send.begin();
			receiver != unique_cells_to_send.end();
			receiver++
		) {
			const int process = receiver->first;

			#ifdef DEBUG
			if ((uint64_t) process == this->rank) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << process << " would send to self"
					<< std::endl;
				abort();
			}
			#endif

			this->cells_to_send[process].reserve(receiver->second.size());

			BOOST_FOREACH(const uint64_t& cell, receiver->second) {
				this->cells_to_send.at(process).push_back(
					std::make_pair(cell, -1)
				);
			}

			if (this->cells_to_send.at(process).size() > 0) {
				std::sort(
					this->cells_to_send.at(process).begin(),
					this->cells_to_send.at(process).end()
				);
			}

			// sequential tags for messages: 1, 2, ...
			for (size_t i = 0; i < this->cells_to_send.at(process).size(); i++) {
				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					const uint64_t cell = this->cells_to_send.at(process)[i].first;
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for sending cell " << cell
						<< " to process " << process
						<< std::endl;
					abort();
				}
				this->cells_to_send.at(process)[i].second = tag;
			}
		}

		// populate final receive list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			sender = unique_cells_to_receive.begin();
			sender != unique_cells_to_receive.end();
			sender++
		) {
			const int process = sender->first;

			#ifdef DEBUG
			if ((uint64_t) process == this->rank) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << process << " would receive from self"
					<< std::endl;
				abort();
			}
			#endif

			this->cells_to_receive[process].reserve(sender->second.size());

			BOOST_FOREACH(const uint64_t& cell, sender->second) {
				this->cells_to_receive.at(process).push_back(
					std::make_pair(cell, -1)
				);
			}

			if (this->cells_to_receive.at(process).size() > 0) {
				std::sort(
					this->cells_to_receive.at(process).begin(),
					this->cells_to_receive.at(process).end()
				);
			}

			// sequential tags for messages: 1, 2, ...
			for (size_t i = 0; i < this->cells_to_receive.at(process).size(); i++) {
				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					const uint64_t cell = this->cells_to_receive.at(process)[i].first;
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for receiving cell " << cell
						<< " from process " << process
						<< std::endl;
					abort();
				}

				this->cells_to_receive.at(process)[i].second = tag;
			}
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
	void recalculate_neighbor_update_send_receive_lists(const int neighborhood_id)
	{
		#ifdef DEBUG
		if (this->user_local_cells_on_process_boundary.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighborhood with id " << neighborhood_id
				<< std::endl;
			abort();
		}
		#endif

		// clear previous lists
		this->user_neigh_cells_to_send[neighborhood_id].clear();
		this->user_neigh_cells_to_receive[neighborhood_id].clear();

		boost::unordered_map<int, boost::unordered_set<uint64_t> >
			user_neigh_unique_sends,
			user_neigh_unique_receives;

		// calculate new lists for neighbor data updates
		BOOST_FOREACH(const uint64_t cell, this->user_local_cells_on_process_boundary.at(neighborhood_id)) {

			#ifdef DEBUG
			if (cell != this->get_child(cell)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Cell " << cell << " has children"
					<< std::endl;
				abort();
			}
			#endif

			// data must be received from neighbors_of
			BOOST_FOREACH(const uint64_t& neighbor, this->user_neigh_of.at(neighborhood_id).at(cell)) {

				if (neighbor == error_cell) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					user_neigh_unique_receives[int(this->cell_process.at(neighbor))].insert(neighbor);
				}
			}

			// data must be sent to neighbors_to
			BOOST_FOREACH(const uint64_t& neighbor, this->user_neigh_to.at(neighborhood_id).at(cell)) {

				if (neighbor == error_cell) {
					continue;
				}

				if (this->cell_process.at(neighbor) != this->rank) {
					user_neigh_unique_sends[int(this->cell_process.at(neighbor))].insert(cell);
				}
			}
		}

		// populate final send list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			receiver = user_neigh_unique_sends.begin();
			receiver != user_neigh_unique_sends.end();
			receiver++
		) {
			const int receiving_process = receiver->first;

			std::vector<std::pair<uint64_t, int> >& current_cells_to_send
				= this->user_neigh_cells_to_send.at(neighborhood_id)[receiving_process];

			current_cells_to_send.reserve(receiver->second.size());

			BOOST_FOREACH(const uint64_t& cell, receiver->second) {
				current_cells_to_send.push_back(std::make_pair(cell, -1));
			}

			std::sort(current_cells_to_send.begin(), current_cells_to_send.end());

			// sequential tags for messages: 1, 2, ...
			for (unsigned int i = 0; i < current_cells_to_send.size(); i++) {
				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					const uint64_t cell = current_cells_to_send[i].first;
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for sending cell " << cell
						<< " to process " << receiving_process
						<< std::endl;
					abort();
				}

				current_cells_to_send[i].second = tag;
			}
		}

		// populate final receive list data structures and sort them
		for (boost::unordered_map<int, boost::unordered_set<uint64_t> >::const_iterator
			sender = user_neigh_unique_receives.begin();
			sender != user_neigh_unique_receives.end();
			sender++
		) {
			const int sending_process = sender->first;

			std::vector<std::pair<uint64_t, int> >& current_cells_to_receive
				= this->user_neigh_cells_to_receive[neighborhood_id][sending_process];

			current_cells_to_receive.reserve(sender->second.size());

			BOOST_FOREACH(const uint64_t& cell, sender->second) {
				current_cells_to_receive.push_back(std::make_pair(cell, -1));
			}

			std::sort(current_cells_to_receive.begin(), current_cells_to_receive.end());

			// sequential tags for messages: 1, 2, ...
			for (unsigned int i = 0; i < current_cells_to_receive.size(); i++) {
				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					const uint64_t cell = current_cells_to_receive[i].first;
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for receiving cell " << cell
						<< " from process " << sending_process
						<< std::endl;
					abort();
				}

				current_cells_to_receive[i].second = tag;
			}
		}
	}


	/*!
	Updates neighbor and neighbor_to lists around given cell's neighborhood.

	Does nothing in the following cases:
		- given cell doesn't exist in the grid
		- given cell has children
	Assumes that the refinement level difference between given cell and
	its neighborhood is no larger than 1.
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
	void update_user_neighbors(const uint64_t cell, const int neighborhood_id)
	{
		if (this->user_hood_of.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No user neighborhood with id " << neighborhood_id
				<< std::endl;
			abort();
		}

		#ifdef DEBUG
		if (this->user_hood_to.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No user neighborhood to with id " << neighborhood_id
				<< std::endl;
			abort();
		}
		#endif

		// find neighbors_of, should be in order given by user
		this->user_neigh_of[neighborhood_id][cell].clear();
		BOOST_FOREACH(const Types<3>::neighborhood_item_t& item, this->user_hood_of[neighborhood_id]) {
			std::vector<uint64_t> cells_at_offset
				= this->get_neighbors_of_at_offset(cell, item[0], item[1], item[2]);
			this->user_neigh_of[neighborhood_id][cell].insert(
				this->user_neigh_of[neighborhood_id][cell].end(),
				cells_at_offset.begin(),
				cells_at_offset.end()
			);
		}

		// find neighbors_to
		this->user_neigh_to[neighborhood_id][cell]
			= this->find_neighbors_to(cell, this->user_hood_to[neighborhood_id]);
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

		// TODO: also update remote_cells_on_process_boundary
		this->local_cells_on_process_boundary.erase(cell);

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

			if (neighbor == error_cell) {
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
				this->local_cells_on_process_boundary.insert(cell);
				this->remote_cells_on_process_boundary.insert(neighbor);
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
				this->local_cells_on_process_boundary.insert(cell);
				this->remote_cells_on_process_boundary.insert(neighbor_to);
			}
		}

		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info(cell)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info for cell " << cell
				<< " is not consistent"
				<< std::endl;
			abort();
		}
		#endif
	}

	/*!
	Updates the remote neighbor info of given cell on this process without children.

	Uses current neighbor lists of neighborhood with given id.
	Does nothing if given cell doesn't exist on this process or has children.
	*/
	void update_user_remote_neighbor_info(const uint64_t cell, const int neighborhood_id)
	{
		if (this->cells.count(cell) == 0) {
			return;
		}

		if (cell != this->get_child(cell)) {
			return;
		}

		#ifdef DEBUG
		if (this->user_hood_of.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No user neighborhood with id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_hood_to.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No user neighborhood to with id " << neighborhood_id
				<< std::endl;
			abort();
		}

		if (this->user_neigh_of.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighborhood with id " << neighborhood_id
				<< " in neighbor lists"
				<< std::endl;
			abort();
		}

		if (this->user_neigh_of.at(neighborhood_id).count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighbor list with neighborhood id " << neighborhood_id
				<< " for cell " << cell
				<< std::endl;
			abort();
		}

		if (this->user_neigh_to.count(neighborhood_id) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighborhood with id " << neighborhood_id
				<< " in neighbor_to lists"
				<< std::endl;
			abort();
		}

		if (this->user_neigh_to.at(neighborhood_id).count(cell) == 0) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " No neighbor_to list with neighborhood id " << neighborhood_id
				<< " for cell " << cell
				<< std::endl;
			abort();
		}
		#endif

		this->user_local_cells_on_process_boundary.at(neighborhood_id).erase(cell);

		// neighbors of given cell
		BOOST_FOREACH(const uint64_t& neighbor, this->user_neigh_of.at(neighborhood_id).at(cell)) {

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
				this->user_local_cells_on_process_boundary.at(neighborhood_id).insert(cell);
				this->user_remote_cells_on_process_boundary.at(neighborhood_id).insert(neighbor);
			}
		}

		// cells with given cell as neighbor
		BOOST_FOREACH(const uint64_t& neighbor_to, this->user_neigh_to.at(neighborhood_id).at(cell)) {

			#ifdef DEBUG
			if (neighbor_to == error_cell) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid cell in neighbor_to list of cell " << cell
					<< std::endl;
				abort();
			}

			if (this->cell_process.count(neighbor_to) == 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Neighbor_to " << neighbor_to
					<< " doesn't exist in cell_process"
					<< std::endl;
				abort();
			}
			#endif

			if (this->cell_process.at(neighbor_to) != this->rank) {
				this->user_local_cells_on_process_boundary.at(neighborhood_id).insert(cell);
				this->user_remote_cells_on_process_boundary.at(neighborhood_id).insert(neighbor_to);
			}
		}

		#ifdef DEBUG
		if (!this->verify_remote_neighbor_info(cell)) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info for cell " << cell
				<< " is not consistent"
				<< std::endl;
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
		// TODO this probably can't be optimized without
		// storing neighbor lists also for remote neighbors
		this->local_cells_on_process_boundary.clear();
		this->remote_cells_on_process_boundary.clear();

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
	Updates the remote neighbor info of all cells on this process without children.

	Uses current neighbor lists of neighborhood with given id.
	*/
	void update_user_remote_neighbor_info(const int neighborhood_id)
	{
		this->user_local_cells_on_process_boundary[neighborhood_id].clear();
		this->user_remote_cells_on_process_boundary[neighborhood_id].clear();

		BOOST_FOREACH(const cell_and_data_pair_t& item, this->cells) {

			if (item.first != this->get_child(item.first)) {
				continue;
			}

			this->update_user_remote_neighbor_info(item.first, neighborhood_id);

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
		const uint64_t cell1_length = this->get_cell_length_in_indices(cell1);
		const uint64_t cell2_length = this->get_cell_length_in_indices(cell2);

		// distance in indices between given cells
		Types<3>::indices_t distance = {{0, 0, 0}};

		const uint64_t grid_length[3] = {
			this->get_length_x() * (uint64_t(1) << this->max_refinement_level),
			this->get_length_y() * (uint64_t(1) << this->max_refinement_level),
			this->get_length_z() * (uint64_t(1) << this->max_refinement_level)
		};

		uint64_t max_distance = 0;

		for (unsigned int i = 0; i < 3; i++) {
			if (indices1[i] <= indices2[i]) {
				if (indices2[i] <= indices1[i] + cell1_length) {
					distance[i] = 0;
				} else {
					distance[i] = indices2[i] - (indices1[i] + cell1_length);
				}

				if (this->is_periodic(i)) {
					const uint64_t distance_to_end = grid_length[i] - (indices2[i] + cell2_length);
					distance[i] = std::min(distance[i], indices1[i] + distance_to_end);
				}
			} else {
				if (indices1[i] <= indices2[i] + cell2_length) {
					distance[i] = 0;
				} else {
					distance[i] = indices1[i] - (indices2[i] + cell2_length);
				}

				if (this->is_periodic(i)) {
					const uint64_t distance_to_end = grid_length[i] - (indices1[i] + cell1_length);
					distance[i] = std::min(distance[i], indices2[i] + distance_to_end);
				}
			}

			max_distance = std::max(max_distance, distance[i]);
		}

		if (this->neighborhood_length == 0) {
			if (max_distance < cell1_length
			&& this->overlapping_indices(cell1, cell2) >= 2) {
				return true;
			// diagonal cell isn't a neighbor
			} else {
				return false;
			}
		}

		if (max_distance < this->neighborhood_length * cell1_length) {
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

					if (this->remote_cells_on_process_boundary.count(refined) == 0) {
						continue;
					}

					// refine all local cells that are too large and neighboring the refined cell
					/*
					TODO: probably faster to search for local neighbors of refined
					cell, even faster would be to also store neighbors lists of
					remote cells with local neighbors
					*/
					BOOST_FOREACH(const uint64_t& local, this->local_cells_on_process_boundary) {

						if (this->is_neighbor(local, refined)
						&& this->get_refinement_level(local) < this->get_refinement_level(refined)
						&& this->cells_to_refine.count(local) == 0) {
							unique_induced_refines.insert(local);
						}
					}
				}
			}
			all_new_refines.clear();

			new_refines.insert(
				new_refines.end(),
				unique_induced_refines.begin(),
				unique_induced_refines.end()
			);
			this->cells_to_refine.insert(
				unique_induced_refines.begin(),
				unique_induced_refines.end()
			);
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
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Invalid refinement level for parent"
					<< std::endl;
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
			this->cells_to_unrefine.insert(
				all_unrefines[process].begin(),
				all_unrefines[process].end()
			);
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
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Remote neighbor info is not consistent"
				<< std::endl;
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
		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
		this->incoming_data.clear();
		this->outgoing_data.clear();
		#endif

		#ifdef DEBUG
		// check that cells_to_refine is identical between processes
		std::vector<uint64_t> ordered_cells_to_refine(this->cells_to_refine.begin(), this->cells_to_refine.end());
		std::sort(ordered_cells_to_refine.begin(), ordered_cells_to_refine.end());

		std::vector<std::vector<uint64_t> > all_ordered_cells_to_refine;
		All_Gather()(ordered_cells_to_refine, all_ordered_cells_to_refine, this->comm);

		for (unsigned int process = 0; process < this->comm_size; process++) {
			if (!std::equal(
				all_ordered_cells_to_refine[process].begin(),
				all_ordered_cells_to_refine[process].end(),
				all_ordered_cells_to_refine[0].begin()
			)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " cells_to_refine differ between processes 0 and " << process
					<< std::endl;
				exit(EXIT_FAILURE);
			}
		}

		// check that cells_to_unrefine is identical between processes
		std::vector<uint64_t> ordered_cells_to_unrefine(this->cells_to_unrefine.begin(), this->cells_to_unrefine.end());
		std::sort(ordered_cells_to_unrefine.begin(), ordered_cells_to_unrefine.end());

		std::vector<std::vector<uint64_t> > all_ordered_cells_to_unrefine;
		All_Gather()(ordered_cells_to_unrefine, all_ordered_cells_to_unrefine, this->comm);

		for (unsigned int process = 0; process < this->comm_size; process++) {
			if (!std::equal(
				all_ordered_cells_to_unrefine[process].begin(),
				all_ordered_cells_to_unrefine[process].end(),
				all_ordered_cells_to_unrefine[0].begin()
			)) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " cells_to_unrefine differ between processes 0 and " << process
					<< std::endl;
				exit(EXIT_FAILURE);
			}
		}
		#endif

		// cells whose neighbor lists have to be updated afterwards
		boost::unordered_set<uint64_t> update_neighbors;

		// a separate neighborhood update function has to be used
		// for cells whose children were removed by unrefining
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

			// without using local neighbor lists figure out rest of the
			// neighbor lists that need updating
			if (this->remote_cells_on_process_boundary.count(refined) > 0) {

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

		// initially only one sibling is recorded per process when unrefining,
		// insert the rest of them now
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
					std::make_pair(unrefined, -1)
				);

			// receive user data of removed cell from its process
			} else if (this->rank == process_of_parent) {
				this->cells_to_receive[process_of_unrefined].push_back(
					std::make_pair(unrefined, -1)
				);
			}
		}

		// receive cells in known order and add message tags
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			sender = this->cells_to_receive.begin();
			sender != this->cells_to_receive.end();
			sender++
		) {
			std::sort(sender->second.begin(), sender->second.end());
			// TODO: merge with identical code in make_new_partition
			for (unsigned int i = 0; i < sender->second.size(); i++) {

				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for receiving cell " << sender->second[i].first
						<< " from process " << sender->first
						<< std::endl;
					abort();
				}

				sender->second[i].second = tag;
			}
		}

		// send cells in known order and add message tags
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::iterator
			receiver = this->cells_to_send.begin();
			receiver != this->cells_to_send.end();
			receiver++
		) {
			std::sort(receiver->second.begin(), receiver->second.end());
			// TODO: check that message tags don't overflow
			for (unsigned int i = 0; i < receiver->second.size(); i++) {

				const int tag = (int) i + 1;
				if (tag > (int) this->max_tag) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Process " << this->rank
						<< ": Message tag would overflow for sending cell " << receiver->second[i].first
						<< " to process " << receiver->first
						<< std::endl;
					abort();
				}

				receiver->second[i].second = tag;
			}
		}


		this->start_user_data_transfers(
			this->unrefined_cell_data,
			this->cells_to_receive,
			this->cells_to_send
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
		// also remote neighbor info of user neighborhoods
		for (boost::unordered_map<int, std::vector<Types<3>::neighborhood_item_t> >::const_iterator
			item = this->user_hood_of.begin();
			item != this->user_hood_of.end();
			item++
		) {
			this->update_user_remote_neighbor_info(item->first);
		}

		#ifdef DEBUG
		if (!this->verify_neighbors()) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " Neighbor lists are inconsistent"
				<< std::endl;
			exit(EXIT_FAILURE);
		}
		#endif

		this->wait_user_data_transfer_receives(
		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
		this->unrefined_cell_data, this->cells_to_receive
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
		for (boost::unordered_set<uint64_t>::const_iterator
			unrefined = all_to_unrefine.begin();
			unrefined != all_to_unrefine.end();
			unrefined++
		) {
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
	data then call this beforehand and use get_remote_neighbors_to() to
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

	User data arriving to this process is saved in given destination except when using boost
	and sending all cells in one message in which case the destination is actually used by
	wait_user_data_transfer_receives(...).
	*/
	bool start_user_data_transfers(
		boost::unordered_map<uint64_t, Cell_Data>& destination,
		const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >& receive_item,
		const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >& send_item
	) {

		#ifndef DCCRG_TRANSFER_USING_BOOST_MPI
		int ret_val = -1;
		#endif

		// post receives
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			sender = receive_item.begin();
			sender != receive_item.end();
			sender++
		) {
			const int sending_process = sender->first;

			#if defined(DEBUG) || !defined(DCCRG_TRANSFER_USING_BOOST_MPI)
			const size_t number_of_receives = sender->second.size();
			#endif

			#ifdef DEBUG
			if (sending_process == (int) this->rank
			&& number_of_receives > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Process " << this->rank
					<< " trying to transfer to self"
					<< std::endl;
				abort();
			}
			#endif

			if (this->send_single_cells) {

				for (std::vector<std::pair<uint64_t, int> >::const_iterator
					item = sender->second.begin();
					item != sender->second.end();
					item++
				) {
					const uint64_t cell = item->first;

					if (destination.count(cell) == 0) {
						destination[cell];
					}

					#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

					this->receive_requests[sending_process].push_back(
						this->boost_comm.irecv(
							sending_process,
							item->second,
							destination.at(cell)
						)
					);

					#else // ifdef DCCRG_TRANSFER_USING_BOOST_MPI

					this->receive_requests[sending_process].push_back(MPI_Request());

					void* address = NULL;
					int count = -1;
					MPI_Datatype user_datatype = MPI_DATATYPE_NULL;
					destination.at(cell).mpi_datatype(
						address,
						count,
						user_datatype,
						cell,
						sending_process,
						(int) this->rank,
						true
					);

					const bool is_named_datatype = Is_Named_Datatype()(user_datatype);

					if (!is_named_datatype) {
						ret_val = MPI_Type_commit(&user_datatype);
						if (ret_val != MPI_SUCCESS) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " MPI_Type_commit failed on process " << this->rank
								<< ", for datatype of cell " << cell
								<< " with returned buffer address " << address
								<< " and returned count " << count
								<< ": " << Error_String()(ret_val)
								<< std::endl;
							abort();
						}
					}

					ret_val = MPI_Irecv(
						address,
						count,
						user_datatype,
						sending_process,
						item->second,
						this->comm,
						&(this->receive_requests[sending_process].back())
					);

					if (ret_val != MPI_SUCCESS) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " MPI_Irecv failed on process " << this->rank
							<< ", for cell " << cell
							<< ", from process " << sending_process
							<< ": " << Error_String()(ret_val)
							<< std::endl;
						abort();
					}

					if (!is_named_datatype) {
						ret_val = MPI_Type_free(&user_datatype);
						if (ret_val != MPI_SUCCESS) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " MPI_Type_free failed on process " << this->rank
								<< ", for a derived datatype of cell " << cell
								<< ": " << Error_String()(ret_val)
								<< std::endl;
							abort();
						}
					}

					#endif
				}

			} else { // if this->send_single_cells

				#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

				this->receive_requests[sending_process].push_back(
					this->boost_comm.irecv(
						sending_process,
						0,
						this->incoming_data[sending_process]
					)
				);

				#else // ifdef DCCRG_TRANSFER_USING_BOOST_MPI

				// reserve space for incoming user data in this end
				// TODO: move into a separate function callable by user
				for (size_t i = 0; i < number_of_receives; i++) {
					const uint64_t cell = sender->second[i].first;
					if (destination.count(cell) == 0) {
						destination[cell];
					}
				}

				// get mpi transfer info from cells
				std::vector<void*> addresses(number_of_receives, NULL);
				std::vector<int> counts(number_of_receives, -1);
				std::vector<MPI_Datatype> datatypes(number_of_receives, MPI_DATATYPE_NULL);

				for (size_t i = 0; i < number_of_receives; i++) {
					const uint64_t cell = sender->second[i].first;
					destination.at(cell).mpi_datatype(
						addresses[i],
						counts[i],
						datatypes[i],
						cell,
						sending_process,
						(int) this->rank,
						true
					);
				}

				// get displacements in bytes for incoming user data
				std::vector<MPI_Aint> displacements(number_of_receives, 0);
				for (size_t i = 0; i < number_of_receives; i++) {
					displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
				}

				MPI_Datatype receive_datatype;
				ret_val = MPI_Type_create_struct(
					(int) number_of_receives,
					&counts[0],
					&displacements[0],
					&datatypes[0],
					&receive_datatype
				);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Type_create_struct failed for process " << this->rank
						<< ": " << Error_String()(ret_val)
						<< std::endl;
					abort();
				}

				ret_val = MPI_Type_commit(&receive_datatype);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Type_commit failed for process " << this->rank
						<< ": " << Error_String()(ret_val)
						<< std::endl;
					abort();
				}

				this->receive_requests[sender->first].push_back(MPI_Request());

				ret_val = MPI_Irecv(
					addresses[0],
					1,
					receive_datatype,
					sender->first,
					0,
					this->comm,
					&(this->receive_requests[sender->first].back())
				);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Irecv failed for process " << this->rank
						<< ", source process " << sender->first
						<< ": " << Error_String()(ret_val)
						<< std::endl;
					abort();
				}

				MPI_Type_free(&receive_datatype);
				BOOST_FOREACH(MPI_Datatype& type, datatypes) {
					if (!Is_Named_Datatype()(type)) {
						ret_val = MPI_Type_free(&type);
						if (ret_val != MPI_SUCCESS) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " MPI_Type_free failed on process " << this->rank
								<< ", for a derived datatype of user data: "
								<< Error_String()(ret_val)
								<< std::endl;
							abort();
						}
					}
				}
				#endif // ifdef DCCRG_TRANSFER_USING_BOOST_MPI
			}
		}

		// post sends
		for (boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >::const_iterator
			receiver = send_item.begin();
			receiver != send_item.end();
			receiver++
		) {
			const int receiving_process = receiver->first;

			#if defined(DEBUG) || !defined(DCCRG_TRANSFER_USING_BOOST_MPI)
			const size_t number_of_sends = receiver->second.size();
			#endif

			#ifdef DEBUG
			if (receiving_process == (int) this->rank
			&& number_of_sends > 0) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " Trying to transfer to self"
					<< std::endl;
				abort();
			}
			#endif

			if (this->send_single_cells) {

				for (std::vector<std::pair<uint64_t, int> >::const_iterator
					item = receiver->second.begin();
					item != receiver->second.end();
					item++
				) {
					const uint64_t cell = item->first;

					#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

					this->send_requests[receiving_process].push_back(
						this->boost_comm.isend(
							receiving_process,
							item->second,
							this->cells.at(cell)
						)
					);

					#else // ifdef DCCRG_TRANSFER_USING_BOOST_MPI

					this->send_requests[receiving_process].push_back(MPI_Request());

					void* address = NULL;
					int count = -1;
					MPI_Datatype user_datatype = MPI_DATATYPE_NULL;

					this->cells.at(cell).mpi_datatype(
						address,
						count,
						user_datatype,
						cell,
						(int) this->rank,
						receiving_process,
						false
					);

					const bool is_named_datatype = Is_Named_Datatype()(user_datatype);

					if (!is_named_datatype) {
						ret_val = MPI_Type_commit(&user_datatype);
						if (ret_val != MPI_SUCCESS) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " MPI_Type_commit failed on process " << this->rank
								<< ", for datatype of cell " << cell
								<< " with returned buffer address " << address
								<< " and returned count " << count
								<< ": " << Error_String()(ret_val)
								<< std::endl;
							abort();
						}
					}

					ret_val = MPI_Isend(
						address,
						count,
						user_datatype,
						receiving_process,
						item->second,
						this->comm,
						&(this->send_requests[receiving_process].back())
					);

					if (ret_val != MPI_SUCCESS) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " MPI_Isend failed on process " << this->rank
							<< ", for cell " << cell
							<< ", to process " << receiving_process
							<< ": " << Error_String()(ret_val)
							<< std::endl;
						abort();
					}

					if (!is_named_datatype) {
						ret_val = MPI_Type_free(&user_datatype);
						if (ret_val != MPI_SUCCESS) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " MPI_Type_free failed on process " << this->rank
								<< ", for a derived datatype of cell " << cell
								<< ": " << Error_String()(ret_val)
								<< std::endl;
							abort();
						}
					}

					#endif
				}

			} else { // if this->send_single_cells


				#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

				// construct the outgoing data vector
				for (std::vector<std::pair<uint64_t, int> >::const_iterator
					item = receiver->second.begin();
					item != receiver->second.end();
					item++
				) {
					const uint64_t cell = item->first;
					Cell_Data* user_data = &(this->cells.at(cell));
					if (user_data == NULL) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " No data for cell " << cell
							<< std::endl;
						abort();
					}
					this->outgoing_data[receiving_process].push_back(*user_data);
				}

				// send all cells
				this->send_requests[receiving_process].push_back(
					this->boost_comm.isend(
						receiving_process,
						0,
						this->outgoing_data[receiving_process]
					)
				);

				#else	// ifdef DCCRG_TRANSFER_USING_BOOST_MPI

				// get mpi transfer info from cells
				std::vector<void*> addresses(number_of_sends, NULL);
				std::vector<int> counts(number_of_sends, -1);
				std::vector<MPI_Datatype> datatypes(number_of_sends, MPI_DATATYPE_NULL);

				for (size_t i = 0; i < number_of_sends; i++) {
					const uint64_t cell = receiver->second[i].first;
					this->cells.at(cell).mpi_datatype(
						addresses[i],
						counts[i],
						datatypes[i],
						cell,
						(int) this->rank,
						receiving_process,
						false
					);
				}

				// get displacements in bytes for outgoing user data
				std::vector<MPI_Aint> displacements(number_of_sends, 0);
				for (size_t i = 0; i < number_of_sends; i++) {
					displacements[i] = (uint8_t*) addresses[i] - (uint8_t*) addresses[0];
				}

				MPI_Datatype send_datatype;

				ret_val = MPI_Type_create_struct(
					(int) number_of_sends,
					&counts[0],
					&displacements[0],
					&datatypes[0],
					&send_datatype
				);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Type_create_struct failed for process " << this->rank
						<< ": " << Error_String()(ret_val)
						<< std::endl;
					abort();
				}

				MPI_Type_commit(&send_datatype);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Type_commit failed for process " << this->rank
						<< ": " << Error_String()(ret_val)
						<< std::endl;
					abort();
				}

				this->send_requests[receiver->first].push_back(MPI_Request());

				ret_val = MPI_Isend(
					addresses[0],
					1,
					send_datatype,
					receiver->first,
					0,
					this->comm,
					&(this->send_requests[receiver->first].back())
				);
				if (ret_val != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Isend failed from process " << this->rank
						<< ", target process " << receiver->first
						<< ": " << Error_String()(ret_val)
						<< std::endl;
					abort();
				}

				MPI_Type_free(&send_datatype);
				BOOST_FOREACH(MPI_Datatype& type, datatypes) {
					if (!Is_Named_Datatype()(type)) {
						ret_val = MPI_Type_free(&type);
						if (ret_val != MPI_SUCCESS) {
							std::cerr << __FILE__ << ":" << __LINE__
								<< " MPI_Type_free failed on process " << this->rank
								<< ", for a derived datatype of user data: "
								<< Error_String()(ret_val)
								<< std::endl;
							abort();
						}
					}
				}

				#endif	// ifdef DCCRG_TRANSFER_USING_BOOST_MPI
			}
		}

		return true;
	}


	/*!
	Waits for the receives of user data transfers between processes to complete.

	User data arriving to this process is saved in given destination.
	*/
	bool wait_user_data_transfer_receives(
	#ifdef DCCRG_TRANSFER_USING_BOOST_MPI
	boost::unordered_map<uint64_t, Cell_Data>& destination,
	const boost::unordered_map<int, std::vector<std::pair<uint64_t, int> > >& receive_item
	#endif
	) {
		bool success = true;

		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		// wait for data to arrive
		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator
			process = this->receive_requests.begin();
			process != this->receive_requests.end();
			process++
		) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}

		// incorporate received data
		if (!this->send_single_cells) {

			for (typename boost::unordered_map<int, std::vector<Cell_Data> >::const_iterator
				sender = this->incoming_data.begin();
				sender != this->incoming_data.end();
				sender++
			) {
				const int sending_process = sender->first;

				size_t i = 0;
				// cells were received in the following order
				for (std::vector<std::pair<uint64_t, int> >::const_iterator
					cell_item = receive_item.at(sending_process).begin();
					cell_item != receive_item.at(sending_process).end();
					cell_item++
				) {
					const uint64_t cell = cell_item->first;
					// TODO move data instead of copying
					destination[cell] = this->incoming_data.at(sending_process)[i];
					i++;
				}
			}
			this->incoming_data.clear();
		}

		#else	// ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		int ret_val = -1;

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator
			process = this->receive_requests.begin();
			process != this->receive_requests.end();
			process++
		) {
			std::vector<MPI_Status> statuses;
			statuses.resize(process->second.size());

			ret_val = MPI_Waitall(process->second.size(), &(process->second[0]), &(statuses[0]));
			if (ret_val != MPI_SUCCESS) {
				BOOST_FOREACH(const MPI_Status& status, statuses) {
					if (status.MPI_ERROR != MPI_SUCCESS) {
						success = false;
						std::cerr << __FILE__ << ":" << __LINE__
							<< " MPI receive failed from process " << status.MPI_SOURCE
							<< " with tag " << status.MPI_TAG
							<< std::endl;
					}
				}
			}
		}

		#endif	// ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		this->receive_requests.clear();

		return success;
	}


	/*!
	Waits for the sends of user data transfers between processes to complete.
	*/
	bool wait_user_data_transfer_sends()
	{
		#ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		for (boost::unordered_map<int, std::vector<boost::mpi::request> >::iterator
			process = this->send_requests.begin();
			process != this->send_requests.end();
			process++
		) {
			boost::mpi::wait_all(process->second.begin(), process->second.end());
		}

		if (!this->send_single_cells) {
			this->outgoing_data.clear();
		}

		#else	// ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		for (boost::unordered_map<int, std::vector<MPI_Request> >::iterator
			process = this->send_requests.begin();
			process != this->send_requests.end();
			process++
		) {
			std::vector<MPI_Status> statuses;
			statuses.resize(process->second.size());

			if (
				MPI_Waitall(process->second.size(), &(process->second[0]), &(statuses[0]))
				!= MPI_SUCCESS
			) {
				BOOST_FOREACH(const MPI_Status& status, statuses) {
					if (status.MPI_ERROR != MPI_SUCCESS) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " MPI receive failed from process " << status.MPI_SOURCE
							<< " with tag " << status.MPI_TAG
							<< std::endl;
						abort();
					}
				}
			}
		}

		#endif	// ifdef DCCRG_TRANSFER_USING_BOOST_MPI

		this->send_requests.clear();

		return true;
	}


	/*!
	Returns true if cells with given index properties overlap.

	Sizes are also given in indices.
	*/
	bool indices_overlap(const uint64_t index1, const uint64_t size1, const uint64_t index2, const uint64_t size2) const
	{
		#ifdef DEBUG
		if (index1 >= this->get_length_x() * (uint64_t(1) << this->max_refinement_level)
		&& index1 >= this->get_length_y() * (uint64_t(1) << this->max_refinement_level)
		&& index1 >= this->get_length_z() * (uint64_t(1) << this->max_refinement_level)) {
			std::cerr << __FILE__ << ":" << __LINE__ << " Invalid index given" << std::endl;
			exit(EXIT_FAILURE);
		}

		if (index2 >= this->get_length_x() * (uint64_t(1) << this->max_refinement_level)
		&& index2 >= this->get_length_y() * (uint64_t(1) << this->max_refinement_level)
		&& index2 >= this->get_length_z() * (uint64_t(1) << this->max_refinement_level)) {
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
		if (cell1 == error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell1 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell2 == error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell2 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		const uint64_t index1 = this->get_x_index(cell1);
		const uint64_t index2 = this->get_x_index(cell2);
		const uint64_t size1 = this->get_cell_length_in_indices(cell1);
		const uint64_t size2 = this->get_cell_length_in_indices(cell2);

		return this->indices_overlap(index1, size1, index2, size2);
	}

	/*!
	Returns true if y indices of given cells overlap, even if they don't exist
	*/
	bool y_indices_overlap(const uint64_t cell1, const uint64_t cell2) const
	{
		if (cell1 == error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell1 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell2 == error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell2 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		const uint64_t index1 = this->get_y_index(cell1);
		const uint64_t index2 = this->get_y_index(cell2);
		const uint64_t size1 = this->get_cell_length_in_indices(cell1);
		const uint64_t size2 = this->get_cell_length_in_indices(cell2);

		return this->indices_overlap(index1, size1, index2, size2);
	}

	/*!
	Returns true if z indices of given cells overlap, even if they don't exist
	*/
	bool z_indices_overlap(const uint64_t cell1, const uint64_t cell2) const
	{
		if (cell1 == error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell1 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell2 == error_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}
		if (cell2 > this->last_cell) {
			std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
			abort();
		}

		const uint64_t index1 = this->get_z_index(cell1);
		const uint64_t index2 = this->get_z_index(cell2);
		const uint64_t size1 = this->get_cell_length_in_indices(cell1);
		const uint64_t size2 = this->get_cell_length_in_indices(cell2);

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

		const uint64_t size1 = this->get_cell_length_in_indices(cell1);
		const uint64_t size2 = this->get_cell_length_in_indices(cell2);

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
		if (indices[0] >= this->length_x * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (indices[1] >= this->length_y * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (indices[2] >= this->length_z * (uint64_t(1) << this->max_refinement_level)) {
			return error_cell;
		}

		if (minimum_refinement_level > maximum_refinement_level) {
			return error_cell;
		}

		int average_refinement_level
			= (maximum_refinement_level + minimum_refinement_level) / 2;

		const uint64_t average_cell
			= this->get_cell_from_indices(indices, average_refinement_level);

		// use binary search recursively (assumes that all cells refine to 8 children)
		if (this->cell_process.count(average_cell) == 0) {

			// search for larger cell
			if (average_refinement_level > minimum_refinement_level) {

				uint64_t larger_cell
					= this->get_existing_cell(
						indices,
						minimum_refinement_level,
						average_refinement_level - 1
					);

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
				uint64_t smaller_cell
					= this->get_existing_cell(
						indices,
						average_refinement_level + 1,
						maximum_refinement_level
					);

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
	static void fill_with_cell_coordinates(
		void *data,
		int /*global_id_size*/,
		int /*local_id_size*/,
		int number_of_cells,
		ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR /*local_ids*/,
		int /*number_of_dimensions*/,
		double *geom_vec,
		int *error
	) {
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
		*error = ZOLTAN_OK;
		return int(dccrg_instance->cells.size());
	}


	/*!
	Writes all cell ids on this process to the global_ids array
	*/
	static void fill_cell_list(
		void* data,
		int /*global_id_size*/,
		int /*local_id_size*/,
		ZOLTAN_ID_PTR global_ids,
		ZOLTAN_ID_PTR /*local_ids*/,
		int number_of_weights_per_object,
		float* object_weights,
		int* error
	) {
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
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
					object_weights[i] = float(dccrg_instance->cell_weights.at(item.first));
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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
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
					= int(dccrg_instance->cell_process.at(neighbor));

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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
		*error = ZOLTAN_OK;

		*number_of_hyperedges = int(dccrg_instance->cells.size());
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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
		*error = ZOLTAN_OK;

		*number_of_edge_weights = int(dccrg_instance->cells.size());
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
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
		*error = ZOLTAN_OK;

		if ((unsigned int) number_of_hyperedges != dccrg_instance->cells.size()) {
			std::cerr
				<< "Zoltan is expecting wrong number of hyperedges: " << number_of_hyperedges
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

				hyperedge_weights[i] = float(1.0 * number_of_hyperedges);
			}

			i++;
		}
	}


	/*!
	Returns the number of hierarchies to use for load balancing.
	*/
	static int get_number_of_load_balancing_hierarchies(void* data, int* error)
	{
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);
		*error = ZOLTAN_OK;
		return int(dccrg_instance->processes_per_part.size());
	}


	/*!
	Returns the part number of this process on given hierarchy level of load balancing.
	*/
	static int get_part_number(void* data, int level, int* error)
	{
		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);

		if (level < 0 || level >= int(dccrg_instance->processes_per_part.size())) {
			std::cerr
				<< "Zoltan wanted a part number for an invalid hierarchy level (should be [0, "
				<< dccrg_instance->processes_per_part.size() - 1
				<< "]): " << level
				<< std::endl;
			*error = ZOLTAN_FATAL;
			return -1;
		} else {
			*error = ZOLTAN_OK;
		}

		int process = int(dccrg_instance->rank);
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
	static void set_partitioning_options(
		void* data,
		int level,
		struct Zoltan_Struct* zz,
		int* error
	) {
		if (zz == NULL) {
			std::cerr << "Zoltan gave a NULL pointer for zz" << std::endl;
			*error = ZOLTAN_FATAL;
			return;
		}

		Dccrg<Cell_Data, Geometry>* dccrg_instance
			= reinterpret_cast<Dccrg<Cell_Data, Geometry> *>(data);

		if (level < 0 || level >= int(dccrg_instance->processes_per_part.size())) {
			std::cerr
				<< "Zoltan wanted partitioning options for an invalid hierarchy "
				<< "level (level should be between 0 and "
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

	Remote neighbor info consists of local_cells_on_process_boundary and
	remote_cells_on_process_boundary.
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

			if (neighbor == error_cell) {
				continue;
			}

			if (this->cell_process.at(neighbor) != this->rank) {

				if (this->local_cells_on_process_boundary.count(cell) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Local cell " << cell
						<< " should be in local_cells_on_process_boundary"
						<< std::endl;
					return false;
				}

				if (this->remote_cells_on_process_boundary.count(neighbor) == 0) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " Remote cell " << neighbor
						<< " should be in remote_cells_on_process_boundary"
						<< std::endl;
					return false;
				}
			}
		}

		return true;
	}


	/*!
	Returns false if remote neighbor info on this process is inconsistent.

	Remote neighbor info consists of local_cells_on_process_boundary
	and remote_cells_on_process_boundary.
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

			// check whether this cell should be in remote_cells_on_process_boundary
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

					if (this->remote_cells_on_process_boundary.count(item->first) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Remote cell " << item->first
							<< " should be in remote_cells_on_process_boundary because:"
							<< std::endl;

						BOOST_FOREACH(const cell_and_data_pair_t& cell, this->cells) {
							if (item->first == cell.first) {
								std::cerr << __FILE__ << ":" << __LINE__
									<< " Same cell."
									<< std::endl;
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

					if (this->remote_cells_on_process_boundary.count(item->first) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Remote cell " << item->first
							<< " should not be in remote_cells_on_process_boundary"
							<< std::endl;
						return false;
					}
				}

			// check whether this cell should be in cells_with_remote_neighbor
			} else {

				bool no_remote_neighbor = true;

				// search in neighbors_of
				const std::vector<uint64_t> neighbors_of
					= this->find_neighbors_of(
						item->first,
						this->neighborhood_of,
						this->max_ref_lvl_diff
					);

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
					if (this->local_cells_on_process_boundary.count(item->first) > 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Local cell " << item->first
							<< " should not be in local_cells_on_process_boundary"
							<< std::endl;
						return false;
					}
				} else {
					if (this->local_cells_on_process_boundary.count(item->first) == 0) {
						std::cerr << __FILE__ << ":" << __LINE__
							<< " Local cell " << item->first
							<< " should be in local_cells_on_process_boundary"
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

