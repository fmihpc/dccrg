/*
Supporting routines for MPI of dccrg.

Copyright 2012 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_MPI_SUPPORT_HPP
#define DCCRG_MPI_SUPPORT_HPP

#include "algorithm"
#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"
#include "iostream"
#include "mpi.h"
#include "stdint.h"
#include "string"
#include "vector"

namespace dccrg {


/*!
Returns the error string of given MPI error.
*/
class Error_String
{
public:

	std::string operator()(int mpi_return_value)
	{
		char mpi_error_string[MPI_MAX_ERROR_STRING + 1];
		int string_length;
		MPI_Error_string(mpi_return_value, mpi_error_string, &string_length);
		mpi_error_string[string_length + 1] = '\0';

		const std::string result(mpi_error_string);
		return result;
	}
};


/*!
Wrapper for MPI_Allgatherv(..., uint64_t, ...).
*/
class All_Gather
{
public:

	void operator()(
	) {
	}
};


/*!
Wrapper for MPI_Allreduce(uint64_t, ..., MPI_SUM).
*/
class All_Reduce
{
public:

	uint64_t operator()(
		uint64_t value,
		MPI_Comm& comm
	) {
		uint64_t result;
		MPI_Allreduce(&value, &result, 1, MPI_UINT64_T, MPI_SUM, comm);
		return result;
	}

};


/*!
Similar to MPI_Allreduce but communicates only with processes given in neighbors.

Given neighbors must not include the process itself.
Any pair of processes in the given communicator either must or must not have each
other in neighbors.
Does not modify the given value or communicator (will be const correct once MPI is).
Example:
MPI_Comm w = MPI_COMM_WORLD;
result = dccrg::Some_Reduce()(3, n, w);
*/
class Some_Reduce
{
public:

	uint64_t operator()(
		uint64_t value,
		const boost::unordered_set<int>& neighbors,
		MPI_Comm& comm
	) {
		// send / receive with asynchronous messages to / from each neighbor
		// TODO: doesn't seem possible to reduce number of messages even if
		//       > 2 processes are all neighbors of each other
		boost::unordered_map<int, std::pair<MPI_Request, uint64_t> > receive_requests;
		boost::unordered_map<int, MPI_Request> send_requests;

		// post receives
		BOOST_FOREACH(const int process, neighbors) {
			receive_requests[process];
			if (
				MPI_Irecv(
					&(receive_requests.at(process).second),
					sizeof(uint64_t),
					MPI_BYTE,
					process,
					1,
					comm,
					&(receive_requests.at(process).first)
				) != MPI_SUCCESS
			) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " MPI_Irecv from process " << process
					<< " failed"
					<< std::endl;
				abort();
			}
		}

		// post sends
		BOOST_FOREACH(const int process, neighbors) {
			send_requests[process];
			if (
				MPI_Isend(
					&value,
					sizeof(uint64_t),
					MPI_BYTE,
					process,
					1,
					comm,
					&(send_requests.at(process))
				) != MPI_SUCCESS
			) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " MPI_Isend to process " << process
					<< " failed"
					<< std::endl;
				abort();
			}
		}

		// wait for receives
		BOOST_FOREACH(const int process, neighbors) {
			if (
				MPI_Wait(&(receive_requests.at(process).first), MPI_STATUS_IGNORE) != MPI_SUCCESS
			) {
				std::cerr << __FILE__ << ":" << __LINE__
					<< " MPI_Wait for receive from process " << process
					<< " failed"
					<< std::endl;
				abort();
			}
		}

		// reduce
		uint64_t result = value;
		BOOST_FOREACH(const int process, neighbors) {
			result += receive_requests.at(process).second;
		}

		// wait for sends
		BOOST_FOREACH(const int process, neighbors) {
			if (MPI_Wait(&(send_requests.at(process)), MPI_STATUS_IGNORE) != MPI_SUCCESS) {
					std::cerr << __FILE__ << ":" << __LINE__
						<< " MPI_Wait for send to process " << process
						<< " failed"
						<< std::endl;
					abort();
			}
		}

		return result;
	}

}; // class

} // namespace

#endif

