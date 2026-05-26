/*
Functions for dccrg to obtain the MPI Datatype from cell data.

Copyright 2014, 2015, 2016 Ilja Honkonen
Copyright 2018 Finnish Meteorological Institute

Dccrg is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

Dccrg is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with dccrg. If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef DCCRG_GET_CELL_DATATYPE_HPP
#define DCCRG_GET_CELL_DATATYPE_HPP


#include "cstdint"
#include "tuple"
#include "type_traits"

#include "mpi.h"

namespace dccrg {
namespace detail {

/*!
Returns the MPI transfer info from given cell.

Version for get_mpi_datatype(const uint64_t, ..., const int).
*/
#define gen_get_cell_mpi_datatype(CONST_ATTR, RESTRAINT, INPUTS...)     \
template<class Cell_Data_T>                                  \
RESTRAINT						     \
std::tuple<void*, int, MPI_Datatype> get_cell_mpi_datatype(  \
	CONST_ATTR Cell_Data_T& cell,                        \
	const uint64_t cell_id,    \
	const int sender,          \
	const int receiver,        \
	const bool receiving,      \
	const int neighborhood_id  \
) {                                                          \
	return cell.get_mpi_datatype(INPUTS);                \
}                                                            
#define CELL_MPI_DATATYPE\
	cell_id,         \
	sender,          \
	receiver,        \
	receiving,       \
	neighborhood_id
gen_get_cell_mpi_datatype(/*non-const*/, /*no constraint*/, CELL_MPI_DATATYPE)
gen_get_cell_mpi_datatype(const        , /*no constraint*/, CELL_MPI_DATATYPE)
#if __cplusplus >= 202002L
template<typename T>
concept has_mem_fn_get_mpi_datatype_const = requires (const T& x, const uint64_t a, const int b, const int c, const bool d, const int e) {
	{x.get_mpi_datatype(a,b,c,d,e)} -> std::same_as<std::tuple<void*, int, MPI_Datatype>>;
};
template<typename T>
concept has_mem_fn_get_mpi_datatype_noArgs_const = requires (const T& x) {
	{x.get_mpi_datatype()} -> std::same_as<std::tuple<void*, int, MPI_Datatype>>;
	!has_mem_fn_get_mpi_datatype_const<T>;
};
template<typename T>
concept has_mem_fn_get_mpi_datatype = requires (T& x, const uint64_t a, const int b, const int c, const bool d, const int e) {
	{x.get_mpi_datatype(a,b,c,d,e)} -> std::same_as<std::tuple<void*, int, MPI_Datatype>>;
	!has_mem_fn_get_mpi_datatype_const<T>;
	!has_mem_fn_get_mpi_datatype_noArgs_const<T>;
};
template<typename T>
concept has_mem_fn_get_mpi_datatype_noArgs = requires (T& x) {
	{x.get_mpi_datatype()} -> std::same_as<std::tuple<void*, int, MPI_Datatype>>;
	!has_mem_fn_get_mpi_datatype<T>;
	!has_mem_fn_get_mpi_datatype_const<T>;
	!has_mem_fn_get_mpi_datatype_noArgs_const<T>;
};

/* Requires C++20
 * For the following to be not a redefinition of above
 */
gen_get_cell_mpi_datatype(/*non-const*/, requires has_mem_fn_get_mpi_datatype_noArgs      <Cell_Data_T>,/*no args*/)
gen_get_cell_mpi_datatype(const        , requires has_mem_fn_get_mpi_datatype_noArgs_const<Cell_Data_T>, /*no args*/)

/*!
Returns the MPI transfer info from given cell.

Version for cell that doesn't have get_mpi_datatype().
*/
template<class Cell_Data>
requires (!(
	//has_mem_fn_get_mpi_datatype_noArgs<Cell_Data> ||
	//has_mem_fn_get_mpi_datatype_noArgs_const<Cell_Data> ||
	//has_mem_fn_get_mpi_datatype_const<Cell_Data> ||
	has_mem_fn_get_mpi_datatype<Cell_Data>))
std::tuple<void*, int, MPI_Datatype> get_cell_mpi_datatype(
	Cell_Data& cell,
	const uint64_t /*cell_id*/,
	const int /*sender*/,
	const int /*receiver*/,
	const bool /*receiving*/,
	const int /*neighborhood_id*/
) {
	// give a human-readable error message
	static_assert(
		not std::is_same<Cell_Data, Cell_Data>::value,
		"Cell_Data given to dccrg is not a supported type and "
			"doesn't have get_mpi_datatype() member function either"
	);
	return std::make_tuple(nullptr, -1, MPI_DATATYPE_NULL);
}
#endif

#define DCCRG_GET_MPI_DATATYPE_BASIC(CPP, MPI) \
	std::tuple< \
		void*, int, MPI_Datatype \
	> inline get_mpi_datatype_basic(CPP& cell) { \
		return std::make_tuple((void*) &cell, 1, MPI); \
	}
DCCRG_GET_MPI_DATATYPE_BASIC(char, MPI_CHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(signed char, MPI_CHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned char, MPI_UNSIGNED_CHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(short int, MPI_SHORT)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned short int, MPI_UNSIGNED_SHORT)
DCCRG_GET_MPI_DATATYPE_BASIC(int, MPI_INT)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned int, MPI_UNSIGNED)
DCCRG_GET_MPI_DATATYPE_BASIC(long int, MPI_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned long int, MPI_UNSIGNED_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(long long int, MPI_LONG_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned long long int, MPI_UNSIGNED_LONG_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(float, MPI_FLOAT)
DCCRG_GET_MPI_DATATYPE_BASIC(double, MPI_DOUBLE)
DCCRG_GET_MPI_DATATYPE_BASIC(long double, MPI_LONG_DOUBLE)
DCCRG_GET_MPI_DATATYPE_BASIC(wchar_t, MPI_WCHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(bool, MPI_CXX_BOOL)
#ifdef DCCRG_USER_COMPLEX
DCCRG_GET_MPI_DATATYPE_BASIC(std::complex<float>, MPI_CXX_FLOAT_COMPLEX)
DCCRG_GET_MPI_DATATYPE_BASIC(std::complex<double>, MPI_CXX_DOUBLE_COMPLEX)
DCCRG_GET_MPI_DATATYPE_BASIC(std::complex<long double>, MPI_CXX_LONG_DOUBLE_COMPLEX)
#endif
#undef DCCRG_GET_MPI_DATATYPE_BASIC

#define DCCRG_GET_MPI_DATATYPE_ARRAY(CPP, MPI) \
	template<std::size_t N> std::tuple< \
		void*, int, MPI_Datatype \
	> inline get_mpi_datatype_basic(std::array<CPP, N>& cell) { \
		return std::make_tuple((void*) cell.data(), cell.size(), MPI); \
	}
DCCRG_GET_MPI_DATATYPE_ARRAY(char, MPI_CHAR)
DCCRG_GET_MPI_DATATYPE_ARRAY(signed char, MPI_CHAR)
DCCRG_GET_MPI_DATATYPE_ARRAY(unsigned char, MPI_UNSIGNED_CHAR)
DCCRG_GET_MPI_DATATYPE_ARRAY(short int, MPI_SHORT)
DCCRG_GET_MPI_DATATYPE_ARRAY(unsigned short int, MPI_UNSIGNED_SHORT)
DCCRG_GET_MPI_DATATYPE_ARRAY(int, MPI_INT)
DCCRG_GET_MPI_DATATYPE_ARRAY(unsigned int, MPI_UNSIGNED)
DCCRG_GET_MPI_DATATYPE_ARRAY(long int, MPI_LONG)
DCCRG_GET_MPI_DATATYPE_ARRAY(unsigned long int, MPI_UNSIGNED_LONG)
DCCRG_GET_MPI_DATATYPE_ARRAY(long long int, MPI_LONG_LONG)
DCCRG_GET_MPI_DATATYPE_ARRAY(unsigned long long int, MPI_UNSIGNED_LONG_LONG)
DCCRG_GET_MPI_DATATYPE_ARRAY(float, MPI_FLOAT)
DCCRG_GET_MPI_DATATYPE_ARRAY(double, MPI_DOUBLE)
DCCRG_GET_MPI_DATATYPE_ARRAY(long double, MPI_LONG_DOUBLE)
DCCRG_GET_MPI_DATATYPE_ARRAY(wchar_t, MPI_WCHAR)
#ifdef DCCRG_USER_COMPLEX
DCCRG_GET_MPI_DATATYPE_ARRAY(std::complex<float>, MPI_CXX_FLOAT_COMPLEX)
DCCRG_GET_MPI_DATATYPE_ARRAY(std::complex<double>, MPI_CXX_DOUBLE_COMPLEX)
DCCRG_GET_MPI_DATATYPE_ARRAY(std::complex<long double>, MPI_CXX_LONG_DOUBLE_COMPLEX)
#endif
#undef DCCRG_GET_MPI_DATATYPE_ARRAY

}} // namespaces

#endif
