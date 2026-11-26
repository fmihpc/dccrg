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
template<class Cell_Data_T>
std::tuple<void*, int, MPI_Datatype> get_cell_mpi_datatype(
	Cell_Data_T& cell,
	const uint64_t cell_id,
	const int sender,
	const int receiver,
	const bool receiving,
	const int neighborhood_id
) {
	return cell.get_mpi_datatype(
		cell_id,
		sender,
		receiver,
		receiving,
		neighborhood_id
	);
}

// Compatibility macros
#define UNUSED_GMD_ARGS \
	const uint64_t cell_id, \
	const int sender,       \
	const int receiver,     \
	const bool receiving,   \
	const int neighborhood_id
#define IGNORE_UNUSED_GMD_ARGS do {\
	(void)cell_id;             \
	(void)sender;              \
	(void)receiver;            \
	(void)receiving;           \
	(void)neighborhood_id;     \
	} while(0)
#define DEFAULT_GMD_ARGS 0,0,0,0,0
// struct types weapper
template<typename T, std::size_t N, MPI_Datatype D>
struct array_wrapper:std::array<T, N> {
	std::tuple<void*, int, MPI_Datatype>
	inline get_mpi_datatype(UNUSED_GMD_ARGS) {
		IGNORE_UNUSED_GMD_ARGS;
		auto cell = (std::array<T, N>)(*this);
		return std::make_tuple((void*) cell.data(), cell.size(), D);
	}
	T operator[] (std::size_t i) const {
		return std::array<T, N>::operator[](i);
	}
	T& operator[] (std::size_t i) {
		return std::array<T, N>::operator[](i);
	}
};
template<typename T, std::size_t N, MPI_Datatype M, typename OS> // OS for std::ostream, w/o including header
OS& operator<< (OS& os, const array_wrapper<T, N, M> x) {
    return os << (std::array<T, N>)x;
}
// primitive types weapper
template<typename T>
struct wrapper {
	T value;
	/*
	std::tuple<void*, int, MPI_Datatype>
	get_mpi_datatype(UNUSED_GMD_ARGS) {
		return std::make_tuple((void*)this, 1, MPI_##type);
	}
	*/
};
template<typename T, typename OS> // OS for std::ostream, w/o including header
OS& operator<< (OS& os, const wrapper<T> x) {
    return os << x.value;
}
#define DCCRG_GET_MPI_DATATYPE_BASIC(T, type)                               \
	struct dccrg_##type:wrapper<T> {                                    \
		T value;						    \
		std::tuple<void*, int, MPI_Datatype>                        \
		get_mpi_datatype(UNUSED_GMD_ARGS) {                         \
			IGNORE_UNUSED_GMD_ARGS;				    \
			return std::make_tuple((void*)this, 1, MPI_##type); \
		}                                                           \
		inline operator T() const { return this->value; } 	    \
		dccrg_##type(const T x): value(x) {} 			    \
		dccrg_##type() {} 					    \
	};
DCCRG_GET_MPI_DATATYPE_BASIC(char, CHAR)
//DCCRG_GET_MPI_DATATYPE_BASIC(signed char, CHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned char, UNSIGNED_CHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(short int, SHORT)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned short int, UNSIGNED_SHORT)
DCCRG_GET_MPI_DATATYPE_BASIC(int, INT)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned int, UNSIGNED)
DCCRG_GET_MPI_DATATYPE_BASIC(long int, LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned long int, UNSIGNED_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(long long int, LONG_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(unsigned long long int, UNSIGNED_LONG_LONG)
DCCRG_GET_MPI_DATATYPE_BASIC(float, FLOAT)
DCCRG_GET_MPI_DATATYPE_BASIC(double, DOUBLE)
DCCRG_GET_MPI_DATATYPE_BASIC(long double, LONG_DOUBLE)
DCCRG_GET_MPI_DATATYPE_BASIC(wchar_t, WCHAR)
DCCRG_GET_MPI_DATATYPE_BASIC(bool, CXX_BOOL)
#ifdef DCCRG_USER_COMPLEX
DCCRG_GET_MPI_DATATYPE_BASIC(std::complex<float>, CXX_FLOAT_COMPLEX)
DCCRG_GET_MPI_DATATYPE_BASIC(std::complex<double>, CXX_DOUBLE_COMPLEX)
DCCRG_GET_MPI_DATATYPE_BASIC(std::complex<long double>, CXX_LONG_DOUBLE_COMPLEX)
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
