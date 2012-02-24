/*
Datatype definitions used by dccrg

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

#ifndef DCCRG_TYPES_HPP
#define DCCRG_TYPES_HPP

#include "boost/array.hpp"
#include "stdint.h"

namespace dccrg {

/*!
\brief Class wrapper around typedefs used in dccrg.

Only C++11 allows typedefs to have templated parameters,
this class wrapper is used until enough compilers support C++11.
*/
template <unsigned int Dimensions> class Types
{
public:

/*!
\brief Defines the indices of a cell in the grid, first value is the cell's index in x direction, second in the y direction, etc.

Indices start from 0 and the index of a cell is the one closest to the starting corner of the grid within the cell.
Indices are of the same size as the smallest possible cell in the grid, e.g. the size in indices of cells of maximum refinement level is 1.
For example in the following 2d grid the indices of cell 2 are (2, 0) and cell 8 are (1, 1) assuming a maximum refinement level of 1.
\verbatim
--------->  x
|3|4|   |
|-|-| 2 |
|7|8|   |
|--------
V

y
\endverbatim

*/
typedef boost::array<uint64_t, Dimensions> indices_t;

/*!
\brief Defines one item of cells' neighborhood in the grid.

Defines an offset (first value is offset in x direction, etc.) of one neighbor to a cell of the same size as the cell itself.
A cell's neighborhood (e.g. stencil) in the grid consists of a list of offsets that define which cells a cell considers as neighbors.
Offsets are given relative to the cell and cells within a volume of equal size as the cell are considered the cell's neighbors.
For example in the case of the game of life in two dimensions (x and y), neighborhood items making up the neighborhood would be:
\verbatim
-1, -1, 0
-1,  0, 0
-1, +1, 0
 0, -1, 0
 0, +1, 0
+1, -1, 0
+1,  0, 0
+1, +1, 0
\endverbatim
*/
typedef boost::array<int, Dimensions> neighborhood_item_t;

};	// class

}	// namespace
#endif

