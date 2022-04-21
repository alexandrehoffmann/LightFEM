/*
 * MpiRange.hpp
 * 
 * Copyright 2022 Alexandre Hoffmann <alexandre.hoffmann.etu@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#ifndef MPI_RANGE_HPP
#define MPI_RANGE_HPP

#include <mpi.h>

class MpiRange
{
public:
	MpiRange(MPI_Comm com, const size_t ntasks);

	inline int getRank() const { return m_rank; }
	
	inline size_t begin() const { return m_begin; } 
	inline size_t end()   const { return m_end; } 
private:
	int m_rank;
	
	size_t m_begin;
	size_t m_end;
};

#endif // MPI_RANGE_HPP
