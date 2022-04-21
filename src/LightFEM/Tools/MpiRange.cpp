/*
 * MpiRange.cpp
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

#include <LightFEM/Tools/MpiRange.hpp>

#include <cmath>
#include <iostream>

MpiRange::MpiRange(MPI_Comm com, const size_t ntasks)
{
	int nProcess;
	
	MPI_Comm_rank(com, &m_rank);
	MPI_Comm_size(com, &nProcess);
	
	const int maxRank = nProcess-1;
	const size_t ntasksPerProcess = ntasks / nProcess;
	
	m_begin = m_rank*ntasksPerProcess;
	if (m_rank == maxRank)
	{
		m_end = ntasks;
	}
	else
	{
		m_end = (m_rank + 1)*ntasksPerProcess;
	}
}
