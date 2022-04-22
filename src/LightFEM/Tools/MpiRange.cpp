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

MpiRange::MpiRange(MPI_Comm com, const size_t ntasks) : 
	m_ntasks(ntasks)
{
	MPI_Comm_rank(com, &m_rank);
	MPI_Comm_size(com, &m_nprocess);
	
	const int maxRank = m_nprocess-1;
	const size_t ntasksPerProcess = size_t(double(m_ntasks) / double(m_nprocess));
	
	m_begin = m_rank*ntasksPerProcess;
	m_end = (m_rank == maxRank) ? m_ntasks : (m_rank + 1)*ntasksPerProcess;
}

int MpiRange::getWorker(const size_t taskId) const
{
	const int maxRank = m_nprocess-1;
	
	const size_t ntasksPerProcess = size_t(double(m_ntasks) / double(m_nprocess));
	for (int rank=0;rank<m_nprocess;++rank)
	{
		const size_t begin = rank*ntasksPerProcess;
		const size_t end = (rank == maxRank) ? m_ntasks : (rank + 1)*ntasksPerProcess;
		if (begin <= taskId and taskId < end) { return rank; }
	}
	return maxRank;
}
