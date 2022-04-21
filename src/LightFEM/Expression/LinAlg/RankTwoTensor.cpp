/*
 * RankTwoTensor.cpp
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

#include <LightFEM/Expression/LinAlg/RankTwoTensor.hpp>

RankTwoTensor::RankTwoTensor(const size_t shape1, const size_t shape2, const double value)
{
	m_shape[0] = shape1;
	m_shape[1] = shape2;
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()), value);
}

////////////////////////////////////////////////////////////////////////

CpxRankTwoTensor::CpxRankTwoTensor(const size_t shape1, const size_t shape2, const std::complex<double> value)
{
	m_shape[0] = shape1;
	m_shape[1] = shape2;
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()), value);
}
