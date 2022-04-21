/*
 * FixedSizeMatrix.hpp
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

#ifndef FIXED_SIZE_MATRIX_HPP
#define FIXED_SIZE_MATRIX_HPP

#include <array>

template<size_t Nrows, size_t Ncols>
class FixedSizeMatrix
{
public:
	FixedSizeMatrix() : m_core({}) {}
	FixedSizeMatrix(const FixedSizeMatrix& other) : m_core(other.m_core) {}
public:
	double  operator() (const size_t i, const size_t j) const { return m_core[j + i*Ncols]; }
	double& operator() (const size_t i, const size_t j)       { return m_core[j + i*Ncols]; }

	constexpr size_t getNRows() const { return Nrows; }
	constexpr size_t getNCols() const { return Ncols; }
private:
	std::array< double, Nrows*Ncols > m_core;
};

#endif // FIXED_SIZE_MATRIX_HPP
