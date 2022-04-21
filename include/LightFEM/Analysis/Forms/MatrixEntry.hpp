/*
 * MatrixEntry.hpp
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

#ifndef MATRIX_ENTRY_HPP
#define MATRIX_ENTRY_HPP

#include <complex>

struct MatrixEntry 
{ 
	MatrixEntry(const int _row=0, const int _col=0, const double _value=0.0) : row(_row), col(_col), value(_value) {}
	MatrixEntry(const MatrixEntry& other) : row(other.row), col(other.col), value(other.value) {}
	
	int row;
	int col;
	double value;
};


////////////////////////////////////////////////////////////////////////

struct CpxMatrixEntry 
{ 
	CpxMatrixEntry(const int _row=0, const int _col=0, const std::complex< double > _value=0.0) : row(_row), col(_col), value(_value) {}
	CpxMatrixEntry(const CpxMatrixEntry& other) : row(other.row), col(other.col), value(other.value) {}
	
	int row;
	int col;
	std::complex< double > value;
};

#endif // MATRIX_ENTRY_HPP
