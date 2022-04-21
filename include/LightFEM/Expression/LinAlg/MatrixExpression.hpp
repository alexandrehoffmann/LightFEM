/*
 * MatrixExpression.hpp
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

#ifndef MATRIX_EXPRESSION_HPP
#define MATRIX_EXPRESSION_HPP

#include <complex>

template <typename Expr>
class MatrixExpression 
{
public:
	// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
	inline size_t getNrows() const { return static_cast<Expr const&>(*this).getNrows(); }
	inline size_t getNcols() const { return static_cast<Expr const&>(*this).getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return static_cast<Expr const&>(*this)(i,j); }
};

template <typename Expr>
class CpxMatrixExpression 
{
public:
	// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
	inline size_t getNrows() const { return static_cast<Expr const&>(*this).getNrows(); }
	inline size_t getNcols() const { return static_cast<Expr const&>(*this).getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return static_cast<Expr const&>(*this)(i,j); }
};

#endif // MATRIX_EXPRESSION_HPP
