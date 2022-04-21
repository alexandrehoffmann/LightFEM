/*
 * VectorExpression.hpp
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

#ifndef VECTOR_EXPRESSION_HPP
#define VECTOR_EXPRESSION_HPP

#include <complex>

template <typename Expr>
class VectorExpression 
{
public:
	// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
	inline size_t getSize() const { return static_cast<Expr const&>(*this).getSize(); }
public:
	inline double operator[](const size_t i) const { return static_cast<Expr const&>(*this)[i]; }
};

template <typename Expr>
class CpxVectorExpression 
{
public:
	// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
	inline size_t getSize() const { return static_cast<Expr const&>(*this).getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return static_cast<Expr const&>(*this)[i]; }
};

#endif // VECTOR_EXPRESSION_HPP
