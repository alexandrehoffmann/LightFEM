/*
 * RankFourTensorExpression.hpp
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

#ifndef RANK_FOUR_TENSOR_EXPRESSION_HPP
#define RANK_FOUR_TENSOR_EXPRESSION_HPP

#include <complex>
#include <tuple>

template <typename Expr>
class RankFourTensorExpression 
{
// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
public:
	inline std::array<size_t, 4> getShape() const { return static_cast<Expr const&>(*this).getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return static_cast<Expr const&>(*this)(i,j,k,l); }
};

template <typename Expr>
class CpxRankFourTensorExpression 
{
// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
public:
	inline std::array<size_t, 4> getShape() const { return static_cast<Expr const&>(*this).getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return static_cast<Expr const&>(*this)(i,j,k,l); }
};

#endif // RANK_FOUR_TENSOR_EXPRESSION_HPP
