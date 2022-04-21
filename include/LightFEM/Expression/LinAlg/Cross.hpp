/*
 * Cross.hpp
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

#ifndef CROSS_HPP
#define CROSS_HPP

template<> struct BinaryOpType<BinaryOp::CROSS, ExprType::VECTOR, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public VectorExpression< BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline size_t getSize() const { return 3; }
public:
	inline double operator[](const size_t i) const;
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline size_t getSize() const { return 3; }
public:
	inline std::complex< double > operator[](const size_t i) const;
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

#include <LightFEM/Expression/LinAlg/Cross.tpp>

#endif // CROSS_HPP
