/*
 * Outer.hpp
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

#ifndef OUTER_HPP
#define OUTER_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/MatrixExpression.hpp>

template<> struct BinaryOpType<BinaryOp::OUTER, ExprType::VECTOR, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public MatrixExpression< BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getSize() != m_rhs.getSize()) { throw std::invalid_argument("the two vectors must share the same size"); } }
public:
	inline size_t getNrows() const { return m_lhs.getSize(); }
	inline size_t getNcols() const { return m_rhs.getSize(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_lhs[i]*m_rhs[j]; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getSize() != m_rhs.getSize()) { throw std::invalid_argument("the two vectors must share the same size"); } }
public:
	inline size_t getNrows() const { return m_lhs.getSize(); }
	inline size_t getNcols() const { return m_rhs.getSize(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_lhs[i]*m_rhs[j]; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

#endif // OUTER_HP
