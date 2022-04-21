/*
 * FunctionBinaryExpression.hpp
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

#ifndef FUNCTION_BINARY_EXPRESSION_HPP
#define FUNCTION_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
class FunctionBinaryExpression : public FunctionExpression< BinaryOpType<Op, LeftType, RightType>::Type, FunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
public:
	typedef typename Traits< FunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::ReturnType ReturnType;
	typedef typename Traits< FunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::LeftExpr_ReturnType LeftExpr_ReturnType;
	typedef typename Traits< FunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::RightExpr_ReturnType RightExpr_ReturnType;
public:
	FunctionBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getMesh() != m_rhs.getMesh()) { throw std::invalid_argument("lhs and rhs must be defined on the same mesh."); } }
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_lhs[e], m_rhs[e]); }
public:
	inline bool containsTrial() const { return m_lhs.containsTrial() or m_rhs.containsTrial(); }
	inline bool containsTest()  const { return m_lhs.containsTest() or m_rhs.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_lhs.getMesh(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

////////////////////////////////////////////////////////////////////////

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
class CpxFunctionBinaryExpression : public CpxFunctionExpression< BinaryOpType<Op, LeftType, RightType>::Type, CpxFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
public:
	typedef typename Traits< CpxFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >::ReturnType ReturnType;
public:
	CpxFunctionBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getMesh() != m_rhs.getMesh()) { throw std::invalid_argument("lhs and rhs must be defined on the same mesh."); } }
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_lhs[e], m_rhs[e]); }
public:
	inline bool containsTrial() const { return m_lhs.containsTrial() or m_rhs.containsTrial(); }
	inline bool containsTest()  const { return m_lhs.containsTest() or m_rhs.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_lhs.getMesh(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

#endif // FUNCTION_BINARY_EXPRESSION_HPP
