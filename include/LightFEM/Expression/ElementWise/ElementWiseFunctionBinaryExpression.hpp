/*
 * ElementWiseFunctionBinaryExpression.hpp
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

#ifndef ELEMENT_WISE_FUNCTION_BINARY_EXPRESSION_HPP
#define ELEMENT_WISE_FUNCTION_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>

#include <LightFEM/Expression/Function/FunctionExpression.hpp>

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
class ElementWiseFunctionBinaryExpression : public ElementWiseFunctionExpression< BinaryOpType<Op, LeftType, RightType>::Type, ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
public:
	typedef typename Traits< ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::ReturnType ReturnType;
public:
	ElementWiseFunctionBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (   m_lhs.getMesh() != m_rhs.getMesh())    { throw std::invalid_argument("lhs and rhs must be defined on the same mesh."); } 
		if (m_lhs.getElement() != m_rhs.getElement()) { throw std::invalid_argument("lhs and rhs must be defined on the same element."); } 
	}
public:
	inline ReturnType operator[] (const size_t k) const { return ReturnType(m_lhs[k], m_rhs[k]); }
public:
	inline bool containsTrial() const { return m_lhs.containsTrial() or m_rhs.containsTrial(); }
	inline bool containsTest()  const { return m_lhs.containsTest() or m_rhs.containsTest(); }
public:
	inline const Mesh*    getMesh()    const { return m_lhs.getMesh();    }
	inline const Element* getElement() const { return m_lhs.getElement(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

////////////////////////////////////////////////////////////////////////

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
class CpxElementWiseFunctionBinaryExpression : public CpxElementWiseFunctionExpression< BinaryOpType<Op, LeftType, RightType>::Type, CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
public:
	typedef typename Traits< CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >::ReturnType ReturnType;
public:
	CpxElementWiseFunctionBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getElement() != m_rhs.getElement()) { throw std::invalid_argument("lhs and rhs must be defined on the same element."); } }
public:
	inline ReturnType operator[] (const size_t k) const { return ReturnType(m_lhs[k], m_rhs[k]); }
public:
	inline bool containsTrial() const { return m_lhs.containsTrial() or m_rhs.containsTrial(); }
	inline bool containsTest()  const { return m_lhs.containsTest() or m_rhs.containsTest(); }
public:
	inline const Mesh*    getMesh()    const { return m_lhs.getMesh();    }
	inline const Element* getElement() const { return m_lhs.getElement(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

#endif // ELEMENT_WISE_FUNCTION_BINARY_EXPRESSION_HPP
