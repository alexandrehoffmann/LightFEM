/*
 * ElementWiseFunctionUnaryExpression.hpp
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

#ifndef ELEMENT_WISE_FUNCTION_UNARY_EXPRESSION_HPP
#define ELEMENT_WISE_FUNCTION_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>

template<UnaryOp Op, ExprType Type, typename Expr>
class ElementWiseFunctionUnaryExpression : public ElementWiseFunctionExpression< UnaryOpType<Op, Type>::Type, ElementWiseFunctionUnaryExpression<Op, Type, Expr> >
{
public:
	typedef typename Traits< ElementWiseFunctionUnaryExpression<Op, Type, Expr> >::ReturnType ReturnType;
public:
	ElementWiseFunctionUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline ReturnType operator[] (const size_t k) const { return ReturnType(m_expr[k]); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Element* getElement() const { return m_expr.getElement(); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

////////////////////////////////////////////////////////////////////////

template<UnaryOp Op, ExprType Type, typename Expr>
class CpxElementWiseFunctionUnaryExpression : public CpxElementWiseFunctionExpression< UnaryOpType<Op, Type>::Type, CpxElementWiseFunctionUnaryExpression<Op, Type, Expr> >
{
public:
	typedef typename Traits< CpxElementWiseFunctionUnaryExpression<Op, Type, Expr> >::ReturnType ReturnType;
public:
	CpxElementWiseFunctionUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline ReturnType operator[] (const size_t k) const { return ReturnType(m_expr[k]); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Element* getElement() const { return m_expr.getElement(); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

#endif // ELEMENT_WISE_FUNCTION_UNARY_EXPRESSION_HPP
