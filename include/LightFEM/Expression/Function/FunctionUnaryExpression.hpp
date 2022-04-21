/*
 * FunctionUnaryExpression.hpp
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

#ifndef FUNCTION_UNARY_EXPRESSION_HPP
#define FUNCTION_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>

template<UnaryOp Op, ExprType Type, typename Expr>
class FunctionUnaryExpression : public FunctionExpression< UnaryOpType<Op, Type>::Type, FunctionUnaryExpression<Op, Type, Expr> >
{
public:
	typedef typename Traits< FunctionUnaryExpression<Op, Type, Expr> >::ReturnType ReturnType;
public:
	FunctionUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_expr[e]); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_expr.getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
};

////////////////////////////////////////////////////////////////////////

template<UnaryOp Op, ExprType Type, typename Expr>
class CpxFunctionUnaryExpression : public CpxFunctionExpression< UnaryOpType<Op, Type>::Type, CpxFunctionUnaryExpression<Op, Type, Expr> >
{
public:
	typedef typename Traits< CpxFunctionUnaryExpression<Op, Type, Expr> >::ReturnType ReturnType;
public:
	CpxFunctionUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_expr[e]); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_expr.getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
};

#endif // FUNCTION_UNARY_EXPRESSION_HPP
