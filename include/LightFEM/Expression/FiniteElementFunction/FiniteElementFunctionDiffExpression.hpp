/*
 * FiniteElementFunctionDiffExpression.hpp
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

#ifndef FINITE_ELEMENT_FUNCTION_DIFF_EXPRESSION_HPP
#define FINITE_ELEMENT_FUNCTION_DIFF_EXPRESSION_HPP

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>

template<ExprType Type> struct GradOpType {};
template<> struct GradOpType<ExprType::SCALAR> { static constexpr ExprType Type = ExprType::VECTOR; };
template<> struct GradOpType<ExprType::VECTOR> { static constexpr ExprType Type = ExprType::MATRIX; };

template<ExprType Type, typename Expr>
class FiniteElementFunctionDiffExpression : public FunctionExpression< Type, FiniteElementFunctionDiffExpression<Type, Expr> >
{
public:
	typedef typename Traits< FiniteElementFunctionDiffExpression<Type, Expr> >::ReturnType ReturnType;
public:
	FiniteElementFunctionDiffExpression(Expr&& expr, const size_t dim) : m_expr(std::forward<Expr>( expr )), m_dim(dim) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_expr.getD(m_dim, e)); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_expr.getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
	size_t m_dim;
};

template<ExprType Type, typename Expr>
class FiniteElementFunctionGradExpression : public FunctionExpression< GradOpType<Type>::Type, FiniteElementFunctionGradExpression<Type, Expr> >
{
public:
	typedef typename Traits< FiniteElementFunctionGradExpression<Type, Expr> >::ReturnType ReturnType;
public:
	FiniteElementFunctionGradExpression(Expr&& expr) : m_expr(std::forward<Expr>( expr )) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_expr.getGrad(e)); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_expr.getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
};

////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr>
class CpxFiniteElementFunctionDiffExpression : public CpxFunctionExpression< Type, CpxFiniteElementFunctionDiffExpression<Type, Expr> >
{
public:
	typedef typename Traits< CpxFiniteElementFunctionDiffExpression<Type, Expr> >::ReturnType ReturnType;
public:
	CpxFiniteElementFunctionDiffExpression(Expr&& expr, const size_t dim) : m_expr(std::forward<Expr>( expr )), m_dim(dim) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_expr.getD(m_dim, e)); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_expr.getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
	size_t m_dim;
};

template<ExprType Type, typename Expr>
class CpxFiniteElementFunctionGradExpression : public CpxFunctionExpression< GradOpType<Type>::Type, CpxFiniteElementFunctionGradExpression<Type, Expr> >
{
public:
	typedef typename Traits< CpxFiniteElementFunctionGradExpression<Type, Expr> >::ReturnType ReturnType;
public:
	CpxFiniteElementFunctionGradExpression(Expr&& expr) : m_expr(std::forward<Expr>( expr )) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_expr.getGrad(e)); }
public:
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
public:
	inline const Mesh* getMesh() const { return m_expr.getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
};

#endif // FINITE_ELEMENT_FUNCTION_DIFF_EXPRESSION_HPP
