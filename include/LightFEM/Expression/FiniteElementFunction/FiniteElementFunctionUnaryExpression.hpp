/*
 * FiniteElementFunctionUnaryExpression.hpp
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

#ifndef FINITE_ELEMENT_FUNCTION_UNARY_EXPRESSION_HPP
#define FINITE_ELEMENT_FUNCTION_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/CoefUnaryOP.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

template<UnaryOp Op, ExprType Type, typename Expr>
class FiniteElementFunctionUnaryExpression : public FiniteElementFunctionExpression< UnaryOpType<Op, Type>::Type, FiniteElementFunctionUnaryExpression<Op, Type, Expr> >
{
public:
	typedef typename Traits< FiniteElementFunctionUnaryExpression<Op, Type, Expr> >::ReturnType ReturnType;
	typedef typename Traits< FiniteElementFunctionUnaryExpression<Op, Type, Expr> >::GradReturnType GradReturnType;
	typedef typename Traits< FiniteElementFunctionUnaryExpression<Op, Type, Expr> >::DiffReturnType DiffReturnType;

	typedef typename FSpace< UnaryOpType<Op, Type>::Type >::Type FSpaceType;
public:
	FiniteElementFunctionUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)), m_coefUnaryOp(std::forward<Expr>(expr)) {}
public:
	inline ReturnType     operator[] (const size_t e)                   const { return ReturnType(m_expr[e]); }
	inline GradReturnType getGrad    (const size_t e)                   const { return GradReturnType(m_expr.getGrad(e)); }
	inline DiffReturnType getD       (const size_t dim, const size_t e) const { return DiffReturnType(getGrad(e), dim); }
	inline double         getCoef    (const size_t globId)              const { return m_coefUnaryOp.getCoef(globId); }
	
	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }
	
	inline const FSpaceType* getFunctionSpace() const { return m_expr.getFunctionSpace(); }
	inline const Mesh*       getMesh()          const { return getFunctionSpace()->getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
	CoefUnaryOp<Op, Expr> m_coefUnaryOp;
};

////////////////////////////////////////////////////////////////////////

template<UnaryOp Op, ExprType Type, typename Expr>
class CpxFiniteElementFunctionUnaryExpression : public CpxFiniteElementFunctionExpression< UnaryOpType<Op, Type>::Type, CpxFiniteElementFunctionUnaryExpression<Op, Type, Expr> >
{
public:
	typedef typename Traits< CpxFiniteElementFunctionUnaryExpression<Op, Type, Expr> >::ReturnType ReturnType;
	typedef typename Traits< CpxFiniteElementFunctionUnaryExpression<Op, Type, Expr> >::GradReturnType GradReturnType;
	typedef typename Traits< CpxFiniteElementFunctionUnaryExpression<Op, Type, Expr> >::DiffReturnType DiffReturnType;

	typedef typename FSpace< UnaryOpType<Op, Type>::Type >::Type FSpaceType;
public:
	CpxFiniteElementFunctionUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)), m_coefUnaryOp(std::forward<Expr>(expr)) {}
public:
	inline ReturnType     operator[] (const size_t e)                   const { return ReturnType(m_expr[e]); }
	inline GradReturnType getGrad    (const size_t e)                   const { return GradReturnType(m_expr.getGrad(e)); }
	inline DiffReturnType getD       (const size_t dim, const size_t e) const { return DiffReturnType(getGrad(e), dim); }
	inline double         getCoef    (const size_t globId)              const { return m_coefUnaryOp.getCoef(globId); }

	inline bool containsTrial() const { return m_expr.containsTrial(); }
	inline bool containsTest()  const { return m_expr.containsTest(); }

	inline const FSpaceType* getFunctionSpace() const { return m_expr.getFunctionSpace(); }
	inline const Mesh*       getMesh()          const { return getFunctionSpace()->getMesh(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
	CpxCoefUnaryOp<Op, Expr> m_coefUnaryOp;
};

#endif // FINITE_ELEMENT_FUNCTION_UNARY_EXPRESSION_HPP
