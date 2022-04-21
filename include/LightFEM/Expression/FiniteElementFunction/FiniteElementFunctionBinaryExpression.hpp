/*
 * FiniteElementFunctionBinaryExpression.hpp
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

#ifndef FINITE_ELEMENT_FUNCTION_BINARY_EXPRESSION_HPP
#define FINITE_ELEMENT_FUNCTION_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/CoefBinaryOp.hpp>
#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
class FiniteElementFunctionBinaryExpression : public FiniteElementFunctionExpression< BinaryOpType<Op, LeftType, RightType>::Type, FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
public:
	typedef typename Traits< FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::ReturnType ReturnType;
	typedef typename Traits< FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::GradReturnType GradReturnType;
	typedef typename Traits< FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::DiffReturnType DiffReturnType;
	
	typedef typename FSpace< BinaryOpType<Op, LeftType, RightType>::Type >::Type FSpaceType;
public:
	FiniteElementFunctionBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)),
		m_coefBinaryOp(std::forward<LeftExpr>(lhs), std::forward<RightExpr>(rhs)) {}
public:
	inline ReturnType     operator[] (const size_t e)                   const { return ReturnType(m_lhs[e], m_rhs[e]); }
	inline GradReturnType getGrad    (const size_t e)                   const { return GradReturnType(m_lhs.getGrad(e), m_rhs.getGrad(e)); }
	inline DiffReturnType getD       (const size_t dim, const size_t e) const { return DiffReturnType(getGrad(e), dim); }
	inline double         getCoef    (const size_t globId)              const { return m_coefBinaryOp.getCoef(globId); }

	inline bool containsTrial() const { return m_lhs.containsTrial() or m_rhs.containsTrial(); }
	inline bool containsTest()  const { return m_lhs.containsTest() or m_rhs.containsTest(); }

	inline const FSpaceType* getFunctionSpace() const { return m_coefBinaryOp.getFunctionSpace(); }
	inline const Mesh*       getMesh()          const { return getFunctionSpace()->getMesh(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
	CoefBinaryOp<Op, LeftExpr, RightExpr> m_coefBinaryOp;
};

////////////////////////////////////////////////////////////////////////

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
class CpxFiniteElementFunctionBinaryExpression : public CpxFiniteElementFunctionExpression< BinaryOpType<Op, LeftType, RightType>::Type, CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
public:
	typedef typename Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::ReturnType ReturnType;
	typedef typename Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::GradReturnType GradReturnType;
	typedef typename Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr > >::DiffReturnType DiffReturnType;
	
	typedef typename FSpace< BinaryOpType<Op, LeftType, RightType>::Type >::Type FSpaceType;
public:
	CpxFiniteElementFunctionBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)),
		m_coefBinaryOp(std::forward<LeftExpr>(lhs), std::forward<RightExpr>(rhs)) {}
public:
	inline ReturnType             operator[] (const size_t e)                   const { return ReturnType(m_lhs[e], m_rhs[e]); }
	inline GradReturnType         getGrad    (const size_t e)                   const { return GradReturnType(m_lhs.getGrad(e), m_rhs.getGrad(e)); }
	inline DiffReturnType         getD       (const size_t dim, const size_t e) const { return DiffReturnType(getGrad(e), dim); }
	inline std::complex< double > getCoef    (const size_t globId)              const { return m_coefBinaryOp.getCoef(globId); }

	inline bool containsTrial() const { return m_lhs.containsTrial() or m_rhs.containsTrial(); }
	inline bool containsTest()  const { return m_lhs.containsTest() or m_rhs.containsTest(); }

	inline const FSpaceType* getFunctionSpace() const { return m_coefBinaryOp.getFunctionSpace(); }
	inline const Mesh*       getMesh()          const { return getFunctionSpace()->getMesh(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
	CpxCoefBinaryOp<Op, LeftExpr, RightExpr> m_coefBinaryOp;
};

#endif // FINITE_ELEMENT_FUNCTION_BINARY_EXPRESSION_HPP
