/*
 * operators.tpp
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

#ifndef ANALYSIS_OPERATORS_TPP
#define ANALYSIS_OPERATORS_TPP

#include <LightFEM/Analysis/operators.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>
#include <LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp>

#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.hpp>

template<ExprType Type, typename Expr> 
FiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, const FunctionExpression<Type, Expr>& expr)
{
	std::vector< Scalar > coefs(Vh->getNBasisFunction());

	for (size_t e=0;e<Vh->getMesh()->getNElem();++e)
	{
		ElementWiseInterpolator interp_expr(expr[e]);
		for (size_t interpNodeId=0;interpNodeId<Vh->getNLocalInterpolationNodes();++interpNodeId)
		{
			const size_t globId = Vh->getGlobalId(e, Vh->getLocalInterpolationFunctionId(interpNodeId));
			coefs[globId] = interp_expr(Vh->getLocalInterpolationNode(interpNodeId));
		}
	}
	return FiniteElementFunction<Type>(Vh, coefs.data());
}

template<ExprType Type>
FiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, typename Functor<Type>::type expr)
{
	std::vector< Scalar > coefs(Vh->getNBasisFunction());

	for (size_t e=0;e<Vh->getMesh()->getNElem();++e)
	{
		for (size_t interpNodeId=0;interpNodeId<Vh->getNLocalInterpolationNodes();++interpNodeId)
		{
			const size_t globId = Vh->getGlobalId(e, Vh->getLocalInterpolationFunctionId(interpNodeId));
			const NodeRef Xi = Vh->getLocalInterpolationNode(interpNodeId);
			const auto [x, y] = Vh->getMesh()->getElem(e)->getXworld(Xi);
			coefs[globId] = expr(x, y);
		}
	}
	return FiniteElementFunction<Type>(Vh, coefs.data());
}

template<ExprType Type, typename Expr>
CpxFiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, const CpxFunctionExpression<Type, Expr>& expr)
{
	std::vector< CpxScalar > coefs(Vh->getNBasisFunction());

	for (size_t e=0;e<Vh->getMesh()->getNElem();++e)
	{
		CpxElementWiseInterpolator interp_expr(expr[e]);
		for (size_t interpNodeId=0;interpNodeId<Vh->getNLocalInterpolationNodes();++interpNodeId)
		{
			const size_t globId = Vh->getGlobalId(e, Vh->getLocalInterpolationFunctionId(interpNodeId));
			coefs[globId] = interp_expr(Vh->getLocalInterpolationNode(interpNodeId));
		}
	}
	return CpxFiniteElementFunction<Type>(Vh, coefs.data());
}

template<ExprType Type>
FiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, typename CpxFunctor<Type>::type expr)
{
	std::vector< CpxScalar > coefs(Vh->getNBasisFunction());

	for (size_t e=0;e<Vh->getMesh()->getNElem();++e)
	{
		for (size_t interpNodeId=0;interpNodeId<Vh->getNLocalInterpolationNodes();++interpNodeId)
		{
			const size_t globId = Vh->getGlobalId(e, Vh->getLocalInterpolationFunctionId(interpNodeId));
			const NodeRef Xi = Vh->getLocalInterpolationNode(interpNodeId);
			const auto [x, y] = Vh->getMesh()->getElem(e)->getXworld(Xi);
			coefs[globId] = expr(x, y);
		}
	}
	return CpxFiniteElementFunction<Type>(Vh, coefs.data());
}

#endif //ANALYSIS_OPERATORS_TPP
