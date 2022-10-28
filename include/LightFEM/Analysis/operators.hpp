/*
 * operators.hpp
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

#ifndef ANALYSIS_OPERATORS_HPP
#define ANALYSIS_OPERATORS_HPP

#include <LightFEM/Analysis/TrialFunction.hpp>
#include <LightFEM/Analysis/TestFunction.hpp>
#include <LightFEM/Analysis/VectorTrialFunction.hpp>
#include <LightFEM/Analysis/VectorTestFunction.hpp>

#include <LightFEM/Expression/ElementWise/ElementWiseDiff.hpp>
#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/Function/Function.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunction.hpp>

typedef ElementWiseDiff<ExprType::SCALAR, GradTrialFunction> DiffTrialFunction;
typedef ElementWiseDiff<ExprType::SCALAR, GradTestFunction> DiffTestFunction;

typedef ElementWiseDiff<ExprType::VECTOR, GradVectorTrialFunction> DiffVectorTrialFunction;
typedef ElementWiseDiff<ExprType::VECTOR, GradVectorTestFunction> DiffVectorTestFunction;

inline DiffTrialFunction dx  (const TrialFunction& u) { return DiffTrialFunction(GradTrialFunction(u), 0); }
inline DiffTrialFunction dy  (const TrialFunction& u) { return DiffTrialFunction(GradTrialFunction(u), 1); }
inline GradTrialFunction grad(const TrialFunction& u) { return GradTrialFunction(u);    }

inline DiffTestFunction dx  (const TestFunction& u) { return DiffTestFunction(GradTestFunction(u), 0); }
inline DiffTestFunction dy  (const TestFunction& u) { return DiffTestFunction(GradTestFunction(u), 1); }
inline GradTestFunction grad(const TestFunction& u) { return GradTestFunction(u);    }

inline DiffVectorTrialFunction dx  (const VectorTrialFunction& u) { return DiffVectorTrialFunction(GradVectorTrialFunction(u), 0); }
inline DiffVectorTrialFunction dy  (const VectorTrialFunction& u) { return DiffVectorTrialFunction(GradVectorTrialFunction(u), 1); }
inline GradVectorTrialFunction grad(const VectorTrialFunction& u) { return GradVectorTrialFunction(u);    }

inline DiffVectorTestFunction dx  (const VectorTestFunction& u) { return DiffVectorTestFunction(GradVectorTestFunction(u), 0); }
inline DiffVectorTestFunction dy  (const VectorTestFunction& u) { return DiffVectorTestFunction(GradVectorTestFunction(u), 1); }
inline GradVectorTestFunction grad(const VectorTestFunction& u) { return GradVectorTestFunction(u);    }

////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr> 
FiniteElementFunctionDiffExpression<Type, Expr> dx(const FiniteElementFunctionExpression<Type, Expr>& expr) { return FiniteElementFunctionDiffExpression<Type, Expr>(expr, 0); }

template<ExprType Type, typename Expr> 
FiniteElementFunctionDiffExpression<Type, Expr> dx(FiniteElementFunctionExpression<Type, Expr>&& expr) { return FiniteElementFunctionDiffExpression<Type, Expr>(expr, 0); }

template<ExprType Type, typename Expr> 
CpxFiniteElementFunctionDiffExpression<Type, Expr> dx(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionDiffExpression<Type, Expr>(expr, 0); }

template<ExprType Type, typename Expr> 
CpxFiniteElementFunctionDiffExpression<Type, Expr> dx(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionDiffExpression<Type, Expr>(expr, 0); }

template<ExprType Type, typename Expr> 
FiniteElementFunctionDiffExpression<Type, Expr> dy(const FiniteElementFunctionExpression<Type, Expr>& expr) { return FiniteElementFunctionDiffExpression<Type, Expr>(expr, 1); }

template<ExprType Type, typename Expr> 
FiniteElementFunctionDiffExpression<Type, Expr> dy(FiniteElementFunctionExpression<Type, Expr>&& expr) { return FiniteElementFunctionDiffExpression<Type, Expr>(expr, 1); }

template<ExprType Type, typename Expr> 
CpxFiniteElementFunctionDiffExpression<Type, Expr> dy(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionDiffExpression<Type, Expr>(expr, 1); }

template<ExprType Type, typename Expr> 
CpxFiniteElementFunctionDiffExpression<Type, Expr> dy(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionDiffExpression<Type, Expr>(expr, 1); }

template<ExprType Type, typename Expr> 
FiniteElementFunctionGradExpression<Type, Expr> grad(const FiniteElementFunctionExpression<Type, Expr>& expr) { return FiniteElementFunctionGradExpression<Type, Expr>(expr); }

template<ExprType Type, typename Expr> 
FiniteElementFunctionGradExpression<Type, Expr> grad(FiniteElementFunctionExpression<Type, Expr>&& expr) { return FiniteElementFunctionGradExpression<Type, Expr>(expr); }

template<ExprType Type, typename Expr> 
CpxFiniteElementFunctionGradExpression<Type, Expr> grad(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionGradExpression<Type, Expr>(expr); }

template<ExprType Type, typename Expr> 
CpxFiniteElementFunctionGradExpression<Type, Expr> grad(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionGradExpression<Type, Expr>(expr); }

////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr> 
FiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, const FunctionExpression<Type, Expr>& expr);

template<ExprType Type>
FiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, typename Functor<Type>::type expr);

template<ExprType Type, typename Expr>
CpxFiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, const CpxFunctionExpression<Type, Expr>& expr);

template<ExprType Type>
CpxFiniteElementFunction<Type> interp(const typename FiniteElementFunction<Type>::FSpaceType* Vh, typename CpxFunctor<Type>::type expr);

#include <LightFEM/Analysis/operators.tpp>

#endif //ANALYSIS_OPERATORS_HPP
