/*
 * Traits.hpp
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

#ifndef FUNCTION_TRAITS_HPP
#define FUNCTION_TRAITS_HPP

#include <type_traits>

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

template<ExprType Type> class ElementWiseFunction; 
template<ExprType Type> class Function; 
template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class FunctionBinaryExpression;
template<UnaryOp Op, ExprType Type, typename Expr> class FunctionUnaryExpression;

class ElementWiseConst;
class Const;

template< ExprType Type > struct Traits< Function<Type> > { typedef const ElementWiseFunction< Type >& ReturnType; };

template<> struct Traits< Const > { typedef const ElementWiseConst ReturnType; };

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
struct Traits< FunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
};

template<UnaryOp Op, ExprType Type, typename Expr>
struct Traits< FunctionUnaryExpression<Op, Type, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef ElementWiseFunctionUnaryExpression<Op, Type, Expr_ReturnType> ReturnType;
};

////////////////////////////////////////////////////////////////////////

template<ExprType Type> class CpxElementWiseFunction; 
template<ExprType Type> class CpxFunction; 
template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class CpxFunctionBinaryExpression;
template<UnaryOp Op, ExprType Type, typename Expr> class CpxFunctionUnaryExpression;

class CpxElementWiseConst;
class CpxConst;

template< ExprType Type > struct Traits< CpxFunction<Type> > { typedef const CpxElementWiseFunction< Type >& ReturnType; };

template<> struct Traits< CpxConst > { typedef const CpxElementWiseConst ReturnType; };

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
struct Traits< CpxFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
};

template<UnaryOp Op, ExprType Type, typename Expr>
struct Traits< CpxFunctionUnaryExpression<Op, Type, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef CpxElementWiseFunctionUnaryExpression<Op, Type, Expr_ReturnType> ReturnType;
};

#endif // FUNCTION_TRAITS_HPP
