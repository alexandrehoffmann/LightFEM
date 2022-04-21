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

#ifndef FINITE_ELEMENT_FUNCTION_TRAITS_HPP
#define FINITE_ELEMENT_FUNCTION_TRAITS_HPP

#include <type_traits>

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

template<ExprType Type> class ElementWiseFunction;
template<ExprType Type> class FiniteElementFunction;
template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class FiniteElementFunctionBinaryExpression;
template<UnaryOp Op, ExprType Type, typename Expr> class FiniteElementFunctionUnaryExpression;
template<ExprType Type, typename Expr> class ElementWiseDiff;
template<ExprType Type, typename Expr> class FiniteElementFunctionDiffExpression;
template<ExprType Type, typename Expr> class FiniteElementFunctionGradExpression;

template<> struct Traits< FiniteElementFunction< ExprType::SCALAR > > 
{ 
	typedef const ElementWiseFunction< ExprType::SCALAR >& ReturnType; 
	typedef const ElementWiseFunction< ExprType::VECTOR >& GradReturnType;
	typedef ElementWiseDiff<ExprType::SCALAR, ElementWiseFunction< ExprType::VECTOR > > DiffReturnType;
	
	typedef ElementWiseFunction< ExprType::SCALAR > ValueType; 
	typedef ElementWiseFunction< ExprType::VECTOR > GradValueType;
};

template<> struct Traits< FiniteElementFunction< ExprType::VECTOR > > 
{ 
	typedef const ElementWiseFunction< ExprType::VECTOR >& ReturnType; 
	typedef const ElementWiseFunction< ExprType::MATRIX >& GradReturnType;
	typedef ElementWiseDiff<ExprType::VECTOR, ElementWiseFunction< ExprType::MATRIX > > DiffReturnType;
	
	typedef ElementWiseFunction< ExprType::VECTOR > ValueType; 
	typedef ElementWiseFunction< ExprType::MATRIX > GradValueType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
struct Traits< FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType RightType, typename RightExpr>
struct Traits< FiniteElementFunctionBinaryExpression<Op, ExprType::SCALAR, const Const&, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<const Const&>::type >::ReturnType  LeftExpr_ReturnType;

	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef ElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType RightType, typename RightExpr>
struct Traits< FiniteElementFunctionBinaryExpression<Op, ExprType::SCALAR, Const, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<Const>::type >::ReturnType  LeftExpr_ReturnType;

	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef ElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr>
struct Traits< FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, ExprType::SCALAR, const Const&> >
{
private:
	typedef typename Traits< typename std::decay<const Const&>::type >::ReturnType  RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
public:
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, ExprType::SCALAR, LeftExpr_ReturnType> ReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, ExprType::SCALAR, LeftExpr_ReturnType> GradReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, ExprType::SCALAR, LeftExpr_ReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr>
struct Traits< FiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, ExprType::SCALAR, Const> >
{
private:
	typedef typename Traits< typename std::decay<Const>::type >::ReturnType  RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
public:
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, ExprType::SCALAR, LeftExpr_ReturnType> ReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, ExprType::SCALAR, LeftExpr_ReturnType> GradReturnType;
	typedef ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, ExprType::SCALAR, LeftExpr_ReturnType> DiffReturnType;
};

template<UnaryOp Op, ExprType Type, typename Expr>
struct Traits< FiniteElementFunctionUnaryExpression<Op, Type, Expr> >
{
private:
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef typename Traits< typename std::decay<Expr>::type >::GradReturnType  Expr_GradReturnType;
	typedef typename Traits< typename std::decay<Expr>::type >::DiffReturnType  Expr_DiffReturnType;
public:
	typedef ElementWiseFunctionUnaryExpression<Op, Type, Expr_ReturnType> ReturnType;
	typedef ElementWiseFunctionUnaryExpression<Op, Type, Expr_GradReturnType> GradReturnType;
	typedef ElementWiseFunctionUnaryExpression<Op, Type, Expr_DiffReturnType> DiffReturnType;
};

template<ExprType Type, typename Expr>
struct Traits< FiniteElementFunctionDiffExpression<Type, Expr> >
{
	typedef typename Traits< Expr >::DiffReturnType ReturnType;
};

template<ExprType Type, typename Expr>
struct Traits< FiniteElementFunctionGradExpression<Type, Expr> >
{
	typedef typename Traits< Expr >::GradReturnType ReturnType;
};

////////////////////////////////////////////////////////////////////////

template<ExprType Type> class CpxElementWiseFunction;
template<ExprType Type> class CpxFiniteElementFunction;
template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class CpxFiniteElementFunctionBinaryExpression;
template<UnaryOp Op, ExprType Type, typename Expr> class CpxFiniteElementFunctionUnaryExpression;
template<ExprType Type, typename Expr> class CpxElementWiseDiff;
template<ExprType Type, typename Expr> class CpxFiniteElementFunctionDiffExpression;
template<ExprType Type, typename Expr> class CpxFiniteElementFunctionGradExpression;

template<> struct Traits< CpxFiniteElementFunction< ExprType::SCALAR > > 
{ 
	typedef const CpxElementWiseFunction< ExprType::SCALAR >& ReturnType; 
	typedef const CpxElementWiseFunction< ExprType::VECTOR >& GradReturnType;
	typedef CpxElementWiseDiff<ExprType::SCALAR, ElementWiseFunction< ExprType::VECTOR > > DiffReturnType;
	
	typedef CpxElementWiseFunction< ExprType::SCALAR > ValueType; 
	typedef CpxElementWiseFunction< ExprType::VECTOR > GradValueType;
};

template<> struct Traits< CpxFiniteElementFunction< ExprType::VECTOR > > 
{ 
	typedef const CpxElementWiseFunction< ExprType::VECTOR >& ReturnType; 
	typedef const CpxElementWiseFunction< ExprType::MATRIX >& GradReturnType;
	typedef CpxElementWiseDiff<ExprType::VECTOR, ElementWiseFunction< ExprType::MATRIX > > DiffReturnType;
	
	typedef CpxElementWiseFunction< ExprType::VECTOR > ValueType; 
	typedef CpxElementWiseFunction< ExprType::MATRIX > GradValueType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType RightType, typename RightExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, ExprType::SCALAR, const Const&, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<const Const&>::type >::ReturnType  LeftExpr_ReturnType;

	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType RightType, typename RightExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, ExprType::SCALAR, Const, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<Const>::type >::ReturnType  LeftExpr_ReturnType;

	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType RightType, typename RightExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, ExprType::SCALAR, const CpxConst&, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<const CpxConst&>::type >::ReturnType  LeftExpr_ReturnType;

	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType RightType, typename RightExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, ExprType::SCALAR, CpxConst, RightType, RightExpr> >
{
private:
	typedef typename Traits< typename std::decay<CpxConst>::type >::ReturnType  LeftExpr_ReturnType;

	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::GradReturnType RightExpr_GradReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::DiffReturnType RightExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_GradReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, ExprType::SCALAR, LeftExpr_ReturnType, RightType, RightExpr_DiffReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, ExprType::SCALAR, const Const&> >
{
private:
	typedef typename Traits< typename std::decay<const Const&>::type >::ReturnType  RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, ExprType::SCALAR, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, ExprType::SCALAR, RightExpr_ReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, ExprType::SCALAR, RightExpr_ReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, ExprType::SCALAR, Const> >
{
private:
	typedef typename Traits< typename std::decay<Const>::type >::ReturnType  RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, ExprType::SCALAR, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, ExprType::SCALAR, RightExpr_ReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, ExprType::SCALAR, RightExpr_ReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&> >
{
private:
	typedef typename Traits< typename std::decay<const CpxConst&>::type >::ReturnType  RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, ExprType::SCALAR, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, ExprType::SCALAR, RightExpr_ReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, ExprType::SCALAR, RightExpr_ReturnType> DiffReturnType;
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr>
struct Traits< CpxFiniteElementFunctionBinaryExpression<Op, LeftType, LeftExpr, ExprType::SCALAR, CpxConst> >
{
private:
	typedef typename Traits< typename std::decay<CpxConst>::type >::ReturnType  RightExpr_ReturnType;

	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::GradReturnType  LeftExpr_GradReturnType;
	typedef typename Traits< typename std::decay<LeftExpr>::type >::DiffReturnType  LeftExpr_DiffReturnType;
public:
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_ReturnType, ExprType::SCALAR, RightExpr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_GradReturnType, ExprType::SCALAR, RightExpr_ReturnType> GradReturnType;
	typedef CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr_DiffReturnType, ExprType::SCALAR, RightExpr_ReturnType> DiffReturnType;
};

template<UnaryOp Op, ExprType Type, typename Expr>
struct Traits< CpxFiniteElementFunctionUnaryExpression<Op, Type, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef typename Traits< typename std::decay<Expr>::type >::GradReturnType  Expr_GradReturnType;
	typedef typename Traits< typename std::decay<Expr>::type >::DiffReturnType  Expr_DiffReturnType;
	
	typedef CpxElementWiseFunctionUnaryExpression<Op, Type, Expr_ReturnType> ReturnType;
	typedef CpxElementWiseFunctionUnaryExpression<Op, Type, Expr_GradReturnType> GradReturnType;
	typedef CpxElementWiseFunctionUnaryExpression<Op, Type, Expr_DiffReturnType> DiffReturnType;
};

template<ExprType Type, typename Expr>
struct Traits< CpxFiniteElementFunctionDiffExpression<Type, Expr> >
{
	typedef typename Traits< Expr >::DiffReturnType ReturnType;
};

template<ExprType Type, typename Expr>
struct Traits< CpxFiniteElementFunctionGradExpression<Type, Expr> >
{
	typedef typename Traits< Expr >::GradReturnType ReturnType;
};


#endif // FINITE_ELEMENT_FUNCTION_TRAITS_HPP
