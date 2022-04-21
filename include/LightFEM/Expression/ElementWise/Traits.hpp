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

#ifndef ELEMENT_WISE_TRAITS_HPP
#define ELEMENT_WISE_TRAITS_HPP

#include <type_traits>

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

class Scalar;
class Vector;
class Matrix;
class RankTwoTensor;
class RankFourTensor;
template<ExprType Type> class ElementWiseFunction;
class ElementWiseConst;
template<ExprType Type, typename Expr> class ElementWiseDiff;

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class ElementWiseFunctionBinaryExpression;
template<UnaryOp Op, ExprType Type, typename Expr> class ElementWiseFunctionUnaryExpression;

template<> struct Traits< ElementWiseFunction<ExprType::SCALAR> >     { typedef Scalar         ValueType; typedef const Scalar&         ReturnType; };
template<> struct Traits< ElementWiseFunction<ExprType::VECTOR> >     { typedef Vector         ValueType; typedef const Vector&         ReturnType; };
template<> struct Traits< ElementWiseFunction<ExprType::MATRIX> >     { typedef Matrix         ValueType; typedef const Matrix&         ReturnType; };
template<> struct Traits< ElementWiseFunction<ExprType::RK2_TENSOR> > { typedef RankTwoTensor  ValueType; typedef const RankTwoTensor&  ReturnType; };
template<> struct Traits< ElementWiseFunction<ExprType::RK4_TENSOR> > {	typedef RankFourTensor ValueType; typedef const RankFourTensor& ReturnType; };

template<> struct Traits< ElementWiseConst > { typedef const Scalar& ReturnType; };

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
struct Traits< ElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef BinaryExpression<Op, LeftType, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
};

template<UnaryOp Op, ExprType Type, typename Expr>
struct Traits< ElementWiseFunctionUnaryExpression<Op, Type, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef UnaryExpression<Op, Type, Expr_ReturnType> ReturnType;
};

template<typename Expr> struct Traits< ElementWiseDiff<ExprType::SCALAR, Expr> > 
{ 
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const Vector&, ExprType::VECTOR, Expr_ReturnType> ReturnType; 
};

template<typename Expr> struct Traits< ElementWiseDiff<ExprType::VECTOR, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr_ReturnType> TransposeType;
	typedef BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, TransposeType, ExprType::VECTOR, const Vector&> ReturnType;
};

////////////////////////////////////////////////////////////////////////

class CpxScalar;
class CpxVector;
class CpxMatrix;
class CpxRankTwoTensor;
class CpxRankFourTensor;
template<ExprType Type> class CpxElementWiseFunction;
class CpxElementWiseConst;
template<ExprType Type, typename Expr> class CpxElementWiseDiff;

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class CpxElementWiseFunctionBinaryExpression;
template<UnaryOp Op, ExprType Type, typename Expr> class CpxElementWiseFunctionUnaryExpression;

template<> struct Traits< CpxElementWiseFunction<ExprType::SCALAR> >     { typedef CpxScalar         ValueType; typedef const CpxScalar&         ReturnType; };
template<> struct Traits< CpxElementWiseFunction<ExprType::VECTOR> >     { typedef CpxVector         ValueType; typedef const CpxVector&         ReturnType; };
template<> struct Traits< CpxElementWiseFunction<ExprType::MATRIX> >     { typedef CpxMatrix         ValueType; typedef const CpxMatrix&         ReturnType; };
template<> struct Traits< CpxElementWiseFunction<ExprType::RK2_TENSOR> > { typedef CpxRankTwoTensor  ValueType; typedef const CpxRankTwoTensor&  ReturnType; };
template<> struct Traits< CpxElementWiseFunction<ExprType::RK4_TENSOR> > { typedef CpxRankFourTensor ValueType; typedef const CpxRankFourTensor& ReturnType; };

template<> struct Traits< CpxElementWiseConst > { typedef const CpxScalar& ReturnType; };

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr>
struct Traits< CpxElementWiseFunctionBinaryExpression<Op, LeftType, LeftExpr, RightType, RightExpr> >
{
	typedef typename Traits< typename std::decay<LeftExpr>::type >::ReturnType  LeftExpr_ReturnType;
	typedef typename Traits< typename std::decay<RightExpr>::type >::ReturnType RightExpr_ReturnType;
	typedef CpxBinaryExpression<Op, LeftType, LeftExpr_ReturnType, RightType, RightExpr_ReturnType> ReturnType;
};

template<UnaryOp Op, ExprType Type, typename Expr>
struct Traits< CpxElementWiseFunctionUnaryExpression<Op, Type, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef CpxUnaryExpression<Op, Type, Expr_ReturnType> ReturnType;
};

template<typename Expr> struct Traits< CpxElementWiseDiff<ExprType::SCALAR, Expr> > 
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const Vector&, ExprType::VECTOR, Expr_ReturnType> ReturnType; 
};

template<typename Expr> struct Traits< CpxElementWiseDiff<ExprType::VECTOR, Expr> >
{
	typedef typename Traits< typename std::decay<Expr>::type >::ReturnType  Expr_ReturnType;
	typedef CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr_ReturnType> TransposeType;
	typedef CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, TransposeType, ExprType::VECTOR, const Vector&> ReturnType;
};

#endif // ELEMENT_WISE_TRAITS_HPP
