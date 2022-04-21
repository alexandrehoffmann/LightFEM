/*
 * BinaryExpression.hpp
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

#ifndef BINARY_EXPRESSION_HPP
#define BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/ExprType.hpp>

enum class BinaryOp
{
	SUM,
	SUB,
	PROD,
	DIV,
	INNER,
	OUTER,
	CROSS,
	DDOT 
};

template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class BinaryExpression {};
template<BinaryOp Op, ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr> class CpxBinaryExpression {};

template<BinaryOp Op, ExprType LeftType, ExprType RightType> struct BinaryOpType { static constexpr bool isDefined = false; };

#include <LightFEM/Expression/LinAlg/RankFourTensorBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/MatrixBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/VectorBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/ScalarBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/Inner.hpp>
#include <LightFEM/Expression/LinAlg/Outer.hpp>
#include <LightFEM/Expression/LinAlg/Cross.hpp>

#endif // BINARY_EXPRESSION_HPP
