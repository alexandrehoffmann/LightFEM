/*
 * UnaryExpression.hpp
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

#ifndef UNARY_EXPRESSION_HPP
#define UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/ExprType.hpp>

enum UnaryOp
{
	MINUS,
	CONJ,
	NORM,
	TRANSPOSE,
	ADJOINT
};

template<UnaryOp Op, ExprType Type, typename Expr> class UnaryExpression {};
template<UnaryOp Op, ExprType Type, typename Expr> class CpxUnaryExpression {};

template<UnaryOp Op, ExprType Type> struct UnaryOpType { static constexpr bool isDefined = false; };

#include <LightFEM/Expression/LinAlg/RankFourTensorUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/MatrixUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/VectorUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/ScalarUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/Norm.hpp>

#endif // UNARY_EXPRESSION_HPP
