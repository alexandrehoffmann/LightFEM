/*
 * Expression.hpp
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

#ifndef EXPRESSION_INCLUDE_HPP
#define EXPRESSION_INCLUDE_HPP

#include <LightFEM/Expression/Traits.hpp>

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/ExprType.hpp>
#include <LightFEM/Expression/LinAlg/Inner.hpp>
#include <LightFEM/Expression/LinAlg/Norm.hpp>
#include <LightFEM/Expression/LinAlg/MatrixBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/MatrixExpression.hpp>
#include <LightFEM/Expression/LinAlg/Matrix.hpp>
#include <LightFEM/Expression/LinAlg/MatrixUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/operators.hpp>
#include <LightFEM/Expression/LinAlg/Outer.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensorBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensorExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensor.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensorUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorExpression.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensor.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/ScalarBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>
#include <LightFEM/Expression/LinAlg/Scalar.hpp>
#include <LightFEM/Expression/LinAlg/ScalarUnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/VectorBinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/VectorExpression.hpp>
#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Expression/LinAlg/VectorUnaryExpression.hpp>

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionBinaryExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunction.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionUnaryExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseConst.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseDiff.hpp>
#include <LightFEM/Expression/ElementWise/Traits.hpp>
#include <LightFEM/Expression/ElementWise/operators.hpp>

#include <LightFEM/Expression/Function/FunctionBinaryExpression.hpp>
#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/Function/Function.hpp>
#include <LightFEM/Expression/Function/FunctionUnaryExpression.hpp>
#include <LightFEM/Expression/Function/Const.hpp>
#include <LightFEM/Expression/Function/Traits.hpp>
#include <LightFEM/Expression/Function/operators.hpp>

#include <LightFEM/Expression/FiniteElementFunction/CoefBinaryOp.hpp>
#include <LightFEM/Expression/FiniteElementFunction/CoefUnaryOP.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionBinaryExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunction.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionUnaryExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionDiffExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/Traits.hpp>
#include <LightFEM/Expression/FiniteElementFunction/operators.hpp>

#endif // EXPRESSION_INCLUDE_HPP
