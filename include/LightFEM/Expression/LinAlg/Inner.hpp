/*
 * Inner.hpp
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

#ifndef INNER_HPP
#define INNER_HPP

#include <complex>
#include <cassert>

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>

template<> struct BinaryOpType<BinaryOp::INNER, ExprType::VECTOR, ExprType::VECTOR>         { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::INNER, ExprType::MATRIX, ExprType::MATRIX>         { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::INNER, ExprType::RK2_TENSOR, ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline std::complex< double > eval() const { return m_value; }
private:
	std::complex< double > m_value;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline std::complex< double > eval() const { return m_value; }
private:
	std::complex< double > m_value;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs);
public:
	inline std::complex< double > eval() const { return m_value; }
private:
	std::complex< double > m_value;
};

#include <LightFEM/Expression/LinAlg/Inner.tpp>

#endif // INNER_HPP
