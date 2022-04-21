/*
 * ScalarUnaryExpression.hpp
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

#ifndef SCALAR_UNARY_EXPRESSION_HPP
#define SCALAR_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>

#include <cmath>

template<> struct UnaryOpType<UnaryOp::MINUS, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::CONJ,  ExprType::SCALAR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };

template <typename Expr>
class UnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr> : public ScalarExpression< UnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline double eval() const { return -m_expr.eval(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
};

////////////////////////////////////////////////////////////////////////

template <typename Expr>
class CpxUnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr> : public CpxScalarExpression< CpxUnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::complex< double > eval() const { return -m_expr.eval(); }
	inline double real() const { return -m_expr.real(); }
	inline double imag() const { return -m_expr.imag(); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::CONJ, ExprType::SCALAR, Expr> : public CpxScalarExpression< CpxUnaryExpression<UnaryOp::CONJ, ExprType::SCALAR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::complex< double > eval() const { return std::conj(m_expr.eval()); }
	inline double real() const { return  m_expr.real(); }
	inline double imag() const { return -m_expr.imag(); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

#endif // SCALAR_UNARY_EXPRESSION_HPP
