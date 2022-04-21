/*
 * VectorUnaryExpression.hpp
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

#ifndef VECTOR_UNARY_EXPRESSION_HPP
#define VECTOR_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/VectorExpression.hpp>

template<> struct UnaryOpType<UnaryOp::MINUS, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::CONJ,  ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };

template<typename Expr>
class UnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr> : public VectorExpression< UnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getSize() const { return m_expr.getSize(); }
public:
	inline double operator[](const size_t i) const { return -m_expr[i]; }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 

////////////////////////////////////////////////////////////////////////

template<typename Expr>
class CpxUnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr> : public CpxVectorExpression< CpxUnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getSize() const { return m_expr.getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return -m_expr[i]; }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 

template<typename Expr>
class CpxUnaryExpression<UnaryOp::CONJ, ExprType::VECTOR, Expr> : public CpxVectorExpression< CpxUnaryExpression<UnaryOp::CONJ, ExprType::VECTOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getSize() const { return m_expr.getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return std::conj(m_expr[i]); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 

#endif // VECTOR_UNARY_EXPRESSION_HPP
