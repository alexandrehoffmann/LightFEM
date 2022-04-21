/*
 * MatrixUnaryExpression.hpp
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

#ifndef MATRIX_UNARY_EXPRESSION_HPP
#define MATRIX_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/MatrixExpression.hpp>

template<> struct UnaryOpType<UnaryOp::MINUS,     ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::TRANSPOSE, ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::CONJ,      ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::ADJOINT,   ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };

template<typename Expr>
class UnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr> : public MatrixExpression< UnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getNrows() const { return m_expr.getNcols(); }
	inline size_t getNcols() const { return m_expr.getNrows(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return -m_expr(i,j); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 

template<typename Expr>
class UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr> : public MatrixExpression< UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getNrows() const { return m_expr.getNcols(); }
	inline size_t getNcols() const { return m_expr.getNrows(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_expr(j,i); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 

////////////////////////////////////////////////////////////////////////

template<typename Expr>
class CpxUnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr> : public CpxMatrixExpression< CpxUnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getNrows() const { return m_expr.getNcols(); }
	inline size_t getNcols() const { return m_expr.getNrows(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return -m_expr(i,j); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 

template<typename Expr>
class CpxUnaryExpression<UnaryOp::CONJ, ExprType::MATRIX, Expr> : public CpxMatrixExpression< CpxUnaryExpression<UnaryOp::CONJ, ExprType::MATRIX, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getNrows() const { return m_expr.getNrows(); }
	inline size_t getNcols() const { return m_expr.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return std::conj(m_expr(i,j)); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};
 

template<typename Expr>
class CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr> : public CpxMatrixExpression< CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getNrows() const { return m_expr.getNcols(); }
	inline size_t getNcols() const { return m_expr.getNrows(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_expr(j,i); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
}; 
 
template<typename Expr>
class CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::MATRIX, Expr> : public CpxMatrixExpression< CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::MATRIX, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline size_t getNrows() const { return m_expr.getNcols(); }
	inline size_t getNcols() const { return m_expr.getNrows(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return std::conj(m_expr(j,i)); }
private:
	typename RefTypeSelector<Expr>::Type m_expr;
}; 

#endif // MATRIX_UNARY_EXPRESSION_HPP
