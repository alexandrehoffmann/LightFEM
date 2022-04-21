/*
 * RankTwoTensorUnaryExpression.hpp
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

#ifndef RANK_TWO_TENSOR_UNARY_EXPRESSION_HPP
#define RANK_TWO_TENSOR_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorExpression.hpp>

template<> struct UnaryOpType<UnaryOp::MINUS,     ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::CONJ,      ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::ADJOINT,   ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };

template<typename Expr>
class UnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr> : public RankTwoTensorExpression< UnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 2> getShape() const { return m_expr.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return -m_expr(i,j); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr> : public RankTwoTensorExpression< UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 2> getShape() const { return m_expr.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_expr(j,i); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

////////////////////////////////////////////////////////////////////////

template<typename Expr>
class CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr> : public CpxRankTwoTensorExpression< CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 2> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return -m_expr(i,j); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK2_TENSOR, Expr> : public CpxRankTwoTensorExpression< CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK2_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 2> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return std::conj(m_expr(i,j)); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr> : public CpxRankTwoTensorExpression< CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 2> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_expr(j,i); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK2_TENSOR, Expr> : public CpxRankTwoTensorExpression< CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK2_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 2> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return std::conj(m_expr(j,i)); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

#endif // RANK_TWO_TENSOR_UNARY_EXPRESSION_HPP
