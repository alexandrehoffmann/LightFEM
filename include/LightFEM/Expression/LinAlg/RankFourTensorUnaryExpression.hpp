/*
 * RankFourTensorUnaryExpression.hpp
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

#ifndef RANK_FOUR_TENSOR_UNARY_EXPRESSION_HPP
#define RANK_FOUR_TENSOR_UNARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensorExpression.hpp>

template<> struct UnaryOpType<UnaryOp::MINUS,     ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::CONJ,      ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::ADJOINT,   ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };

template<typename Expr>
class UnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr> : public RankFourTensorExpression< UnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 4> getShape() const { return m_expr.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return -m_expr(i,j,k,l); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr> : public RankFourTensorExpression< UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 4> getShape() const { return m_expr.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_expr(k,l,i,j); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

////////////////////////////////////////////////////////////////////////

template<typename Expr>
class CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr> : public CpxRankFourTensorExpression< CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 4> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return -m_expr(i,j,k,l); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK4_TENSOR, Expr> : public CpxRankFourTensorExpression< CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK4_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 4> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return std::conj(m_expr(i,j,k,l)); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr> : public CpxRankFourTensorExpression< CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 4> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_expr(k,l,i,j); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr>
class CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK4_TENSOR, Expr> : public CpxRankFourTensorExpression< CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK4_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) : m_expr(std::forward<Expr>(expr)) {}
public:
	inline std::array<size_t, 4> getShape() const { return m_expr.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return std::conj(m_expr(k,l,i,j)); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

#endif // RANK_FOUR_TENSOR_UNARY_EXPRESSION_HPP
