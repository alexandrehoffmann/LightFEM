/*
 * operators.hpp
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
 
#ifndef LINALG_OPERATORS_HPP
#define LINALG_OPERATORS_HPP

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

////////////////////////////////////////////////////////////////////////
////                      rk4_tensor operators                      ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, const Expr&> operator-(const RankFourTensorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr> operator-(RankFourTensorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, const Expr&> operator-(const CpxRankFourTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr> operator-(CpxRankFourTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK4_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK4_TENSOR, const Expr&> conj(const CpxRankFourTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK4_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK4_TENSOR, Expr> conj(CpxRankFourTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK4_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, const Expr&> transpose(const RankFourTensorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr> transpose(RankFourTensorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, const Expr&> transpose(const CpxRankFourTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr> transpose(CpxRankFourTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK4_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK4_TENSOR, const Expr&> adjoint(const CpxRankFourTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK4_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK4_TENSOR, Expr> adjoint(CpxRankFourTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK4_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> ddot(const RankFourTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> ddot(const RankFourTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> ddot(RankFourTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> ddot(RankFourTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> ddot(const RankFourTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> ddot(const RankFourTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> ddot(RankFourTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> ddot(RankFourTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> ddot(const RankFourTensorExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr> ddot(const RankFourTensorExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&> ddot(RankFourTensorExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr> ddot(RankFourTensorExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> ddot(const RankFourTensorExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr> ddot(const RankFourTensorExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&> ddot(RankFourTensorExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr> ddot(RankFourTensorExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr> ddot(const CpxRankFourTensorExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr> ddot(CpxRankFourTensorExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator+(const RankFourTensorExpression<LeftExpr>& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator+(const RankFourTensorExpression<LeftExpr>& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator+(RankFourTensorExpression<LeftExpr>&& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator+(RankFourTensorExpression<LeftExpr>&& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator+(const RankFourTensorExpression<LeftExpr>& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator+(const RankFourTensorExpression<LeftExpr>& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator+(RankFourTensorExpression<LeftExpr>&& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator+(RankFourTensorExpression<LeftExpr>&& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator+(const CpxRankFourTensorExpression<LeftExpr>& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator+(const CpxRankFourTensorExpression<LeftExpr>& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator+(CpxRankFourTensorExpression<LeftExpr>&& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator+(CpxRankFourTensorExpression<LeftExpr>&& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator+(const CpxRankFourTensorExpression<LeftExpr>& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator+(const CpxRankFourTensorExpression<LeftExpr>& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator+(CpxRankFourTensorExpression<LeftExpr>&& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator+(CpxRankFourTensorExpression<LeftExpr>&& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator-(const RankFourTensorExpression<LeftExpr>& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator-(const RankFourTensorExpression<LeftExpr>& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator-(RankFourTensorExpression<LeftExpr>&& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator-(RankFourTensorExpression<LeftExpr>&& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator-(const RankFourTensorExpression<LeftExpr>& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator-(const RankFourTensorExpression<LeftExpr>& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator-(RankFourTensorExpression<LeftExpr>&& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator-(RankFourTensorExpression<LeftExpr>&& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator-(const CpxRankFourTensorExpression<LeftExpr>& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator-(const CpxRankFourTensorExpression<LeftExpr>& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator-(CpxRankFourTensorExpression<LeftExpr>&& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator-(CpxRankFourTensorExpression<LeftExpr>&& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator-(const CpxRankFourTensorExpression<LeftExpr>& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator-(const CpxRankFourTensorExpression<LeftExpr>& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator-(CpxRankFourTensorExpression<LeftExpr>&& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator-(CpxRankFourTensorExpression<LeftExpr>&& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const RankFourTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const RankFourTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(RankFourTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(RankFourTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const RankFourTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const RankFourTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(RankFourTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(RankFourTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxRankFourTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxRankFourTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxRankFourTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxRankFourTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxRankFourTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxRankFourTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxRankFourTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxRankFourTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const RankFourTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const RankFourTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(RankFourTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(RankFourTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const RankFourTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const RankFourTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(RankFourTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(RankFourTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxRankFourTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxRankFourTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxRankFourTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxRankFourTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxRankFourTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxRankFourTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxRankFourTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxRankFourTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                      rk2_tensor operators                      ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, const Expr&> operator-(const RankTwoTensorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr> operator-(RankTwoTensorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, const Expr&> operator-(const CpxRankTwoTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr> operator-(CpxRankTwoTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK2_TENSOR, const Expr&> conj(const CpxRankTwoTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK2_TENSOR, Expr> conj(CpxRankTwoTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, const Expr&> norm(const RankTwoTensorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr> norm(RankTwoTensorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, const Expr&> norm(const CpxRankTwoTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr> norm(CpxRankTwoTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, const Expr&> transpose(const RankTwoTensorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr> transpose(RankTwoTensorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, const Expr&> transpose(const CpxRankTwoTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr> transpose(CpxRankTwoTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK2_TENSOR, const Expr&> adjoint(const CpxRankTwoTensorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK2_TENSOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK2_TENSOR, Expr> adjoint(CpxRankTwoTensorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::RK2_TENSOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> inner(const RankTwoTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> inner(const RankTwoTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> inner(RankTwoTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> inner(RankTwoTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> inner(const RankTwoTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> inner(const RankTwoTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> inner(RankTwoTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> inner(RankTwoTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> inner(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> inner(const CpxRankTwoTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> inner(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> inner(CpxRankTwoTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> inner(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> inner(const CpxRankTwoTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> inner(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> inner(CpxRankTwoTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator+(const RankTwoTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator+(const RankTwoTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator+(RankTwoTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator+(RankTwoTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator+(const RankTwoTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator+(const RankTwoTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator+(RankTwoTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator+(RankTwoTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator+(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator+(const CpxRankTwoTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator+(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator+(CpxRankTwoTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator+(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator+(const CpxRankTwoTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator+(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator+(CpxRankTwoTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator-(const RankTwoTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator-(const RankTwoTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator-(RankTwoTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator-(RankTwoTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator-(const RankTwoTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator-(const RankTwoTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator-(RankTwoTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator-(RankTwoTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator-(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator-(const CpxRankTwoTensorExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator-(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator-(CpxRankTwoTensorExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator-(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator-(const CpxRankTwoTensorExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator-(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator-(CpxRankTwoTensorExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const RankTwoTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const RankTwoTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(RankTwoTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(RankTwoTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const RankTwoTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const RankTwoTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(RankTwoTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(RankTwoTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxRankTwoTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxRankTwoTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxRankTwoTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxRankTwoTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const RankTwoTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const RankTwoTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(RankTwoTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(RankTwoTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const RankTwoTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const RankTwoTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(RankTwoTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(RankTwoTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxRankTwoTensorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxRankTwoTensorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxRankTwoTensorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxRankTwoTensorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxRankTwoTensorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxRankTwoTensorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        matrix operators                        ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, const Expr&> operator-(const MatrixExpression<Expr>& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr> operator-(MatrixExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, const Expr&> operator-(const CpxMatrixExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr> operator-(CpxMatrixExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::MATRIX, const Expr&> conj(const CpxMatrixExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::MATRIX, Expr> conj(CpxMatrixExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, const Expr&> norm(const MatrixExpression<Expr>& expr) { return UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr> norm(MatrixExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, const Expr&> norm(const CpxMatrixExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr> norm(CpxMatrixExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Expr&> transpose(const MatrixExpression<Expr>& expr) { return UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr> transpose(MatrixExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Expr&> transpose(const CpxMatrixExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr> transpose(CpxMatrixExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::MATRIX, const Expr&> adjoint(const CpxMatrixExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::MATRIX, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::MATRIX, Expr> adjoint(CpxMatrixExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::ADJOINT, ExprType::MATRIX, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> inner(const MatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> inner(const MatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> inner(MatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> inner(MatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> inner(const MatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> inner(const MatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> inner(MatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> inner(MatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> inner(const CpxMatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> inner(const CpxMatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> inner(CpxMatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> inner(CpxMatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> inner(const CpxMatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> inner(const CpxMatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> inner(CpxMatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> inner(CpxMatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator+(const MatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator+(const MatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator+(MatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator+(MatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator+(const MatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator+(const MatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator+(MatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator+(MatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator+(const CpxMatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator+(const CpxMatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator+(CpxMatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator+(CpxMatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator+(const CpxMatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator+(const CpxMatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator+(CpxMatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator+(CpxMatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator-(const MatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator-(const MatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator-(MatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator-(MatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator-(const MatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator-(const MatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator-(MatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator-(MatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator-(const CpxMatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator-(const CpxMatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator-(CpxMatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator-(CpxMatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator-(const CpxMatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator-(const CpxMatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator-(CpxMatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator-(CpxMatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const MatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const MatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(MatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator*(MatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const MatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const MatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(MatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator*(MatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const CpxMatrixExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const CpxMatrixExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(CpxMatrixExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator*(CpxMatrixExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const CpxMatrixExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const CpxMatrixExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(CpxMatrixExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> operator*(CpxMatrixExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const MatrixExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const MatrixExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(MatrixExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr> operator*(MatrixExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const MatrixExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const MatrixExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(MatrixExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr> operator*(MatrixExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const CpxMatrixExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const CpxMatrixExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(CpxMatrixExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr> operator*(CpxMatrixExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const CpxMatrixExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const CpxMatrixExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(CpxMatrixExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr> operator*(CpxMatrixExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const MatrixExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const MatrixExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(MatrixExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator*(MatrixExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const MatrixExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const MatrixExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(MatrixExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator*(MatrixExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxMatrixExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxMatrixExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxMatrixExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxMatrixExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxMatrixExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxMatrixExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxMatrixExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxMatrixExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const MatrixExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const MatrixExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(MatrixExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator/(MatrixExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const MatrixExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const MatrixExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(MatrixExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator/(MatrixExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxMatrixExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxMatrixExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxMatrixExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxMatrixExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxMatrixExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxMatrixExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxMatrixExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxMatrixExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        vector operators                        ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, const Expr&> operator-(const VectorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr> operator-(VectorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, const Expr&> operator-(const CpxVectorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr> operator-(CpxVectorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::VECTOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::VECTOR, const Expr&> conj(const CpxVectorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::VECTOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::VECTOR, Expr> conj(CpxVectorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::VECTOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, const Expr&> norm(const VectorExpression<Expr>& expr) { return UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr> norm(VectorExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, const Expr&> norm(const CpxVectorExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr> norm(CpxVectorExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> inner(const VectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> inner(const VectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> inner(VectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> inner(VectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> inner(const VectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> inner(const VectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> inner(VectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> inner(VectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> inner(const CpxVectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> inner(const CpxVectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> inner(CpxVectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> inner(CpxVectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> inner(const CpxVectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> inner(const CpxVectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> inner(CpxVectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> inner(CpxVectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> outer(const VectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> outer(const VectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> outer(VectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> outer(VectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> outer(const VectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> outer(const VectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> outer(VectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> outer(VectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> outer(const CpxVectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> outer(const CpxVectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> outer(CpxVectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> outer(CpxVectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> outer(const CpxVectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> outer(const CpxVectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> outer(CpxVectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> outer(CpxVectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::OUTER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> cross(const VectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> cross(const VectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> cross(VectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> cross(VectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> cross(const VectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> cross(const VectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> cross(VectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> cross(VectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> cross(const CpxVectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> cross(const CpxVectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> cross(CpxVectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> cross(CpxVectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> cross(const CpxVectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> cross(const CpxVectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> cross(CpxVectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> cross(CpxVectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator+(const VectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator+(const VectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator+(VectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator+(VectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator+(const VectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator+(const VectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator+(VectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator+(VectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator+(const CpxVectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator+(const CpxVectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator+(CpxVectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator+(CpxVectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator+(const CpxVectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator+(const CpxVectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator+(CpxVectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator+(CpxVectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator-(const VectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator-(const VectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator-(VectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator-(VectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator-(const VectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator-(const VectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator-(VectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator-(VectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator-(const CpxVectorExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator-(const CpxVectorExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator-(CpxVectorExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator-(CpxVectorExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator-(const CpxVectorExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator-(const CpxVectorExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator-(CpxVectorExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> operator-(CpxVectorExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const VectorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const VectorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(VectorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(VectorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const VectorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const VectorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(VectorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(VectorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxVectorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxVectorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxVectorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxVectorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxVectorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxVectorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxVectorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxVectorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const VectorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const VectorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(VectorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(VectorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const VectorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const VectorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(VectorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(VectorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxVectorExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxVectorExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxVectorExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxVectorExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxVectorExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxVectorExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxVectorExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxVectorExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        scalar operators                        ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, const Expr&> operator-(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr> operator-(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, const Expr&> operator-(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr> operator-(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::MINUS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, const Expr&> abs(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, Expr> abs(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, const Expr&> abs(const CpxScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, Expr> abs(CpxScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::ABS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::SCALAR, const Expr&> conj(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::CONJ, ExprType::SCALAR, Expr> conj(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::CONJ, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::EXP, ExprType::SCALAR, const Expr&> exp(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::EXP, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::EXP, ExprType::SCALAR, Expr> exp(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::EXP, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::EXP, ExprType::SCALAR, const Expr&> exp(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::EXP, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::EXP, ExprType::SCALAR, Expr> exp(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::EXP, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::LOG, ExprType::SCALAR, const Expr&> log(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::LOG, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::LOG, ExprType::SCALAR, Expr> log(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::LOG, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::LOG, ExprType::SCALAR, const Expr&> log(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::LOG, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::LOG, ExprType::SCALAR, Expr> log(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::LOG, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, const Expr&> sqrt(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, Expr> sqrt(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, const Expr&> sqrt(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, Expr> sqrt(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::SQRT, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::SIN, ExprType::SCALAR, const Expr&> sin(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::SIN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::SIN, ExprType::SCALAR, Expr> sin(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::SIN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::SIN, ExprType::SCALAR, const Expr&> sin(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::SIN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::SIN, ExprType::SCALAR, Expr> sin(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::SIN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::COS, ExprType::SCALAR, const Expr&> cos(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::COS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::COS, ExprType::SCALAR, Expr> cos(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::COS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::COS, ExprType::SCALAR, const Expr&> cos(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::COS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::COS, ExprType::SCALAR, Expr> cos(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::COS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TAN, ExprType::SCALAR, const Expr&> tan(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::TAN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::TAN, ExprType::SCALAR, Expr> tan(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::TAN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TAN, ExprType::SCALAR, const Expr&> tan(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::TAN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::TAN, ExprType::SCALAR, Expr> tan(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::TAN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, const Expr&> asin(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, Expr> asin(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, const Expr&> asin(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, Expr> asin(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::ASIN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, const Expr&> acos(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, Expr> acos(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, const Expr&> acos(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, Expr> acos(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::ACOS, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, const Expr&> atan(const ScalarExpression<Expr>& expr) { return UnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline UnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, Expr> atan(ScalarExpression<Expr>&& expr) { return UnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, const Expr&> atan(const CpxScalarExpression<Expr>& expr) { return CpxUnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<typename Expr>
inline CpxUnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, Expr> atan(CpxScalarExpression<Expr>&& expr) { return CpxUnaryExpression<UnaryOp::ATAN, ExprType::SCALAR, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator+(const ScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator+(const ScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator+(ScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator+(ScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator+(const ScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator+(const ScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator+(ScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator+(ScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator+(const CpxScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator+(const CpxScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator+(CpxScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator+(CpxScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator+(const CpxScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator+(const CpxScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator+(CpxScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator+(CpxScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator-(const ScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator-(const ScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator-(ScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator-(ScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator-(const ScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator-(const ScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator-(ScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator-(ScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator-(const CpxScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator-(const CpxScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator-(CpxScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator-(CpxScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator-(const CpxScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator-(const CpxScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator-(CpxScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator-(CpxScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const RankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, RankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK4_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const CpxRankFourTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, CpxRankFourTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK4_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const RankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, RankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::RK2_TENSOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const CpxRankTwoTensorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, CpxRankTwoTensorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const MatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, MatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::MATRIX, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const CpxMatrixExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, CpxMatrixExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::MATRIX, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const VectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, VectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::VECTOR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const CpxVectorExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, CpxVectorExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::VECTOR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const ScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const ScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(ScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(ScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator*(const CpxScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator*(const CpxScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator*(CpxScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator*(CpxScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const ScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const ScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(ScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(ScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const ScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const ScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(ScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(ScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxScalarExpression<LeftExpr>& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxScalarExpression<LeftExpr>& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxScalarExpression<LeftExpr>&& lhs, const ScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxScalarExpression<LeftExpr>&& lhs, ScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&> operator/(const CpxScalarExpression<LeftExpr>& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr> operator/(const CpxScalarExpression<LeftExpr>& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, const LeftExpr&, ExprType::SCALAR, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&> operator/(CpxScalarExpression<LeftExpr>&& lhs, const CpxScalarExpression<RightExpr>& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<typename LeftExpr, typename RightExpr>
inline CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> operator/(CpxScalarExpression<LeftExpr>&& lhs, CpxScalarExpression<RightExpr>&& rhs) { return CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

#endif // LINALG_OPERATORS_HPP
