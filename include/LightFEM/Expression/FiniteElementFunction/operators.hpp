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

#ifndef FINITE_ELEMENT_FUNCTION_OPERATORS_HPP
#define FINITE_ELEMENT_FUNCTION_OPERATORS_HPP

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionBinaryExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionUnaryExpression.hpp>

////////////////////////////////////////////////////////////////////////
////                         conj operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&> conj(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr> conj(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                      transpose operators                       ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline FiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&> transpose(const FiniteElementFunctionExpression<Type, Expr>& expr) { return FiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&> transpose(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline FiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr> transpose(FiniteElementFunctionExpression<Type, Expr>&& expr) { return FiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr> transpose(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                       adjoint operators                        ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ADJOINT, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::ADJOINT, Type, const Expr&> adjoint(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::ADJOINT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ADJOINT, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::ADJOINT, Type, Expr> adjoint(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::ADJOINT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                        minus operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline FiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&> operator-(const FiniteElementFunctionExpression<Type, Expr>& expr) { return FiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&> operator-(const CpxFiniteElementFunctionExpression<Type, Expr>& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline FiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr> operator-(FiniteElementFunctionExpression<Type, Expr>&& expr) { return FiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline CpxFiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr> operator-(CpxFiniteElementFunctionExpression<Type, Expr>&& expr) { return CpxFiniteElementFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         sum operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         sub operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         prod operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, const RightExpr&> operator*(const Const& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, const RightExpr&>(std::forward<const Const>( static_cast<const Const&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, RightExpr> operator*(const Const& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, RightExpr>(std::forward<const Const>( static_cast<const Const&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, const RightExpr&> operator*(Const&& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, const RightExpr&>(std::forward<Const>( static_cast<Const&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, RightExpr> operator*(Const&& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, RightExpr>(std::forward<Const>( static_cast<Const&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, const RightExpr&> operator*(const Const& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, const RightExpr&>(std::forward<const Const>( static_cast<const Const&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, RightExpr> operator*(const Const& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const Const&, RightType, RightExpr>(std::forward<const Const>( static_cast<const Const&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, const RightExpr&> operator*(Const&& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, const RightExpr&>(std::forward<Const>( static_cast<Const&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, RightExpr> operator*(Const&& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, Const, RightType, RightExpr>(std::forward<Const>( static_cast<Const&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, const RightExpr&> operator*(const CpxConst& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, const RightExpr&>(std::forward<const CpxConst>( static_cast<const CpxConst&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, RightExpr> operator*(const CpxConst& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, RightExpr>(std::forward<const CpxConst>( static_cast<const CpxConst&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, const RightExpr&> operator*(CpxConst&& lhs, const FiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, const RightExpr&>(std::forward<CpxConst>( static_cast<CpxConst&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, RightExpr> operator*(CpxConst&& lhs, FiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, RightExpr>(std::forward<CpxConst>( static_cast<CpxConst&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, const RightExpr&> operator*(const CpxConst& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, const RightExpr&>(std::forward<const CpxConst>( static_cast<const CpxConst&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, RightExpr> operator*(const CpxConst& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, const CpxConst&, RightType, RightExpr>(std::forward<const CpxConst>( static_cast<const CpxConst&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, const RightExpr&> operator*(CpxConst&& lhs, const CpxFiniteElementFunctionExpression<RightType, RightExpr>& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, const RightExpr&>(std::forward<CpxConst>( static_cast<CpxConst&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, RightType>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, RightExpr> operator*(CpxConst&& lhs, CpxFiniteElementFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, CpxConst, RightType, RightExpr>(std::forward<CpxConst>( static_cast<CpxConst&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&> operator*(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const Const& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, Const> operator*(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, Const&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, Const>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const Const&> operator*(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const Const& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const Const&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, Const> operator*(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, Const&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, Const>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&> operator*(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst> operator*(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&> operator*(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, CpxConst> operator*(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, CpxConst>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&> operator*(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const Const& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, Const> operator*(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, Const&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, Const>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const Const&> operator*(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const Const& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const Const&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, Const> operator*(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, Const&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, Const>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&> operator*(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst> operator*(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&> operator*(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, CpxConst> operator*(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, ExprType::SCALAR, CpxConst>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&> operator/(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const Const& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, Const> operator/(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, Const&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, Const>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const Const&> operator/(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const Const& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const Const&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, Const> operator/(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, Const&& rhs) { return FiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, Const>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&> operator/(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst> operator/(const FiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&> operator/(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, CpxConst> operator/(FiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, CpxConst>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&> operator/(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const Const& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const Const&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, Const> operator/(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, Const&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, Const>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const Const&> operator/(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const Const& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const Const&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const Const>( static_cast<const Const&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, Const> operator/(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, Const&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, Const>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<Const>( static_cast<Const&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&> operator/(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, const CpxConst&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst> operator/(const CpxFiniteElementFunctionExpression<LeftType, LeftExpr>& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, ExprType::SCALAR, CpxConst>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&> operator/(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxConst& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, const CpxConst&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const CpxConst>( static_cast<const CpxConst&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, ExprType::SCALAR>::isDefined, bool> = true>
inline CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, CpxConst> operator/(CpxFiniteElementFunctionExpression<LeftType, LeftExpr>&& lhs, CpxConst&& rhs) { return CpxFiniteElementFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, ExprType::SCALAR, CpxConst>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<CpxConst>( static_cast<CpxConst&&>(rhs) )); }

#endif // FINITE_ELEMENT_FUNCTION_OPERATORS_HPP
