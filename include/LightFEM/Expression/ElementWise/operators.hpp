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

#ifndef ELEMENT_WISE_FUNCTION_OPERATORS_HPP
#define ELEMENT_WISE_FUNCTION_OPERATORS_HPP

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionBinaryExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionUnaryExpression.hpp>

////////////////////////////////////////////////////////////////////////
////                         conj operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&> conj(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr> conj(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         norm operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&> norm(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&> norm(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, Expr> norm(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, Expr> norm(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::NORM, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                      transpose operators                       ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&> transpose(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&> transpose(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr> transpose(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr> transpose(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                       adjoint operators                        ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ADJOINT, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ADJOINT, Type, const Expr&> adjoint(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ADJOINT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ADJOINT, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ADJOINT, Type, Expr> adjoint(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ADJOINT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                        minus operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&> operator-(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&> operator-(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr> operator-(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr> operator-(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         ddot operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        inner operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        outer operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        cross operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         sum operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         sub operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         prod operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         div operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

#endif // ELEMENT_WISE_FUNCTION_OPERATORS_HPP
