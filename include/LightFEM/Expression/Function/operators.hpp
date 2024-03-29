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

#ifndef FUNCTION_OPERATORS_HPP
#define FUNCTION_OPERATORS_HPP

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionBinaryExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionUnaryExpression.hpp>

#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/Function/FunctionBinaryExpression.hpp>
#include <LightFEM/Expression/Function/FunctionUnaryExpression.hpp>

////////////////////////////////////////////////////////////////////////
////                       min/max operators                        ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double min(const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > min(const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double max(const FunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > max(const CpxFunctionExpression<ExprType::SCALAR, Expr>& expr);

////////////////////////////////////////////////////////////////////////
////                        minus operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&> operator-(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&> operator-(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::MINUS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::MINUS, Type, Expr> operator-(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::MINUS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::MINUS, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr> operator-(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::MINUS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         abs operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&> abs(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ABS, Type, Expr> abs(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::ABS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&> abs(const CpxFunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ABS, Type, Expr> abs(CpxFunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::ABS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         conj operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&> conj(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr> conj(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         exp operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&> exp(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&> exp(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::EXP, Type, Expr> exp(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::EXP, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::EXP, Type, Expr> exp(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::EXP, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         log operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&> log(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&> log(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::LOG, Type, Expr> log(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::LOG, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::LOG, Type, Expr> log(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::LOG, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         sqrt operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&> sqrt(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&> sqrt(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::SQRT, Type, Expr> sqrt(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::SQRT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::SQRT, Type, Expr> sqrt(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::SQRT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         sin operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&> sin(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&> sin(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::SIN, Type, Expr> sin(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::SIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::SIN, Type, Expr> sin(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::SIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         cos operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::COS, Type, const Expr&> cos(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::COS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::COS, Type, const Expr&> cos(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::COS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::COS, Type, Expr> cos(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::COS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::COS, Type, Expr> cos(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::COS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         tan operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&> tan(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&> tan(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::TAN, Type, Expr> tan(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::TAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::TAN, Type, Expr> tan(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::TAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         asin operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&> asin(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&> asin(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ASIN, Type, Expr> asin(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::ASIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ASIN, Type, Expr> asin(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::ASIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         acos operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&> acos(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&> acos(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ACOS, Type, Expr> acos(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::ACOS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ACOS, Type, Expr> acos(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::ACOS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         atan operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&> atan(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&> atan(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::ATAN, Type, Expr> atan(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::ATAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ATAN, Type, Expr> atan(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::ATAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         norm operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&> norm(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&> norm(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::NORM, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::NORM, Type, Expr> norm(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::NORM, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::NORM, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::NORM, Type, Expr> norm(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::NORM, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                      transpose operators                       ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&> transpose(const FunctionExpression<Type, Expr>& expr) { return FunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&> transpose(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline FunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr> transpose(FunctionExpression<Type, Expr>&& expr) { return FunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TRANSPOSE, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr> transpose(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::TRANSPOSE, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                       adjoint operators                        ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ADJOINT, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ADJOINT, Type, const Expr&> adjoint(const CpxFunctionExpression<Type, Expr>& expr) { return CpxFunctionUnaryExpression<UnaryOp::ADJOINT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ADJOINT, Type>::isDefined, bool> = true>
inline CpxFunctionUnaryExpression<UnaryOp::ADJOINT, Type, Expr> adjoint(CpxFunctionExpression<Type, Expr>&& expr) { return CpxFunctionUnaryExpression<UnaryOp::ADJOINT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         ddot operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> ddot(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> ddot(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DDOT, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> ddot(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DDOT, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        inner operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> inner(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> inner(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::INNER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> inner(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::INNER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        outer operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> outer(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> outer(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::OUTER, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> outer(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::OUTER, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                        cross operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> cross(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> cross(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::CROSS, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> cross(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::CROSS, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         sum operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator+(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator+(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUM, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator+(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUM, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         sub operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator-(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator-(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::SUB, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator-(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::SUB, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         prod operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator*(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator*(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::PROD, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator*(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::PROD, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

////////////////////////////////////////////////////////////////////////
////                         div operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return FunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline FunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return FunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return ElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const ElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(const FunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(ElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(FunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const FunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const ElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, FunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, ElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, const RightExpr&>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, RightExpr>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType> operator/(const CpxElementWiseFunctionExpression<LeftType, LeftExpr>& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, const LeftExpr&, RightType, typename RightExpr::ReturnType>(std::forward<const LeftExpr>( static_cast<const LeftExpr&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(const CpxFunctionExpression<LeftType, LeftExpr>& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, const RightExpr&>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, const CpxElementWiseFunctionExpression<RightType, RightExpr>& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, const RightExpr&>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<const RightExpr>( static_cast<const RightExpr&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, RightExpr>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType> operator/(CpxElementWiseFunctionExpression<LeftType, LeftExpr>&& lhs, CpxFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, LeftExpr, RightType, typename RightExpr::ReturnType>(std::forward<LeftExpr>( static_cast<LeftExpr&&>(lhs) ), rhs[rhs.getMesh()->getElemId(lhs.getElement())]); }

template<ExprType LeftType, typename LeftExpr, ExprType RightType, typename RightExpr, std::enable_if_t<BinaryOpType<BinaryOp::DIV, LeftType, RightType>::isDefined, bool> = true>
inline CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr> operator/(CpxFunctionExpression<LeftType, LeftExpr>&& lhs, CpxElementWiseFunctionExpression<RightType, RightExpr>&& rhs) { return CpxElementWiseFunctionBinaryExpression<BinaryOp::DIV, LeftType, typename LeftExpr::ReturnType, RightType, RightExpr>(lhs[lhs.getMesh()->getElemId(rhs.getElement())], std::forward<RightExpr>( static_cast<RightExpr&&>(rhs) )); }

#include <LightFEM/Expression/Function/operators.tpp>

#endif // FUNCTION_OPERATORS_HPP
