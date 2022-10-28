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
////                       min/max operators                        ////
////////////////////////////////////////////////////////////////////////

template<typename Expr>
double min(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > min(const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
double max(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

template<typename Expr>
std::complex< double > max(const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& expr);

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
////                         abs operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&> abs(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, Expr> abs(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&> abs(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ABS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, Expr> abs(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ABS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         conj operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&> conj(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::CONJ, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr> conj(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::CONJ, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         exp operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&> exp(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&> exp(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, Expr> exp(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::EXP, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, Expr> exp(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::EXP, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         log operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&> log(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&> log(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, Expr> log(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::LOG, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, Expr> log(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::LOG, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         sqrt operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&> sqrt(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&> sqrt(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, Expr> sqrt(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SQRT, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, Expr> sqrt(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::SQRT, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         sin operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&> sin(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&> sin(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, Expr> sin(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::SIN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, Expr> sin(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::SIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         cos operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, const Expr&> cos(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, const Expr&> cos(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, Expr> cos(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::COS, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, Expr> cos(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::COS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         tan operators                          ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&> tan(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&> tan(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, Expr> tan(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::TAN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, Expr> tan(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::TAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         asin operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&> asin(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&> asin(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, Expr> asin(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ASIN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, Expr> asin(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ASIN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         acos operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&> acos(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&> acos(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, Expr> acos(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ACOS, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, Expr> acos(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ACOS, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

////////////////////////////////////////////////////////////////////////
////                         atan operators                         ////
////////////////////////////////////////////////////////////////////////

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&> atan(const ElementWiseFunctionExpression<Type, Expr>& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&> atan(const CpxElementWiseFunctionExpression<Type, Expr>& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, const Expr&>(std::forward<const Expr>( static_cast<const Expr&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline ElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, Expr> atan(ElementWiseFunctionExpression<Type, Expr>&& expr) { return ElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

template<ExprType Type, typename Expr, std::enable_if_t<UnaryOpType<UnaryOp::ATAN, Type>::isDefined, bool> = true>
inline CpxElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, Expr> atan(CpxElementWiseFunctionExpression<Type, Expr>&& expr) { return CpxElementWiseFunctionUnaryExpression<UnaryOp::ATAN, Type, Expr>(std::forward<Expr>( static_cast<Expr&&>(expr) )); }

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

#include <LightFEM/Expression/ElementWise/operators.tpp>

#endif // ELEMENT_WISE_FUNCTION_OPERATORS_HPP
