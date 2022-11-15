/*
 * ScalarBinaryExpression.hpp
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

#ifndef SCALAR_BINARY_EXPRESSION_HPP
#define SCALAR_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>

template<> struct BinaryOpType<BinaryOp::SUM,  ExprType::SCALAR, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::SUB,  ExprType::SCALAR, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DIV,  ExprType::SCALAR, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline double eval() const { return m_lhs.eval() + m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline double eval() const { return m_lhs.eval() - m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline double eval() const { return m_lhs.eval()*m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public ScalarExpression< BinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline double eval() const { return m_lhs.eval() / m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::SUM, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline std::complex< double > eval() const { return m_lhs.eval() + m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::SUB, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline std::complex< double > eval() const { return m_lhs.eval() - m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline std::complex< double > eval() const { return m_lhs.eval()*m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> : public CpxScalarExpression< CpxBinaryExpression<BinaryOp::DIV, ExprType::SCALAR, LeftExpr, ExprType::SCALAR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::SCALAR;
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}
public:
	inline std::complex< double > eval() const { return m_lhs.eval() / m_rhs.eval(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

#endif // SCALAR_EXPRESSION_HPP
