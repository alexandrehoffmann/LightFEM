/*
 * Norm.hpp
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

#ifndef NORM_HPP
#define NORM_HPP

#include <complex>
#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>

template<> struct UnaryOpType<UnaryOp::NORM, ExprType::SCALAR>     { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::NORM, ExprType::VECTOR>     { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::NORM, ExprType::MATRIX>     { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };
template<> struct UnaryOpType<UnaryOp::NORM, ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::SCALAR; static constexpr bool isDefined = true; };

template<typename Expr>
class UnaryExpression<UnaryOp::NORM, ExprType::SCALAR, Expr> : public ScalarExpression< UnaryExpression<UnaryOp::NORM, ExprType::SCALAR, Expr> >
{
public:
	UnaryExpression(Expr&& expr) { const double value = expr.eval(); m_value = std::sqrt(value*value); }
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

template<typename Expr>
class UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr> : public ScalarExpression< UnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

template<typename Expr>
class UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr> : public ScalarExpression< UnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr> >
{
public:
	UnaryExpression(Expr&& expr);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

template<typename Expr>
class UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr> : public ScalarExpression< UnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr> >
{
public:
	UnaryExpression(Expr&& expr);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

////////////////////////////////////////////////////////////////////////

template<typename Expr>
class CpxUnaryExpression<UnaryOp::NORM, ExprType::SCALAR, Expr> : public ScalarExpression< CpxUnaryExpression<UnaryOp::NORM, ExprType::SCALAR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr) { const double value = expr.eval(); m_value = std::sqrt(value*std::conj(value)); }
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

template<typename Expr>
class CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr> : public ScalarExpression< CpxUnaryExpression<UnaryOp::NORM, ExprType::VECTOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

template<typename Expr>
class CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr> : public ScalarExpression< CpxUnaryExpression<UnaryOp::NORM, ExprType::MATRIX, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

template<typename Expr>
class CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr> : public ScalarExpression< CpxUnaryExpression<UnaryOp::NORM, ExprType::RK2_TENSOR, Expr> >
{
public:
	CpxUnaryExpression(Expr&& expr);
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
}; 

#include <LightFEM/Expression/LinAlg/Norm.tpp>

#endif // NORM_HPP
