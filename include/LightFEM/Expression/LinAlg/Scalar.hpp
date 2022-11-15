/*
 * Scalar.hpp
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

#ifndef SCALAR_HPP
#define SCALAR_HPP

#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>

class Scalar : public ScalarExpression<Scalar>
{
public:
	Scalar() : m_value(0) {}
	Scalar(const double value) : m_value(value) {}
	template<typename Expr>
	Scalar(const ScalarExpression< Expr >& expr) : m_value(expr.eval()) {}
public:
	template<typename Expr> Scalar& operator=  (const ScalarExpression< Expr >& expr) { m_value  = expr.eval(); return *this; }
	template<typename Expr> Scalar& operator+= (const ScalarExpression< Expr >& expr) { m_value += expr.eval(); return *this; }
	template<typename Expr> Scalar& operator-= (const ScalarExpression< Expr >& expr) { m_value -= expr.eval(); return *this; }
	template<typename Expr> Scalar& operator*= (const ScalarExpression< Expr >& expr) { m_value *= expr.eval(); return *this; }
	template<typename Expr> Scalar& operator/= (const ScalarExpression< Expr >& expr) { m_value /= expr.eval(); return *this; }
public:
	Scalar& operator=  (const double rhs) { m_value = rhs;  return *this; }
	Scalar& operator+= (const double rhs) { m_value += rhs; return *this; }
	Scalar& operator-= (const double rhs) { m_value -= rhs; return *this; }
	Scalar& operator*= (const double rhs) { m_value *= rhs; return *this; }
	Scalar& operator/= (const double rhs) { m_value /= rhs; return *this; }
public:
	inline double eval() const { return m_value; }
private:
	double m_value;
};

////////////////////////////////////////////////////////////////////////

class CpxScalar : public CpxScalarExpression<CpxScalar>
{
public:
	CpxScalar() : m_value(0) {}
	CpxScalar(const double value) : m_value(value) {}
	CpxScalar(const double real_value, const double imag_value) : m_value(real_value, imag_value) {}
	CpxScalar(const std::complex< double > value) : m_value(value) {}
	template<typename Expr>
	CpxScalar(const ScalarExpression< Expr >& expr) : m_value(expr.eval()) {}
	template<typename Expr>
	CpxScalar(const CpxScalarExpression< Expr >& expr) : m_value(expr.eval()) {}
public:
	template<typename Expr> CpxScalar& operator=  (const ScalarExpression< Expr >& expr)    { m_value  = expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator=  (const CpxScalarExpression< Expr >& expr) { m_value  = expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator+= (const ScalarExpression< Expr >& expr)    { m_value += expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator+= (const CpxScalarExpression< Expr >& expr) { m_value += expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator-= (const ScalarExpression< Expr >& expr)    { m_value -= expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator-= (const CpxScalarExpression< Expr >& expr) { m_value -= expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator*= (const ScalarExpression< Expr >& expr)    { m_value *= expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator*= (const CpxScalarExpression< Expr >& expr) { m_value *= expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator/= (const ScalarExpression< Expr >& expr)    { m_value /= expr.eval(); return *this; }
	template<typename Expr> CpxScalar& operator/= (const CpxScalarExpression< Expr >& expr) { m_value /= expr.eval(); return *this; }
public:
	CpxScalar& operator=  (const double rhs)                  { m_value  = rhs; return *this; }
	CpxScalar& operator=  (const std::complex< double >& rhs) { m_value  = rhs; return *this; }
	CpxScalar& operator+= (const double rhs)                  { m_value += rhs; return *this; }
	CpxScalar& operator+= (const std::complex< double >& rhs) { m_value += rhs; return *this; }
	CpxScalar& operator-= (const double rhs)                  { m_value -= rhs; return *this; }
	CpxScalar& operator-= (const std::complex< double >& rhs) { m_value -= rhs; return *this; }
	CpxScalar& operator*= (const double rhs)                  { m_value *= rhs; return *this; }
	CpxScalar& operator*= (const std::complex< double >& rhs) { m_value *= rhs; return *this; }
	CpxScalar& operator/= (const double rhs)                  { m_value /= rhs; return *this; }
	CpxScalar& operator/= (const std::complex< double >& rhs) { m_value /= rhs; return *this; }
public:
	inline std::complex< double > eval() const { return m_value; }
private:
	std::complex< double > m_value;
};

#endif // SCALAR_HPP
