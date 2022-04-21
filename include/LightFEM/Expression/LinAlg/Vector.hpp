/*
 * Vector.hpp
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

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <valarray>
#include <complex>

#include <LightFEM/Expression/LinAlg/VectorExpression.hpp>
#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>

class Vector : public VectorExpression<Vector>
{
public:
	Vector(const size_t size=1, const double value = 0.0) : m_core(value, size) {}
	template<typename Expr>
	Vector(const VectorExpression< Expr >& expr) : m_core(expr.getSize()) { for (size_t i=0; i<expr.getSize();++i) { m_core[i] = expr[i]; } }
	Vector(std::initializer_list< double > values) : m_core(values) {}
public:
	Vector& operator=(std::initializer_list< double > values) { m_core = values; return *this; }
	template<typename Expr> Vector& operator=  (const VectorExpression< Expr >& expr) { m_core.resize(expr.getSize()); for (size_t i=0; i<expr.getSize();++i) { m_core[i]  = expr[i]; } return *this; }
	template<typename Expr> Vector& operator+= (const VectorExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] += expr[i]; } return *this; }
	template<typename Expr> Vector& operator-= (const VectorExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] -= expr[i]; } return *this; }
	template<typename Expr> Vector& operator*= (const ScalarExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] /= expr.eval(); } return *this; }
	template<typename Expr> Vector& operator/= (const ScalarExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] *= expr.eval(); } return *this; }
public:
	inline size_t getSize() const { return m_core.size(); }
public:
	inline double  operator[](const size_t i) const { return m_core[i]; }
	inline double& operator[](const size_t i)       { return m_core[i]; }
private:
	std::valarray< double > m_core;
};

class CpxVector : public CpxVectorExpression<CpxVector>
{
public:
	CpxVector(const size_t size=1, std::complex< double > value = 0.0) : m_core(value, size) {}
	template<typename Expr>
	CpxVector(const VectorExpression< Expr >& expr) : m_core(expr.getSize()) { for (size_t i=0; i<expr.getSize();++i) { m_core[i] = expr[i]; } }
	template<typename Expr>
	CpxVector(const CpxVectorExpression< Expr >& expr) : m_core(expr.getSize()) { for (size_t i=0; i<expr.getSize();++i) { m_core[i] = expr[i]; } }
	CpxVector(std::initializer_list< std::complex< double > > values) : m_core(values) {}
public:
	template<typename Expr> CpxVector& operator=  (const VectorExpression< Expr >& expr)    { m_core.resize(expr.getSize()); for (size_t i=0; i<expr.getSize();++i) { m_core[i]  = expr[i]; } return *this; }
	template<typename Expr> CpxVector& operator=  (const CpxVectorExpression< Expr >& expr) { m_core.resize(expr.getSize()); for (size_t i=0; i<expr.getSize();++i) { m_core[i]  = expr[i]; } return *this; }
	template<typename Expr> CpxVector& operator+= (const VectorExpression< Expr >& expr)    { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] += expr[i]; } return *this; }
	template<typename Expr> CpxVector& operator+= (const CpxVectorExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] += expr[i]; } return *this; }
	template<typename Expr> CpxVector& operator-= (const VectorExpression< Expr >& expr)    { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] -= expr[i]; } return *this; }
	template<typename Expr> CpxVector& operator-= (const CpxVectorExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] -= expr[i]; } return *this; }
	template<typename Expr> CpxVector& operator*= (const ScalarExpression< Expr >& expr)    { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] /= expr.eval(); } return *this; }
	template<typename Expr> CpxVector& operator*= (const CpxScalarExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] /= expr.eval(); } return *this; }
	template<typename Expr> CpxVector& operator/= (const ScalarExpression< Expr >& expr)    { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] *= expr.eval(); } return *this; }
	template<typename Expr> CpxVector& operator/= (const CpxScalarExpression< Expr >& expr) { if (m_core.size() != expr.getSize()) { throw std::invalid_argument("the two vector must share the same size"); } for (size_t i=0; i<expr.getSize();++i) { m_core[i] *= expr.eval(); } return *this; }
public:
	inline size_t getSize() const { return m_core.size(); }
public:
	inline std::complex< double >  operator[](const size_t i) const { return m_core[i]; }
	inline std::complex< double >& operator[](const size_t i)       { return m_core[i]; }
private:
	std::valarray< std::complex< double > > m_core;
};

#endif // VECTOR_HPP
