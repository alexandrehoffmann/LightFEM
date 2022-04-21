/*
 * VectorBinaryExpression.hpp
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

#ifndef VECTOR_BINARY_EXPRESSION_HPP
#define VECTOR_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/VectorExpression.hpp>

template<> struct BinaryOpType<BinaryOp::SUM,  ExprType::VECTOR, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::SUB,  ExprType::VECTOR, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::VECTOR, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::MATRIX, ExprType::VECTOR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DIV,  ExprType::VECTOR, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::VECTOR; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public VectorExpression< BinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getSize() != m_rhs.getSize()) { throw std::invalid_argument("lhs and rhs must share the same size"); }}
public:
	inline size_t getSize() const { return m_lhs.getSize(); }
public:
	inline double operator[](const size_t i) const { return m_lhs[i] + m_rhs[i]; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public VectorExpression< BinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getSize() != m_rhs.getSize()) { throw std::invalid_argument("lhs and rhs must share the same size"); }}
public:
	inline size_t getSize() const { return m_lhs.getSize(); }
public:
	inline double operator[](const size_t i) const { return m_lhs[i] - m_rhs[i]; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename VectorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::VECTOR, VectorExpr> : public VectorExpression< BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::VECTOR, VectorExpr> >
{
public:
	BinaryExpression(ScalarExpr&& scalar, VectorExpr&& vector) : m_scalar(std::forward<ScalarExpr>(scalar)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_vector.getSize(); }
public:
	inline double operator[](const size_t i) const { return m_scalar.eval()*m_vector[i]; }
private:
	typename RefTypeSelector<ScalarExpr>::Type  m_scalar;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

template<typename ScalarExpr, typename VectorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> : public VectorExpression< BinaryExpression<BinaryOp::PROD, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(VectorExpr&& vector, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_vector.getSize(); }
public:
	inline double operator[](const size_t i) const { return m_scalar.eval()*m_vector[i]; }
private:
	typename RefTypeSelector<ScalarExpr>::Type  m_scalar;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

template<typename MatrixExpr, typename VectorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, MatrixExpr, ExprType::VECTOR, VectorExpr>  : public VectorExpression< BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, MatrixExpr, ExprType::VECTOR, VectorExpr> >
{
public:
	BinaryExpression(MatrixExpr&& matrix, VectorExpr&& vector) : m_matrix(std::forward<MatrixExpr>(matrix)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_matrix.getNrows(); }
public:
	inline double operator[](const size_t i) const { double value = 0.0; for (size_t j=0;j<m_vector.getSize();++j) { value += m_matrix(i,j)*m_vector[j]; } return value; }
private:
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

template<typename ScalarExpr, typename VectorExpr>
class BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> : public VectorExpression< BinaryExpression<BinaryOp::DIV, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(VectorExpr&& vector, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_vector.getSize(); }
public:
	inline double operator[](const size_t i) const { return m_vector[i] / m_scalar.eval(); }
public:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::SUM, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getSize() != m_rhs.getSize()) { throw std::invalid_argument("lhs and rhs must share the same size"); }}
public:
	inline size_t getSize() const { return m_lhs.getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return m_lhs[i] + m_rhs[i]; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::SUB, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getSize() != m_rhs.getSize()) { throw std::invalid_argument("lhs and rhs must share the same size"); }}
public:
	inline size_t getSize() const { return m_lhs.getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return m_lhs[i] - m_rhs[i]; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename VectorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::VECTOR, VectorExpr> : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::VECTOR, VectorExpr> >
{
public:
	CpxBinaryExpression(ScalarExpr&& scalar, VectorExpr&& vector) : m_scalar(std::forward<ScalarExpr>(scalar)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_vector.getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return m_scalar.eval()*m_vector[i]; }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

template<typename ScalarExpr, typename VectorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(VectorExpr&& vector, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_vector.getSize(); }
public:
	inline std::complex< double > operator[](const size_t i) const { return m_scalar.eval()*m_vector[i]; }
private:
	typename RefTypeSelector<ScalarExpr>::Type  m_scalar;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

template<typename MatrixExpr, typename VectorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, MatrixExpr, ExprType::VECTOR, VectorExpr>  : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, MatrixExpr, ExprType::VECTOR, VectorExpr> >
{
public:
	CpxBinaryExpression(MatrixExpr&& matrix, VectorExpr&& vector) : m_matrix(std::forward<MatrixExpr>(matrix)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_matrix.getNrows(); }
public:
	inline std::complex< double > operator[](const size_t i) const { std::complex< double > value = 0.0; for (size_t j=0;j<m_vector.getSize();++j) { value += m_matrix(i,j)*m_vector[j]; } return value; }
private:
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

template<typename ScalarExpr, typename VectorExpr>
class CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> : public CpxVectorExpression< CpxBinaryExpression<BinaryOp::DIV, ExprType::VECTOR, VectorExpr, ExprType::SCALAR, ScalarExpr> > 
{
public:
	CpxBinaryExpression(VectorExpr&& vector, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_vector(std::forward<VectorExpr>(vector)) {}
public:
	inline size_t getSize() const { return m_vector.getSize(); }
public:
	inline double operator[](const size_t i) const { return m_vector[i] / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<VectorExpr>::Type m_vector;
};

#endif // VECTOR_BINARY_EXPRESSION_HPP
