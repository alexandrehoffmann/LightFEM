/*
 * MatrixBinaryExpression.hpp
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

#ifndef MATRIX_BINARY_EXPRESSION_HPP
#define MATRIX_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/MatrixExpression.hpp>

template<> struct BinaryOpType<BinaryOp::SUM,  ExprType::MATRIX, ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::SUB,  ExprType::MATRIX, ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::MATRIX, ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::SCALAR, ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::MATRIX, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DIV,  ExprType::MATRIX, ExprType::SCALAR> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DDOT, ExprType::RK4_TENSOR, ExprType::MATRIX> { static constexpr ExprType Type = ExprType::MATRIX; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public MatrixExpression< BinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getNrows() != m_rhs.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
		if (m_lhs.getNcols() != m_rhs.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); } }
public:
	inline size_t getNrows() const { return m_lhs.getNrows(); }
	inline size_t getNcols() const { return m_lhs.getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_lhs(i,j) + m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public MatrixExpression< BinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getNrows() != m_rhs.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
		if (m_lhs.getNcols() != m_rhs.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); } }
public:
	inline size_t getNrows() const { return m_lhs.getNrows(); }
	inline size_t getNcols() const { return m_lhs.getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_lhs(i,j) - m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public MatrixExpression< BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getNcols() != m_rhs.getNrows()) { throw std::invalid_argument("invalid matrix size."); } }
public:
	inline size_t getNrows() const { return m_lhs.getNrows(); }
	inline size_t getNcols() const { return m_rhs.getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { double value = 0.0; for (size_t k=0;k<m_lhs.getNcols();++k) { value += m_lhs(i,k)*m_rhs(k,j); } return value; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename MatrixExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::MATRIX, MatrixExpr> : public MatrixExpression< BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::MATRIX, MatrixExpr> >
{
public:
	BinaryExpression(ScalarExpr&& scalar, MatrixExpr&& matrix) : m_scalar(std::forward<ScalarExpr>(scalar)), m_matrix(std::forward<MatrixExpr>(matrix)) {}
public:
	inline size_t getNrows() const { return m_matrix.getNrows(); }
	inline size_t getNcols() const { return m_matrix.getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_matrix(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
};

template<typename ScalarExpr, typename MatrixExpr>
class BinaryExpression<BinaryOp::PROD,  ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> : public MatrixExpression< BinaryExpression<BinaryOp::PROD,  ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(MatrixExpr&& matrix, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_matrix(std::forward<MatrixExpr>(matrix)) {}
public:
	inline size_t getNrows() const { return m_matrix.getNrows(); }
	inline size_t getNcols() const { return m_matrix.getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_matrix(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
};

#include <iostream>

template<typename ScalarExpr, typename MatrixExpr>
class BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> : public MatrixExpression< BinaryExpression<BinaryOp::DIV, ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(MatrixExpr&& matrix, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_matrix(std::forward<MatrixExpr>(matrix)) {}
public:
	inline size_t getNrows() const { return m_matrix.getNrows(); }
	inline size_t getNcols() const { return m_matrix.getNcols(); }
public:
	inline double operator()(const size_t i, const size_t j) const { std::cout << m_matrix(i,j) << " / " <<  m_scalar.eval() << " = " << m_matrix(i,j) / m_scalar.eval() << std::endl; return m_matrix(i,j) / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
};

template<typename Rk4TensorExpr, typename MatrixExpr>
class BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr> : public MatrixExpression< BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr> >
{
public:
	BinaryExpression(Rk4TensorExpr&& rk4tensor, MatrixExpr&& matrix);
public:
	inline size_t getNrows() const { return m_nRows; }
	inline size_t getNcols() const { return m_nCols; }
public:
	inline double operator()(const size_t i, const size_t j) const;
private:
	typename RefTypeSelector<Rk4TensorExpr>::Type m_rk4tensor;
	typename RefTypeSelector<MatrixExpr>::Type    m_matrix;
	const size_t m_nRows; 
	const size_t m_nCols; 
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::SUM, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getNrows() != m_rhs.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
		if (m_lhs.getNcols() != m_rhs.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); } }
public:
	inline size_t getNrows() const { return m_lhs.getNrows(); }
	inline size_t getNcols() const { return m_lhs.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_lhs(i,j) + m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::SUB, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getNrows() != m_rhs.getNrows()) { throw std::invalid_argument("lhs and rhs must share the same number of rows"); }
		if (m_lhs.getNcols() != m_rhs.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same number of cols"); } }
public:
	inline size_t getNrows() const { return m_lhs.getNrows(); }
	inline size_t getNcols() const { return m_lhs.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_lhs(i,j) - m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getNcols() != m_rhs.getNrows()) { throw std::invalid_argument("invalid matrix size."); } }
public:
	inline size_t getNrows() const { return m_lhs.getNrows(); }
	inline size_t getNcols() const { return m_rhs.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { std::complex< double > value = 0.0; for (size_t k=0;k<m_lhs.getNcols();++k) { value += m_lhs(i,k)*m_rhs(k,j); } return value; }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename MatrixExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::MATRIX, MatrixExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::MATRIX, MatrixExpr> >
{
public:
	CpxBinaryExpression(ScalarExpr&& scalar, MatrixExpr&& matrix) : m_scalar(std::forward<ScalarExpr>(scalar)), m_matrix(std::forward<MatrixExpr>(matrix)) {}
public:
	inline size_t getNrows() const { return m_matrix.getNrows(); }
	inline size_t getNcols() const { return m_matrix.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_matrix(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
};

template<typename ScalarExpr, typename MatrixExpr>
class CpxBinaryExpression<BinaryOp::PROD,  ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::PROD,  ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(MatrixExpr&& matrix, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_matrix(std::forward<MatrixExpr>(matrix)) {}
public:
	inline size_t getNrows() const { return m_matrix.getNrows(); }
	inline size_t getNcols() const { return m_matrix.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_matrix(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
};

template<typename ScalarExpr, typename MatrixExpr>
class CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::DIV, ExprType::MATRIX, MatrixExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(MatrixExpr&& matrix, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_matrix(std::forward<MatrixExpr>(matrix)) {}
public:
	inline size_t getNrows() const { return m_matrix.getNrows(); }
	inline size_t getNcols() const { return m_matrix.getNcols(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_matrix(i,j) / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<MatrixExpr>::Type m_matrix;
};

template<typename Rk4TensorExpr, typename MatrixExpr>
class CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr> : public CpxMatrixExpression< CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr> >
{
public:
	CpxBinaryExpression(Rk4TensorExpr&& rk4tensor, MatrixExpr&& matrix);
public:
	inline size_t getNrows() const { return m_nRows; }
	inline size_t getNcols() const { return m_nCols; }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const;
private:
	typename RefTypeSelector<Rk4TensorExpr>::Type m_rk4tensor;
	typename RefTypeSelector<MatrixExpr>::Type    m_matrix;
	const size_t m_nRows; 
	const size_t m_nCols; 
};

#include <LightFEM/Expression/LinAlg/MatrixBinaryExpression.tpp>

#endif // MATRIX_BINARY_EXPRESSION_HPP
