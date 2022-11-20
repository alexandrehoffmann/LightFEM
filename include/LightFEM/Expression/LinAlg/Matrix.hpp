/*
 * Matrix.hpp
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

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <valarray>
#include <complex>

#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>
#include <LightFEM/Expression/LinAlg/MatrixExpression.hpp>

class Matrix : public MatrixExpression< Matrix >
{
public:
	Matrix(const size_t nrows = 1, const size_t ncols = 1, const double value = 0.0) : m_nrows(nrows), m_ncols(ncols), m_core(value, nrows*ncols) {}
	Matrix(std::initializer_list< std::initializer_list< double > > values);
	template<typename Expr>
	Matrix(const MatrixExpression< Expr >& expr);
public:
	template<typename Expr> Matrix& operator=  (const MatrixExpression< Expr >& expr);
	template<typename Expr> Matrix& operator+= (const MatrixExpression< Expr >& expr);
	template<typename Expr> Matrix& operator-= (const MatrixExpression< Expr >& expr);
	template<typename Expr> Matrix& operator*= (const ScalarExpression< Expr >& expr);
	template<typename Expr> Matrix& operator/= (const ScalarExpression< Expr >& expr);
public:
	inline size_t getNrows() const { return m_nrows; }
	inline size_t getNcols() const { return m_ncols; }
	inline void resize(const size_t nrows, const size_t ncols) { m_nrows = nrows; m_ncols = ncols; m_core.resize(nrows*ncols); }
public:
	inline double  operator()(const size_t i, const size_t j) const { return m_core[flatIndex(i,j)]; }
	inline double& operator()(const size_t i, const size_t j)       { return m_core[flatIndex(i,j)]; }
private:
	inline size_t flatIndex(const size_t i, const size_t j) const { return i*m_ncols + j; } 
private:
	size_t m_nrows;
	size_t m_ncols;
	std::valarray< double > m_core;
};

class CpxMatrix : public CpxMatrixExpression< CpxMatrix >
{
public:
	CpxMatrix(const size_t nrows = 1, const size_t ncols = 1, std::complex< double > value = 0.0) : m_nrows(nrows), m_ncols(ncols), m_core(value, nrows*ncols) {}
	CpxMatrix(std::initializer_list< std::initializer_list< std::complex< double > > > values);
	template<typename Expr>
	CpxMatrix(const MatrixExpression< Expr >& expr) : m_nrows(expr.getNrows()), m_ncols(expr.getNcols()), m_core(m_nrows*m_ncols) { for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) { m_core[flatIndex(i,j)] =  expr(i,j); } } }
	template<typename Expr>
	CpxMatrix(const CpxMatrixExpression< Expr >& expr) : m_nrows(expr.getNrows()), m_ncols(expr.getNcols()), m_core(m_nrows*m_ncols) { for (size_t i=0; i<expr.getNrows();++i) { for (size_t j=0; j<expr.getNcols();++j) { m_core[flatIndex(i,j)] =  expr(i,j); } } }
public:
	template<typename Expr> CpxMatrix& operator=  (const MatrixExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator=  (const CpxMatrixExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator+= (const MatrixExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator+= (const CpxMatrixExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator-= (const MatrixExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator-= (const CpxMatrixExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator*= (const ScalarExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator*= (const CpxScalarExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator/= (const ScalarExpression< Expr >& expr);
	template<typename Expr> CpxMatrix& operator/= (const CpxScalarExpression< Expr >& expr);
public:
	inline size_t getNrows() const { return m_nrows; }
	inline size_t getNcols() const { return m_ncols; }
public:
	inline std::complex< double >  operator()(const size_t i, const size_t j) const { return m_core[i*m_ncols + j]; }
	inline std::complex< double >& operator()(const size_t i, const size_t j)       { return m_core[i*m_ncols + j]; }
private:
	inline size_t flatIndex(const size_t i, const size_t j) const { return i*m_ncols + j; } 
private:
	size_t m_nrows;
	size_t m_ncols;
	std::valarray< std::complex< double > > m_core;
};

#include <LightFEM/Expression/LinAlg/Matrix.tpp>

#endif // MATRIX_EXPRESSION_HPP
