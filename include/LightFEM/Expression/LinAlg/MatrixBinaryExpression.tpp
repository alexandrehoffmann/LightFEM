/*
 * MatrixBinaryExpression.tpp
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

#ifndef MATRIX_BINARY_EXPRESSION_TPP
#define MATRIX_BINARY_EXPRESSION_TPP

#include <LightFEM/Expression/LinAlg/MatrixBinaryExpression.hpp>

template<typename Rk4TensorExpr, typename MatrixExpr>
BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr>::BinaryExpression(Rk4TensorExpr&& rk4tensor, MatrixExpr&& matrix) : 
	m_rk4tensor(std::forward<Rk4TensorExpr>(rk4tensor)), 
	m_matrix(std::forward<MatrixExpr>(matrix)) 
{ 
	std::array<size_t, 4> rk4shape = m_rk4tensor.getShape(); 
	
	if (rk4shape[2] != m_matrix.getNrows() or rk4shape[3] != m_matrix.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } 
	
	m_nRows = rk4shape[0]; 
	m_nCols = rk4shape[1];  
}

template<typename Rk4TensorExpr, typename MatrixExpr>
double BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr>::operator()(const size_t i, const size_t j) const
{
	double value = 0; 
	for (size_t k=0;k<m_matrix.getNrows();++k) { for (size_t l=0;l<m_matrix.getNcols();++l)
	{ 
		value += m_rk4tensor(i,j,k,l)*m_matrix(k,l); 
	}} 
	return value; 
}

//////////////////////////////////////////////////////////////////

template<typename Rk4TensorExpr, typename MatrixExpr>
CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr>::CpxBinaryExpression(Rk4TensorExpr&& rk4tensor, MatrixExpr&& matrix) : 
	m_rk4tensor(std::forward<Rk4TensorExpr>(rk4tensor)), 
	m_matrix(std::forward<MatrixExpr>(matrix)) 
{ 
	std::array<size_t, 4> rk4shape = m_rk4tensor.getShape(); 
	
	if (rk4shape[2] != m_matrix.getNrows() or rk4shape[3] != m_matrix.getNcols()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } 
	
	m_nRows = rk4shape[0]; 
	m_nCols = rk4shape[1];  
}

template<typename Rk4TensorExpr, typename MatrixExpr>
std::complex< double > CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::MATRIX, MatrixExpr>::operator()(const size_t i, const size_t j) const
{
	std::complex< double > value = 0; 
	for (size_t k=0;k<m_matrix.getNrows();++k) { for (size_t l=0;l<m_matrix.getNcols();++l)
	{ 
		value += m_rk4tensor(i,j,k,l)*m_matrix(k,l); 
	}} 
	return value; 
}

#endif // MATRIX_BINARY_EXPRESSION_TPP
