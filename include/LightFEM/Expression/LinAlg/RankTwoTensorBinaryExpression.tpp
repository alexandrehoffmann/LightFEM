/*
 * RankTwoTensorBinaryExpression.tpp
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

#ifndef RANK_TWO_TENSOR_BINARY_EXPRESSION_TPP
#define RANK_TWO_TENSOR_BINARY_EXPRESSION_TPP

#include <LightFEM/Expression/LinAlg/RankTwoTensorBinaryExpression.hpp>

template<typename Rk4TensorExpr, typename Rk2TensorExpr>
BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr>::BinaryExpression(Rk4TensorExpr&& rk4tensor, Rk2TensorExpr&& rk2tensor) :
	m_rk4tensor(std::forward<Rk4TensorExpr>(rk4tensor)), 
	m_rk2tensor(std::forward<Rk2TensorExpr>(rk2tensor)) 
{ 
	std::array<size_t, 4> rk4shape = m_rk4tensor.getShape();
	std::array<size_t, 2> rk2shape = m_rk2tensor.getShape();
	
	if (rk4shape[2] != rk2shape[0] or rk4shape[3] != rk2shape[1]) { throw std::invalid_argument("lhs and rhs must share the same shape"); } 
	
	m_shape[0] = rk4shape[0]; 
	m_shape[1] = rk4shape[1]; 
}

template<typename Rk4TensorExpr, typename Rk2TensorExpr>
inline double BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr>::BinaryExpression::operator()(const size_t i, const size_t j) const
{ 
	std::array<size_t, 2> shape = m_rk2tensor.getShape();
	
	double value = 0; 
	for (size_t k=0;k<shape[0];++k) { for (size_t l=0;l<shape[1];++l)
	{ 
		value += m_rk4tensor(i,j,k,l)*m_rk2tensor(k,l); 
	}} 
	return value; 
}

//////////////////////////////////////////////////////////////////

template<typename Rk4TensorExpr, typename Rk2TensorExpr>
CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr>::CpxBinaryExpression(Rk4TensorExpr&& rk4tensor, Rk2TensorExpr&& rk2tensor) :
	m_rk4tensor(std::forward<Rk4TensorExpr>(rk4tensor)),
	m_rk2tensor(std::forward<Rk2TensorExpr>(rk2tensor))
{
	std::array<size_t, 4> rk4shape = m_rk4tensor.getShape();
	std::array<size_t, 2> rk2shape = m_rk2tensor.getShape();

	if (rk4shape[2] != rk2shape[0] or rk4shape[3] != rk2shape[1]) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	m_shape[0] = rk4shape[0];
	m_shape[1] = rk4shape[1];
}

template<typename Rk4TensorExpr, typename Rk2TensorExpr>
inline std::complex< double > CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr>::operator()(const size_t i, const size_t j) const
{
	std::array<size_t, 2> shape = m_rk2tensor.getShape();

	std::complex< double > value = 0;
	for (size_t k=0;k<shape[0];++k) { for (size_t l=0;l<shape[1];++l)
	{
		value += m_rk4tensor(i,j,k,l)*m_rk2tensor(k,l);
	}}
	return value;
}

#endif // RANK_TWO_TENSOR_BINARY_EXPRESSION_TPP
