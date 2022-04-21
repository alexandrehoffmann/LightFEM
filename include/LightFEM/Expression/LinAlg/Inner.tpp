/*
 * Inner.tpp
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

#ifndef INNER_TPP
#define INNER_TPP

#include <LightFEM/Expression/LinAlg/Inner.hpp>

template<typename LeftExpr, typename RightExpr>
BinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>::BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs)
{
	if (lhs.getSize() != rhs.getSize()) { throw std::invalid_argument("the two vectors must share the same size"); }

	m_value = 0.0;
	for (size_t i=0;i<lhs.getSize();++i) { m_value += lhs[i]*rhs[i]; }
}

template<typename LeftExpr, typename RightExpr>
BinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>::BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs)
{
	if (lhs.getNrows() != rhs.getNrows()) { throw std::invalid_argument("the two matrices must share the same size"); }
	if (lhs.getNcols() != rhs.getNcols()) { throw std::invalid_argument("the two matrices must share the same size"); }
	
	m_value = 0.0;
	for (size_t i=0;i<lhs.getNrows();++i) { for (size_t j=0;j<lhs.getNcols();++j)
	{
		m_value += lhs(i,j)*rhs(i,j);
	}}
}

template<typename LeftExpr, typename RightExpr>
BinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>::BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs)
{
	if (lhs.getShape() != rhs.getShape()) { throw std::invalid_argument("the two tensors must share the same shape"); }

	const std::array<size_t, 2> shape = lhs.getShape();

	m_value = 0.0;
	for (size_t i=0;i<shape[0];++i) { for (size_t j=0;j<shape[1];++j)
	{
		m_value += lhs(i,j)*rhs(i,j);
	}}
}

///////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
CpxBinaryExpression<BinaryOp::INNER, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>::CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs)
{
	if (lhs.getSize() != rhs.getSize()) { throw std::invalid_argument("the two vectors must share the same size"); }

	m_value = 0.0;
	for (size_t i=0;i<lhs.getSize();++i) { m_value += lhs[i]*std::conj(rhs[i]); }
}

template<typename LeftExpr, typename RightExpr>
CpxBinaryExpression<BinaryOp::INNER, ExprType::MATRIX, LeftExpr, ExprType::MATRIX, RightExpr>::CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs)
{
	if (lhs.getNrows() != rhs.getNrows()) { throw std::invalid_argument("the two matrices must share the same size"); }
	if (lhs.getNcols() != rhs.getNcols()) { throw std::invalid_argument("the two matrices must share the same size"); }
	
	m_value = 0.0;
	for (size_t i=0;i<lhs.getNrows();++i) { for (size_t j=0;j<lhs.getNcols();++j)
	{
		m_value += lhs(i,j)*std::conj(rhs(i,j));
	}}
}

template<typename LeftExpr, typename RightExpr>
CpxBinaryExpression<BinaryOp::INNER, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr>::CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs)
{
	if (lhs.getShape() != rhs.getShape()) { throw std::invalid_argument("the two tensors must share the same shape"); }

	const std::array<size_t, 2> shape = lhs.getShape();

	m_value = 0.0;
	for (size_t i=0;i<shape[0];++i) { for (size_t j=0;j<shape[1];++j)
	{
		m_value += lhs(i,j)*std::conj(rhs(i,j));
	}}
}

#endif // INNER_TPP
