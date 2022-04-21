/*
 * Cross.tpp
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

#ifndef CROSS_TPP
#define CROSS_TPP

#include <LightFEM/Expression/LinAlg/Outer.hpp>

template<typename LeftExpr, typename RightExpr>
BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>::BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : 
	m_lhs(std::forward<LeftExpr>(lhs)), 
	m_rhs(std::forward<RightExpr>(rhs)) 
{ 
	if (m_lhs.getSize() != m_rhs.getSize()) 
	{ 
		throw std::invalid_argument("the two vectors must share the same size"); 
	} 
	if (m_lhs.getSize() != 3) 
	{ 
		throw std::invalid_argument("the two vectors must be of size 3"); 
	}
}

template<typename LeftExpr, typename RightExpr>
double BinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>::operator[](const size_t i) const 
{ 
	if (i == 0) { return m_lhs[1]*m_rhs[2] - m_lhs[2]*m_rhs[1]; }; 
	if (i == 1) { return m_lhs[2]*m_rhs[0] - m_lhs[0]*m_rhs[2]; }; 
	if (i == 2) { return m_lhs[0]*m_rhs[1] - m_lhs[1]*m_rhs[0]; }; 
	
	throw std::invalid_argument("vector of size 3"); 
}

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>::CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : 
	m_lhs(std::forward<LeftExpr>(lhs)), 
	m_rhs(std::forward<RightExpr>(rhs)) 
{ 
	if (m_lhs.getSize() != m_rhs.getSize()) 
	{ 
		throw std::invalid_argument("the two vectors must share the same size"); 
	} 
	if (m_lhs.getSize() != 3) 
	{ 
		throw std::invalid_argument("the two vectors must be of size 3"); 
	}
}

template<typename LeftExpr, typename RightExpr>
std::complex< double > CpxBinaryExpression<BinaryOp::CROSS, ExprType::VECTOR, LeftExpr, ExprType::VECTOR, RightExpr>::operator[](const size_t i) const 
{ 
	if (i == 0) { return m_lhs[1]*m_rhs[2] - m_lhs[2]*m_rhs[1]; }; 
	if (i == 1) { return m_lhs[2]*m_rhs[0] - m_lhs[0]*m_rhs[2]; }; 
	if (i == 2) { return m_lhs[0]*m_rhs[1] - m_lhs[1]*m_rhs[0]; }; 
	
	throw std::invalid_argument("vector of size 3"); 
}

#endif // CROSS_TPP
