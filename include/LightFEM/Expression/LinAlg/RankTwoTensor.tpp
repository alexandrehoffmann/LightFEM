/*
 * RankTwoTensor.tpp
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

#ifndef RANK_TWO_TENSOR_TPP
#define RANK_TWO_TENSOR_TPP

#include <LightFEM/Expression/LinAlg/RankTwoTensor.hpp>

#include <numeric>
#include <functional>

template<typename Expr>
RankTwoTensor::RankTwoTensor(const RankTwoTensorExpression< Expr >& expr) : 
	m_shape(expr.getShape()),
	m_core(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()))
{ 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] =  expr(i,j); 
		} 
	} 
}

template<typename Expr>
RankTwoTensor& RankTwoTensor::operator= (const RankTwoTensorExpression< Expr >& expr)
{ 
	m_shape = expr.getShape();
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()));

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] = expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankTwoTensor& RankTwoTensor::operator+= (const RankTwoTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] += expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankTwoTensor& RankTwoTensor::operator-= (const RankTwoTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] -= expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankTwoTensor& RankTwoTensor::operator*= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] *= expr.eval(); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
RankTwoTensor& RankTwoTensor::operator/= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] /= expr.eval(); 
		} 
	} 
	
	return *this;
}

////////////////////////////////////////////////////////////////////////

template<typename Expr>
CpxRankTwoTensor::CpxRankTwoTensor(const RankTwoTensorExpression< Expr >& expr) : 
	m_shape(expr.getShape()),
	m_core(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()))
{ 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] =  expr(i,j); 
		} 
	} 
}

template<typename Expr>
CpxRankTwoTensor::CpxRankTwoTensor(const CpxRankTwoTensorExpression< Expr >& expr) : 
	m_shape(expr.getShape()),
	m_core(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()))
{ 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] =  expr(i,j); 
		} 
	} 
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator= (const RankTwoTensorExpression< Expr >& expr)
{
	m_shape = expr.getShape();
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()));
 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] =  expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator= (const CpxRankTwoTensorExpression< Expr >& expr)
{
	m_shape = expr.getShape();
	m_core.resize(std::accumulate(m_shape.begin(), m_shape.end(), 1, std::multiplies<size_t>()));
	 
	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] =  expr(i,j); 
		} 
	}
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator+= (const RankTwoTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] += expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator+= (const CpxRankTwoTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] += expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator-= (const RankTwoTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] -= expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator-= (const CpxRankTwoTensorExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] -= expr(i,j); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator*= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] *= expr.eval(); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator*= (const CpxScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] *= expr.eval(); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator/= (const ScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] /= expr.eval(); 
		} 
	} 
	
	return *this;
}

template<typename Expr>
CpxRankTwoTensor& CpxRankTwoTensor::operator/= (const CpxScalarExpression< Expr >& expr)
{ 
	if (getShape() != expr.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); }

	for (size_t i=0; i<m_shape[0];++i) 
	{ 
		for (size_t j=0; j<m_shape[1];++j) 
		{ 
			m_core[flatIndex(i,j)] /= expr.eval(); 
		} 
	} 
	
	return *this;
}

#endif // RANK_TWO_TENSOR_TPP
