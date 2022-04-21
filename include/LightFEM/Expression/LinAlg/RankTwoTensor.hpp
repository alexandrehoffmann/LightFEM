/*
 * RankTwoTensor.hpp
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

#ifndef RANK_TWO_TENSOR_HPP
#define RANK_TWO_TENSOR_HPP

#include <valarray>
#include <complex>

#include <LightFEM/Expression/LinAlg/RankTwoTensorExpression.hpp>
#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>

class RankTwoTensor : public RankTwoTensorExpression< RankTwoTensor >
{
public:
	RankTwoTensor(const size_t shape1=1, const size_t shape2=1, const double value=0);
	template<typename Expr>
	RankTwoTensor(const RankTwoTensorExpression< Expr >& expr);
public:
	template<typename Expr> RankTwoTensor& operator=  (const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr> RankTwoTensor& operator+= (const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr> RankTwoTensor& operator-= (const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr> RankTwoTensor& operator*= (const ScalarExpression< Expr >& expr);
	template<typename Expr> RankTwoTensor& operator/= (const ScalarExpression< Expr >& expr);
public:
	inline std::array<size_t, 2> getShape() const { return m_shape; }
public:
	inline double  operator()(const size_t i, const size_t j) const { return m_core[flatIndex(i,j)]; }
	inline double& operator()(const size_t i, const size_t j)       { return m_core[flatIndex(i,j)]; }
private:
	inline size_t flatIndex(const size_t i, const size_t j) const { return i*m_shape[1] + j; }
private:
	std::array< size_t, 2 > m_shape;
	std::valarray< double > m_core;
};

////////////////////////////////////////////////////////////////////////

class CpxRankTwoTensor : public CpxRankTwoTensorExpression< CpxRankTwoTensor >
{
public:
	CpxRankTwoTensor(const size_t shape1=1, const size_t shape2=1, const std::complex< double > value = 0.0);
	template<typename Expr>
	CpxRankTwoTensor(const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr>
	CpxRankTwoTensor(const CpxRankTwoTensorExpression< Expr >& expr);
public:
public:
	template<typename Expr> CpxRankTwoTensor& operator=  (const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator=  (const CpxRankTwoTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator+= (const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator+= (const CpxRankTwoTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator-= (const RankTwoTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator-= (const CpxRankTwoTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator*= (const ScalarExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator*= (const CpxScalarExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator/= (const ScalarExpression< Expr >& expr);
	template<typename Expr> CpxRankTwoTensor& operator/= (const CpxScalarExpression< Expr >& expr);
public:
	inline std::array<size_t, 2> getShape() const { return m_shape; }
public:
	inline std::complex< double >  operator()(const size_t i, const size_t j) const { return m_core[flatIndex(i,j)]; }
	inline std::complex< double >& operator()(const size_t i, const size_t j)       { return m_core[flatIndex(i,j)]; }
private:
	inline size_t flatIndex(const size_t i, const size_t j) const { return i*m_shape[1] + j; }
private:
	std::array< size_t, 2 > m_shape;
	std::valarray< std::complex< double > > m_core;
};

#include <LightFEM/Expression/LinAlg/RankTwoTensor.tpp>

#endif // RANK_TWO_TENSOR_HPP
