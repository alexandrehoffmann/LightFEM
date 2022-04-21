/*
 * RankFourTensor.hpp
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

#ifndef RANK_FOUR_TENSOR_HPP
#define RANK_FOUR_TENSOR_HPP

#include <valarray>
#include <complex>

#include <LightFEM/Expression/LinAlg/RankFourTensorExpression.hpp>
#include <LightFEM/Expression/LinAlg/ScalarExpression.hpp>

class RankFourTensor : public RankFourTensorExpression< RankFourTensor >
{
public:
	RankFourTensor(const size_t shape1=1, const size_t shape2=1, const size_t shape3=1, const size_t shape4=1, const double value=0.0);
	template<typename Expr>
	RankFourTensor(const RankFourTensorExpression< Expr >& expr);
public:
	template<typename Expr> RankFourTensor& operator=  (const RankFourTensorExpression< Expr >& expr);
	template<typename Expr> RankFourTensor& operator+= (const RankFourTensorExpression< Expr >& expr);
	template<typename Expr> RankFourTensor& operator-= (const RankFourTensorExpression< Expr >& expr);
	template<typename Expr> RankFourTensor& operator*= (const ScalarExpression< Expr >& expr);
	template<typename Expr> RankFourTensor& operator/= (const ScalarExpression< Expr >& expr);
public:
	inline std::array<size_t, 4> getShape() const { return m_shape; }
public:
	inline double  operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_core[flatIndex(i,j,k,l)]; }
	inline double& operator()(const size_t i, const size_t j, const size_t k, const size_t l)       { return m_core[flatIndex(i,j,k,l)]; }
private:
	inline size_t flatIndex(const size_t i, const size_t j, const size_t k, const size_t l) const { return l + m_shape[3]*(k + m_shape[2]*(j + i*m_shape[1])); }
public:
	std::array< size_t, 4 > m_shape;
	std::valarray< double > m_core;
};

////////////////////////////////////////////////////////////////////////

class CpxRankFourTensor : public CpxRankFourTensorExpression< CpxRankFourTensor >
{
public:
	CpxRankFourTensor(const size_t shape1=1, const size_t shape2=1, const size_t shape3=1, const size_t shape4=1, const std::complex< double > value = 0.0);
	template<typename Expr>
	CpxRankFourTensor(const RankFourTensorExpression< Expr >& expr);
	template<typename Expr>
	CpxRankFourTensor(const CpxRankFourTensorExpression< Expr >& expr);
public:
	template<typename Expr> CpxRankFourTensor& operator=  (const RankFourTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator=  (const CpxRankFourTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator+= (const RankFourTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator+= (const CpxRankFourTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator-= (const RankFourTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator-= (const CpxRankFourTensorExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator*= (const ScalarExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator*= (const CpxScalarExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator/= (const ScalarExpression< Expr >& expr);
	template<typename Expr> CpxRankFourTensor& operator/= (const CpxScalarExpression< Expr >& expr);
public:
	inline std::array<size_t, 4> getShape() const { return m_shape; }
public:
	inline std::complex< double >  operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_core[flatIndex(i,j,k,l)]; }
	inline std::complex< double >& operator()(const size_t i, const size_t j, const size_t k, const size_t l)       { return m_core[flatIndex(i,j,k,l)]; }
private:
	inline size_t flatIndex(const size_t i, const size_t j, const size_t k, const size_t l) const { return l + m_shape[3]*(k + m_shape[2]*(j + i*m_shape[1])); }
private:
	std::array< size_t, 4 > m_shape;
	std::valarray< std::complex< double > > m_core;
};

#include <LightFEM/Expression/LinAlg/RankFourTensor.tpp>

#endif // RANK_FOUR_TENSOR_HPP
