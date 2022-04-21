/*
 * RankFourTensorBinaryExpression.hpp
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

#ifndef RANK_FOUR_TENSOR_BINARY_EXPRESSION_HPP
#define RANK_FOUR_TENSOR_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensorExpression.hpp>

template<> struct BinaryOpType<BinaryOp::SUM,  ExprType::RK4_TENSOR, ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::SUB,  ExprType::RK4_TENSOR, ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::SCALAR,     ExprType::RK4_TENSOR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::RK4_TENSOR, ExprType::SCALAR    > { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DIV,  ExprType::RK4_TENSOR,     ExprType::SCALAR> { static constexpr ExprType Type = ExprType::RK4_TENSOR; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> : public RankFourTensorExpression< BinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 4> getShape() const { return m_lhs.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_lhs(i,j,k,l) + m_rhs(i,j,k,l); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> : public RankFourTensorExpression< BinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 4> getShape() const { return m_lhs.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_lhs(i,j,k,l) - m_rhs(i,j,k,l); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename TensorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK4_TENSOR, TensorExpr> : public RankFourTensorExpression< BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK4_TENSOR, TensorExpr> >
{
public:
	BinaryExpression(ScalarExpr&& scalar, TensorExpr&& tensor) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 4> getShape() const { return m_tensor.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_scalar.eval()*m_tensor(i,j,k,l); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename ScalarExpr, typename TensorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public RankFourTensorExpression< BinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 4> getShape() const { return m_tensor.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_scalar.eval()*m_tensor(i,j,k,l); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename ScalarExpr, typename TensorExpr>
class BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public RankFourTensorExpression< BinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 4> getShape() const { return m_tensor.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_tensor(i,j,k,l) / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> : public CpxRankFourTensorExpression< CpxBinaryExpression<BinaryOp::SUM, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::RK4_TENSOR;
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 4> getShape() const { return m_lhs.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_lhs(i,j,k,l) + m_rhs(i,j,k,l); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> : public CpxRankFourTensorExpression< CpxBinaryExpression<BinaryOp::SUB, ExprType::RK4_TENSOR, LeftExpr, ExprType::RK4_TENSOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 4> getShape() const { return m_lhs.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_lhs(i,j,k,l) - m_rhs(i,j,k,l); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename TensorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK4_TENSOR, TensorExpr> : public CpxRankFourTensorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK4_TENSOR, TensorExpr> >
{
public:
	CpxBinaryExpression(ScalarExpr&& scalar, TensorExpr&& tensor) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 4> getShape() const { return m_tensor.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_scalar.eval()*m_tensor(i,j,k,l); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename ScalarExpr, typename TensorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public CpxRankFourTensorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 4> getShape() const { return m_tensor.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_scalar.eval()*m_tensor(i,j,k,l); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename ScalarExpr, typename TensorExpr>
class CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public CpxRankFourTensorExpression< CpxBinaryExpression<BinaryOp::DIV, ExprType::RK4_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 4> getShape() const { return m_tensor.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j, const size_t k, const size_t l) const { return m_tensor(i,j,k,l) / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

#endif // RANK_FOUR_TENSOR_BINARY_EXPRESSION_HPP
