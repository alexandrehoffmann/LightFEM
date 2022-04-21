/*
 * RankTwoTensorBinaryExpression.hpp
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

#ifndef RANK_TWO_TENSOR_BINARY_EXPRESSION_HPP
#define RANK_TWO_TENSOR_BINARY_EXPRESSION_HPP

#include <LightFEM/Expression/LinAlg/RefTypeSelector.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensorExpression.hpp>

template<> struct BinaryOpType<BinaryOp::SUM,  ExprType::RK2_TENSOR, ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::SUB,  ExprType::RK2_TENSOR, ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::SCALAR,     ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::PROD, ExprType::RK2_TENSOR, ExprType::SCALAR    > { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DDOT, ExprType::RK4_TENSOR, ExprType::RK2_TENSOR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };
template<> struct BinaryOpType<BinaryOp::DIV,  ExprType::RK2_TENSOR,     ExprType::SCALAR> { static constexpr ExprType Type = ExprType::RK2_TENSOR; static constexpr bool isDefined = true; };

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> : public RankTwoTensorExpression< BinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 2> getShape() const { return m_lhs.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_lhs(i,j) + m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> : public RankTwoTensorExpression< BinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> >
{
public:
	BinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 2> getShape() const { return m_lhs.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_lhs(i,j) - m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename TensorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK2_TENSOR, TensorExpr> : public RankTwoTensorExpression< BinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK2_TENSOR, TensorExpr> >
{
public:
	BinaryExpression(ScalarExpr&& scalar, TensorExpr&& tensor) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 2> getShape() const { return m_tensor.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_tensor(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename ScalarExpr, typename TensorExpr>
class BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public RankTwoTensorExpression< BinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 2> getShape() const { return m_tensor.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_tensor(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename Rk4TensorExpr, typename Rk2TensorExpr>
class BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr> : public RankTwoTensorExpression< BinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr> >
{
public:
	BinaryExpression(Rk4TensorExpr&& rk4tensor, Rk2TensorExpr&& rk2tensor);
public:
	inline std::array<size_t, 2> getShape() const { return m_shape; }
public:
	inline double operator()(const size_t i, const size_t j) const;
private:
	typename RefTypeSelector<Rk4TensorExpr>::Type m_rk4tensor;
	typename RefTypeSelector<Rk2TensorExpr>::Type m_rk2tensor;
	std::array<size_t, 2> m_shape;
};

template<typename ScalarExpr, typename TensorExpr>
class BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public RankTwoTensorExpression< BinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	BinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 2> getShape() const { return m_tensor.getShape(); }
public:
	inline double operator()(const size_t i, const size_t j) const { return m_tensor(i,j) / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

////////////////////////////////////////////////////////////////////////

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> : public CpxRankTwoTensorExpression< CpxBinaryExpression<BinaryOp::SUM, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> >
{
public:
	static constexpr ExprType Type = ExprType::RK2_TENSOR;
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 2> getShape() const { return m_lhs.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_lhs(i,j) + m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr>
class CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> : public CpxRankTwoTensorExpression< CpxBinaryExpression<BinaryOp::SUB, ExprType::RK2_TENSOR, LeftExpr, ExprType::RK2_TENSOR, RightExpr> >
{
public:
	CpxBinaryExpression(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { if (m_lhs.getShape() != m_rhs.getShape()) { throw std::invalid_argument("lhs and rhs must share the same shape"); } }
public:
	inline std::array<size_t, 2> getShape() const { return m_lhs.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_lhs(i,j) - m_rhs(i,j); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename ScalarExpr, typename TensorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK2_TENSOR, TensorExpr> : public CpxRankTwoTensorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::SCALAR, ScalarExpr, ExprType::RK2_TENSOR, TensorExpr> >
{
public:
	CpxBinaryExpression(ScalarExpr&& scalar, TensorExpr&& tensor) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 2> getShape() const { return m_tensor.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_tensor(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename ScalarExpr, typename TensorExpr>
class CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public CpxRankTwoTensorExpression< CpxBinaryExpression<BinaryOp::PROD, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 2> getShape() const { return m_tensor.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_scalar.eval()*m_tensor(i,j); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

template<typename Rk4TensorExpr, typename Rk2TensorExpr>
class CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr> : public CpxRankTwoTensorExpression< CpxBinaryExpression<BinaryOp::DDOT, ExprType::RK4_TENSOR, Rk4TensorExpr, ExprType::RK2_TENSOR, Rk2TensorExpr> >
{
public:
	CpxBinaryExpression(Rk4TensorExpr&& rk4tensor, Rk2TensorExpr&& rk2tensor);
public:
	inline std::array<size_t, 2> getShape() const { return m_shape; }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const;
private:
	typename RefTypeSelector<Rk4TensorExpr>::Type m_rk4tensor;
	typename RefTypeSelector<Rk2TensorExpr>::Type m_rk2tensor;
	std::array<size_t, 2> m_shape;
};

template<typename ScalarExpr, typename TensorExpr>
class CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> : public CpxRankTwoTensorExpression< CpxBinaryExpression<BinaryOp::DIV, ExprType::RK2_TENSOR, TensorExpr, ExprType::SCALAR, ScalarExpr> >
{
public:
	CpxBinaryExpression(TensorExpr&& tensor, ScalarExpr&& scalar) : m_scalar(std::forward<ScalarExpr>(scalar)), m_tensor(std::forward<TensorExpr>(tensor)) { }
public:
	inline std::array<size_t, 2> getShape() const { return m_tensor.getShape(); }
public:
	inline std::complex< double > operator()(const size_t i, const size_t j) const { return m_tensor(i,j) / m_scalar.eval(); }
private:
	typename RefTypeSelector<ScalarExpr>::Type m_scalar;
	typename RefTypeSelector<TensorExpr>::Type m_tensor;
};

#include <LightFEM/Expression/LinAlg/RankTwoTensorBinaryExpression.tpp>

#endif // RANK_TWO_TENSOR_BINARY_EXPRESSION_HPP
