/*
 * CoefBinaryOp.hpp
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

#ifndef COEF_BINARY_OP_HPP
#define COEF_BINARY_OP_HPP

#include <type_traits>

#include <LightFEM/Expression/LinAlg/BinaryExpression.hpp>
#include <LightFEM/Expression/FiniteElementFunction/FiniteElementFunctionExpression.hpp>

#include <LightFEM/Expression/Function/Const.hpp>

template<BinaryOp Op, typename LeftExpr, typename RightExpr> class CoefBinaryOp {};

template<typename LeftExpr, typename RightExpr> class CoefBinaryOp<BinaryOp::SUM, LeftExpr, RightExpr>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { 
		if (m_lhs.getFunctionSpace() != m_rhs.getFunctionSpace()) { throw std::invalid_argument("lhs and rhs must belong to the same function space."); } }
		
	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) + m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr> class CoefBinaryOp<BinaryOp::SUB, LeftExpr, RightExpr>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) { 
		if (m_lhs.getFunctionSpace() != m_rhs.getFunctionSpace()) { throw std::invalid_argument("lhs and rhs must belong to the same function space."); } }
		
	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) - m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename RightExpr> class CoefBinaryOp<BinaryOp::PROD, const Const&, RightExpr>
{
public:
	typedef typename std::decay_t<RightExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(const Const& lhs, RightExpr&& rhs) : m_lhs(lhs), m_rhs(std::forward<RightExpr>(rhs)) {}
		
	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getValue() * m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_rhs.getFunctionSpace(); }
private:
	const Const& m_lhs;
	typename RefTypeSelector<RightExpr>::Type   m_rhs;
};

template<typename RightExpr> class CoefBinaryOp<BinaryOp::PROD, Const, RightExpr>
{
public:
	typedef typename std::decay_t<RightExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(Const&& lhs, RightExpr&& rhs) : m_lhs(std::forward<Const>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}

	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getValue() * m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_rhs.getFunctionSpace(); }
private:
	Const m_lhs;
	typename RefTypeSelector<RightExpr>::Type   m_rhs;
};

template<typename LeftExpr> class CoefBinaryOp<BinaryOp::PROD, LeftExpr, const Const&>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(LeftExpr&& lhs, const Const& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(rhs) {}
		
	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) * m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	const Const& m_rhs;
};

template<typename LeftExpr> class CoefBinaryOp<BinaryOp::PROD, LeftExpr, Const>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(LeftExpr&& lhs, Const&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<Const>(rhs)) {}

	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) * m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	Const m_rhs;
};

template<typename LeftExpr> class CoefBinaryOp<BinaryOp::DIV, LeftExpr, const Const&>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(LeftExpr&& lhs, const Const& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(rhs) {}
		
	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) / m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	const Const& m_rhs;
};

template<typename LeftExpr> class CoefBinaryOp<BinaryOp::DIV, LeftExpr, Const>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CoefBinaryOp(LeftExpr&& lhs, Const&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<Const>(rhs)) {}

	inline const Scalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) / m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	Const m_rhs;
};

////////////////////////////////////////////////////////////////////////

template<BinaryOp Op, typename LeftExpr, typename RightExpr> class CpxCoefBinaryOp {};

template<typename LeftExpr, typename RightExpr> class CpxCoefBinaryOp<BinaryOp::SUM, LeftExpr, RightExpr>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getFunctionSpace() != m_rhs.getFunctionSpace()) { throw std::invalid_argument("lhs and rhs must belong to the same function space."); } }

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) + m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename LeftExpr, typename RightExpr> class CpxCoefBinaryOp<BinaryOp::SUB, LeftExpr, RightExpr>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, RightExpr&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {
		if (m_lhs.getFunctionSpace() != m_rhs.getFunctionSpace()) { throw std::invalid_argument("lhs and rhs must belong to the same function space."); } }

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) - m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type  m_lhs;
	typename RefTypeSelector<RightExpr>::Type m_rhs;
};

template<typename RightExpr> class CpxCoefBinaryOp<BinaryOp::PROD, const Const&, RightExpr>
{
public:
	typedef typename std::decay_t<RightExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(const Const& lhs, RightExpr&& rhs) : m_lhs(lhs), m_rhs(std::forward<RightExpr>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getValue() * m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_rhs.getFunctionSpace(); }
private:
	const Const& m_lhs;
	typename RefTypeSelector<RightExpr>::Type   m_rhs;
};

template<typename RightExpr> class CpxCoefBinaryOp<BinaryOp::PROD, const CpxConst&, RightExpr>
{
public:
	typedef typename std::decay_t<RightExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(const CpxConst& lhs, RightExpr&& rhs) : m_lhs(lhs), m_rhs(std::forward<RightExpr>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getValue() * m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_rhs.getFunctionSpace(); }
private:
	const CpxConst& m_lhs;
	typename RefTypeSelector<RightExpr>::Type   m_rhs;
};

template<typename RightExpr> class CpxCoefBinaryOp<BinaryOp::PROD, Const, RightExpr>
{
public:
	typedef typename std::decay_t<RightExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(Const&& lhs, RightExpr&& rhs) : m_lhs(std::forward<Const>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getValue() * m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_rhs.getFunctionSpace(); }
private:
	Const m_lhs;
	typename RefTypeSelector<RightExpr>::Type   m_rhs;
};

template<typename RightExpr> class CpxCoefBinaryOp<BinaryOp::PROD, CpxConst, RightExpr>
{
public:
	typedef typename std::decay_t<RightExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(CpxConst&& lhs, RightExpr&& rhs) : m_lhs(std::forward<CpxConst>(lhs)), m_rhs(std::forward<RightExpr>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getValue() * m_rhs.getCoef(globId); }

	inline const FSpaceType* getFunctionSpace() const { return m_rhs.getFunctionSpace(); }
private:
	CpxConst m_lhs;
	typename RefTypeSelector<RightExpr>::Type   m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::PROD, LeftExpr, const Const&>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, const Const& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(rhs) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) * m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	const Const& m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::PROD, LeftExpr, const CpxConst&>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, const CpxConst& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(rhs) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) * m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	const CpxConst& m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::PROD, LeftExpr, Const>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, Const&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<Const>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) * m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	Const m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::PROD, LeftExpr, CpxConst>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, CpxConst&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<CpxConst>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) * m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	CpxConst m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::DIV, LeftExpr, const Const&>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, const Const& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(rhs) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) / m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	const Const& m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::DIV, LeftExpr, const CpxConst&>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, const CpxConst& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(rhs) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) / m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	const CpxConst& m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::DIV, LeftExpr, Const>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, Const&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<Const>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) / m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	Const m_rhs;
};

template<typename LeftExpr> class CpxCoefBinaryOp<BinaryOp::DIV, LeftExpr, CpxConst>
{
public:
	typedef typename std::decay_t<LeftExpr>::FSpaceType FSpaceType;
public:
	CpxCoefBinaryOp(LeftExpr&& lhs, CpxConst&& rhs) : m_lhs(std::forward<LeftExpr>(lhs)), m_rhs(std::forward<CpxConst>(rhs)) {}

	inline const CpxScalar& getCoef(const size_t globId) const { return m_lhs.getCoef(globId) / m_rhs.getValue(); }

	inline const FSpaceType* getFunctionSpace() const { return m_lhs.getFunctionSpace(); }
private:
	typename RefTypeSelector<LeftExpr>::Type    m_lhs;
	CpxConst m_rhs;
};

#endif // COEF_BINARY_OP_HPP
