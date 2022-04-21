/*
 * CoefUnaryOP.hpp
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

#ifndef COEF_UNARY_OP_HPP
#define COEF_UNARY_OP_HPP

#include <LightFEM/Expression/LinAlg/UnaryExpression.hpp>

template<UnaryOp Op, typename Expr> class CoefUnaryOp {};

template<typename Expr> class CoefUnaryOp<UnaryOp::MINUS, Expr>
{
public:
	CoefUnaryOp(Expr&& expr) : m_expr(expr) {}
	
	inline double getCoef(const size_t globId) const { return -m_expr.getCoef(globId); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

////////////////////////////////////////////////////////////////////////

template<UnaryOp Op, typename Expr> class CpxCoefUnaryOp {};

template<typename Expr> class CpxCoefUnaryOp<UnaryOp::MINUS, Expr>
{
public:
	CpxCoefUnaryOp(Expr&& expr) : m_expr(expr) {}

	inline std::complex< double > getCoef(const size_t globId) const { return -m_expr.getCoef(globId); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

template<typename Expr> class CpxCoefUnaryOp<UnaryOp::CONJ, Expr>
{
public:
	CpxCoefUnaryOp(Expr&& expr) : m_expr(expr) {}

	inline std::complex< double > getCoef(const size_t globId) const { return std::conj(m_expr.getCoef(globId)); }
private:
	typename RefTypeSelector<Expr>::Type  m_expr;
};

#endif // COEF_UNARY_OP_HPP
