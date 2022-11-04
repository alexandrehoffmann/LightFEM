/*
 * ScalarExpression.hpp
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

#ifndef SCALAR_EXPRESSION_HPP
#define SCALAR_EXPRESSION_HPP

#include <complex>

template <typename Expr>
class ScalarExpression 
{
public:
	inline double eval() const { return static_cast<Expr const&>(*this).eval(); }
	
	inline operator double() const { return static_cast<Expr const&>(*this).eval(); }
};

////////////////////////////////////////////////////////////////////////

template <typename Expr>
class CpxScalarExpression 
{
public:
	inline std::complex< double > eval() const { return static_cast<Expr const&>(*this).eval(); }
	
	inline operator std::complex< double >() const { return static_cast<Expr const&>(*this).eval(); }
	inline double real() const { return static_cast<Expr const&>(*this).real(); }
	inline double imag() const { return static_cast<Expr const&>(*this).imag(); }
};

////////////////////////////////////////////////////////////////////////

template<typename T, typename Expr>  constexpr decltype(auto) operator+=(std::complex<T>& lhs, const CpxScalarExpression<Expr>& rhs) { return lhs += std::complex<double>(rhs); }
template<typename T, typename Expr>  constexpr decltype(auto) operator-=(std::complex<T>& lhs, const CpxScalarExpression<Expr>& rhs) { return lhs -= std::complex<double>(rhs); }
template<typename T, typename Expr>  constexpr decltype(auto) operator*=(std::complex<T>& lhs, const CpxScalarExpression<Expr>& rhs) { return lhs *= std::complex<double>(rhs); }
template<typename T, typename Expr>  constexpr decltype(auto) operator/=(std::complex<T>& lhs, const CpxScalarExpression<Expr>& rhs) { return lhs /= std::complex<double>(rhs); }

template<typename Expr> constexpr decltype(auto) operator+(const CpxScalarExpression<Expr>& lhs, const double rhs) { return std::complex<double>(lhs) + rhs; }
template<typename Expr> constexpr decltype(auto) operator-(const CpxScalarExpression<Expr>& lhs, const double rhs) { return std::complex<double>(lhs) - rhs; }
template<typename Expr> constexpr decltype(auto) operator*(const CpxScalarExpression<Expr>& lhs, const double rhs) { return std::complex<double>(lhs) * rhs; }
template<typename Expr> constexpr decltype(auto) operator/(const CpxScalarExpression<Expr>& lhs, const double rhs) { return std::complex<double>(lhs) / rhs; }

template<typename Expr> constexpr decltype(auto) operator+(const double lhs, const CpxScalarExpression<Expr>& rhs) { return lhs + std::complex<double>(rhs); }
template<typename Expr> constexpr decltype(auto) operator-(const double lhs, const CpxScalarExpression<Expr>& rhs) { return lhs - std::complex<double>(rhs); }
template<typename Expr> constexpr decltype(auto) operator*(const double lhs, const CpxScalarExpression<Expr>& rhs) { return lhs * std::complex<double>(rhs); }
template<typename Expr> constexpr decltype(auto) operator/(const double lhs, const CpxScalarExpression<Expr>& rhs) { return lhs / std::complex<double>(rhs); }

#endif // SCALAR_EXPRESSION_HPP
