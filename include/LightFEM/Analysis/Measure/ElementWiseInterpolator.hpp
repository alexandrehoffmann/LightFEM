/*
 * ElementWiseInterpolator.hpp
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

#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

#include <valarray>
#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>

class ElementWiseInterpolator
{
public:
	template<typename Expr>
	ElementWiseInterpolator(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr >& f);
public:
	inline double operator() (const NodeRef& Xi) const { return eval(Xi.xi1, Xi.xi2); }
	inline double operator() (const double xi1, const double xi2) const { return eval(xi1, xi2); }
private:
	double eval(const double xi1, const double xi2) const;
	
	static double lagrange_1d(const double xi, const size_t i);
	
	std::valarray< double > m_fi;
};

////////////////////////////////////////////////////////////////////////

class CpxElementWiseInterpolator
{
public:
	template<typename Expr>
	CpxElementWiseInterpolator(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr >& f);
	template<typename Expr>
	CpxElementWiseInterpolator(const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr >& f);
public:
	inline std::complex< double > operator() (const NodeRef& Xi) const { return eval(Xi.xi1, Xi.xi2); }
	inline std::complex< double > operator() (const double xi1, const double xi2) const { return eval(xi1, xi2); }
private:
	std::complex< double > eval(const double xi1, const double xi2) const;
	
	static double lagrange_1d(const double xi1, const size_t i);
	
	std::valarray< std::complex< double > > m_fi;
};

#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.tpp>

#endif // INTERPOLATOR_HPP
