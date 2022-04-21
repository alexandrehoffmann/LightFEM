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
	inline double operator() (const NodeRef& Xi) { return eval(Xi.xi1, Xi.xi2); }
	inline double operator() (const double xi1, const double xi2) { return eval(xi1, xi2); }
private:
	inline void initDists(const double xi1, const double xi2) { for (size_t i=0;i<Element::getNxi();++i) {m_dist_xi1[i] = dist(xi1, Element::get_xi(i)); m_dist_xi2[i] = dist(xi2, Element::get_xi(i)); } }
	std::pair<double, double > l();
	double eval(const double xi1, const double xi2);
private:
	std::valarray< double > m_fi;
	std::valarray< double > m_dist_xi1;
	std::valarray< double > m_dist_xi2;
	std::valarray< double > m_wi1;
	std::valarray< double > m_wi2;
private:
	static void initWeights();

	inline static double dist(const double t, const double ti) { return (fabs(t - ti) > 1.0e-16) ? t - ti : 1.0e-16; }
private:
	static std::valarray< double > s_wi;
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
	inline std::complex< double > operator() (const NodeRef& Xi) { return eval(Xi.xi1, Xi.xi2); }
	inline std::complex< double > operator() (const double xi1, const double xi2) { return eval(xi1, xi2); }
private:
	inline void initDists(const double xi1, const double xi2) { for (size_t i=0;i<Element::getNxi();++i) {m_dist_xi1[i] = dist(xi1, Element::get_xi(i)); m_dist_xi2[i] = dist(xi2, Element::get_xi(i)); } }
	std::pair<double, double > l();
	std::complex< double > eval(const double xi1, const double xi2);
private:
	std::valarray< std::complex< double > > m_fi;
	std::valarray< double > m_dist_xi1;
	std::valarray< double > m_dist_xi2;
	std::valarray< double > m_wi1;
	std::valarray< double > m_wi2;
private:
	static void initWeights();

	inline static double dist(const double t, const double ti) { return (fabs(t - ti) > 1.0e-16) ? t - ti : 1.0e-16; }
private:
	static std::valarray< double > s_wi;
};

#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.tpp>

#endif // INTERPOLATOR_HPP
