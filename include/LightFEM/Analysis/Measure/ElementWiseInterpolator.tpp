/*
 * ElementWiseInterpolator.tpp
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

#ifndef INTERPOLATOR_TPP
#define INTERPOLATOR_TPP

#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.hpp>

template<typename Expr>
ElementWiseInterpolator::ElementWiseInterpolator(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& f) :
	m_fi(Element::getNxiNd()),
	m_dist_xi1(Element::getNxiNd()),
	m_dist_xi2(Element::getNxiNd()),
	m_wi1(Element::getNxiNd()),
	m_wi2(Element::getNxiNd())
{
	if (s_wi.size() != Element::getNxi()) { initWeights(); }

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_fi[k] = f[k];
	}
}

////////////////////////////////////////////////////////////////////////

template<typename Expr>
CpxElementWiseInterpolator::CpxElementWiseInterpolator(const ElementWiseFunctionExpression<ExprType::SCALAR, Expr>& f) :
	m_fi(Element::getNxiNd()),
	m_dist_xi1(Element::getNxiNd()),
	m_dist_xi2(Element::getNxiNd()),
	m_wi1(Element::getNxiNd()),
	m_wi2(Element::getNxiNd())
{
	if (s_wi.size() != Element::getNxi()) { initWeights(); }

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_fi[k] = f[k];
	}
}

template<typename Expr>
CpxElementWiseInterpolator::CpxElementWiseInterpolator(const CpxElementWiseFunctionExpression<ExprType::SCALAR, Expr>& f) :
	m_fi(Element::getNxiNd()),
	m_dist_xi1(Element::getNxiNd()),
	m_dist_xi2(Element::getNxiNd()),
	m_wi1(Element::getNxiNd()),
	m_wi2(Element::getNxiNd())
{
	if (s_wi.size() != Element::getNxi()) { initWeights(); }

	for (size_t k=0;k<Element::getNxiNd();++k)
	{
		m_fi[k] = f[k];
	}
}

#endif // INTERPOLATOR_TPP
