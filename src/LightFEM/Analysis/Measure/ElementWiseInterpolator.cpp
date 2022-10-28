/*
 * ElementWiseInterpolator.cpp
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

#include <LightFEM/Analysis/Measure/ElementWiseInterpolator.hpp>

double ElementWiseInterpolator::eval(const double xi1, const double xi2) const
{
	double value = 0.0;
	
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			const double lx = lagrange_1d(xi1, i);
			const double ly = lagrange_1d(xi2, j);
			
			value += m_fi[Element::index2d(i,j)]*lx*ly;
		}
	}
	
	return value;
}

double ElementWiseInterpolator::lagrange_1d(const double xi, const size_t i)
{
	double li = 1.0;
	
	for (size_t j=0;j<Element::getNxi();++j) { if (i != j)
	{
		li *= (xi - Element::get_xi(j)) / (Element::get_xi(i) - Element::get_xi(j));
	}}
	
	return li;
}

////////////////////////////////////////////////////////////////////////

std::complex< double > CpxElementWiseInterpolator::eval(const double xi1, const double xi2) const
{
	std::complex< double > value = 0.0;
	
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			const double lx = lagrange_1d(xi1, i);
			const double ly = lagrange_1d(xi2, j);
			
			value += m_fi[Element::index2d(i,j)]*lx*ly;
		}
	}
	
	return value;
}

double CpxElementWiseInterpolator::lagrange_1d(const double xi, const size_t i)
{
	double li = 1.0;
	
	for (size_t j=0;j<Element::getNxi();++j) { if (i != j)
	{
		li *= (xi - Element::get_xi(j)) / (Element::get_xi(i) - Element::get_xi(j));
	}}
	
	return li;
}
