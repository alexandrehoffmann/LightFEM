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

//// https://github.com/gregvw/barycentric-lagrange/blob/master/lagrange.cpp

std::valarray< double > ElementWiseInterpolator::s_wi;

double ElementWiseInterpolator::eval(const double xi1, const double xi2)
{
	initDists(xi1, xi2);

	for (size_t i=0;i<Element::getNxi();++i)
	{
		m_wi1[i] = s_wi[i] / m_dist_xi1[i];
		m_wi2[i] = s_wi[i] / m_dist_xi2[i];
	}

	const auto [l1, l2] = l();

	double value(0.0);
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			value += m_fi[Element::index2d(i,j)]*m_wi1[i]*m_wi2[j];
		}
	}
	return value*l1*l2;
}

std::pair<double, double> ElementWiseInterpolator::l()
{
	std::pair< double, double > l = std::make_pair(1.0, 1.0);
	for(size_t i=0;i<Element::getNxi();i++)
	{
		std::get<0>(l) *= m_dist_xi1[i];
		std::get<1>(l) *= m_dist_xi2[i];
	}
	return l;
}

void ElementWiseInterpolator::initWeights()
{
	s_wi.resize(Element::getNxi());
	for (size_t i=0;i<Element::getNxi();++i)
	{
		double prod = 1.0;
		for (size_t j=0;j<Element::getNxi();++j)
		{
			if (i != j)
			{
				prod *= (Element::get_xi(i) - Element::get_xi(j));
			}
		}
		s_wi[i] = 1.0 / prod;
	}
}

////////////////////////////////////////////////////////////////////////

//// https://github.com/gregvw/barycentric-lagrange/blob/master/lagrange.cpp

std::valarray< double > CpxElementWiseInterpolator::s_wi;

std::complex< double > CpxElementWiseInterpolator::eval(const double xi1, const double xi2)
{
	initDists(xi1, xi2);

	for (size_t i=0;i<Element::getNxi();++i)
	{
		m_wi1[i] = s_wi[i] / m_dist_xi1[i];
		m_wi2[i] = s_wi[i] / m_dist_xi2[i];
	}

	const auto [l1, l2] = l();

	std::complex< double > value(0.0, 0.0);
	for (size_t i=0;i<Element::getNxi();++i)
	{
		for (size_t j=0;j<Element::getNxi();++j)
		{
			value += m_fi[Element::index2d(i,j)]*m_wi1[i]*m_wi2[j];
		}
	}
	return value*l1*l2;
}

std::pair<double, double> CpxElementWiseInterpolator::l()
{
	std::pair< double, double > l = std::make_pair(1.0, 1.0);
	for(size_t i=0;i<Element::getNxi();i++)
	{
		std::get<0>(l) *= m_dist_xi1[i];
		std::get<1>(l) *= m_dist_xi2[i];
	}
	return l;
}

void CpxElementWiseInterpolator::initWeights()
{
	s_wi.resize(Element::getNxi());
	for (size_t i=0;i<Element::getNxi();++i)
	{
		double prod = 1.0;
		for (size_t j=0;j<Element::getNxi();++j)
		{
			if (i != j)
			{
				prod *= (Element::get_xi(i) - Element::get_xi(j));
			}
		}
		s_wi[i] = 1.0 / prod;
	}
}
