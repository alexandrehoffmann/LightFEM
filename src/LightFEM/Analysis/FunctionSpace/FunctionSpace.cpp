/*
 * FunctionSpace.cpp
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

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>

#include <LightFEM/Mesh/Mesh.hpp>
#include <LightFEM/Analysis/TrialFunction.hpp>
#include <LightFEM/Analysis/TestFunction.hpp>

FunctionSpace::FunctionSpace(const Mesh* mesh) : 
	m_mesh(mesh),
	m_zero(Element::getNxiNd(), Scalar(0.0)),
	m_gradZero(Element::getNxiNd(), Vector(2, 0.0))
{
	
}

TrialFunction FunctionSpace::getTrialFunction(const size_t e, const int locId) const
{
	return TrialFunction(this, m_mesh->getElem(e), locId);
}

TestFunction FunctionSpace::getTestFunction(const size_t e, const int locId) const
{
	return TestFunction(this, m_mesh->getElem(e), locId);
}

void FunctionSpace::initDiscretization()
{
	const size_t N = getNBasisFunctionPerElement();

	m_discBasisFunctions.resize( N );
	m_discGradXiBasisFunctions.resize( N );

	for (size_t locId=0;locId<N;++locId)
	{
		m_discBasisFunctions[locId].resize(Element::getNxiNd());
		m_discGradXiBasisFunctions[locId].resize(Element::getNxiNd());

		for (size_t i=0;i<Element::getNxi();++i)
		{
			for (size_t j=0;j<Element::getNxi();++j)
			{
				m_discBasisFunctions[locId][Element::index2d(i, j)]       = evalDiscBasisFunction(locId, Element::get_xi(i), Element::get_xi(j));
				m_discGradXiBasisFunctions[locId][Element::index2d(i, j)] = evalDiscGradXiBasisFunction(locId, Element::get_xi(i), Element::get_xi(j));
			}
		}
	}
}

void FunctionSpace::initIsIdOnBoundary()
{
	m_isIdOnBoundary.resize( getNBasisFunction() );

	for (size_t be=0;be<m_mesh->getNBoundaryElem();++be)
	{
		const size_t e = m_mesh->getElemIdFromBoundaryElemId(be);
		const Element::Boundary b = m_mesh->getBoundaryElem(be)->getBoundary();

		for (size_t locId=0;locId<getNBasisFunctionPerElement();++locId)
		{
			m_isIdOnBoundary[getGlobalId(e, locId)] = isLocIdOnBoundary(locId, b);
		}
	}
}

