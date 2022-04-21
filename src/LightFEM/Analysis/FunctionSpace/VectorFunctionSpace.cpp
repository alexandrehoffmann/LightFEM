/*
 * VectorFunctionSpace.cpp
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

#include <LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp>

#include <LightFEM/Mesh/Mesh.hpp>
#include <LightFEM/Analysis/VectorTestFunction.hpp>
#include <LightFEM/Analysis/VectorTrialFunction.hpp>

VectorFunctionSpace::VectorFunctionSpace(const std::initializer_list<FunctionSpace *> functionSpaces) :
	m_functionSpaces(functionSpaces),
	m_zero(Element::getNxiNd(), Vector(2, 0.0)),
	m_gradZero(Element::getNxiNd(), Matrix(2, functionSpaces.size(), 0.0))
{
	for (size_t i=1;i<m_functionSpaces.size();++i)
	{
		if (m_functionSpaces[i]->getMesh() != m_functionSpaces[i-1]->getMesh()) { throw std::invalid_argument("all function space must be defined for the same mesh"); }
	}

	m_nBasisFunctionPerElement = std::accumulate(std::begin(m_functionSpaces), std::end(m_functionSpaces), 0, [](int accumulator, const FunctionSpace * functionSpace)
	{
		return accumulator + functionSpace->getNBasisFunctionPerElement();
	});
	
	const size_t NInterpNodes = std::accumulate(std::begin(m_functionSpaces), std::end(m_functionSpaces), 0, [](int accumulator, const FunctionSpace * functionSpace)
	{
		return accumulator + functionSpace->getNLocalInterpolationNodes();
	});

	/////// init locId to fspaceId and fspace's locId
	
	m_interpNodes.resize(NInterpNodes);
	m_interpFuncLocId.resize(NInterpNodes);
	
	std::vector< std::pair< size_t, size_t > > locIdTofspaceAndLocId(m_nBasisFunctionPerElement);
	std::vector< std::pair< size_t, size_t > > interpFuncLocIdTofspaceAndinterpFuncLocId(NInterpNodes);
	std::vector< std::vector< size_t > > invLocIdTofspaceAndLocId(m_functionSpaces.size());
	
	size_t k=0;
	size_t kp=0;
	for (size_t fspaceId=0;fspaceId<m_functionSpaces.size();++fspaceId)
	{
		invLocIdTofspaceAndLocId[fspaceId].resize(m_functionSpaces[fspaceId]->getNBasisFunctionPerElement());
		for (size_t locId=0;locId<m_functionSpaces[fspaceId]->getNBasisFunctionPerElement();++locId)
		{
			locIdTofspaceAndLocId[k] = std::make_pair(fspaceId, locId);
			invLocIdTofspaceAndLocId[fspaceId][locId] = k;
			++k;
		}
		for (size_t interpFuncLocId=0;interpFuncLocId<m_functionSpaces[fspaceId]->getNLocalInterpolationNodes();++interpFuncLocId)
		{
			interpFuncLocIdTofspaceAndinterpFuncLocId[kp] = std::make_pair(fspaceId, interpFuncLocId);
			++kp;
		}
	}

	/////// init local interpolation function

	m_discBasisFunctions.resize(m_nBasisFunctionPerElement);
	m_discGradXiBasisFunctions.resize(m_nBasisFunctionPerElement);

	Vector zero(m_functionSpaces.size(), 0.0);
	Matrix gradZero(2, m_functionSpaces.size(), 0.0);

	for (size_t locId=0;locId<m_nBasisFunctionPerElement;++locId)
	{
		m_discBasisFunctions[locId].resize(Element::getNxiNd(), zero);
		m_discGradXiBasisFunctions[locId].resize(Element::getNxiNd(), gradZero);
	}

	for (size_t locId=0;locId<m_nBasisFunctionPerElement;++locId)
	{
		const auto [fspaceId, fspaceLocId] = locIdTofspaceAndLocId[locId];

		const std::vector< Scalar >& u = m_functionSpaces[fspaceId]->getDiscBasisFunction(fspaceLocId);
		const std::vector< Vector >& grad_u = m_functionSpaces[fspaceId]->getDiscGradXiBasisFunction(fspaceLocId);
		for (size_t i=0;i<Element::getNxi();++i)
		{
			for (size_t j=0;j<Element::getNxi();++j)
			{
				m_discBasisFunctions[locId][Element::index2d(i, j)][fspaceId] = u[Element::index2d(i, j)];
				m_discGradXiBasisFunctions[locId][Element::index2d(i, j)](0, fspaceId) = grad_u[Element::index2d(i, j)][0];
				m_discGradXiBasisFunctions[locId][Element::index2d(i, j)](1, fspaceId) = grad_u[Element::index2d(i, j)][1];
			}
		}
	}
	
	for (size_t interpFuncLocId=0;interpFuncLocId<m_interpFuncLocId.size();++interpFuncLocId)
	{
		const auto [fspaceId, fspaceInterpFuncLocId] = interpFuncLocIdTofspaceAndinterpFuncLocId[interpFuncLocId];
		const size_t fspaceLocId = m_functionSpaces[fspaceId]->getLocalInterpolationFunctionId(fspaceInterpFuncLocId);
		m_interpNodes[interpFuncLocId] = m_functionSpaces[fspaceId]->getLocalInterpolationNode(fspaceInterpFuncLocId);
		m_interpFuncLocId[interpFuncLocId] = invLocIdTofspaceAndLocId[fspaceId][fspaceLocId];
	}

	/////// init global numbering

	std::vector< size_t > globalNumberingPadding(m_functionSpaces.size());
	globalNumberingPadding[0] = 0;
	for (size_t fspaceId=1;fspaceId<m_functionSpaces.size();++fspaceId)
	{
		globalNumberingPadding[fspaceId] = globalNumberingPadding[fspaceId - 1] + m_functionSpaces[fspaceId-1]->getNBasisFunction();
	}

	std::vector< std::pair<size_t, size_t> > globalIdToFspaceIdAndGlobalIde(getMesh()->getNElem()*m_nBasisFunctionPerElement);
	m_globalIds.resize(getMesh()->getNElem()*m_nBasisFunctionPerElement);

	for (size_t e=0;e<getMesh()->getNElem();++e)
	{
		for (size_t locId=0;locId<m_nBasisFunctionPerElement;++locId)
		{
			const auto [fspaceId, fspaceLocId] = locIdTofspaceAndLocId[locId];
			m_globalIds[index2d(e, locId)] = globalNumberingPadding[fspaceId] + m_functionSpaces[fspaceId]->getGlobalId(e, fspaceLocId);
			globalIdToFspaceIdAndGlobalIde[m_globalIds[index2d(e, locId)]] = std::make_pair(fspaceId, m_functionSpaces[fspaceId]->getGlobalId(e, fspaceLocId));
		}
	}
	m_nDofs = m_globalIds.max() + 1;

	/////// init boundaries

	m_isIdOnBoundary.resize( getNBasisFunction() );
	for (size_t globId=0;globId<getNBasisFunction();++globId)
	{
		const auto [fspaceId, fspaceGlobId] = globalIdToFspaceIdAndGlobalIde[globId];
		m_isIdOnBoundary[globId] = m_functionSpaces[fspaceId]->isIdOnBoundary(fspaceGlobId);
	}
}

VectorTrialFunction VectorFunctionSpace::getTrialFunction(const size_t e, const int locId) const
{
	return VectorTrialFunction(this, getMesh()->getElem(e), locId);
}

VectorTestFunction VectorFunctionSpace::getTestFunction(const size_t e, const int locId) const
{
	return VectorTestFunction(this, getMesh()->getElem(e), locId);
}
