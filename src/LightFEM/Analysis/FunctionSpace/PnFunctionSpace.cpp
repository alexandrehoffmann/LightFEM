/*
 * PnFunctionSpace.cpp
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

#include <LightFEM/Analysis/FunctionSpace/PnFunctionSpace.hpp>

#include <LightFEM/Analysis/Measure/GaussLobattoQuadrature.hpp>

#include <LightFEM/Mesh/Node.hpp>

#include <numeric>

PnFunctionSpace::PnFunctionSpace(const Mesh *mesh, const size_t order) :
	FunctionSpace(mesh)
{
	/////// init 1d local interpolation function

	m_interpNodes1D = GaussLobattoQuadrature::getNodes(order+1);

	const size_t NInterpNodes = m_interpNodes1D.size();

	/////// init local interpolation function

	m_basisFunctions.resize(NInterpNodes*NInterpNodes);
	m_gradXibasisFunctions.resize(NInterpNodes*NInterpNodes);

	m_isLocIdOnBoundary.resize(4*NInterpNodes*NInterpNodes);

	for (size_t i=0;i<NInterpNodes;++i)
	{
		for (size_t j=0;j<NInterpNodes;++j)
		{
			m_basisFunctions[j+i*NInterpNodes] = [this,i,j](const double xi1, const double xi2) -> Scalar
			{
				return Scalar(l(i,xi1)*l(j,xi2));
			};
			m_gradXibasisFunctions[j+i*NInterpNodes] = [this,i,j](const double xi1, const double xi2) -> Vector
			{
				return Vector({ dl(i,xi1)*l(j,xi2), l(i,xi1)*dl(j,xi2) });
			};
			m_isLocIdOnBoundary[index2d(j+i*NInterpNodes, Element::Boundary::TOP)] = (j == NInterpNodes-1);
			m_isLocIdOnBoundary[index2d(j+i*NInterpNodes, Element::Boundary::BOTTTOM)] = (j == 0);
			m_isLocIdOnBoundary[index2d(j+i*NInterpNodes, Element::Boundary::LEFT)] = (i == 0);
			m_isLocIdOnBoundary[index2d(j+i*NInterpNodes, Element::Boundary::RIGHT)] = (i == NInterpNodes-1);
		}
	}

	/////// init global numbering

	const size_t N = NInterpNodes;

	std::vector< NodeWorld > nodesWorld(m_mesh->getNElem()*NInterpNodes*NInterpNodes);
	std::vector< size_t > idx(m_mesh->getNElem()*NInterpNodes*NInterpNodes);
	size_t maxIdx = 0;

	m_interpNodes.resize(N*N);

	for (size_t e=0;e<m_mesh->getNElem();++e)
	{
		for (size_t i=0;i<N;++i)
		{
			for (size_t j=0;j<N;++j)
			{
				m_interpNodes[j+i*N] = NodeRef(m_interpNodes1D[i], m_interpNodes1D[j]);
				nodesWorld[index2d(e,j+i*N)] = m_mesh->getElem(e)->getXworld(m_interpNodes[j+i*N]);

				if (fabs(nodesWorld[index2d(e,j+i*N)].x) < epsilon ) { nodesWorld[index2d(e,j+i*N)].x = 0.0; }
				if (fabs(nodesWorld[index2d(e,j+i*N)].y) < epsilon ) { nodesWorld[index2d(e,j+i*N)].y = 0.0; }

				idx[index2d(e,j+i*N)] = maxIdx;
				++maxIdx;
			}
		}
	}
	std::stable_sort(std::begin(idx), std::end(idx), [&nodesWorld](const size_t i1, const size_t i2) -> bool
	{
		if (fabs(nodesWorld[i1].x - nodesWorld[i2].x) > epsilon)
		{
			return nodesWorld[i1].x < nodesWorld[i2].x;
		}
		return nodesWorld[i1].y < nodesWorld[i2].y;
	});

	m_globalIds.resize(nodesWorld.size());
	m_globalIds[idx[0]] = 0;
	
	for (size_t i=1;i<idx.size();++i)
	{
		if ( eq(nodesWorld[idx[i]], nodesWorld[idx[i-1]], epsilon) ) { m_globalIds[idx[i]] = m_globalIds[idx[i-1]]; }
		else { m_globalIds[idx[i]] = m_globalIds[idx[i-1]] + 1; }
	}

	m_nDofs = m_globalIds.max() + 1;

	initDiscretization();
	initIsIdOnBoundary();
}

double PnFunctionSpace::getHMin() const
{
	const size_t N = m_interpNodes1D.size();
	
	std::vector< NodeWorld > nodesWorld(getNBasisFunction());
	std::valarray< bool > isComputed(false, getNBasisFunction());
	
	for (size_t e=0;e<m_mesh->getNElem();++e)
	{
		for (size_t i=0;i<N;++i)
		{
			for (size_t j=0;j<N;++j)
			{
				const size_t globId = getGlobalId(e,j+i*N);
				
				if (not isComputed[globId])
				{
					nodesWorld[globId] = m_mesh->getElem(e)->getXworld(m_interpNodes[j+i*N]);
					if (fabs(nodesWorld[globId].x) < epsilon ) { nodesWorld[globId].x = 0.0; }
					if (fabs(nodesWorld[globId].y) < epsilon ) { nodesWorld[globId].y = 0.0; }
					
					isComputed[globId] = true;
				}
			}
		}
	}
	
	double hmin = dist(nodesWorld[0], nodesWorld[1]);
	
	for (size_t i=0;i<getNBasisFunction();++i)
	{
		for (size_t j=i+1;j<getNBasisFunction();++j)
		{
			hmin = std::min(hmin, dist(nodesWorld[i], nodesWorld[j]));
		}
	}
	
	return hmin;
}

std::vector< double > PnFunctionSpace::getHMinPerElement() const
{	
	const size_t N = m_interpNodes1D.size();
	
	std::vector< double > hminPerElement(m_mesh->getNElem());
	
	std::vector< NodeWorld > nodesWorldPerElement(N*N);
	
	for (size_t e=0;e<m_mesh->getNElem();++e)
	{
		for (size_t i=0;i<N;++i)
		{
			for (size_t j=0;j<N;++j)
			{
				const size_t globId = getGlobalId(e,j+i*N);
				
				nodesWorldPerElement[j+i*N] = m_mesh->getElem(e)->getXworld(m_interpNodes[j+i*N]);
				if (fabs(nodesWorldPerElement[j+i*N].x) < epsilon ) { nodesWorldPerElement[globId].x = 0.0; }
				if (fabs(nodesWorldPerElement[j+i*N].y) < epsilon ) { nodesWorldPerElement[globId].y = 0.0; }
			}
		}
		
		hminPerElement[e] = dist(nodesWorldPerElement[0], nodesWorldPerElement[1]);
		
		for (size_t i=0;i<nodesWorldPerElement.size();++i)
		{
			for (size_t j=i+1;j<nodesWorldPerElement.size();++j)
			{
				hminPerElement[e] = std::min(hminPerElement[e], dist(nodesWorldPerElement[i], nodesWorldPerElement[j]));
			}
		}
		
	}
	
	return hminPerElement;
}

double PnFunctionSpace::l(const size_t i, const double t) const
{
	double lt = 1.0;

	for (size_t j=0;j<m_interpNodes1D.size();++j)
	{
		if (j != i) { lt *= (t - m_interpNodes1D[j]) / (m_interpNodes1D[i] - m_interpNodes1D[j]); }
	}

	return lt;
}

double PnFunctionSpace::dl(const size_t i, const double t) const
{
	double dlt = 0.0;

	for (size_t k=0;k<m_interpNodes1D.size();++k)
	{
		if (k != i)
		{
			double prod = 1.0 / (m_interpNodes1D[i] - m_interpNodes1D[k]);
			for (size_t j=0;j<m_interpNodes1D.size();++j)
			{
				if (j != i and j != k) { prod *= (t - m_interpNodes1D[j]) / (m_interpNodes1D[i] - m_interpNodes1D[j]);  }
			}
			dlt += prod;
		}
	}

	return dlt;
}
