/*
 * SamplingFunction.cpp
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

#include <LightFEM/Analysis/Measure/SamplingFunction.hpp>

#include <list>

double SamplingFunction::epsilon = 1.0e-11;
double CpxSamplingFunction::epsilon = 1.0e-11;

SamplingFunction::SamplingFunction(const Mesh* mesh, const NodeWorld& X) :
	Measure(mesh)
{
	std::vector< double > _weights{1.0};
	std::vector< NodeWorld > _X{X};
	init(mesh, _X, _weights);
}

SamplingFunction::SamplingFunction(const Mesh* mesh, const NodeWorld& X, const double weight) :
	Measure(mesh)
{
	std::vector< double > _weights{weight};
	std::vector< NodeWorld > _X{X};
	init(mesh, _X, _weights);
}

SamplingFunction::SamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X) :
	Measure(mesh)
{
	std::vector< double > weights(X.size(), 1.0);
	init(mesh, X, weights);
}

SamplingFunction::SamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< double >& weights) :
	Measure(mesh)
{
	if (X.size() != weights.size()) { throw std::invalid_argument("X and weights must share the same size"); }
	
	init(mesh, X, weights);
}

//void SamplingFunction::init(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< double >& weights)
//{
	//for (size_t i=0;i<X.size();++i)
	//{
		//bool elemFound = false;
		//for (size_t e=0;e<mesh->getNElem();++e)
		//{
			//const auto [isInElem, Xi] = mesh->getElem(e)->getXRef(X[i], epsilon);
			//if (isInElem)
			//{
				//m_nodesAndWeightsPerElement[e].push_back( NodeAndWeight(Xi, weights[i]) );
				//elemFound = true;
				//break;
			//}
		//}
		//if (not elemFound) { throw std::invalid_argument("X[" + std::to_string(i) + "] (" + std::to_string(X[i].x) + ", " + std::to_string(X[i].y) + ") cannot be found on mesh."); }
	//}
//}
void SamplingFunction::init(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< double >& weights)
{
	for (size_t i=0;i<X.size();++i)
	{
		std::list< std::pair<size_t, NodeRef> > elemIdAndNode;
		for (size_t e=0;e<mesh->getNElem();++e)
		{
			const auto [isInElem, Xi] = mesh->getElem(e)->getXRef(X[i], epsilon);
			if (isInElem)
			{
				elemIdAndNode.push_back(std::make_pair(e, Xi));
			}
		}
		if (elemIdAndNode.size() == 0) { throw std::invalid_argument("X[" + std::to_string(i) + "] (" + std::to_string(X[i].x) + ", " + std::to_string(X[i].y) + ") cannot be found on mesh."); }
		for (const auto& [e, Xi] : elemIdAndNode)
		{
			m_nodesAndWeightsPerElement[e].push_back( NodeAndWeight(Xi, weights[i]) );
		}
	}
}

////////////////////////////////////////////////////////////////////////

CpxSamplingFunction::CpxSamplingFunction(const Mesh* mesh, const NodeWorld& X) :
	CpxMeasure(mesh)
{
	std::vector< std::complex< double > > _weights{std::complex<double>(1.0,0.0)};
	std::vector< NodeWorld > _X{X};
	init(mesh, _X, _weights);
}

CpxSamplingFunction::CpxSamplingFunction(const Mesh* mesh, const NodeWorld& X, const std::complex< double >  weight) :
	CpxMeasure(mesh)
{
	std::vector< std::complex< double >  > _weights{weight};
	std::vector< NodeWorld > _X{X};
	init(mesh, _X, _weights);
}

CpxSamplingFunction::CpxSamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X) :
	CpxMeasure(mesh)
{
	std::vector< std::complex< double >  > weights(X.size(), std::complex< double > (1.0, 0.0));
	init(mesh, X, weights);
}

CpxSamplingFunction::CpxSamplingFunction(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< std::complex< double >  >& weights) :
	CpxMeasure(mesh)
{
	if (X.size() != weights.size()) { throw std::invalid_argument("X and weights must share the same size"); }

	init(mesh, X, weights);
}

void CpxSamplingFunction::init(const Mesh* mesh, const std::vector< NodeWorld >& X, const std::vector< std::complex< double > >& weights)
{
	for (size_t i=0;i<X.size();++i)
	{
		bool elemFound = false;
		for (size_t e=0;e<mesh->getNElem();++e)
		{
			const auto [isInElem, Xi] = mesh->getElem(e)->getXRef(X[i], epsilon);
			if (isInElem)
			{
				m_nodesAndWeightsPerElement[e].push_back( CpxNodeAndWeight(Xi, weights[i]) );
				elemFound = true;
				break;
			}
		}
		if (not elemFound) { throw std::invalid_argument("X[" + std::to_string(i) + "] (" + std::to_string(X[i].x) + ", " + std::to_string(X[i].y) + ") cannot be found on mesh."); }
	}
}
