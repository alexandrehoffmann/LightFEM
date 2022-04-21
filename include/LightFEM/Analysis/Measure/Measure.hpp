/*
 * Measure.hpp
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

#ifndef MEASURE_HPP
#define MEASURE_HPP

#include <vector>
#include <LightFEM/Mesh/Mesh.hpp>
#include <LightFEM/Mesh/Node.hpp>

struct NodeAndWeight
{
	NodeAndWeight(const double _xi1=0.0, const double _xi2=0.0, const double _wi=0.0) : xi1(_xi1), xi2(_xi2), wi(_wi) {}
	NodeAndWeight(const NodeRef& Xi, const double _wi) : xi1(Xi.xi1), xi2(Xi.xi2), wi(_wi) {}
	NodeAndWeight(const NodeAndWeight& other) : xi1(other.xi1), xi2(other.xi2), wi(other.wi) {}

	double xi1;
	double xi2;
	double wi;
};

struct BoundaryNodeAndWeight
{
	BoundaryNodeAndWeight(const double _t=0.0, const double _topWi=0.0, const double _bottomWi=0.0, const double _leftWi=0.0, const double _rightWi=0.0) 
		: t(_t), topWi(_topWi), bottomWi(_bottomWi), leftWi(_leftWi), rightWi(_rightWi) {}
	BoundaryNodeAndWeight(const BoundaryNodeAndWeight& other) : t(other.t), topWi(other.topWi), bottomWi(other.bottomWi), leftWi(other.leftWi), rightWi(other.rightWi) {}

	double t;
	double topWi;
	double bottomWi;
	double leftWi;
	double rightWi;
};

class Measure
{
public:
	Measure(const Mesh* mesh) : m_mesh(mesh), m_nodesAndWeightsPerElement(mesh->getNElem()), m_boundaryNodesAndWeightsPerElement(mesh->getNElem()) {}
public:
	inline const std::vector< NodeAndWeight >&         getNodesAndWeights         (const size_t e) const { return m_nodesAndWeightsPerElement[e]; }
	inline const std::vector< BoundaryNodeAndWeight >& getBoundaryNodesAndWeights (const size_t e) const { return m_boundaryNodesAndWeightsPerElement[e]; }
	inline const Mesh* getMesh() const { return m_mesh; }
protected:
	const Mesh* m_mesh;
	std::vector< std::vector< NodeAndWeight > > m_nodesAndWeightsPerElement;
	std::vector< std::vector< BoundaryNodeAndWeight > > m_boundaryNodesAndWeightsPerElement;
};

////////////////////////////////////////////////////////////////////////

struct CpxNodeAndWeight
{
	CpxNodeAndWeight(const double _xi1=0.0, const double _xi2=0.0, const std::complex< double > _wi=0.0) : xi1(_xi1), xi2(_xi2), wi(_wi) {}
	CpxNodeAndWeight(const NodeRef& Xi, const std::complex< double > _wi) : xi1(Xi.xi1), xi2(Xi.xi2), wi(_wi) {}
	CpxNodeAndWeight(const CpxNodeAndWeight& other) : xi1(other.xi1), xi2(other.xi2), wi(other.wi) {}

	double xi1;
	double xi2;
	std::complex< double > wi;
};

struct CpxBoundaryNodeAndWeight
{
	CpxBoundaryNodeAndWeight(const double _t=0.0, const std::complex< double > _topWi=0.0, const std::complex< double > _bottomWi=0.0, const std::complex< double > _leftWi=0.0, const std::complex< double > _rightWi=0.0) 
		: t(_t), topWi(_topWi), bottomWi(_bottomWi), leftWi(_leftWi), rightWi(_rightWi) {}
	CpxBoundaryNodeAndWeight(const CpxBoundaryNodeAndWeight& other) : t(other.t), topWi(other.topWi), bottomWi(other.bottomWi), leftWi(other.leftWi), rightWi(other.rightWi) {}

	double t;
	std::complex< double > topWi;
	std::complex< double > bottomWi;
	std::complex< double > leftWi;
	std::complex< double > rightWi;
};

class CpxMeasure
{
public:
	CpxMeasure(const Mesh* mesh) : m_mesh(mesh), m_nodesAndWeightsPerElement(mesh->getNElem()), m_boundaryNodesAndWeightsPerElement(mesh->getNElem()) {}
public:
	inline const std::vector< CpxNodeAndWeight >&         getNodesAndWeights         (const size_t e) const { return m_nodesAndWeightsPerElement[e]; }
	inline const std::vector< CpxBoundaryNodeAndWeight >& getBoundaryNodesAndWeights (const size_t e) const { return m_boundaryNodesAndWeightsPerElement[e]; }
	inline const Mesh* getMesh() const { return m_mesh; }
protected:
	const Mesh* m_mesh;
	std::vector< std::vector< CpxNodeAndWeight > > m_nodesAndWeightsPerElement;
	std::vector< std::vector< CpxBoundaryNodeAndWeight > > m_boundaryNodesAndWeightsPerElement;
};

#endif // MEASURE_HPP
