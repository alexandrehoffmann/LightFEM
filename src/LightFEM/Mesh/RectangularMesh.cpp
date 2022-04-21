/*
 * RectangularMesh.cpp
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

#include <LightFEM/Mesh/RectangularMesh.hpp>

#include <iostream>

RectangularMesh::RectangularMesh(const double xmin, const double xmax, const double ymin, const double ymax, const size_t Nx, const size_t Ny) :
	m_xmin(xmin),
	m_xmax(xmax),
	m_ymin(ymin),
	m_ymax(ymax),
	m_Nx(Nx),
	m_Ny(Ny)
{
	std::vector< double > x = linspace(m_xmin, m_xmax, m_Nx);
	std::vector< double > y = linspace(m_ymin, m_ymax, m_Ny);

	m_domainsName.push_back("Omega");

	m_boundariesName.push_back("Bottom");
	m_boundariesName.push_back("Left");
	m_boundariesName.push_back("Top");
	m_boundariesName.push_back("Right");
	m_boundariesName.push_back("dOmega");

	const int omegaId = getDomainId("Omega");
	const int bottomId = getBoundaryId("Bottom");
	const int leftId = getBoundaryId("Left");
	const int topId = getBoundaryId("Top");
	const int rightId = getBoundaryId("Right");
	const int domegaId = getBoundaryId("dOmega");

	m_nodes.resize(x.size()*y.size());

	for (size_t i=0;i<m_Nx;++i) { for (size_t j=0;j<m_Ny;++j)
	{
		m_nodes[index2d(i,j)].x = x[i];
		m_nodes[index2d(i,j)].y = y[j];
	}}

	for (size_t i=0;i<m_Nx-1;++i) { for (size_t j=0;j<m_Ny-1;++j)
	{
		std::array<std::array<NodeWorld*,2>, 2> X;
		for (size_t ii=0;ii<2;++ii) { for (size_t jj=0;jj<2;++jj)
		{
			X[ii][jj] = &m_nodes[index2d(ii+i,jj+j)];
		}}
		m_elements.push_back( new Element(X, {omegaId}) );

		if (j == 0)      { m_boundaryElements.push_back( new BoundaryElement(m_elements.back(), Element::Boundary::BOTTTOM, {bottomId, domegaId} ) ); m_boundaryElementIdToElementId.push_back(m_elements.size()-1);}
		if (i == 0)      { m_boundaryElements.push_back( new BoundaryElement(m_elements.back(), Element::Boundary::LEFT   , {leftId, domegaId} ) ); m_boundaryElementIdToElementId.push_back(m_elements.size()-1);}
		if (j == m_Ny-2) { m_boundaryElements.push_back( new BoundaryElement(m_elements.back(), Element::Boundary::TOP    , {topId, domegaId} ) ); m_boundaryElementIdToElementId.push_back(m_elements.size()-1);}
		if (i == m_Nx-2) { m_boundaryElements.push_back( new BoundaryElement(m_elements.back(), Element::Boundary::RIGHT  , {rightId, domegaId} ) ); m_boundaryElementIdToElementId.push_back(m_elements.size()-1);}

	}}
}
