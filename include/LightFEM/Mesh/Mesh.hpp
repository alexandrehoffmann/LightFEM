/*
 * Mesh.hpp
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

#ifndef MESH_HPP
#define MESH_HPP

#include <LightFEM/Mesh/Element.hpp>
#include <LightFEM/Mesh/BoundaryElement.hpp>

#include <vector>
#include <map>

class Mesh
{
public:
	virtual ~Mesh() { for (size_t e=0;e<getNElem();++e) { delete m_elements[e]; } for (size_t be=0;be<getNBoundaryElem();++be) { delete m_boundaryElements[be]; } }
public:
	const Element* getElem(const size_t e) const { return m_elements[e]; }
	size_t getNElem() const { return m_elements.size(); }
	const BoundaryElement* getBoundaryElem(const size_t be) const { return m_boundaryElements[be]; }
	size_t getNBoundaryElem() const { return m_boundaryElements.size(); }

	size_t getElemId(const Element* element) const;

	inline size_t getElemIdFromBoundaryElemId(const size_t be) const {
		if(getBoundaryElem(be)->getParentElement() != getElem(m_boundaryElementIdToElementId[be])) { throw std::logic_error("m_boundaryElementIdToElementId isn't properly filled."); }
		return m_boundaryElementIdToElementId[be]; } // getBoundaryElem(be).getParent() == getElem(getElemIdFromBoundaryElemId(be))

	int getDomainId(const std::string& domaineName) const; // -1 means the domain do not exists
	int getBoundaryId(const std::string& boundaryName) const; // -1 means the domain do not exists
protected:
	std::vector< Element* >          m_elements;
	std::vector< BoundaryElement* >  m_boundaryElements;
	std::vector< NodeWorld >         m_nodes;

	std::vector< size_t > m_boundaryElementIdToElementId; // getBoundaryElem(be).getParent() == getElem(m_boundaryElementIdToElementId[be])

	std::vector<std::string> m_domainsName;
	std::vector<std::string> m_boundariesName;
};

#endif // MESH_HPP
