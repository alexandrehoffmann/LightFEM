/*
 * BoundaryElement.hpp
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

#ifndef BOUNDARY_ELEMENT_HPP
#define BOUNDARY_ELEMENT_HPP

#include <LightFEM/Mesh/Element.hpp>
#include <vector>

class BoundaryElement
{
public:
	BoundaryElement(const Element* parent_element = nullptr, Element::Boundary boundary = Element::Boundary::BOTTTOM) : m_parent_element(parent_element), m_boundary(boundary) {}
	BoundaryElement(const Element* parent_element, Element::Boundary boundary, std::initializer_list< int > ids) : m_parent_element(parent_element), m_boundary(boundary), m_ids(ids) {}
	BoundaryElement(const Element* parent_element, Element::Boundary boundary, const std::vector< int >& ids) : m_parent_element(parent_element), m_boundary(boundary), m_ids(ids) {}
public:
	inline const Element*    getParentElement() const { return m_parent_element; }
	inline Element::Boundary getBoundary()      const { return m_boundary; }
	inline bool isInDomain(const int domainId)  const { return std::find(m_ids.cbegin(), m_ids.cend(), domainId) !=  m_ids.cend(); }
private:
	const Element* m_parent_element;
	Element::Boundary m_boundary;
	std::vector< int > m_ids;
};

#endif // BOUNDARY_ELEMENT_HPP
