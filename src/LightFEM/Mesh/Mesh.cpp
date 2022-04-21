/*
 * Mesh.cpp
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
 
#include <LightFEM/Mesh/Mesh.hpp>

#include <iostream>

size_t Mesh::getElemId(const Element* element) const
{
	for (size_t e=0;e<m_elements.size();++e)
	{
		if (m_elements[e] == element) { return e; }
	}
	return m_elements.size();
}

int Mesh::getDomainId(const std::string &domaineName) const
{
	std::vector< std::string >::const_iterator itr = std::find(m_domainsName.cbegin(), m_domainsName.cend(), domaineName);
	if (itr != m_domainsName.cend()) { return std::distance(m_domainsName.cbegin(), itr); }

	return -1;
}

int Mesh::getBoundaryId(const std::string &boundaryName) const
{
	std::vector< std::string >::const_iterator itr = std::find(m_boundariesName.cbegin(), m_boundariesName.cend(), boundaryName);
	if (itr != m_boundariesName.cend()) { return std::distance(m_boundariesName.cbegin(), itr); }

	return -1;
}
