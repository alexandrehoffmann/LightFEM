/*
 * GmshMesh.cpp
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

#include <LightFEM/Mesh/GmshMesh.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

std::string read_line(std::ifstream& filestream)
{
    std::string line;
    std::getline(filestream, line);
    return line;
} 

bool findSection(std::ifstream& filestream, const std::string& searchname)
{
	filestream.clear();
	filestream.seekg(0, std::ios::beg);

	for(std::string line; std::getline( filestream, line ); )
	{
		if (line.find(searchname) != std::string::npos)
		{
			return true;
		}
	}
	return false;
} 

GmshMesh::GmshMesh(const std::string& filename)
{
	std::ifstream in(filename);
	if(!in.is_open())
	{
		std::cerr << "can't open " << filename << std::endl;
		exit(EXIT_FAILURE);
	}
	readPhysicalNames(in);
	readEntities(in);
	readNodes(in);
	readElements(in);

	initElements();
	initBoundaryElements();

	m_points.clear();
	m_curves.clear();
	m_surfaces.clear();
	m_entities.clear();
}

void GmshMesh::initElements()
{
	for (const Entity& entity : m_entities)
	{
		for (size_t e=0;e<entity.elements.size();++e)
		{
			if (entity.elements[e].size() == 4)
			{
				std::array<std::array<NodeWorld*,2>, 2> X;

				const int node00Tag = entity.elements[e][0];
				const int node10Tag = entity.elements[e][1];
				const int node11Tag = entity.elements[e][2];
				const int node01Tag = entity.elements[e][3];

				const size_t node00Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node00Tag));
				const size_t node10Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node10Tag));
				const size_t node11Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node11Tag));
				const size_t node01Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node01Tag));

				X[0][0] = &m_nodes[node00Index];
				X[1][0] = &m_nodes[node10Index];
				X[1][1] = &m_nodes[node11Index];
				X[0][1] = &m_nodes[node01Index];

				const size_t surfaceIdx = std::distance(m_surfaces.cbegin(), std::find_if(m_surfaces.cbegin(), m_surfaces.cend(), [&entity](const Surface& surface) { return surface.tag == entity.tag; }));

				std::vector< int > domainIds;
				for (const int physicalTag : m_surfaces[surfaceIdx].physicalTags)
				{
					domainIds.push_back( std::distance(m_domainsTag.cbegin(), std::find(m_domainsTag.cbegin(), m_domainsTag.cend(), physicalTag)) );
				}

				m_elements.push_back( new Element(X, domainIds));
			}
			else if (entity.elements[e].size() == 9)
			{
				std::array<std::array<NodeWorld*,3>, 3> X;

				const int node00Tag = entity.elements[e][0];
				const int node20Tag = entity.elements[e][1];
				const int node22Tag = entity.elements[e][2];
				const int node02Tag = entity.elements[e][3];

				const int node10Tag = entity.elements[e][4];
				const int node21Tag = entity.elements[e][5];
				const int node12Tag = entity.elements[e][6];
				const int node01Tag = entity.elements[e][7];

				const int node11Tag = entity.elements[e][8];

				const size_t node00Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node00Tag));
				const size_t node10Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node10Tag));
				const size_t node20Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node20Tag));

				const size_t node01Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node01Tag));
				const size_t node11Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node11Tag));
				const size_t node21Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node21Tag));

				const size_t node02Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node02Tag));
				const size_t node12Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node12Tag));
				const size_t node22Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node22Tag));

				X[0][0] = &m_nodes[node00Index];
				X[1][0] = &m_nodes[node10Index];
				X[2][0] = &m_nodes[node20Index];

				X[0][1] = &m_nodes[node01Index];
				X[1][1] = &m_nodes[node11Index];
				X[2][1] = &m_nodes[node21Index];

				X[0][2] = &m_nodes[node02Index];
				X[1][2] = &m_nodes[node12Index];
				X[2][2] = &m_nodes[node22Index];

				const size_t surfaceIdx = std::distance(m_surfaces.cbegin(), std::find_if(m_surfaces.cbegin(), m_surfaces.cend(), [&entity](const Surface& surface) { return surface.tag == entity.tag; }));

				std::vector< int > domainIds;
				for (const int physicalTag : m_surfaces[surfaceIdx].physicalTags)
				{
					domainIds.push_back( std::distance(m_domainsTag.cbegin(), std::find(m_domainsTag.cbegin(), m_domainsTag.cend(), physicalTag)) );
				}

				m_elements.push_back( new Element(X, domainIds));
			}
		}
	}
}

void GmshMesh::initBoundaryElements()
{
	for (const Entity& entity : m_entities)
	{
		for (size_t e=0;e<entity.elements.size();++e)
		{
			if (entity.elements[e].size() == 2)
			{
				const int node0Tag = entity.elements[e][0];
				const int node1Tag = entity.elements[e][1];

				const size_t node0Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node0Tag));
				const size_t node1Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node1Tag));

				std::vector< Element* >::const_iterator itTop = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index](const Element* elem)
				{
					return elem->getNode(0,2) == &m_nodes[node1Index] and elem->getNode(2,2) == &m_nodes[node0Index];
				});
				std::vector< Element* >::const_iterator itBottom = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index](const Element* elem)
				{
					return elem->getNode(0,0) == &m_nodes[node0Index] and elem->getNode(2,0) == &m_nodes[node1Index];
				});
				std::vector< Element* >::const_iterator itLeft = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index](const Element* elem)
				{
					return elem->getNode(0,0) == &m_nodes[node1Index] and elem->getNode(0,2) == &m_nodes[node0Index];
				});
				std::vector< Element* >::const_iterator itRight = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index](const Element* elem)
				{
					return elem->getNode(2,0) == &m_nodes[node0Index] and elem->getNode(2,2) == &m_nodes[node1Index];
				});

				const size_t curveIdx = std::distance(m_curves.cbegin(), std::find_if(m_curves.cbegin(), m_curves.cend(), [&entity](const Curve& curve) { return curve.tag == entity.tag; }));

				std::vector< int > boundaryIds;
				for (const int physicalTag : m_curves[curveIdx].physicalTags)
				{
					boundaryIds.push_back( std::distance(m_boundariesTag.cbegin(), std::find(m_boundariesTag.cbegin(), m_boundariesTag.cend(), physicalTag)) );
				}

				if (itTop != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itTop, Element::Boundary::TOP, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itTop));
				}
				if (itBottom != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itBottom, Element::Boundary::BOTTTOM, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itBottom));
				}
				if (itLeft != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itLeft, Element::Boundary::LEFT, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itLeft));
				}
				if (itRight != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itRight, Element::Boundary::RIGHT, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itRight));
				}
			}
			else if (entity.elements[e].size() == 3)
			{
				const int node0Tag = entity.elements[e][0];
				const int node1Tag = entity.elements[e][2];
				const int node2Tag = entity.elements[e][1];

				const size_t node0Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node0Tag));
				const size_t node1Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node1Tag));
				const size_t node2Index = std::distance(m_nodeTags.cbegin(), std::find(m_nodeTags.cbegin(), m_nodeTags.cend(), node2Tag));

				std::vector< Element* >::const_iterator itTop = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index, &node2Index](const Element* elem)
				{
					return elem->getNode(0,2) == &m_nodes[node2Index] and elem->getNode(1,2) == &m_nodes[node1Index] and elem->getNode(2,2) == &m_nodes[node0Index];
				});
				std::vector< Element* >::const_iterator itBottom = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index, &node2Index](const Element* elem)
				{
					return elem->getNode(0,0) == &m_nodes[node0Index] and elem->getNode(1,0) == &m_nodes[node1Index] and elem->getNode(2,0) == &m_nodes[node2Index];
				});
				std::vector< Element* >::const_iterator itLeft = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index, &node2Index](const Element* elem)
				{
					return elem->getNode(0,0) == &m_nodes[node2Index] and elem->getNode(0,1) == &m_nodes[node1Index] and elem->getNode(0,2) == &m_nodes[node0Index];
				});
				std::vector< Element* >::const_iterator itRight = std::find_if(m_elements.cbegin(), m_elements.cend(), [this, &node0Index, &node1Index, &node2Index](const Element* elem)
				{
					return elem->getNode(2,0) == &m_nodes[node0Index] and elem->getNode(2,1) == &m_nodes[node1Index] and elem->getNode(2,2) == &m_nodes[node2Index];
				});

				const size_t curveIdx = std::distance(m_curves.cbegin(), std::find_if(m_curves.cbegin(), m_curves.cend(), [&entity](const Curve& curve) { return curve.tag == entity.tag; }));

				std::vector< int > boundaryIds;
				for (const int physicalTag : m_curves[curveIdx].physicalTags)
				{
					boundaryIds.push_back( std::distance(m_boundariesTag.cbegin(), std::find(m_boundariesTag.cbegin(), m_boundariesTag.cend(), physicalTag)) );
				}

				if (itTop != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itTop, Element::Boundary::TOP, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itTop));
				}
				if (itBottom != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itBottom, Element::Boundary::BOTTTOM, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itBottom));
				}
				if (itLeft != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itLeft, Element::Boundary::LEFT, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itLeft));
				}
				if (itRight != m_elements.cend())
				{
					m_boundaryElements.push_back( new BoundaryElement(*itRight, Element::Boundary::RIGHT, boundaryIds) );
					m_boundaryElementIdToElementId.push_back(std::distance(m_elements.cbegin(), itRight));
				}
			}
		}
	}
}

void GmshMesh::readPhysicalNames(std::ifstream& in)
{
	if (findSection(in, "PhysicalNames"))
	{
		std::stringstream firstLinestream(read_line(in));
		
		size_t nPhysicalNames;
		firstLinestream >> nPhysicalNames;
		for (size_t i=0;i<nPhysicalNames;++i)
		{
			size_t dimension, physicalTag;
			std::string name;
			
			std::stringstream linestream(read_line(in));
			
			linestream >> dimension >> physicalTag >> name;

			name.erase(std::remove(name.begin(), name.end(), '"'), name.end());

			if (dimension == 1)      { m_boundariesName.push_back(name); m_boundariesTag.push_back(physicalTag); }
			else if (dimension == 2) { m_domainsName.push_back(name); m_domainsTag.push_back(physicalTag); }
		}
	}
}

void GmshMesh::readEntities(std::ifstream& in)
{
	if (findSection(in, "Entities"))
	{
		std::stringstream firstLinestream(read_line(in));
		size_t numPoints, numCurves, numSurfaces, numVolumes;
		firstLinestream >> numPoints >> numCurves >> numSurfaces >> numVolumes;
		
		m_points.resize(numPoints);	
		m_curves.resize(numCurves);
		m_surfaces.resize(numSurfaces);
		
		for (size_t i=0;i<numPoints;++i) 
		{ 
			std::stringstream linestream(read_line(in));
			size_t numPhysicalTags;
			linestream >> m_points[i].tag >> m_points[i].x >> m_points[i].y >> m_points[i].z >> numPhysicalTags;
			
			m_points[i].physicalTags.resize(numPhysicalTags);
			for (size_t tagId=0;i<numPhysicalTags;++tagId) { linestream >> m_points[i].physicalTags[tagId]; }
		} 
		for (size_t i=0;i<numCurves;++i)
		{
			std::stringstream linestream(read_line(in));
			size_t numPhysicalTags;
			size_t numBoundingPoints;
			linestream >> m_curves[i].tag >> m_curves[i].minX >> m_curves[i].minY >> m_curves[i].minZ >> m_curves[i].maxX >> m_curves[i].maxY >> m_curves[i].maxZ >> numPhysicalTags;

			m_curves[i].physicalTags.resize(numPhysicalTags);
			for (size_t tagId=0;tagId<numPhysicalTags;++tagId) { linestream >> m_curves[i].physicalTags[tagId]; } linestream >> numBoundingPoints;

			m_curves[i].pointTags.resize(numBoundingPoints);
			for (size_t tagId=0;tagId<numBoundingPoints;++tagId) { linestream >> m_curves[i].pointTags[tagId]; }
		}
		for (size_t i=0;i<numSurfaces;++i)
		{
			std::stringstream linestream(read_line(in));
			size_t numPhysicalTags;
			size_t numBoundingCurves;
			linestream >> m_surfaces[i].tag >> m_surfaces[i].minX >> m_surfaces[i].minY >> m_surfaces[i].minZ >> m_surfaces[i].maxX >> m_surfaces[i].maxY >> m_surfaces[i].maxZ >> numPhysicalTags;

			m_surfaces[i].physicalTags.resize(numPhysicalTags);
			for (size_t tagId=0;tagId<numPhysicalTags;++tagId) { linestream >> m_surfaces[i].physicalTags[tagId]; } linestream >> numBoundingCurves;

			m_surfaces[i].curveTags.resize(numBoundingCurves);
			for (size_t tagId=0;tagId<numBoundingCurves;++tagId) { linestream >> m_surfaces[i].curveTags[tagId]; }
		}
	}
}

void GmshMesh::readNodes(std::ifstream& in)
{
	if (findSection(in, "Nodes"))
	{
		std::stringstream firstLinestream(read_line(in));
		size_t numEntityBlocks, numNodes, minNodeTag, maxNodeTag;
		firstLinestream >> numEntityBlocks >> numNodes >> minNodeTag >> maxNodeTag;
		m_entities.resize(numEntityBlocks);
		for (size_t i=0;i<numEntityBlocks;++i)
		{
			std::stringstream linestream(read_line(in));
			size_t numNodesInBlock;
			linestream >> m_entities[i].dim >> m_entities[i].tag >> m_entities[i].parametric >> numNodesInBlock;

			m_entities[i].nodeTags.resize(numNodesInBlock);
			for (size_t j=0;j<numNodesInBlock;++j)
			{
				std::stringstream _linestream(read_line(in));
				_linestream >> m_entities[i].nodeTags[j];
			}
			m_entities[i].nodes.resize(numNodesInBlock);
			for (size_t j=0;j<numNodesInBlock;++j)
			{
				std::stringstream _linestream(read_line(in));
				double x,y,z;
				_linestream >> x >> y >> z;
				m_nodes.push_back( NodeWorld(x,y) );
				m_nodeTags.push_back(m_entities[i].nodeTags[j]);
				m_entities[i].nodes[j] = &m_nodes.back(); // m_nodes[entity[i].nodeIndex[j]] jth node of the ith entity
			}
		}
	}
}

void GmshMesh::readElements(std::ifstream& in)
{
	if (findSection(in, "Elements"))
	{
		std::stringstream firstLinestream(read_line(in));
		size_t numEntityBlocks, numElements, minElementTag, maxElementTag;
		firstLinestream >> numEntityBlocks >> numElements >> minElementTag >> maxElementTag;
		for (size_t i=0;i<numEntityBlocks;++i)
		{
			std::stringstream linestream(read_line(in));
			int entityDim, entityTag;
			linestream >> entityDim >> entityTag;

			std::vector< Entity >::iterator it  = std::find_if(m_entities.begin(), m_entities.end(), [entityDim, entityTag](const Entity& entity) -> bool
			{
				return entity.dim == entityDim and entity.tag == entityTag;
			});
			if (it != m_entities.end())
			{
				Entity& entity = *it;
				size_t numElementsInBlock;
				linestream >> entity.elementType >> numElementsInBlock;
				entity.elements.resize(numElementsInBlock);
				entity.elementTags.resize(numElementsInBlock);
				for (size_t j=0;j<numElementsInBlock;++j)
				{
					std::stringstream _linestream(read_line(in));
					_linestream >> entity.elementTags[j];
					if (entity.elementType == 3) // 4-node quadrangle
					{
						entity.elements[j].resize(4);
					}
					else if (entity.elementType == 10) // 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
					{
						entity.elements[j].resize(9);
					}
					else if (entity.elementType == 1) // 2-node line (used for boundary elements)
					{
						entity.elements[j].resize(2);
					}
					else if (entity.elementType == 8) // 3-node second order line (2 nodes associated with the vertices and 1 with the edge) (used for boundary elements)
					{
						entity.elements[j].resize(3);
					}
					for (size_t k=0;k<entity.elements[j].size();++k) { _linestream >> entity.elements[j][k]; }
				}
			}
		}
	}
}
