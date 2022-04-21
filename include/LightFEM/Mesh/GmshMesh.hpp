/*
 * GmshMesh.hpp
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

#ifndef GMSH_MESH_HPP
#define GMSH_MESH_HPP

#include <LightFEM/Mesh/Mesh.hpp>

#include <vector>
#include <map>

std::string read_line(std::ifstream& filestream);
bool findSection(std::ifstream& filestream, const std::string& searchname);

class GmshMesh : public Mesh
{
public:
	struct Point
	{
		int tag;
		double x;
		double y;
		double z;
		std::vector< int > physicalTags;
	};
	struct Curve
	{
		int tag;
		double minX;
		double minY;
		double minZ;
		double maxX;
		double maxY;
		double maxZ;
		std::vector< int > physicalTags;
		std::vector< int > pointTags;
	};
	struct Surface
	{
		int tag;
		double minX;
		double minY;
		double minZ;
		double maxX;
		double maxY;
		double maxZ;
		std::vector< int > physicalTags;
		std::vector< int > curveTags;
	};
	struct Entity
	{
		int dim;
		int tag;
		bool parametric;
		std::vector< int > nodeTags;
		std::vector< NodeWorld* > nodes;
		int elementType;
		std::vector< std::vector< size_t > > elements;
		std::vector< int > elementTags;
	};

public:
	GmshMesh(const std::string& filename);
private:
	void readPhysicalNames(std::ifstream& in);
	void readEntities(std::ifstream& in);
	void readNodes(std::ifstream& in);
	void readElements(std::ifstream& in);
private:
	void initElements();
	void initBoundaryElements();
private:
	std::vector< int > m_domainsTag;
	std::vector< int > m_nodeTags;

	std::vector< int > m_boundariesTag;
	
	std::vector< Point > m_points;
	std::vector< Curve > m_curves;
	std::vector< Surface > m_surfaces;
	std::vector< Entity > m_entities;
};

#endif // GMSH_MESH_HPP
