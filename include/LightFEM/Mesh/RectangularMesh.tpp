/*
 * RectangularMesh.tpp
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

#ifndef RECTANGULAR_MESH_TPP
#define RECTANGULAR_MESH_TPP

#include <LightFEM/Mesh/RectangularMesh.hpp>

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in, bool endpoint)
{
	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	double delta = (endpoint) ? (end - start) / (num - 1) : (end - start) / (num);

	for(int i=0; i < num; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	return linspaced;
}

#endif // RECTANGULAR_MESH_TPP
