/*
 * RectangularMesh.hpp
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

#ifndef RECTANGULAR_MESH_HPP
#define RECTANGULAR_MESH_HPP

#include <LightFEM/Mesh/Mesh.hpp>

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in, bool endpoint=true);

class RectangularMesh : public Mesh
{
public:
	RectangularMesh(
			const double xmin, const double xmax,
			const double ymin, const double ymax,
			const size_t Nx,
			const size_t Ny);
private:
	inline size_t index2d(const size_t i, const size_t j) const { return i*m_Ny + j; }
private:
	double m_xmin;
	double m_xmax;
	double m_ymin;
	double m_ymax;

	size_t m_Nx;
	size_t m_Ny;
};

#include <LightFEM/Mesh/RectangularMesh.tpp>

#endif // RECTANGULAR_MESH_HPP
