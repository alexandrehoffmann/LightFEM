/*
 * GaussLobattoQuadrature.hpp
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

#ifndef GAUSS_LOBATTO_QUADRATURE_HPP
#define GAUSS_LOBATTO_QUADRATURE_HPP

#include <LightFEM/Analysis/Measure/Measure.hpp>
#include <cmath>

class GaussLobattoQuadrature : public Measure
{
public:
	GaussLobattoQuadrature(const Mesh* mesh, const size_t N);
public:
	static std::valarray<long double> getNodes(const size_t N);
private:
	static long double P(const size_t N, const long double x) { return std::legendre(N, x); }
	static long double dP(const size_t N, const long double x);
	static long double d2P(const size_t N, const long double x);
};

#endif // GAUSS_LOBATTO_QUADRATURE_HPP
