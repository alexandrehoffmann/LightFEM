/*
 * Node.hpp
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

#ifndef NODE_HPP
#define NODE_HPP

#include <LightFEM/Expression/LinAlg/Vector.hpp>

struct NodeWorld
{
	NodeWorld(double xx = 0, double yy = 0) : x(xx), y(yy) {}
	NodeWorld(const NodeWorld& other) : x(other.x), y(other.y) {}
	double x;
	double y;

	NodeWorld& operator= (const NodeWorld& other) { x = other.x; y = other.y; return *this; }

	operator Vector() const { return Vector({x, y}); }
};

struct NodeRef
{
	NodeRef(double xxi1 = 0, double xxi2 = 0) : xi1(xxi1), xi2(xxi2) {}
	NodeRef(const NodeRef& other) : xi1(other.xi1), xi2(other.xi2) {}
	double xi1;
	double xi2;

	NodeRef& operator= (const NodeRef& other) { xi1 = other.xi1; xi2 = other.xi2; return *this; }

	operator Vector() const { return Vector({xi1, xi2}); }
};

inline double norm(const NodeWorld& node) { return sqrt( node.x*node.x + node.y*node.y ); }
inline double norm(const NodeRef& node) { return sqrt( node.xi1*node.xi1 + node.xi2*node.xi2 ); }

inline double dist(const NodeWorld& lhs, const NodeWorld& rhs) { return sqrt( (lhs.x - rhs.x)*(lhs.x - rhs.x) + (lhs.y - rhs.y)*(lhs.y - rhs.y) ); }
inline double dist(const NodeRef& lhs, const NodeRef& rhs) { return sqrt( (lhs.xi1 - rhs.xi1)*(lhs.xi1 - rhs.xi1) + (lhs.xi2 - rhs.xi2)*(lhs.xi2 - rhs.xi2) ); }

inline bool eq(const NodeWorld& lhs, const NodeWorld& rhs, const double eps = 1.0e-11) { return dist(lhs, rhs) <= eps*std::min(norm(lhs), norm(rhs)); }
inline bool eq(const NodeRef& lhs, const NodeRef& rhs, const double eps = 1.0e-11)     { return dist(lhs, rhs) <= eps*std::min(norm(lhs), norm(rhs)); }

#endif // NODE_HPP
