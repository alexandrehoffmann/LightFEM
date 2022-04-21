/*
 * Element.hpp
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

#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <vector>
#include <LightFEM/Mesh/Node.hpp>

#include <LightFEM/Expression/LinAlg/Matrix.hpp>
#include <LightFEM/Expression/LinAlg/Vector.hpp>

class Element
{
public:
	enum class Boundary 
	{
		TOP,
		BOTTTOM,
		LEFT,
		RIGHT
	};
	static constexpr std::array<Boundary, 4> AllBoundaries {Boundary::TOP, Boundary::BOTTTOM, Boundary::LEFT, Boundary::RIGHT};
public:
	static constexpr size_t ORDER = 2;
public:
	Element(const std::array<std::array<NodeWorld*, 2>, 2>& X);
	Element(const std::array<std::array<NodeWorld*, 2>, 2>& X, std::initializer_list< int > domainIds);
	Element(const std::array<std::array<NodeWorld*, 2>, 2>& X, const std::vector< int >& domainIds);
	Element(const std::array<std::array<NodeWorld*, 3>, 3>& X);
	Element(const std::array<std::array<NodeWorld*, 3>, 3>& X, std::initializer_list< int > domainIds);
	Element(const std::array<std::array<NodeWorld*, 3>, 3>& X, const std::vector< int >& domainIds);
	~Element();
public:
	NodeWorld getXworld(const NodeRef& Xi) const;
	Matrix    getJacobian(const NodeRef& Xi) const;
	double    getAbsDetJacobian(const NodeRef& Xi) const; // allows to compute integral on real-world element
	Matrix    getInvJacobian(const NodeRef& Xi) const;  // allows to compute derivatives in real-world coordinates from derivatives on ref element
public:
	NodeWorld getXworld(const double xi1, const double xi2) const { return getXworld(NodeRef(xi1, xi2)); }
	Matrix    getJacobian(const double xi1, const double xi2) const { return getJacobian(NodeRef(xi1, xi2)); }
	Vector    getNormal(const Boundary b, const double t) const;
	
	double    getAbsDetJacobian(const double xi1, const double xi2) const { return getAbsDetJacobian(NodeRef(xi1, xi2)); } // allows to compute integral on real-world element
	double    getDs(const Boundary boundary, const double t) const; // allows to compute boundary integral on real-world element
	Matrix    getInvJacobian(const double xi1, const double xi2) const { return getInvJacobian(NodeRef(xi1, xi2)); }  // allows to compute derivatives in real-world coordinates from derivatives on ref element
public:
	inline        Vector getNormalDisc         (const Boundary b, const size_t k) const { return m_normal[index2d(b, k)]; }
	
	inline        double getAbsDetJacobianDisc (const size_t i, const size_t j)   const { return m_abs_dejJ[index2d(i, j)]; } // allows to compute integral on real-world element
	inline        double getDsDisc             (const Boundary b, const size_t k) const { return m_ds[index2d(b, k)]; } // allows to compute boundary integral on real-world element
	inline const Matrix& getInvJacobianDisc    (const size_t i, const size_t j)   const { return m_invJ[index2d(i, j)]; }  // allows to compute derivatives in real-world coordinates from derivatives on ref element

	inline const Matrix& getInvJacobianDisc (const size_t k)   const { return m_invJ[k]; }  // allows to compute derivatives in real-world coordinates from derivatives on ref element

	inline double getWeight(const size_t i, const size_t j)   const { return get_w(i)*get_w(j)*getAbsDetJacobianDisc(i, j); }
	inline double getWeight(const Boundary b, const size_t k) const { return get_w(k)*getDsDisc(b, k); }
public:
	inline NodeWorld* getNode(const size_t i, const size_t j) const { return m_X[i][j]; }
public:
	inline bool isInDomain(const int domainId) const { return std::find(m_domainIds.cbegin(), m_domainIds.cend(), domainId) !=  m_domainIds.cend(); }
private:
	void initTransport();
	void init(const std::array<std::array<NodeWorld*, 2>, 2>& X);
	void init();
public:
	std::pair< bool, NodeRef > getXRef(const double x, const double y, const double epsilon = 1.0e-6) const { return getXRef(NodeWorld(x, y), epsilon); }
	std::pair< bool, NodeRef > getXRef(const NodeWorld& X, const double epsilon = 1.0e-6) const;
private:
	std::array<std::array<NodeWorld*, ORDER+1>, ORDER+1> m_X;
	std::array<std::array<NodeWorld, ORDER+1>, ORDER+1>  m_c; // T(Xi) = sum_kl c_kl Xi_1^k Xi_2^l

	std::valarray<double> m_abs_dejJ;
	std::valarray<double> m_ds;
	std::vector< Matrix > m_invJ;
	
	std::vector< Vector > m_normal;

	std::vector< int > m_domainIds;

	std::array<std::array<bool, ORDER+1>, ORDER+1>  m_deleteX; // T(Xi) = sum_kl c_kl Xi_1^k Xi_2^l
public:
	static size_t getNxi() { return xi.size(); }
	static size_t getNxiNd() { return getNxi()*getNxi(); }
	static double get_xi (const size_t i) { return xi[i]; }
	static double get_w  (const size_t i) { return weights[i]; }
public:
	static size_t index2d(const size_t i, const size_t j) { return i*xi.size() + j; }
	static size_t index2d(const Boundary b, const size_t k) { return int(b)*xi.size() + k; }
private:
	static const std::valarray<double> xi;
	static const std::valarray<double> weights;
private:
	static const std::array<std::array<NodeRef, ORDER+1>, ORDER+1> Xi;
};

#endif // ELEMENT_HPP
