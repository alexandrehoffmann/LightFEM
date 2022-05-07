/*
 * PnFunctionSpace.hpp
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

#ifndef PN_FUNCTION_SPACE_HPP
#define PN_FUNCTION_SPACE_HPP

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>

class PnFunctionSpace: public FunctionSpace
{
public:
	PnFunctionSpace(const Mesh* mesh, const size_t order);

	virtual size_t getGlobalId(const size_t e, const size_t locId) const { return m_globalIds[index2d(e, locId)]; }
	virtual size_t getNBasisFunctionPerElement()                   const { return m_basisFunctions.size(); }
	virtual size_t getNBasisFunction()                             const { return m_nDofs; }

	virtual size_t         getNLocalInterpolationNodes     ()               const { return m_interpNodes.size(); }
	virtual const NodeRef& getLocalInterpolationNode       (const size_t i) const { return m_interpNodes[i]; }
	virtual size_t         getLocalInterpolationFunctionId (const size_t i) const { return i; } // all our local functions are interpolation functions
	
	virtual double getHMin() const;
protected:
	virtual bool isLocIdOnBoundary(const size_t locId, const Element::Boundary b) const { return m_isLocIdOnBoundary[index2d(locId, b)]; }

	virtual Scalar evalDiscBasisFunction      (const size_t locId, const double xi1, const double xi2) const { return m_basisFunctions[locId](xi1, xi2); }
	virtual Vector evalDiscGradXiBasisFunction(const size_t locId, const double xi1, const double xi2) const { return m_gradXibasisFunctions[locId](xi1, xi2); }
private:
	inline size_t index2d(const size_t e, const size_t locId) const { return locId + m_basisFunctions.size()*e; }
	inline size_t index2d(const size_t locId, const Element::Boundary b) const { return 4*locId + int(b); }

	double l(const size_t i, const double t) const;
	double dl(const size_t i, const double t) const;
private:
	size_t m_order;

	std::valarray< long double > m_interpNodes1D;
	std::vector< NodeRef > m_interpNodes;

	std::vector< std::function<Scalar(const double, const double)> > m_basisFunctions;
	std::vector< std::function<Vector(const double, const double)> > m_gradXibasisFunctions;

	std::valarray< size_t > m_globalIds;
	size_t m_nDofs;

	std::valarray< bool > m_isLocIdOnBoundary;
};

#endif // PN_FUNCTION_SPACE_HPP
