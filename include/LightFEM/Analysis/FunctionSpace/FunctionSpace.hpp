/*
 * FunctionSpace.hpp
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

#ifndef FUNCTION_SPACE_HPP
#define FUNCTION_SPACE_HPP

#include <vector>

class TrialFunction;
class TestFunction;
class Mesh;

#include <LightFEM/Expression/LinAlg/Scalar.hpp>
#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Mesh/Element.hpp>

class FunctionSpace
{
public:
	FunctionSpace(const Mesh* mesh);

	TrialFunction getTrialFunction (const size_t e, const int locId) const; // if id = -1 returns 0
	TestFunction  getTestFunction  (const size_t e, const int locId) const; // if id = -1 returns 0

	virtual size_t getGlobalId(const size_t e, const size_t locId) const = 0;
	virtual size_t getNBasisFunctionPerElement()                   const = 0;
	virtual size_t getNBasisFunction()                             const = 0;

	virtual size_t         getNLocalInterpolationNodes     ()               const = 0; // some functions may not be interpolating
	virtual const NodeRef& getLocalInterpolationNode       (const size_t i) const = 0;
	virtual size_t         getLocalInterpolationFunctionId (const size_t i) const = 0;
	
	virtual double getHMin() const = 0;

	inline bool isIdOnBoundary(const size_t globId) const { return m_isIdOnBoundary[globId]; }
	
	inline const std::vector< Scalar >& getDiscBasisFunction       (const int locId) const { return (locId == -1) ? m_zero     : m_discBasisFunctions[locId]; }
	inline const std::vector< Vector >& getDiscGradXiBasisFunction (const int locId) const { return (locId == -1) ? m_gradZero : m_discGradXiBasisFunctions[locId]; }
	
	inline const Mesh* getMesh() const { return m_mesh; }
protected:
	virtual bool isLocIdOnBoundary(const size_t locId, const Element::Boundary b) const = 0;

	virtual Scalar evalDiscBasisFunction     (const size_t locId, const double xi1, const double xi2) const = 0;
	virtual Vector evalDiscGradXiBasisFunction (const size_t locId, const double xi1, const double xi2) const = 0;
protected:
	void initDiscretization();
	void initIsIdOnBoundary();
protected:
	const Mesh* m_mesh;

	std::vector< Scalar > m_zero;
	std::vector< Vector > m_gradZero;
	std::vector< std::vector< Scalar > > m_discBasisFunctions;
	std::vector< std::vector< Vector > > m_discGradXiBasisFunctions;

	std::vector< bool > m_isIdOnBoundary;
	size_t m_nDofs;
};

#endif // FUNCTION_SPACE_HPP
