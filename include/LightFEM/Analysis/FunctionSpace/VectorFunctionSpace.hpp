/*
 * VectorFunctionSpace.hpp
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

#ifndef VECTOR_FUNCTION_SPACE_HPP
#define VECTOR_FUNCTION_SPACE_HPP

class VectorTrialFunction;
class VectorTestFunction;
class Mesh;

#include <initializer_list>		

#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensor.hpp>
#include <LightFEM/Mesh/Element.hpp>
#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>

class VectorFunctionSpace
{
public:
	VectorFunctionSpace(const std::initializer_list< FunctionSpace* > functionSpaces);
	
	VectorTrialFunction getTrialFunction (const size_t e, const int locId) const; // if id = -1 returns 0
	VectorTestFunction  getTestFunction  (const size_t e, const int locId) const; // if id = -1 returns 0

	inline size_t getGlobalId(const size_t e, const size_t locId) const { return m_globalIds[index2d(e, locId)]; }
	inline size_t getNBasisFunctionPerElement()                   const { return m_nBasisFunctionPerElement; }
	inline size_t getNBasisFunction()                             const { return m_nDofs; }
	
	virtual size_t         getNLocalInterpolationNodes     ()               const = 0; // some functions may not be interpolating
	virtual const NodeRef& getLocalInterpolationNode       (const size_t i) const = 0;
	virtual size_t         getLocalInterpolationFunctionId (const size_t i) const = 0;
	
	inline bool isIdOnBoundary(const size_t globId) const { return m_isIdOnBoundary[globId]; }
	
	inline const std::vector< Vector >& getDiscBasisFunction       (const int locId) const { return (locId == -1) ? m_zero     : m_discBasisFunctions[locId]; }
	inline const std::vector< Matrix >& getDiscGradXiBasisFunction (const int locId) const { return (locId == -1) ? m_gradZero : m_discGradXiBasisFunctions[locId]; }
	
	inline const Mesh* getMesh() const { return m_functionSpaces[0]->getMesh(); }
private:
	inline std::size_t index2d(const std::size_t e, const std::size_t locId) const { return locId + m_nBasisFunctionPerElement*e; }
private:
	std::vector< FunctionSpace* > m_functionSpaces;

	std::vector< NodeRef > m_interpNodes;
	std::vector< size_t >  m_interpFuncLocId;

	std::vector< Vector > m_zero;
	std::vector< Matrix > m_gradZero;
	std::vector< std::vector< Vector > > m_discBasisFunctions;
	std::vector< std::vector< Matrix > > m_discGradXiBasisFunctions;

	std::valarray< size_t > m_globalIds;
	size_t m_nDofs;
	size_t m_nBasisFunctionPerElement;
	
	std::vector< bool > m_isIdOnBoundary;
};

#endif // VECTOR_FUNCTION_SPACE_HPP
