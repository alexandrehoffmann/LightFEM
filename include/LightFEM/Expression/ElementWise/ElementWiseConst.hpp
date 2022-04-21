/*
 * ElementWiseConst.hpp
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

#ifndef ELEMENT_WISE_CONST_HPP
#define ELEMENT_WISE_CONST_HPP

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/ExprType.hpp>

class ElementWiseConst : public ElementWiseFunctionExpression< ExprType::SCALAR, ElementWiseConst >
{
public:
	ElementWiseConst(const Element* element, const Scalar& value) : m_element(element), m_value(value) {}
public:
	inline const Scalar& operator[] (const size_t ) const { return m_value; }
public:
	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return false; }
public:
	inline const Element* getElement() const { return m_element; }
private:
	const Element* m_element;
	const Scalar&  m_value;
};

////////////////////////////////////////////////////////////////////////

class CpxElementWiseConst : public CpxElementWiseFunctionExpression< ExprType::SCALAR, CpxElementWiseConst >
{
public:
	CpxElementWiseConst(const Element* element, const CpxScalar& value) : m_element(element), m_value(value) {}
public:
	inline const CpxScalar& operator[] (const size_t ) const { return m_value; }
public:
	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return false; }
public:
	inline const Element* getElement() const { return m_element; }
private:
	const Element* m_element;
	const CpxScalar&  m_value;
};

#endif // ELEMENT_WISE_CONST_HPP
