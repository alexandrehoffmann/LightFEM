/*
 * Const.hpp
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

#ifndef FUNCTION_CONST_HPP
#define FUNCTION_CONST_HPP

#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/ExprType.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseConst.hpp>

class Const : public FunctionExpression< ExprType::SCALAR, Const >
{
public:
	typedef typename Traits< Const >::ReturnType ReturnType;
public:
	Const(const Mesh* mesh, const double value) : m_value(value), m_mesh(mesh) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_mesh->getElem(e), m_value); }

	inline const Scalar& getValue() const { return m_value; }
public:
	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return false; }
public:
	inline const Mesh* getMesh() const { return m_mesh; }
private:
	const Scalar m_value;
	const Mesh*  m_mesh;
};

////////////////////////////////////////////////////////////////////////

class CpxConst : public CpxFunctionExpression< ExprType::SCALAR, CpxConst >
{
public:
	typedef typename Traits< CpxConst >::ReturnType ReturnType;
public:
	CpxConst(const Mesh* mesh, const std::complex< double > value) : m_value(value), m_mesh(mesh) {}
public:
	inline ReturnType operator[] (const size_t e) const { return ReturnType(m_mesh->getElem(e), m_value); }

	inline const CpxScalar& getValue() const { return m_value; }
public:
	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return false; }
public:
	inline const Mesh* getMesh() const { return m_mesh; }
private:
	const CpxScalar m_value;
	const Mesh*     m_mesh;
};

////////////////////////////////////////////////////////////////////////

class ConstFactory
{
public:
	ConstFactory(const Mesh* mesh) : m_mesh(mesh) {}
public:
	Const    operator() (const double value) const { return Const(m_mesh, value); }
	CpxConst operator() (const std::complex< double > value) const { return CpxConst(m_mesh, value); }
	CpxConst operator() (const double real_value, const double imag_value) const { return CpxConst(m_mesh, std::complex< double >(real_value, imag_value)); }
private:
	const Mesh* m_mesh; 
};

#endif // FUNCTION_CONST_HPP
