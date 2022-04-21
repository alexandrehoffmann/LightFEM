/*
 * TestFunction.hpp
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

#ifndef TEST_FUNCTION_HPP
#define TEST_FUNCTION_HPP

#include <LightFEM/Expression/Traits.hpp>
#include <LightFEM/Expression/LinAlg/Scalar.hpp>
#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Expression/LinAlg/Matrix.hpp>
#include <LightFEM/Expression/LinAlg/operators.hpp>

class FunctionSpace;
class TestFunction;
class GradTestFunction;
template<> struct Traits< TestFunction > { typedef const Scalar& ReturnType; };
template<> struct Traits< GradTestFunction > 
{ 
	typedef UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Matrix&> TransposeJacobianType;
	typedef BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, TransposeJacobianType, ExprType::VECTOR, const Vector&> ReturnType;
};

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>

class TestFunction : public ElementWiseFunctionExpression<ExprType::SCALAR, TestFunction>
{
public:
	typedef Traits<TestFunction>::ReturnType ReturnType;
public:
	TestFunction(const FunctionSpace* functionSpace, const Element* element, const int locId);
	TestFunction(const TestFunction& other);

	inline ReturnType operator[] (const size_t k) const { return m_discBasisFunction[k]; }

	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return true; }

	inline const FunctionSpace* getFunctionSpace() const { return m_functionSpace; }
	inline const Element*       getElement()       const { return m_elem; }
	inline       int            getLocId()         const { return m_locId; }
private:
	const FunctionSpace* m_functionSpace;
	const Element*       m_elem;
	const int            m_locId;

	const std::vector< Scalar >& m_discBasisFunction;
};

////////////////////////////////////////////////////////////////////////

class GradTestFunction : public ElementWiseFunctionExpression<ExprType::VECTOR, GradTestFunction>
{
public:
	typedef Traits<GradTestFunction>::ReturnType ReturnType;
public:
	GradTestFunction(const TestFunction& testFunction);

	inline ReturnType operator[] (const size_t k) const { return transpose(m_elem->getInvJacobianDisc(k))*m_discGradXiBasisFunctions[k]; }

	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return true; }

	inline const FunctionSpace* getFunctionSpace() const { return m_functionSpace; }
	inline const Element*       getElement()       const { return m_elem; }
	inline       int            getLocId()         const { return m_locId; }
private:
	const FunctionSpace* m_functionSpace;
	const Element*       m_elem;
	const int            m_locId;

	const std::vector< Vector >& m_discGradXiBasisFunctions;
};

#endif // TEST_FUNCTION_HPP
