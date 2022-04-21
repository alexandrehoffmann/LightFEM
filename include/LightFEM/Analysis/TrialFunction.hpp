/*
 * TrialFunction.hpp
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

#ifndef TRIAL_FUNCTION_HPP
#define TRIAL_FUNCTION_HPP

#include <LightFEM/Expression/Traits.hpp>
#include <LightFEM/Expression/LinAlg/Scalar.hpp>
#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Expression/LinAlg/Matrix.hpp>
#include <LightFEM/Expression/LinAlg/operators.hpp>

class FunctionSpace;
class TrialFunction;
class GradTrialFunction;
template<> struct Traits< TrialFunction > { typedef const Scalar& ReturnType; };
template<> struct Traits< GradTrialFunction > 
{ 
	typedef UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Matrix&> TransposeJacobianType;
	typedef BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, TransposeJacobianType, ExprType::VECTOR, const Vector&> ReturnType;
};

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>

class TrialFunction : public ElementWiseFunctionExpression<ExprType::SCALAR, TrialFunction>
{
public:
	typedef Traits<TrialFunction>::ReturnType ReturnType;
public:
	TrialFunction(const FunctionSpace* functionSpace, const Element* element, const int locId);
	TrialFunction(const TrialFunction& other);

	inline ReturnType operator[] (const size_t k) const { return m_discBasisFunction[k]; }

	inline bool containsTrial() const { return true; }
	inline bool containsTest()  const { return false; }

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

class GradTrialFunction : public ElementWiseFunctionExpression<ExprType::VECTOR, GradTrialFunction>
{
public:
	typedef Traits<GradTrialFunction>::ReturnType ReturnType;
public:
	GradTrialFunction(const TrialFunction& trialFunction);

	inline ReturnType operator[] (const size_t k) const { return transpose(m_elem->getInvJacobianDisc(k))*m_discGradXiBasisFunctions[k]; }

	inline bool containsTrial() const { return true; }
	inline bool containsTest()  const { return false; }

	inline const FunctionSpace* getFunctionSpace() const { return m_functionSpace; }
	inline const Element*       getElement()       const { return m_elem; }
	inline       int            getLocId()         const { return m_locId; }
private:
	const FunctionSpace* m_functionSpace;
	const Element*       m_elem;
	const int            m_locId;

	const std::vector< Vector >& m_discGradXiBasisFunctions;
};

#endif // TRIAL_FUNCTION_HPP
