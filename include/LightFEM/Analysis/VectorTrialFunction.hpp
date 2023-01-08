/*
 * VectorTrialFunction.hpp
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

#ifndef VECTOR_TRIAL_FUNCTION_HPP
#define VECTOR_TRIAL_FUNCTION_HPP

#include <LightFEM/Expression/Traits.hpp>

#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Expression/LinAlg/Matrix.hpp>
#include <LightFEM/Expression/LinAlg/operators.hpp>

class VectorFunctionSpace;
class VectorTrialFunction;
class GradVectorTrialFunction;
template<> struct Traits< VectorTrialFunction > { typedef const Vector& ReturnType; };
template<> struct Traits< GradVectorTrialFunction > 
{
	typedef UnaryExpression<UnaryOp::TRANSPOSE, ExprType::MATRIX, const Matrix&> TransposeJacobianType;
	typedef BinaryExpression<BinaryOp::PROD, ExprType::MATRIX, TransposeJacobianType, ExprType::MATRIX, const Matrix&> ReturnType;
};

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>

class VectorTrialFunction : public ElementWiseFunctionExpression<ExprType::VECTOR, VectorTrialFunction>
{
public:
	typedef Traits<VectorTrialFunction>::ReturnType ReturnType;
public:
	VectorTrialFunction(const VectorFunctionSpace* functionSpace, const Element* element, const int locId);
	VectorTrialFunction(const VectorTrialFunction& other);

	inline ReturnType operator[] (const size_t k) const { return m_discBasisFunction[k]; }

	inline bool containsTrial() const { return true; }
	inline bool containsTest()  const { return false; }

	inline const VectorFunctionSpace* getFunctionSpace() const { return m_functionSpace; }
	inline const Mesh*                getMesh()          const { return m_mesh; }
	inline const Element*             getElement()       const { return m_elem; }
	inline       int                  getLocId()         const { return m_locId; }
private:
	const VectorFunctionSpace* m_functionSpace;
	const Mesh*                m_mesh;
	const Element*             m_elem;
	const int                  m_locId;

	const std::vector< Vector >& m_discBasisFunction;
};

////////////////////////////////////////////////////////////////////////

class GradVectorTrialFunction : public ElementWiseFunctionExpression<ExprType::MATRIX, GradVectorTrialFunction>
{
public:
	typedef Traits<GradVectorTrialFunction>::ReturnType ReturnType;
public:
	GradVectorTrialFunction(const VectorTrialFunction& testFunction);

	inline ReturnType operator[] (const size_t k) const { return transpose(m_elem->getInvJacobianDisc(k))*m_discGradXiBasisFunctions[k]; }
	
	inline bool containsTrial() const { return true; }
	inline bool containsTest()  const { return false; }

	inline const VectorFunctionSpace* getFunctionSpace() const { return m_functionSpace; }
	inline const Mesh*                getMesh()          const { return m_mesh; }
	inline const Element*             getElement()       const { return m_elem; }
	inline       int                  getLocId()         const { return m_locId; }
private:
	const VectorFunctionSpace* m_functionSpace;
	const Mesh*                m_mesh;
	const Element*             m_elem;
	const int                  m_locId;

	const std::vector< Matrix >& m_discGradXiBasisFunctions;
};

#endif // VECTOR_TRIAL_FUNCTION_HPP
