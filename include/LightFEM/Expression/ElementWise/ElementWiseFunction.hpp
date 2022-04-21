/*
 * ElementWiseFunction.hpp
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

#ifndef ELEMENT_WISE_FUNCTION_HPP
#define ELEMENT_WISE_FUNCTION_HPP

#include <LightFEM/Expression/ElementWise/ElementWiseFunctionExpression.hpp>
#include <LightFEM/Expression/LinAlg/Scalar.hpp>
#include <LightFEM/Expression/LinAlg/Vector.hpp>
#include <LightFEM/Expression/LinAlg/Matrix.hpp>
#include <LightFEM/Expression/LinAlg/RankTwoTensor.hpp>
#include <LightFEM/Expression/LinAlg/RankFourTensor.hpp>
#include <LightFEM/Expression/LinAlg/ExprType.hpp>

#include <LightFEM/Mesh/Element.hpp>

#include <vector>

template<ExprType Type>
class ElementWiseFunction : public ElementWiseFunctionExpression< Type, ElementWiseFunction<Type> >
{
public:
	typedef typename Traits< ElementWiseFunction<Type> >::ValueType ValueType;
public:
	ElementWiseFunction(const Element* element = nullptr) : m_values(Element::getNxiNd()), m_containsTrial(false), m_containsTest(false), m_element(element)  {}
	template<typename Expr>	ElementWiseFunction(const ElementWiseFunctionExpression< Type, Expr >& expr);
public:
	template<typename Expr> ElementWiseFunction<Type>& operator=  (const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> ElementWiseFunction<Type>& operator+= (const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> ElementWiseFunction<Type>& operator-= (const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr>	ElementWiseFunction<Type>& operator*= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr>	ElementWiseFunction<Type>& operator/= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr);
public:
	void setElement(const Element* element) { m_element = element; }
public:
	inline const ValueType& operator[] (const size_t k) const { return m_values[k]; }
	inline       ValueType& operator[] (const size_t k)       { return m_values[k]; }
public:
	inline bool containsTrial() const { return m_containsTrial; }
	inline bool containsTest()  const { return m_containsTest; }
public:
	inline const Element* getElement() const { return m_element; }
private:
	std::vector< ValueType > m_values;
	bool m_containsTrial;
	bool m_containsTest;
	const Element* m_element;
};

typedef ElementWiseFunction<ExprType::SCALAR>     ElementWiseScalarField;
typedef ElementWiseFunction<ExprType::VECTOR>     ElementWiseVectorField;
typedef ElementWiseFunction<ExprType::MATRIX>     ElementWiseMatrixField;
typedef ElementWiseFunction<ExprType::RK2_TENSOR> ElementWiseRankTwoTensorField;
typedef ElementWiseFunction<ExprType::RK4_TENSOR> ElementWiseRankFourTensorField;

////////////////////////////////////////////////////////////////////////

template<ExprType Type>
class CpxElementWiseFunction : public CpxElementWiseFunctionExpression< Type, CpxElementWiseFunction<Type> >
{
public:
	typedef typename Traits< CpxElementWiseFunction<Type> >::ValueType ValueType;
public:
	CpxElementWiseFunction(const Element* element = nullptr) : m_values(Element::getNxiNd()), m_containsTrial(false), m_containsTest(false), m_element(element)  {}
	template<typename Expr>	CpxElementWiseFunction(const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr>	CpxElementWiseFunction(const CpxElementWiseFunctionExpression< Type, Expr >& expr);
public:
	template<typename Expr> CpxElementWiseFunction<Type>& operator=  (const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator=  (const CpxElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator+= (const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator+= (const CpxElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator-= (const ElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator-= (const CpxElementWiseFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator*= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator*= (const CpxElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator/= (const ElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr> CpxElementWiseFunction<Type>& operator/= (const CpxElementWiseFunctionExpression< ExprType::SCALAR, Expr >& expr);
public:
	void setElement(const Element* element) { m_element = element; }
public:
	inline const ValueType& operator[] (const size_t k) const { return m_values[k]; }
	inline       ValueType& operator[] (const size_t k)       { return m_values[k]; }
public:
	inline bool containsTrial() const { return m_containsTrial; }
	inline bool containsTest()  const { return m_containsTest; }
public:
	inline const Element* getElement() const { return m_element; }
private:
	std::vector< ValueType > m_values;
	bool m_containsTrial;
	bool m_containsTest;
	const Element* m_element;
};

#include <LightFEM/Expression/ElementWise/ElementWiseFunction.tpp>

typedef CpxElementWiseFunction<ExprType::SCALAR>     CpxElementWiseScalarField;
typedef CpxElementWiseFunction<ExprType::VECTOR>     CpxElementWiseVectorField;
typedef CpxElementWiseFunction<ExprType::MATRIX>     CpxElementWiseMatrixField;
typedef CpxElementWiseFunction<ExprType::RK2_TENSOR> CpxElementWiseRankTwoTensorField;
typedef CpxElementWiseFunction<ExprType::RK4_TENSOR> CpxElementWiseRankFourTensorField;

#endif // ELEMENT_WISE_FUNCTION_HPP
