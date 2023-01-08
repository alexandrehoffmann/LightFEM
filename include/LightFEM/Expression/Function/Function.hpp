/*
 * Function.hpp
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

#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <LightFEM/Expression/Function/FunctionExpression.hpp>
#include <LightFEM/Expression/ElementWise/ElementWiseFunction.hpp>
#include <LightFEM/Expression/LinAlg/ExprType.hpp>
#include <LightFEM/Mesh/Mesh.hpp>

#include <vector>

template<ExprType Type>
class Function : public FunctionExpression< Type, Function<Type> >
{
public:
	typedef typename Traits< Function<Type> >::ReturnType ReturnType;
public:
	Function(const Mesh* mesh);
	template<typename FunctorType> Function(const Mesh* mesh, FunctorType f);
	template<typename Expr>	       Function(const FunctionExpression< Type, Expr >& expr);
public:
	template<typename Expr> Function<Type>& operator=  (const FunctionExpression< Type, Expr >& expr);
	template<typename Expr> Function<Type>& operator+= (const FunctionExpression< Type, Expr >& expr);
	template<typename Expr> Function<Type>& operator-= (const FunctionExpression< Type, Expr >& expr);
	template<typename Expr>	Function<Type>& operator*= (const FunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr>	Function<Type>& operator/= (const FunctionExpression< ExprType::SCALAR, Expr >& expr);
public:
	inline const ElementWiseFunction<Type>& operator[] (const size_t e) const { return m_values[e]; }
	inline       ElementWiseFunction<Type>& operator[] (const size_t e)       { return m_values[e]; }
	
	inline const ElementWiseFunction<Type>* data() const { return m_values.data(); }
	inline       ElementWiseFunction<Type>* data()       { return m_values.data(); }
public:
	inline bool containsTrial() const { return m_containsTrial; }
	inline bool containsTest()  const { return m_containsTest; }
public:
	inline const Mesh* getMesh() const { return m_mesh; }
private:
	std::vector< ElementWiseFunction<Type> > m_values;
	
	bool m_containsTrial;
	bool m_containsTest;
	
	const Mesh* m_mesh;
};

typedef Function<ExprType::SCALAR>     ScalarField;
typedef Function<ExprType::VECTOR>     VectorField;
typedef Function<ExprType::MATRIX>     MatrixField;
typedef Function<ExprType::RK2_TENSOR> RankTwoTensorField;
typedef Function<ExprType::RK4_TENSOR> RankFourTensorField;

////////////////////////////////////////////////////////////////////////

template<ExprType Type>
class CpxFunction : public CpxFunctionExpression< Type, CpxFunction<Type> >
{
public:
	typedef typename Traits< CpxFunction<Type> >::ReturnType ReturnType;
public:
	CpxFunction(const Mesh* mesh);
	template<typename CpxFunctorType> CpxFunction(const Mesh* mesh, CpxFunctorType f);
	template<typename Expr>	          CpxFunction(const FunctionExpression< Type, Expr >& expr);
	template<typename Expr>	          CpxFunction(const CpxFunctionExpression< Type, Expr >& expr);
public:
	template<typename Expr> CpxFunction<Type>& operator=  (const FunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator=  (const CpxFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator+= (const FunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator+= (const CpxFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator-= (const FunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator-= (const CpxFunctionExpression< Type, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator*= (const FunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator*= (const CpxFunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator/= (const FunctionExpression< ExprType::SCALAR, Expr >& expr);
	template<typename Expr> CpxFunction<Type>& operator/= (const CpxFunctionExpression< ExprType::SCALAR, Expr >& expr);
public:
	inline const CpxElementWiseFunction<Type>& operator[] (const size_t e) const { return m_values[e]; }
	inline       CpxElementWiseFunction<Type>& operator[] (const size_t e)       { return m_values[e]; }
	
	inline const CpxElementWiseFunction<Type>* data() const { return m_values.data(); }
	inline       CpxElementWiseFunction<Type>* data()       { return m_values.data(); }
public:
	inline bool containsTrial() const { return m_containsTrial; }
	inline bool containsTest()  const { return m_containsTest; }
public:
	inline const Mesh* getMesh() const { return m_mesh; }
private:
	std::vector< CpxElementWiseFunction<Type> > m_values;
	bool m_containsTrial;
	bool m_containsTest;
	const Mesh* m_mesh;
};

#include <LightFEM/Expression/Function/Function.tpp>

typedef CpxFunction<ExprType::SCALAR>     CpxScalarField;
typedef CpxFunction<ExprType::VECTOR>     CpxVectorField;
typedef CpxFunction<ExprType::MATRIX>     CpxMatrixField;
typedef CpxFunction<ExprType::RK2_TENSOR> CpxRankTwoTensorField;
typedef CpxFunction<ExprType::RK4_TENSOR> CpxRankFourTensorField;

#endif // FUNCTION_HPP
