/*
 * FiniteElementFunctionExpression.hpp
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

#ifndef FINITE_ELEMENT_FUNCTION_EXPRESSION_HPP
#define FINITE_ELEMENT_FUNCTION_EXPRESSION_HPP

#include <LightFEM/Expression/Traits.hpp>
#include <LightFEM/Expression/Function/FunctionExpression.hpp>

#include <LightFEM/Analysis/FunctionSpace/FunctionSpace.hpp>
#include <LightFEM/Analysis/FunctionSpace/VectorFunctionSpace.hpp>

template<ExprType Type> struct FSpace {};
template<> struct FSpace< ExprType::SCALAR > { typedef const FunctionSpace Type; typedef Scalar ValueType; typedef Vector GradValueType; };
template<> struct FSpace< ExprType::VECTOR > { typedef const VectorFunctionSpace Type; typedef Vector ValueType; typedef Matrix GradValueType; };

template <ExprType Type, typename Expr>
class FiniteElementFunctionExpression : public FunctionExpression<Type, Expr>
{
public:
	typedef typename Traits< Expr >::ReturnType ReturnType;
	typedef typename Traits< Expr >::GradReturnType GradReturnType;
	typedef typename Traits< Expr >::DiffReturnType DiffReturnType;
	
	typedef typename FSpace< Type >::Type FSpaceType;
// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
public:
	inline ReturnType     operator[] (const size_t e)                   const { return static_cast<Expr const&>(*this)[e]; }
	inline GradReturnType getGrad    (const size_t e)                   const { return static_cast<Expr const&>(*this).getGrad(e); }
	inline DiffReturnType getD       (const size_t dim, const size_t e) const { return static_cast<Expr const&>(*this).getD(dim, e); }
	inline double         getCoef    (const size_t globId)              const { return static_cast<Expr const&>(*this).getCoef(globId); }
public:
	inline bool containsTrial() const { return static_cast<Expr const&>(*this).containsTrial(); }
	inline bool containsTest()  const { return static_cast<Expr const&>(*this).containsTest(); }
public:
	inline const FSpaceType* getFunctionSpace() const { return static_cast<Expr const&>(*this).getFunctionSpace(); }
	inline const Mesh*       getMesh()          const { return static_cast<Expr const&>(*this).getMesh(); }
};

////////////////////////////////////////////////////////////////////////

template <ExprType Type, typename Expr>
class CpxFiniteElementFunctionExpression : public CpxFunctionExpression<Type, Expr>
{
public:
	typedef typename Traits< Expr >::ReturnType ReturnType;
	typedef typename Traits< Expr >::GradReturnType GradReturnType;
	typedef typename Traits< Expr >::GradReturnType DiffReturnType;
	typedef typename FSpace< Type >::Type FSpaceType;
// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
public:
	inline ReturnType            operator[] (const size_t e)                   const { return static_cast<Expr const&>(*this)[e]; }
	inline GradReturnType        getGrad    (const size_t e)                   const { return static_cast<Expr const&>(*this).getGrad(e); }
	inline DiffReturnType        getD       (const size_t dim, const size_t e) const { return static_cast<Expr const&>(*this).getD(dim, e); }
	inline std::complex< double >getCoef    (const size_t globId)              const { return static_cast<Expr const&>(*this).getCoef(globId); }
public:
	inline bool containsTrial() const { return false; }
	inline bool containsTest()  const { return false; }
public:
	inline const FSpaceType* getFunctionSpace() const { return static_cast<Expr const&>(*this).getFunctionSpace(); }
	inline const Mesh*       getMesh()          const { return static_cast<Expr const&>(*this).getMesh(); }
};

#endif // FINITE_ELEMENT_FUNCTION_EXPRESSION_HPP
