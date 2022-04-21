/*
 * ElementWiseFunctionExpression.hpp
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

#ifndef ELEMENT_WISE_FUNCTION_EXPRESSION_HPP
#define ELEMENT_WISE_FUNCTION_EXPRESSION_HPP

#include <LightFEM/Expression/Traits.hpp>
#include <LightFEM/Mesh/Element.hpp>

template <ExprType Type, typename Expr>
class ElementWiseFunctionExpression
{
public:
	typedef typename Traits< Expr >::ReturnType ReturnType;
// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
public:
	inline ReturnType operator[] (const size_t k) const { return static_cast<Expr const&>(*this)[k]; }
public:
	inline bool containsTrial() const { return static_cast<Expr const&>(*this).containsTrial(); }
	inline bool containsTest()  const { return static_cast<Expr const&>(*this).containsTest(); }
public:
	inline const Element* getElement() const { return static_cast<Expr const&>(*this).getElement(); }
};

////////////////////////////////////////////////////////////////////////

template <ExprType Type, typename Expr>
class CpxElementWiseFunctionExpression
{
public:
	typedef typename Traits< Expr >::ReturnType ReturnType;
// Delegation to the actual expression type. This avoids dynamic polymorphism (a.k.a. virtual functions in C++)
public:
	inline ReturnType operator[] (const size_t k) const { return static_cast<Expr const&>(*this)[k]; }
public:
	inline bool containsTrial() const { return static_cast<Expr const&>(*this).containsTrial(); }
	inline bool containsTest()  const { return static_cast<Expr const&>(*this).containsTest(); }
public:
	inline const Element* getElement() const { return static_cast<Expr const&>(*this).getElement(); }
};

#endif // ELEMENT_WISE_FUNCTION_EXPRESSION_HPP
